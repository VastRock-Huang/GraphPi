#include "../include/graphmpi.h"
#include "../include/common.h"
#include <mpi.h>
#include <cstring>
#include <omp.h>
#include <cstdio>

// 工厂模式获取单例Graphmpi对象
Graphmpi& Graphmpi::getinstance() {
    static Graphmpi gm;
    return gm;
}

//! 初始化函数
void Graphmpi::init(int _threadcnt, Graph* _graph, const Schedule& schedule) {
    initialized = true;
    threadcnt = _threadcnt;
    graph = _graph;
    int provided;
    // 初始化MPI执行环境
    //MPI_THREAD_FUNNELED:若进程为多线程,则只有调用MPI_Init_thread的线程才会进行MPI调用
    //arg4:记录程序可用的线程支持级别
    MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
    // 在主从模型框架中,MPI_Comm_size和MPI_Comm_rank可用于确定通信器的各个进程的角色
    // 返回通信器(communicator)中涉及的进程数
    //MPI_COMM_WORLD:表示可用进程的总数
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    // 确定通信器中调用进程的级别
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    // 进行MPI进程之间的同步
    MPI_Barrier(MPI_COMM_WORLD);
    starttime = get_wall_time();
    idlethreadcnt = 0;
    for (int i = 0; i < threadcnt; i++) {
        data[i] = new unsigned int[MESSAGE_SIZE];
        // 初始时每个线程都被上锁
        lock[i].test_and_set(); // 将锁原子地设置标志为true并获得其先前值
    }
    const int CHUNK_CONST = 70;
    omp_chunk_size = std::max(int((long long)(graph->v_cnt) * CHUNK_CONST / graph->e_cnt), 8);
    mpi_chunk_size = (threadcnt - 1) * omp_chunk_size;
    skip_flag = ~schedule.get_restrict_last(1);
    printf("mpi_csize = %d, omp_csize = %d\n", mpi_chunk_size, omp_chunk_size);
    fflush(stdout);
}

//! 主线程运行
long long Graphmpi::runmajor() {
    long long tot_ans = 0;      // 总匹配结果
    const int IDLE = 2, END = 3, OVERWORK = 4, REPORT = 5, SERVER = 0;
    const int ROLL_SIZE = 327680;   // 缓冲区数组大小
    static unsigned int recv[MESSAGE_SIZE];     // 接收数据
    static unsigned int local_data[MESSAGE_SIZE];   // 进程局部数据,记录着分发的结点数据
    MPI_Request sendrqst, recvrqst;
    MPI_Status status;
    // MPI_Irecv():非阻塞调用分配通信请求对象并将其与请求句柄(参数请求)关联
    //arg1:接收数据缓冲区
    //arg3:接收数据的数据类型 MPI_UNSIGNED:无符号整型
    //arg4:来源等级
    //arg6:通信器句柄 所有线程
    //arg7:记录接受请求句柄
    MPI_Irecv(recv, sizeof(recv), MPI_UNSIGNED, MPI_ANY_SOURCE, 0,
              MPI_COMM_WORLD, &recvrqst);
    int idlenodecnt = 0;    // 空闲进程(结点)数 - 只对主进程有效
    int cur = 0;    // 当前数据图结点
    unsigned int edgel, edger;  // 用于记录某个结点的边的CSR左右边界
    unsigned int *send;
    std::queue<int> workq;
    // 获取0结点的边的范围
    graph->get_edge_index(0, edgel, edger);


    // 分发结点,设置发送数据
    auto get_send = [&](unsigned int *send) {
        auto next_cur = [&]() {
            cur++;  // 结点数加1,即遍历下个结点
            // 将所有结点分为100份,输出当前已经完成的部分
            int k = std::max(graph->v_cnt / 100, 1);
            if (cur % k == 0) {
                printf("nearly %d out of 100 task assigned, time = %f\n", cur / k, get_wall_time() - starttime);
                fflush(stdout);
            }
            if (cur < graph->v_cnt) {
                // 设置结点的边的范围
                graph->get_edge_index(cur, edgel, edger);
            }
        };
        send[0] = OVERWORK;
        // 跳过最外层的限制条件不满足的结点
        //skip_flag为true表示有(0,1)限制条件
        //graph->edge[edgel]>=cur表示该结点的邻结点都更大,不满足限制条件
        if (cur < graph->v_cnt && skip_flag && graph->edge[edgel] >= cur) next_cur();

        if (cur == graph->v_cnt) {  // 遍历完所有结点
            send[1] = -1;
        }
        else {  // 未遍历完所有结点
            send[1] = cur;  // 当前结点
            send[2] = edgel;    // 该结点的边的左边界
            // 若该结点的边数超过MPI块大小
            if (edger - edgel > mpi_chunk_size) {
                send[3] = edgel += mpi_chunk_size;  // 设置部分边的范围
            }
            else {
                send[3] = edger;
                next_cur(); // 若所有边都发送完遍历下个结点
            }
        }
    };
    int buft = 0;   // 缓冲区发送数据指针
    int bufw = 0;   // 缓冲区写数据指针
    // 一个循环二维数组作为缓冲区数组
    static unsigned int buf[ROLL_SIZE][MESSAGE_SIZE + 1];
    // 设置send缓冲区,移动写数据指针
    auto roll_send = [&]() {
        send = buf[bufw] + 1;
        bufw = (bufw + 1) % ROLL_SIZE;
    };
    // 发送数据或获取发送数据状态,并移动发送指针
    auto update_send = [&]() {
        // 发送指针不等于写入指针,即有数据未发送
        if (buft != bufw) {
            static bool ini_flag = false;   // 初始化标识,标志是否发送数据
            bool flag;
            if (ini_flag) {
                int testflag;
                // MPI_Test():完成特定发送或接收的测试
                //arg1:发送请求句柄
                //arg2:记录操作是否完成
                //arg3:状态对象 忽略
                MPI_Test(&sendrqst, &testflag, MPI_STATUS_IGNORE);
                flag = ini_flag;
            }
            else {
                ini_flag = flag = true;
            }
            if (flag) {
                // 发送数据到0号进程
                // MPI_Isend():启动标准模式、非阻塞发送
                //arg1:发送数据的缓冲区地址
                //arg3:发送数据的数据类型
                //arg4:目的地的级别
                //arg6:通信器句柄
                MPI_Isend(buf[buft] + 1, MESSAGE_SIZE, MPI_UNSIGNED, buf[buft][0], 0, MPI_COMM_WORLD, &sendrqst);
                buft = (buft + 1) % ROLL_SIZE;      // 发送指针后移
            }
        }
    };
    // 非主进程设置为空闲用于获取新数据
    auto get_new_local_data = [&]() {
        if (my_rank) {  // 非0号进程,数据状态设置为IDLE
            roll_send();    // 设置send数组,缓冲区写指针后移
            send[-1] = SERVER;  // 等价于buf[bufw]=SERVER
            send[0] = IDLE;
        }
        else {  // 0号进程(主进程)直接设置分发数据,将数据状态设置为OVERWOKR
            get_send(local_data);
        }
    };


    get_new_local_data();
    for (;;) {
        // 发送数据
        //初始:对于0号进程,其buft==bufw不会发送数据;其它进程发送空闲状态
        //后来:进程为OVERWORK后也会持续请求数据
        update_send();
        int testflag = 0;
        // 接受数据判断
        MPI_Test(&recvrqst, &testflag, &status);
        // 若接收到数据
        if (testflag) {
            if (recv[0] == IDLE) {  // 若收到数据状态为IDLE--主进程接收
                roll_send();    // 设置send数组,缓冲区写指针后移
                send[-1] = status.MPI_SOURCE;   // 记录接收数据的发送进程号
                get_send(send); // 设置分发结点数据,数据状态设置为OVERWORK
            }
            else if (recv[0] == OVERWORK) { // 若收到数据状态为OVERWORK--非主进程接收
                // 拷贝接收数据到局部数据
                memcpy(local_data, recv, sizeof(recv));
            }
            else if (recv[0] == END) {  // 若收到数据状态为END(全部进程空闲)-非主进程接收
                // 将非主进程的本次匹配结果更新到本进程的总结果,并退出匹配
                tot_ans = (((long long)(recv[1]) << 32) | (unsigned)recv[2]);
                break;
            }
            else if (recv[0] == REPORT) {   // 若收到数据为REPORT-主进程接收
                // 更新本次匹配结果到主进程的总结果
                tot_ans += (((long long)(recv[1]) << 32) | (unsigned)recv[2]);
                idlenodecnt++;  // 空闲结点+1
            }
            // 非阻塞接收数据
            MPI_Irecv(recv, sizeof(recv), MPI_UNSIGNED, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recvrqst);
        }
        // 若队列不为空
        if (!idleq.empty()) {
            if (local_data[1] == -1) {  // 若遍历完所有的结点
                int tmpthread = idleq.front_and_pop();  // 出队一个线程号
                data[tmpthread][1] = -1;
                lock[tmpthread].clear();
            }
            else if (local_data[2] != local_data[3]) {  // 若局部数据有结点边的数据
                int tmpthread = idleq.front_and_pop();  // 从本进程获取一个空闲线程
                memcpy(data[tmpthread], local_data, sizeof(local_data));    // 将数据拷贝给线程
                // 这里改变local_data[2]是防止local_data被重复分发给两个线程
                data[tmpthread][3] = local_data[2] = std::min(local_data[2] + omp_chunk_size, local_data[3]);
                // 释放线程锁
                //线程锁在初始化或者图匹配时被锁定,直到有数据分发到该线程的数据数组时释放
                //由于通过队列缓冲区记录空闲线程,因此不会出现线程在匹配时被分配新数据的情况
                lock[tmpthread].clear();
                // 将数据状态设置为空闲,用于获取新数据
                if (local_data[2] == local_data[3]) get_new_local_data();
            }
        }
        // 当前进程的工作线程全部空闲(主线程不用于匹配),即该结点的匹配完成
        if (idlethreadcnt == threadcnt - 1) {
            idlethreadcnt = -1;
            if (my_rank) {  // 若为非主进程
                roll_send();    // 移动写数据指针
                // 设置本次任务的匹配结果
                send[-1] = SERVER;
                send[0] = REPORT;
                send[1] = node_ans >> 32;
                send[2] = node_ans;
            }
            else {  // 若为主进程
                idlenodecnt++;  // 空闲结点+1
                tot_ans += node_ans;    // 更新结点结果到总结果
            }
        }
        // 若所有进程都空闲--主进程
        if (idlenodecnt == comm_sz) {
            // 遍历所有非主进程
            for (int i = 1; i < comm_sz; i++) {
                roll_send();        // 移动写数据指针
                // 设置发送数据
                send[-1] = i;   // 发送给本进程
                send[0] = END;
                send[1] = tot_ans >> 32;
                send[2] = tot_ans;
            }
            // 主进程会发送数据给所有非主线程
            for (; buft != bufw; update_send());

            if (comm_sz > 1) {  // 若有非主进程
                // 等待发送数据完成
                MPI_Wait(&sendrqst, MPI_STATUS_IGNORE);
            }
            break;
        }
    }
    return tot_ans;
}

Graphmpi::Graphmpi() {
    initialized = false;    
}

Graphmpi::~Graphmpi() {
    if (initialized) {
        for (int i = 0; i < threadcnt; i++) delete[] data[i];
        MPI_Finalize();
    }
}

// 放入队列线程编号,并获取当前线程的数据
unsigned int* Graphmpi::get_edge_range() {
    int thread_num = omp_get_thread_num();
    idleq.push(thread_num); // 记录线程编号
    // 阻塞到获得锁,获得锁时该线程已经得到了数据
    for (;lock[thread_num].test_and_set(););
    // 若当前线程有数据则返回
    return ~data[thread_num][1] ? data[thread_num] : nullptr;
}

//! 更新结点匹配结果,并增加空闲线程数
void Graphmpi::report(long long local_ans) {
// omp 原子方式执行下一语句
#pragma omp atomic
    node_ans += local_ans;
#pragma omp atomic
    idlethreadcnt++;    // 空闲线程+1
    printf("node = %d, thread = %d, local_ans = %lld, time = %f\n",
           my_rank, omp_get_thread_num(), local_ans, get_wall_time() - starttime);
    fflush(stdout);
}

void Graphmpi::set_loop_flag() {
    loop_flag = true;
}

void Graphmpi::set_loop(int *data, int size) {
    int k = omp_get_thread_num();
    loop_data[k] = data;
    loop_size[k] = size;
}

void Graphmpi::get_loop(int *&data, int &size) {
    if (loop_flag) {
        int k = omp_get_thread_num();
        data = loop_data[k];
        size = loop_size[k];
    }
}

Bx2k256Queue::Bx2k256Queue() {
    memset(q, -1, sizeof(q));
    h = t = 0;
    lock.clear();
}

bool Bx2k256Queue::empty() {
    bool ans;
    for (;lock.test_and_set(););
    ans = (h == t);
    lock.clear();
    return ans;
}

void Bx2k256Queue::push(int k) {
    for (;lock.test_and_set(););
    q[t] = k;
    t++;
    lock.clear();
}

int Bx2k256Queue::front_and_pop() {
    int ans;
    for (;lock.test_and_set(););
    ans = q[h++];
    lock.clear();
    return ans;
}
