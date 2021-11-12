#pragma once
#include "graph.h"
#include <queue>
#include <atomic>

class Bx2k256Queue {
public:
    Bx2k256Queue();
    bool empty();
    void push(int);
    int front_and_pop();// can only be called by main thread

private:
    std::atomic_flag lock;  // 队列的自旋锁
    unsigned char h, t;     // 头尾指针
    int q[256];
};

class Graphmpi {
public:
    static Graphmpi& getinstance();
    // 初始化
    void init(int thread_count, Graph *graph, const Schedule& schedule); // get node range
    long long runmajor(); // mpi uses on major thread
    unsigned int* get_edge_range();
    void report(long long local_ans);
    void set_loop_flag();
    void set_loop(int*, int);
    void get_loop(int*&, int&);

private:
    static const int MAXTHREAD = 24;    // 最大线程数
    static const int MESSAGE_SIZE = 5;
    Graph* graph;   // 数据图
    int *loop_data[MAXTHREAD];
    int comm_sz;    // 通信器中的进程数
    int my_rank;    // 当前进程的序号(0号进程为主进程)
    int idlethreadcnt;  // 空闲线程数
    int threadcnt;  // 线程数
    int mpi_chunk_size; // MPI块大小
    int omp_chunk_size; // OpenMP块大小
    unsigned int loop_size[MAXTHREAD];
    unsigned int *data[MAXTHREAD];  // 每个线程的数据,data[t]表示线程t的图匹配数据
    long long node_ans;     // 进程(结点)一次任务的匹配结果
    double starttime;
    // 使用MPI的标识,即分布式时为true
    bool loop_flag = false; // loop_flag is set when using mpi;
    // 有(0,1)限制条件的标识
    bool skip_flag;     // skip_flag is set when there is a restriction on the first pattern edge
    // 进程的空闲线程队列
    Bx2k256Queue idleq;
    //std::queue<int> idleq;
    std::atomic_flag lock[MAXTHREAD];   // 每个线程的锁
    std::atomic_flag qlock;
    // 初始化标识
    bool initialized;
    Graphmpi();
    Graphmpi(const Graphmpi&&) = delete;
    Graphmpi(const Graphmpi&) = delete;
    Graphmpi& operator = (const Graphmpi&) = delete;
    ~Graphmpi();
};
