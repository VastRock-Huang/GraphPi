#include "../include/graph.h"
#include "../include/graphmpi.h"
#include "../include/vertex_set.h"
#include "../include/common.h"
#include <cstdio>
#include <sys/time.h>
#include <unistd.h>
#include <cstdlib>
#include <omp.h>
#include <algorithm>
#include <cstring>
#include <mpi.h>
#include <atomic>
#include <queue>
#include <iostream>

int Graph::intersection_size(int v1, int v2) {
    unsigned int l1, r1;
    get_edge_index(v1, l1, r1);
    unsigned int l2, r2;
    get_edge_index(v2, l2, r2);
    int ans = 0;
    while (l1 < r1 && l2 < r2) {
        if (edge[l1] < edge[l2]) {
            ++l1;
        } else {
            if (edge[l2] < edge[l1]) {
                ++l2;
            } else {
                ++l1;
                ++l2;
                ++ans;
            }
        }
    }
    return ans;
}

/*int Graph::intersection_size_mpi(int v1, int v2) {
    Graphmpi &gm = Graphmpi::getinstance();
    int ans = 0;
    if (gm.include(v2))
        return intersection_size(v1, v2);
    unsigned int l1, r1;
    get_edge_index(v1, l1, r1);
    int *data = gm.getneighbor(v2);
    for (int l2 = 0; l1 < r1 && ~data[l2];) {
        if(edge[l1] < data[l2]) {
            ++l1;
        }
        else if(edge[l1] > data[l2]) {
            ++l2;
        }
        else {
            ++l1;
            ++l2;
            ++ans;
        }
    }
    return ans;
}
*/

int Graph::intersection_size_clique(int v1, int v2) {
    unsigned int l1, r1;
    get_edge_index(v1, l1, r1);
    unsigned int l2, r2;
    get_edge_index(v2, l2, r2);
    int min_vertex = v2;
    int ans = 0;
    if (edge[l1] >= min_vertex || edge[l2] >= min_vertex)
        return 0;
    while (l1 < r1 && l2 < r2) {
        if (edge[l1] < edge[l2]) {
            if (edge[++l1] >= min_vertex)
                break;
        } else {
            if (edge[l2] < edge[l1]) {
                if (edge[++l2] >= min_vertex)
                    break;
            } else {
                ++ans;
                if (edge[++l1] >= min_vertex)
                    break;
                if (edge[++l2] >= min_vertex)
                    break;
            }
        }
    }
    return ans;
}

long long Graph::triangle_counting() {
    long long ans = 0;
    for (int v = 0; v < v_cnt; ++v) {
        // for v in G
        unsigned int l, r;
        get_edge_index(v, l, r);
        for (unsigned int v1 = l; v1 < r; ++v1) {
            //for v1 in N(v)
            ans += intersection_size(v, edge[v1]);
        }
    }
    ans /= 6;
    return ans;
}

long long Graph::triangle_counting_mt(int thread_count) {
    long long ans = 0;
#pragma omp parallel num_threads(thread_count)
    {
        tc_mt(&ans);
    }
    return ans;
}

void Graph::tc_mt(long long *global_ans) {
    long long my_ans = 0;
//    #pragma omp for schedule(dynamic)
    for (int v = 0; v < v_cnt; ++v) {
        // for v in G
        unsigned int l, r;
        get_edge_index(v, l, r);
        for (unsigned int v1 = l; v1 < r; ++v1) {
            if (v <= edge[v1])
                break;
            //for v1 in N(v)
            my_ans += intersection_size_clique(v, edge[v1]);
        }
    }
#pragma omp critical
    {
        *global_ans += my_ans;
    }
}

//! 获取v结点的边在edge数组中的范围,范围为[l,r-1]
void Graph::get_edge_index(int v, unsigned int &l, unsigned int &r) const {
    l = vertex[v];
    r = vertex[v + 1];
}

void Graph::pattern_matching_func(const Schedule &schedule, VertexSet *vertex_set, VertexSet &subtraction_set,
                                  long long &local_ans, int depth, bool clique) {
    int loop_set_prefix_id = schedule.get_loop_set_prefix_id(depth);
    int loop_size = vertex_set[loop_set_prefix_id].get_size();
    if (loop_size <= 0)
        return;
    int *loop_data_ptr = vertex_set[loop_set_prefix_id].get_data_ptr();
    /*if (clique == true)
      {
      int last_vertex = subtraction_set.get_last();
    // The number of this vertex must be greater than the number of last vertex.
    loop_start = std::upper_bound(loop_data_ptr, loop_data_ptr + loop_size, last_vertex) - loop_data_ptr;
    }*/
    if (depth == schedule.get_size() - 1) {
        // TODO : try more kinds of calculation.
        // For example, we can maintain an ordered set, but it will cost more to maintain itself when entering or exiting recursion.
        if (clique == true)
            local_ans += loop_size;
        else if (loop_size > 0)
            local_ans += VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set);
        return;
    }

    int last_vertex = subtraction_set.get_last();
    for (int i = 0; i < loop_size; ++i) {
        if (last_vertex <= loop_data_ptr[i] && clique == true)
            break;
        int vertex = loop_data_ptr[i];
        if (!clique)
            if (subtraction_set.has_data(vertex))
                continue;
        unsigned int l, r;
        get_edge_index(vertex, l, r);
        bool is_zero = false;
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
            vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, &edge[l], (int) r - l, prefix_id, vertex,
                                                   clique);
            if (vertex_set[prefix_id].get_size() == 0) {
                is_zero = true;
                break;
            }
        }
        if (is_zero) continue;
        //subtraction_set.insert_ans_sort(vertex);
        subtraction_set.push_back(vertex);
        pattern_matching_func(schedule, vertex_set, subtraction_set, local_ans, depth + 1, clique);
        subtraction_set.pop_back();
    }
}

long long Graph::pattern_matching(const Schedule &schedule, int thread_count, bool clique) {
    long long global_ans = 0;
// omp parallel: OpenMP定义并行区域
// num_threads(): 设置线程数
// reduction(+: <var>): 表示该变量var最后多个线程求和(+)
//#pragma omp parallel num_threads(thread_count) reduction(+: global_ans)
    {
        double start_time = get_wall_time();
        double current_time;
        VertexSet *vertex_set = new VertexSet[schedule.get_total_prefix_num()];
        VertexSet subtraction_set;
        VertexSet tmp_set;
        subtraction_set.init();
        // 本线程满足匹配的数量
        long long local_ans = 0;
        // TODO : try different chunksize
// omp for: 分配循环到多线程
// schedule(dynamic): 动态分配到多线程
// nowait: 取消for循环后的隐含屏障
//#pragma omp for schedule(dynamic) nowait
        // 遍历输入数据图的每个结点
        for (int vertex = 0; vertex < v_cnt; ++vertex) {
            unsigned int l, r;  // 结点vertex在edge数组中的范围[l,r-1]
            get_edge_index(vertex, l, r);
            // 遍历以模式图中0结点结尾的前缀prefix
            // 经过该循环后,所有以0结点为结尾的前缀的结点集都被设置为了当前在数据图在匹配的结点vertex的邻域
            //理论上以0结点为结尾的前缀只有一个,但该前缀对应的结点不止一个,但至少1结点满足
            for (int prefix_id = schedule.get_last(0); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
                // 求以0结点结尾前缀的父前缀的结点集vertex_set[father_id]与结点vertex的邻域(input_data)的交集
                // 理论上讲结点0是第一个结点,以0为结尾的前缀应该没有父结点,
                //因此vertex_set[prefix_id]=input_data

                vertex_set[prefix_id].build_vertex_set(stats, schedule, vertex_set,
                                                       &edge[l], (int) r - l, prefix_id);
                /*
                if(vertex_set[prefix_id].get_size()==r-l){
                    for (int j = 0; j < vertex_set[prefix_id].get_size(); ++j) {
                        if(vertex_set[prefix_id].get_data(j)!=edge[l+j]){
                            printf("false\n");
                            break;
                        }
                    }
                }else{
                    printf("false\n");
                }*/
            }
            //subtraction_set.insert_ans_sort(vertex);
            // 存入当前数据图结点
            subtraction_set.push_back(vertex);
            //if (schedule.get_total_restrict_num() > 0 && clique == false)
            pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1);
            /*
            if(true)
                pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1);
            else
                pattern_matching_func(schedule, vertex_set, subtraction_set, local_ans, 1, clique); */
            subtraction_set.pop_back();
            /*
            if( (vertex & (-vertex)) == (1<<15) ) {
                current_time = get_wall_time();
                if( current_time - start_time > max_running_time) {
                    printf("TIMEOUT!\n");
                    fflush(stdout);
                    assert(0);
                }
            }*/
        }
/*
        for(int i=0; i<schedule.get_total_prefix_num(); ++i) {
            printf("vertex_set %d: ",i);
            for(int j=0;j<vertex_set[i].get_size();++j){
                printf("%d ", vertex_set[i].get_data(j));
            }
            putchar('\n');
        }*/
        delete[] vertex_set;
        // TODO : Computing multiplicity for a pattern
        global_ans += local_ans;    // 将线程局部变量的计算结果加到总结过中
//        printf("time + 1\n");

    }
    // 除以使用容斥定理引入的冗余倍数,得到正确结果
    //若未使用容斥定理in_exclusion_optimize_redundancy为1
    return global_ans / schedule.get_in_exclusion_optimize_redundancy();
}

//! 深度优先回溯图匹配
//! \param[in] depth 匹配深度,也可以视为第depth个结点,因为结点序号就是匹配的顺序, from 1 to size-1(size为模式图结点数)
//! \param[in] subtraction_set 当前已匹配的结点集,subtraction_set[i]对应着模式图中的第i个结点
//! \param[in] vertex_set 每个前缀的候选结点集,vertex_set[prefix_id]即以prefix_id号前缀为前缀的数据图中的所有结点
void Graph::pattern_matching_aggressive_func(const Schedule &schedule, VertexSet *vertex_set,
                                             VertexSet &subtraction_set, VertexSet &tmp_set,
                                             long long &local_ans,
                                             int depth) // 3 same # or @ in comment are useful in code generation ###
{
    // depth结点的前缀索引
    //第1次调用时获取的为第1个结点的前缀索引,理论上该前缀为[0],即只有第0个结点
    int loop_set_prefix_id = schedule.get_loop_set_prefix_id(depth);// @@@
    // depth结点的前缀的结点集的大小
    //
    int loop_size = vertex_set[loop_set_prefix_id].get_size();

    if (loop_size <= 0)
        return;
    // depth结点的前缀的结点集
    //depth结点的前缀一定以x(0<=x<depth)结点结尾,因此在上层递归调用中已经求出
    int *loop_data_ptr = vertex_set[loop_set_prefix_id].get_data_ptr();
/* @@@ 
    //Case: in_exclusion_optimize_num = 2
    if (depth == schedule.get_size() - 2 && schedule.get_in_exclusion_optimize_num() == 2) { 
        int loop_set_prefix_id_nxt = schedule.get_loop_set_prefix_id( depth + 1);
        int loop_size_nxt = vertex_set[loop_set_prefix_id_nxt].get_size();
        int size1 = VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set);
        int size2 = VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id_nxt], subtraction_set);
        VertexSet tmp_set;
        tmp_set.init();
        tmp_set.intersection(vertex_set[loop_set_prefix_id], vertex_set[loop_set_prefix_id_nxt]);
        int size3 = VertexSet::unorderd_subtraction_size(tmp_set, subtraction_set);
        local_ans += 1ll * size1 * size2 - size3;
        return;
    }
*/
/*
    //Case: in_exclusion_optimize_num = 3
    if( depth == schedule.get_size() - 3 && schedule.get_in_exclusion_optimize_num() == 3) { 
        int in_exclusion_optimize_num = 3;
        int loop_set_prefix_ids[ in_exclusion_optimize_num];
        for(int i = 0; i < in_exclusion_optimize_num; ++i)
            loop_set_prefix_ids[i] = schedule.get_loop_set_prefix_id( depth + i );
        
        int loop_sizes[ in_exclusion_optimize_num ];
        for(int i = 0; i < in_exclusion_optimize_num; ++i)
            loop_sizes[i] = VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_ids[i]], subtraction_set);
        
        local_ans += 1ll * loop_sizes[0] * loop_sizes[1] * loop_sizes[2];

        for(int i = 1; i < 3; ++i) 
            for(int j = 0; j < i; ++j){
                VertexSet tmp_set;
                tmp_set.init();
                tmp_set.intersection(vertex_set[loop_set_prefix_ids[i]], vertex_set[loop_set_prefix_ids[j]]);
                int tmp_size = VertexSet::unorderd_subtraction_size(tmp_set, subtraction_set);
                int size2;
                for(int k = 0; k < 3; ++k)
                    if( i != k && j != k) size2 = loop_sizes[k];
                local_ans -= 1ll * tmp_size * size2;
            }
        VertexSet tmp1;
        tmp1.init();
        tmp1.intersection(vertex_set[loop_set_prefix_ids[0]], vertex_set[loop_set_prefix_ids[1]]);
        VertexSet tmp2;
        tmp2.init();
        tmp2.intersection(vertex_set[loop_set_prefix_ids[2]], tmp1);
        local_ans += 1ll * 2 * VertexSet::unorderd_subtraction_size(tmp2, subtraction_set);
        return;
    }
*/
    //Case: in_exclusion_optimize_num > 1
    // 若匹配的深度达到了倒数第k个
    //若未使用容斥定理,k=0,不会递归到该深度,则以下代码不会生效
    if (depth == schedule.get_size() - schedule.get_in_exclusion_optimize_num()) {
        // 容斥定理的k值
        int in_exclusion_optimize_num = schedule.get_in_exclusion_optimize_num();// @@@
        // 倒数k个结点的前缀索引的数组
        int loop_set_prefix_ids[in_exclusion_optimize_num];
        loop_set_prefix_ids[0] = loop_set_prefix_id;
        for (int i = 1; i < in_exclusion_optimize_num; ++i)
            loop_set_prefix_ids[i] = schedule.get_loop_set_prefix_id(depth + i);
        // 遍历每种容斥定理的情况
        for (int optimize_rank = 0; optimize_rank < schedule.in_exclusion_optimize_group.size(); ++optimize_rank) {
            const std::vector<std::vector<int> > &cur_graph = schedule.in_exclusion_optimize_group[optimize_rank];
            long long val = schedule.in_exclusion_optimize_val[optimize_rank];
            // 遍历一种容斥定理分组情况的每个分组,并求交集
            for (int cur_graph_rank = 0; cur_graph_rank < cur_graph.size(); ++cur_graph_rank) {
                //                VertexSet tmp_set;

                //if size == 1 , we will not call intersection(...)
                //so we will not allocate memory for data
                //otherwise, we need to copy the data to do intersection(...)
                if (cur_graph[cur_graph_rank].size() == 1) {
                    int id = loop_set_prefix_ids[cur_graph[cur_graph_rank][0]];
                    val = val * VertexSet::unorderd_subtraction_size(vertex_set[id], subtraction_set);
                } else {
                    // 求多个结点集的交集
                    int id0 = loop_set_prefix_ids[cur_graph[cur_graph_rank][0]];
                    int id1 = loop_set_prefix_ids[cur_graph[cur_graph_rank][1]];
                    tmp_set.init(this->max_degree);
                    tmp_set.intersection(vertex_set[id0], vertex_set[id1]);

                    for (int i = 2; i < cur_graph[cur_graph_rank].size(); ++i) {
                        int id = loop_set_prefix_ids[cur_graph[cur_graph_rank][i]];
                        tmp_set.intersection_with(vertex_set[id]);
                    }
                    val = val * VertexSet::unorderd_subtraction_size(tmp_set, subtraction_set);
                }
                if (val == 0) break;

            }
            local_ans += val;
        }
        return;// @@@

    }
    //Case: in_exclusion_optimize_num <= 1
    // 匹配深度到达最后一个结点
    if (depth == schedule.get_size() - 1) {
        // TODO : try more kinds of calculation. @@@
        // For example, we can maintain an ordered set, but it will cost more to maintain itself when entering or exiting recursion.
        if (schedule.get_total_restrict_num() > 0)  // 限制条件数量>0
        {
            // 以depth结点为终结点,(在数据图中)序号最小的起始结点
            int min_vertex = v_cnt;
            // 遍历以depth为终结点的限制条件
            for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i)) {
                if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
                    min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));
            }
            // depth结点的前缀的结点集
            const VertexSet &vset = vertex_set[loop_set_prefix_id];
            // lower_bound():得到首个不小于某数的迭代器
            // 小于min_vertex结点的数量
            int size_after_restrict =
                    std::lower_bound(vset.get_data_ptr(), vset.get_data_ptr() + vset.get_size(), min_vertex)
                    - vset.get_data_ptr();
            if (size_after_restrict > 0) {
                // 将候选结点的前size_after_restrict结点排除已匹配的结点,剩余的阶段即为满足的结点
                //其个数加到本线程的满足匹配的个数中
                local_ans += VertexSet::unorderd_subtraction_size(vset, subtraction_set, size_after_restrict);
//                printf("embeding:");
//                for (auto i = 0; i < subtraction_set.get_size(); ++i) {
//                    printf("%d, ", subtraction_set.get_data(i));
//                }
//                printf("\n");
            }
        } else {  // 若无限制条件
            // 则候选结点集减去已匹配的结点,剩余结点的数量就是满足模式图的数量
            local_ans += VertexSet::unorderd_subtraction_size(vertex_set[loop_set_prefix_id], subtraction_set);
        }

//        printf("vertex_set:\n");
//        for(int i=0;i<schedule.get_total_prefix_num(); ++i){
//            printf("No.%d: ", i);
//            for(int j=0;j<vertex_set[i].get_size();++j){
//                printf("%d ", vertex_set[i].get_data(j));
//            }
//            putchar('\n');
//        }

        return;// @@@
    }

    // TODO : min_vertex is also a loop invariant @@@
    // 以depth结点为终结点,(在数据图中)序号最小的起始结点
    int min_vertex = v_cnt;
    // 遍历以depth为终结点的限制条件找到序号最小的结点
    //对于同终结点的限制条件(a,c),(b,c),在模式图中:vertex_a<vertex_c,vertex_b<vertex_c
    //则在数据图中则满足:vertex_a>vertex_c,vertex_b>vertex_c => min(vertex_a,vertex_b)>vertex_c
    //所有下面的循环中找的下一个结点必须小于min_vertex,否则跳出循环
    for (int i = schedule.get_restrict_last(depth); i != -1; i = schedule.get_restrict_next(i))
        if (min_vertex > subtraction_set.get_data(schedule.get_restrict_index(i)))
            min_vertex = subtraction_set.get_data(schedule.get_restrict_index(i));
    if (depth == 1) Graphmpi::getinstance().get_loop(loop_data_ptr, loop_size);
    int ii = 0;
    // 遍历depth结点的前缀的结点集
    //由于求交集intersection()算法的特性,以序号从小到大进行遍历
    for (int &i = ii; i < loop_size; ++i) {
        // 若min_vertex不大于depth结点的前缀的结点集中的第i个结点,则退出
        //对于限制条件,数据图结点满足的大小关系和模式图中结点满足的大小关系正好相反
        //如限制条件(a,b),在模式图中是vertex_a<vertex_b,而在匹配的数据图中是vertex_a>vertex_b
        //由于loop_data_ptr数组是有序的,该结点不满足则后面的都不满足因此跳过
        if (min_vertex <= loop_data_ptr[i])
            break;
        // depth结点的前缀的结点集中的第i个结点
        int vertex = loop_data_ptr[i];
        // 若该结点在已搜索的结点中则跳过
        if (subtraction_set.has_data(vertex))
            continue;
        unsigned int l, r;
        // 获取结点的边数据
        get_edge_index(vertex, l, r);
        bool is_zero = false;
        // 遍历以depth结点为结尾的前缀
        for (int prefix_id = schedule.get_last(depth); prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
            // 求以depth结点结尾前缀的父前缀的结点集vertex_set[father_id]与结点vertex的领域(input_data)的交集
            //以depth结点结尾的前缀的父前缀实际上就是上一个匹配结点结尾的前缀,该函数实际上就是求交集操作,
            //求出了以当前已匹配结点为前缀的所有候选结点,因为相同的前缀的候选结点的范围是相同的
            vertex_set[prefix_id].build_vertex_set(stats, schedule, vertex_set, &edge[l],
                                                   (int) r - l, prefix_id);
            // 若交集为空则退出循环
            if (vertex_set[prefix_id].get_size() == 0) {
                is_zero = true;
                break;
            }
        }
        // 有以depth结点为结尾的前缀, 但交集为零
        //则下一深度的候选结点不存在,则跳过本次循环判断另一个本深度的候选结点
        if (is_zero) continue;
        //subtraction_set.insert_ans_sort(vertex);
        // 置入本次匹配的结点
        subtraction_set.push_back(vertex);
        // 递归匹配下一深度的结点
        pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set,
                                         tmp_set, local_ans, depth + 1);// @@@
        subtraction_set.pop_back(); // @@@
    }
    //if (depth == 1 && ii < loop_size) Graphmpi::getinstance().set_cur(subtraction_set.get_data(0));// @@@
}

// ###
long long Graph::pattern_matching_mpi(const Schedule &schedule, int thread_count, bool clique) {
    Graphmpi &gm = Graphmpi::getinstance();     // 获取Graphmpi对象
    long long global_ans = 0;
// 并行
#pragma omp parallel num_threads(thread_count)
    {
// 仅主线程执行
#pragma omp master
        {
            gm.init(thread_count, this, schedule);
        }
// 同步屏障
#pragma omp barrier //mynodel have to be calculated before running other threads
#pragma omp master
        {
            global_ans = gm.runmajor();
        }

        // 非主线程运行
        if (omp_get_thread_num()) {
            VertexSet *vertex_set = new VertexSet[schedule.get_total_prefix_num()];
            long long local_ans = 0;    // 线程局部匹配结果
            VertexSet subtraction_set;  // 已匹配结点集
            VertexSet tmp_set;
            subtraction_set.init();
            int last = -1;
            // 设置使用MPI分布式
            gm.set_loop_flag();
            auto match_edge = [&](int vertex, int *data, int size) {
                if (vertex != last) {
                    if (~last) subtraction_set.pop_back();
                    unsigned int l, r;
                    get_edge_index(vertex, l, r);
                    // 设置0结点为结尾的前缀的结点集为当前匹配结点vertex的邻域
                    for (int prefix_id = schedule.get_last(0);
                         prefix_id != -1; prefix_id = schedule.get_next(prefix_id)) {
                        vertex_set[prefix_id].build_vertex_set(schedule, vertex_set, edge + l, r - l, prefix_id);
                    }
                    subtraction_set.push_back(vertex);
                    last = vertex;
                }
                gm.set_loop(data, size);
                // 图匹配
                pattern_matching_aggressive_func(schedule, vertex_set, subtraction_set, tmp_set, local_ans, 1);
            };
            for (unsigned int *data; (data = gm.get_edge_range());) {
                // 进行图匹配
                match_edge(data[1], edge + data[2], data[3] - data[2]);
                /*for (int i = 1; i <= data[4]; i++) {
                    int l, r;
                    get_edge_index(data[1] + i, l, r);
                    match_edge(data[1] + i, edge + l, r - l);
                }*/
            }
            delete[] vertex_set;
            gm.report(local_ans);   // 更新结果
        }
    }
    return global_ans;
}
