#include "../include/dataloader.h"
#include "../include/graph.h"
#include "../include/vertex_set.h"
#include "../include/common.h"
#include <cstdlib>
#include <cstring>
#include <algorithm>

// 加载图数据
bool DataLoader::load_data(Graph* &g, DataType type, const char* path, int oriented_type) {
    if(type == Patents || type == Orkut || type == complete8
    || type == LiveJournal || type == MiCo || type == CiteSeer
    || type == Wiki_Vote || type == Test) {
        return general_load_data(g, type, path, oriented_type);
    }

    if( type == Twitter) {
        return twitter_load_data(g, type, path, oriented_type);
    }
    printf("invalid DataType!\n");
    return false;
}

// 通用加载图数据函数
bool DataLoader::general_load_data(Graph *&g, DataType type, const char* path, int oriented_type) {
    // 读取数据集到标准输入
    if (freopen(path, "r", stdin) == NULL)
    {
        printf("File not found. %s\n", path);
        return false;
    }
    printf("Load begin in %s\n",path);
    g = new Graph();

    //load triangle counting information
    switch(type) {
        case DataType::Patents : {
            g->tri_cnt = Patents_tri_cnt;
            break;
        }
        case DataType::LiveJournal : {
            g->tri_cnt = LiveJournal_tri_cnt;
            break;
        }
        case DataType::MiCo : {
            g->tri_cnt = MiCo_tri_cnt;
            break;
        }
        case DataType::CiteSeer : {
            g->tri_cnt = CiteSeer_tri_cnt;
            break;
        }
        case DataType::Wiki_Vote : {
            g->tri_cnt = Wiki_Vote_tri_cnt;
            break;
        }
        case DataType::Orkut : {
            g->tri_cnt = Orkut_tri_cnt;
            break;
        }
        case DataType::Test: {
            g->tri_cnt = Test_tri_cnt;
            break;
        }
        default : {
            g->tri_cnt = -1;
            break;
        }
    }
    // 读取图的结点数和边数
    scanf("%d%u",&g->v_cnt,&g->e_cnt);
    // 记录结点的度
    int* degree = new int[g->v_cnt];
    memset(degree, 0, g->v_cnt * sizeof(int));
    g->e_cnt *= 2;  // QUE:边数*2的原因 ANS:CSR存图法存的相当于有向边,1条无向边会存两次,边数*2用于表示CSR数组的大小
    auto *e = new std::pair<int,int>[g->e_cnt];
    id.clear();
    int x,y;
    int tmp_v;
    unsigned int tmp_e;
    tmp_v = 0;
    tmp_e = 0;
    // 加载图数据
    while(scanf("%d%d",&x,&y)!=EOF) {
        // 读取时过滤自环边
        if(x == y) {
            printf("find self circle\n");
            g->e_cnt -=2;
            continue;
            //return false;
        }
        // 使用map和tmp_v来压缩结点集
        if(!id.count(x)) id[x] = tmp_v ++;
        if(!id.count(y)) id[y] = tmp_v ++;
        x = id[x];
        y = id[y];
        // 添加边
        e[tmp_e++] = std::make_pair(x,y);
        e[tmp_e++] = std::make_pair(y,x);
        ++degree[x];
        ++degree[y];
        //if(tmp_e % 1000000u == 0u) {
        //    printf("load %u edges\n",tmp_e);
        //    fflush(stdout);
        //}
    }

    // oriented_type == 0 do nothing
    //               == 1 high degree first
    //               == 2 low degree first
    if ( oriented_type != 0 ) {
        // 结点-结点度数 对
        std::pair<int,int> *rank = new std::pair<int,int>[g->v_cnt];
        int *new_id = new int[g->v_cnt];
        for(int i = 0; i < g->v_cnt; ++i) rank[i] = std::make_pair(i,degree[i]);
        // 按结点度数排序
        if( oriented_type == 1) std::sort(rank, rank + g->v_cnt, cmp_degree_gt);
        if( oriented_type == 2) std::sort(rank, rank + g->v_cnt, cmp_degree_lt);
        // newid数组用于将原id替换为排名
        for(int i = 0; i < g->v_cnt; ++i) new_id[rank[i].first] = i;
        // 替换所有边
        for(unsigned int i = 0; i < g->e_cnt; ++i) {
            e[i].first = new_id[e[i].first];
            e[i].second = new_id[e[i].second];
        }
        delete[] rank;
        delete[] new_id;
    }
    // 对结点度数组从小到大排序
    std::sort(degree, degree + g->v_cnt);

    // 两个结点邻域交集的大小,是度数第二大的结点的度数
    // The max size of intersections is the second largest degree.
    //TODO VertexSet::max_intersection_size has different value with
    // different dataset, but we use a static variable now.
    VertexSet::max_intersection_size = std::max(VertexSet::max_intersection_size,
                                                degree[g->v_cnt - 2]);
    g->max_degree = degree[g->v_cnt - 1];   // 最大结点度数
    delete[] degree;
    if(tmp_v != g->v_cnt) {
        printf("vertex number error!\n");
    }
    if(tmp_e != g->e_cnt) {
        printf("edge number error!\n");
    }
    if(tmp_v != g->v_cnt || tmp_e != g->e_cnt) {
        fclose(stdin);
        delete g;
        delete[] e;
        return false;
    }
    // 对边按结点升序排序去重
    std::sort(e,e+tmp_e,cmp_pair);
    g->e_cnt = std::unique(e,e+tmp_e) - e;  //  去重后边的实际数量
    for(unsigned int i = 0; i < g->e_cnt - 1; ++i) {
        if (e[i] == e[i + 1]) {
            printf("have same edge\n");
            fclose(stdin);
            delete g;
            delete[] e;
            return false;
        }
    }
    g->edge = new int[g->e_cnt];
    g->vertex = new unsigned int[g->v_cnt + 1];
    // CSR 存图法
    bool* have_edge = new bool[g->v_cnt];   // 记录结点是否右边
    int lst_v = -1;
    for(int i = 0; i < g->v_cnt; ++i) have_edge[i] = false;
    for(unsigned int i = 0; i < g->e_cnt; ++i) {
        if(e[i].first != lst_v) {
            have_edge[e[i].first] = true;
            g->vertex[e[i].first] = i;  // vertex数组记录结点起始边索引, 下标为边的起点
        }
        lst_v = e[i].first;
        g->edge[i] = e[i].second;   // edge数组只记录边的终点
    }
    delete[] e;
    printf("Success! There are %d nodes and %u edges.\n",g->v_cnt,g->e_cnt);
    fflush(stdout);
    g->vertex[g->v_cnt] = g->e_cnt;
    // 无边的情况
    for(int i = g->v_cnt - 1; i >= 0; --i)
        if(!have_edge[i]) {
            g->vertex[i] = g->vertex[i+1];
        }
    delete[] have_edge;

    return true;
}

bool DataLoader::twitter_load_data(Graph *&g, DataType type, const char* path, int oriented_type) {
    if (freopen(path, "r", stdin) == NULL)
    {
        printf("File not found. %s\n", path);
        return false;
    }
    printf("Load begin in %s\n",path);
    g = new Graph();
    g->tri_cnt = Twitter_tri_cnt;
    unsigned int* buffer = new unsigned int[41652230u + 2936729768u + 10u];
    FILE* file = fopen(path, "r");
    fread(buffer, sizeof(unsigned int), 41652230u + 2936729768u + 4, file);
    g->v_cnt = buffer[0];
    g->e_cnt = buffer[1];
    int mx_degree = buffer[2];
    VertexSet::max_intersection_size = std::max( VertexSet::max_intersection_size, mx_degree);
    g->max_degree = mx_degree;
    g->edge = new int [g->e_cnt];
    g->vertex = new unsigned int [g->v_cnt + 1];
    for(int i = 0; i < g->v_cnt + 1; ++i)
        g->vertex[i] = buffer[ 3 + i];
    for(unsigned int i = 0; i < g->e_cnt; ++i)
        g->edge[i] = buffer[4 + g->v_cnt + i];
    delete[] buffer;
    return true;
}

// 构造一个结点数为clique_size的完全图到g中
bool DataLoader::load_complete(Graph* &g, int clique_size) {
    g = new Graph();

    g->v_cnt = clique_size;
    g->e_cnt = clique_size * (clique_size - 1) / 2;

    int* degree = new int[g->v_cnt];
    memset(degree, 0, g->v_cnt * sizeof(int));
    g->e_cnt *= 2;
    std::pair<int,int> *e = new std::pair<int,int>[g->e_cnt];
    id.clear();
    int tmp_v;
    unsigned int tmp_e;
    tmp_v = 0;
    tmp_e = 0;
    for(int i = 0; i < clique_size; ++i)
        for(int j = 0; j < i; ++j) {
            int x = i, y = j;
            if(!id.count(x)) id[x] = tmp_v ++;
            if(!id.count(y)) id[y] = tmp_v ++;
            x = id[x];
            y = id[y];
            e[tmp_e++] = std::make_pair(x,y);
            e[tmp_e++] = std::make_pair(y,x);
            ++degree[x];
            ++degree[y];
        }

    std::sort(degree, degree + g->v_cnt);

    // The max size of intersections is the second largest degree.
    //TODO VertexSet::max_intersection_size has different value with different dataset, but we use a static variable now.
    VertexSet::max_intersection_size = std::max( VertexSet::max_intersection_size, degree[g->v_cnt - 2]);
    g->max_degree = degree[g->v_cnt - 1];
    delete[] degree;
    if(tmp_v != g->v_cnt) {
        printf("vertex number error!\n");
    }
    if(tmp_e != g->e_cnt) {
        printf("edge number error!\n");
    }
    if(tmp_v != g->v_cnt || tmp_e != g->e_cnt) {
        fclose(stdin);
        delete g;
        delete[] e;
        return false;
    }
    std::sort(e,e+tmp_e,cmp_pair);
    g->e_cnt = unique(e,e+tmp_e) - e;
    g->edge = new int[g->e_cnt];
    g->vertex = new unsigned int[g->v_cnt + 1];
    // CSR存图
    bool* have_edge = new bool[g->v_cnt];
    int lst_v = -1;
    for(int i = 0; i < g->v_cnt; ++i) have_edge[i] = false;
    for(unsigned int i = 0; i < g->e_cnt; ++i) {
        if(e[i].first != lst_v) {
            have_edge[e[i].first] = true;
            g->vertex[e[i].first] = i;
        }
        lst_v = e[i].first;
        g->edge[i] = e[i].second;
    }
    delete[] e;
    g->vertex[g->v_cnt] = g->e_cnt;
    for(int i = g->v_cnt - 1; i >= 0; --i)
        if(!have_edge[i]) {
            g->vertex[i] = g->vertex[i+1];
        }
    delete[] have_edge;
    return true;
}

bool DataLoader::cmp_pair(std::pair<int,int>a, std::pair<int,int>b) {
    return a.first < b.first || (a.first == b.first && a.second < b.second);
}

bool DataLoader::cmp_degree_gt(std::pair<int,int> a,std::pair<int,int> b) {
    return a.second > b.second;
}

bool DataLoader::cmp_degree_lt(std::pair<int,int> a,std::pair<int,int> b) {
    return a.second < b.second;
}

long long DataLoader::comb(int n, int k) {
    long long ans = 1;
    for(int i = n; i > n - k; --i)
        ans = ans * i;
    for(int i = 1; i <= k; ++i)
        ans = ans / k;
    return ans;
}
