#include <../include/graph.h>
#include <../include/dataloader.h>
#include "../include/pattern.h"
#include "../include/schedule.h"
#include "../include/common.h"
#include "../include/motif_generator.h"

#include <assert.h>
#include <iostream>
#include <string>
#include <algorithm>

void test_pattern(Graph* g, const Pattern &pattern, int performance_modeling_type,
                  int restricts_type, bool use_in_exclusion_optimize = false) {
    int thread_num = 24;
    double t1,t2;

    bool is_pattern_valid;
    Schedule schedule(pattern, is_pattern_valid, performance_modeling_type,
                      restricts_type, use_in_exclusion_optimize, g->v_cnt, g->e_cnt, g->tri_cnt);
    assert(is_pattern_valid);
    schedule.print_schedule();
    const auto& pairs = schedule.restrict_pair;
    printf("restricts: size=%llu ",pairs.size());
    for(auto& p : pairs)
        printf("(%d,%d)",p.first,p.second);
    puts("\n");

    double tt{};
    t1 = get_wall_time();
    g->stats.start(tt);
    long long ans = g->pattern_matching(schedule, thread_num);
    g->stats.end(TOTAL, tt);
    t2 = get_wall_time();
    printf("---------------------\noutput:\n");
    printf("ans %lld\n", ans);
    printf("time %.6lf\n", t2 - t1);
    std::cout << g->stats.ToString() << std::endl;
    fflush(stdout);

}

int main(int argc,char *argv[]) {

    if(argc < 5) {
        printf("Usage: %s dataset_name graph_file pattern_size pattern_adjacency_matrix\n", argv[0]);
        printf("Example(Triangle counting on dataset WikiVote) : \n");
        printf("%s Wiki-Vote ../../dataset/wiki-vote_input 3 011101110\n", argv[0]);
        return 0;
    }

    Graph *g;
    DataLoader D;

    const std::string type = argv[1];
    const std::string path = argv[2];

    int size = atoi(argv[3]);
    char* adj_mat = argv[4];

    // comments in include/schedule.h explain the meaning of these parameters.
    int test_type = 1; // performance_modeling_type = restricts_type = use_in_exclusion_optimize = 1

    DataType my_type;

    GetDataType(my_type, type);

    if(my_type == DataType::Invalid) {
        printf("Dataset not found!\n");
        return 0;
    }

    assert(D.load_data(g,my_type,path.c_str())==true);

    printf("Load data success!\n");
    fflush(stdout);

    Pattern p(size, adj_mat);   // 模式图
    test_pattern(g, p, test_type, test_type, false);

    delete g;
}
