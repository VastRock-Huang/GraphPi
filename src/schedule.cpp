#include "../include/schedule.h"
#include "../include/graph.h"
#include "../include/dataloader.h"
#include <cstdio>
#include <cstring>
#include <assert.h>
#include <algorithm>
#include <iostream>

Schedule::Schedule(const Pattern& pattern, bool &is_pattern_valid,
                   int performance_modeling_type, int restricts_type,
                   bool use_in_exclusion_optimize ,int v_cnt,
                   unsigned int e_cnt, long long tri_cnt)
{
    if( performance_modeling_type != 0 && tri_cnt == -1) {
        printf("Fatal: Can not use performance modeling if not "
               "have triangle number of this dataset.\n");
        fflush(stdout);
        assert(0);
    }

    is_pattern_valid = true;
    size = pattern.get_size();  // 结点数
    adj_mat = new int[size * size];

    // 拷贝邻接矩阵
    // not use performance_modeling, simply copy the adj_mat from pattern
    memcpy(adj_mat, pattern.get_adj_mat_ptr(),
           size * size * sizeof(int));

    std::vector< std::pair<int,int> > best_pairs;
    best_pairs.clear();
    //Initialize adj_mat
    //If we use performance_modeling, we may change the order of vertex,
    //the best order produced by performance_modeling(...) is saved in
    // best_order[]
    //Finally, we use best_order[] to relocate adj_mat
    if(performance_modeling_type != 0) {    // use performance_modeling
        unsigned int pow = 1;
        for (int i = 2; i <= size; ++i) pow *= i;
        // schedules的候选排列
        std::vector< std::vector<int> > candidate_permutations;
        candidate_permutations.clear();

        bool use[size];
        for (int i = 0; i < size; ++i) use[i] = false;
        std::vector<int> tmp_vec;
        // 获取schedules的n!可能的全部排列
        get_full_permutation(candidate_permutations, use,
                             tmp_vec, 0);
        assert(candidate_permutations.size() == pow);
        // 按Phase1规则移除无效排列
        remove_invalid_permutation(candidate_permutations);
        // use GraphPi's modeling
        if(performance_modeling_type == 1) {
            //reduce candidates
            // 该模式Pattern在Phase2中定义的k值,即所有排列中最大的k值
            int max_val = 0;
            for(const auto &vec : candidate_permutations) {
                max_val = std::max(max_val, get_vec_optimize_num(vec));
            }
            std::vector< std::vector<int> > tmp;
            tmp.clear();
            // 仅保留满足模式k值的排列
            for(const auto &vec : candidate_permutations)
                if(get_vec_optimize_num(vec) == max_val) {
                    tmp.push_back(vec);
                }
            candidate_permutations = tmp;
        }

        // 最佳schedule
        int *best_order = new int[size];
        double min_val;
        bool have_best = false;

        // 遍历每个候选schedule
        for(const auto &vec : candidate_permutations) {
            int rank[size];     // 记录排列的序号
            for(int i = 0; i < size; ++i) rank[vec[i]] = i;

            // Pattern新的邻接矩阵:按照当前排列的序号重新构造
            int* cur_adj_mat;
            cur_adj_mat = new int[size*size];
            for(int i = 0; i < size; ++i)
                for(int j = 0; j < size; ++j)
                    cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

            // 限制条件集合
            std::vector< std::vector< std::pair<int,int> > > restricts_vector;
            restricts_vector.clear();

            // use GraphPi's restricts
            if(restricts_type == 1) {
                restricts_generate(cur_adj_mat, restricts_vector);
            }
            else {  // not use restricts or use GraphZero's restricts
                Schedule schedule(cur_adj_mat, size);

                std::vector< std::pair<int,int> > pairs;
                schedule.GraphZero_aggressive_optimize(pairs);

                restricts_vector.clear();
                restricts_vector.push_back(pairs);
            }
            // 若无限制条件
            if( restricts_vector.size() == 0) {
                std::vector< std::pair<int,int> > Empty;    //空限制条件集合
                Empty.clear();

                // 性能开销参数,代表当前configuration的性能开销,越小越好
                double val;
                // use GraphPi's modeling
                if(performance_modeling_type == 1) {
                    val = our_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt, tri_cnt);
                }
                else {
                    // use GraphZero's modeling
                    if(performance_modeling_type == 2) {
                        val = GraphZero_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt);
                    }
                    else {  // use naive modeling
                        val = Naive_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt);
                    }
                }
                // 记录当前空限制条件的schedule为最佳
                if(have_best == false || val < min_val) {
                    have_best = true;
                    min_val = val;
                    for(int i = 0; i < size; ++i) best_order[i] = vec[i];
                    best_pairs = Empty;
                }
            }
//            std::cout <<"cur_adj:[";
//            for(unsigned i = 0; i < size; ++i) {
//                for(unsigned j =0; j< size; ++j) {
//                    std::cout << cur_adj_mat[INDEX(i,j,size)];
//                }
//                std::cout <<" ";
//            }
//            std::cout <<"]\n";
//
//            std::cout <<"permut:[";
//            for(auto x: vec) {
//                std::cout << x << " ";
//            }
//            std::cout << "]\n";
//
//            std::cout << "restricts_vec:[\n";
//            for (const auto &rs: restricts_vector) {
//                std::cout << " [";
//                for (const auto &p: rs) {
//                    std::cout << "(" << p.first << ", " << p.second << ") ";
//                }
//                std::cout << "]\n";
//            }
//            std::cout << "]\n";

            // 遍历每组限制条件
            for(const auto& pairs : restricts_vector) {
                double val;
                if(performance_modeling_type == 1) {
                    val = our_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt, tri_cnt);
                }
                else {
                    if(performance_modeling_type == 2) {
                        val = GraphZero_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt);
                    }
                    else {
                        val = Naive_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt);
                    }
                }

                if(have_best == false || val < min_val) {
                    have_best = true;
                    min_val = val;
                    for(int i = 0; i < size; ++i) best_order[i] = vec[i];
                    best_pairs = pairs;
                }
            }

        }

        // 根据最佳configuration的顺序重构模式图
        int rank[size];
        for(int i = 0; i < size; ++i) rank[best_order[i]] = i;

        const int* pattern_adj_mat = pattern.get_adj_mat_ptr();
        for(int i = 0; i < size; ++i)
            for(int j = 0; j < size; ++j)
                adj_mat[INDEX(rank[i], rank[j], size)] = pattern_adj_mat[INDEX(i, j, size)];
        delete[] best_order;
    }
    else {  // not use performance modeling
        std::vector< int > I;
        I.clear();
        for(int i = 0; i < size; ++i) I.push_back(i);

        std::vector< std::vector< std::pair<int,int> > > restricts_vector;
        restricts_vector.clear();

        // use restricts
        if(restricts_type != 0) {
            // use GraphPi's restricts
            if(restricts_type == 1) {
                restricts_generate(adj_mat, restricts_vector);
            }
            else {  // use GraphZero's restricts
                std::vector< std::pair<int,int> > pairs;
                GraphZero_aggressive_optimize(pairs);

                restricts_vector.clear();
                restricts_vector.push_back(pairs);
            }
        }

        bool have_best = false;
        double min_val;

        for(const auto& pairs : restricts_vector) {
            double val;
            if(restricts_type == 1) {
                val = our_estimate_schedule_restrict(I, pairs, v_cnt, e_cnt, tri_cnt);
            }
            else {
                val = GraphZero_estimate_schedule_restrict(I, pairs, v_cnt, e_cnt);
            }
            if(have_best == false || val < min_val) {
                have_best = true;
                min_val = val;
                best_pairs = pairs;
            }
        }

    }

    // 是否使用容斥规则
    if(use_in_exclusion_optimize) {
        std::vector<int> I;
        I.clear();
        for(int i = 0; i < size; ++i) I.push_back(i);
        in_exclusion_optimize_num = get_vec_optimize_num(I);
        if( in_exclusion_optimize_num <= 1) {
            printf("Can not use in_exclusion_optimize with this schedule\n");
            in_exclusion_optimize_num = 0;
        }
        else {
            printf("use in_exclusion_optimize with size %d\n", in_exclusion_optimize_num);
            init_in_exclusion_optimize();
        }
    }
    else {
            in_exclusion_optimize_num = 0;
    }

    // 第i个循环涉及的候选结点至多为前i-1个结点的邻结点的交集,是结点是全连通图,
    //第1个结点邻结点有(size-1)个,第2个结点的临界点(不包括1结点)有(size-2)个,...
    //全部则为 size*(size-1)/2 个
    // The I-th loop consists of at most the intersection of i-1 VertexSet.
    // So the max number of prefix = 0 + 1 + ... + size-1 = size * (size-1) / 2
    int max_prefix_num = size * (size - 1) / 2;
    father_prefix_id = new int[max_prefix_num];
    last = new int[size];
    next = new int[max_prefix_num];
    loop_set_prefix_id = new int[size];
    prefix = new Prefix[max_prefix_num];
    restrict_last = new int[size];
    restrict_next = new int[max_prefix_num];
    restrict_index = new int[max_prefix_num];
    memset(father_prefix_id, -1, max_prefix_num * sizeof(int));
    memset(last, -1, size * sizeof(int));
    memset(next, -1, max_prefix_num * sizeof(int));
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));

    total_prefix_num = 0;
    total_restrict_num = 0;

    // Phase1规则检查
    // The I-th vertex must connect with at least one vertex from 0 to i-1.
    for (int i = 1; i < size; ++i)
    {
        bool valid = false;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)])
            {
                valid = true;
                break;
            }
        if (valid == false)
        {
            //Invalid Schedule
            is_pattern_valid = false;
            return;
        }
    }

    build_loop_invariant();
    if( restricts_type != 0) add_restrict(best_pairs);

    set_in_exclusion_optimize_redundancy();

    printf("prefix:\n");
    for(int i=0;i<total_prefix_num;++i){
        printf("\tNo.%d: ", i);
        for(int j=0;j<prefix[i].get_size();++j){
            printf("%d ",prefix[i].get_data(j));
        }
        putchar('\n');
    }

    printf("restrict list:\n");
    for(int i=0;i<size;++i){
        printf("\tvertex%d: ",i);
        for(int id=restrict_last[i];id!=-1;id=restrict_next[id]){
            printf("%d ", restrict_index[id]);
        }
        putchar('\n');
    }

    printf("in_exclusion_optimize_group:\n");
    for(int i=0;i<in_exclusion_optimize_group.size();++i){
        printf("\t%d:\n",i);
        for(int j=0;j<in_exclusion_optimize_group[i].size();++j){
            printf("\t %d: ", j);
            for(int k=0;k<in_exclusion_optimize_group[i][j].size();++k){
                printf("%d ", in_exclusion_optimize_group[i][j][k]);
            }
            putchar('\n');
        }
    }
    printf("in_exclusion_optimize_val: ");
    for(int i=0;i<in_exclusion_optimize_val.size();++i){
        printf("%d ",in_exclusion_optimize_val[i]);
    }
    putchar('\n');

    printf("last:");
    for(int i=0;i<size;++i){
        printf("%d ", last[i]);
    }
    putchar('\n');
    printf("next:");
    for(int i=0;i<total_prefix_num;++i){
        printf("%d ",next[i]);
    }
    putchar('\n');
}

Schedule::Schedule(const int* _adj_mat, int _size)
{
    size = _size;
    adj_mat = new int[size * size];

    memcpy(adj_mat, _adj_mat, size * size * sizeof(int));

    // QUE:为什么这么计算max_prefix_num,是么时候取到最大值? ANS:根据前缀的生成规则,
    //第i个结点可能与前i-1个结点连接,因此其前缀为(0,1,...,i-1),由于前缀的父前缀也会同时生成,
    //因此共有i个前缀:(0),(0,1),(0,1,2)...(0,1,...,i-1),对于size个结点,最多的前缀数
    //即为0+1+2+...size-1 = size * (size-1) / 2,但应该有重复不会达到该值
    // The I-th loop consists of at most the intersection of i-1 VertexSet.
    // So the max number of prefix = 0 + 1 + ... + size-1 = size * (size-1) / 2
    int max_prefix_num = size * (size - 1) / 2;
    father_prefix_id = new int[max_prefix_num];
    last = new int[size];
    next = new int[max_prefix_num];
    loop_set_prefix_id = new int[size];
    prefix = new Prefix[max_prefix_num];
    restrict_last = new int[size];
    restrict_next = new int[max_prefix_num];
    restrict_index = new int[max_prefix_num];
    memset(father_prefix_id, -1, max_prefix_num * sizeof(int));
    memset(last, -1, size * sizeof(int));
    memset(next, -1, max_prefix_num * sizeof(int));
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));

    total_prefix_num = 0;
    total_restrict_num = 0;
    in_exclusion_optimize_num = 0;

    // Phase1规则检查
    // The I-th vertex must connect with at least one vertex from 0 to i-1.
    for (int i = 1; i < size; ++i)
    {
        bool valid = false;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)])
            {
                valid = true;
                break;
            }
        if (valid == false)
        {
            printf("invalid schedule!\n");
            assert(0);
        }
    }

    build_loop_invariant();

    set_in_exclusion_optimize_redundancy();
}

Schedule::~Schedule()
{
    delete[] adj_mat;
    delete[] father_prefix_id;
    delete[] last;
    delete[] next;
    delete[] loop_set_prefix_id;
    delete[] prefix;
    delete[] restrict_last;
    delete[] restrict_next;
    delete[] restrict_index;
}

// 获取模式的最大度数
int Schedule::get_max_degree() const{
    int mx = 0;
    for(int i = 0; i < size; ++i) {
        int cnt = 0;
        for(int j = 0; j < size; ++j)
            cnt += adj_mat[INDEX(i,j,size)];
        if(cnt > mx) mx = cnt;
    }
    return mx;
}

// 构造前缀数组以及链式结构
void Schedule:: build_loop_invariant()
{
    // 临时CSR数组:以当前遍历结点为起点,存边的终点
    int* tmp_data = new int[size];
    loop_set_prefix_id[0] = -1;
    for (int i = 1; i < size; ++i)
    {
        int data_size = 0;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)])
                tmp_data[data_size++] = j;
        loop_set_prefix_id[i] = find_father_prefix(data_size, tmp_data);
    }
    assert(total_prefix_num <= size * (size - 1) / 2);
    delete[] tmp_data;
}

// 寻找当前数组data的前缀数组,返回前缀数组在prefix数组中的索引
int Schedule::find_father_prefix(int data_size, const int* data)
{
    if (data_size == 0)
        return -1;
    int num = data[data_size - 1];  // 前缀数组最后一个元素的值,表示一个结点的序号
    for (int prefix_id = last[num]; prefix_id != -1; prefix_id = next[prefix_id])
        if (prefix[prefix_id].equal(data_size, data))
            return prefix_id;
    // 没有在prefix数组中找到和当前数组data相同的数组,
    //则递归的寻找当前数组data的父前缀的在prefix数组中的索引
    // not found, create new prefix and find its father prefix id recursively
    int father = find_father_prefix(data_size - 1, data);
    // 记录当前前缀数组的父前缀索引
    father_prefix_id[total_prefix_num] = father;
    // 更新前缀链表结构
    next[total_prefix_num] = last[num];
    last[num] = total_prefix_num;
    // 拷贝前缀数组到prefix数组
    prefix[total_prefix_num].init(data_size, data);
    ++total_prefix_num;
    return total_prefix_num - 1;
}

// 设置(替换非添加)限制条件
void Schedule::add_restrict(const std::vector< std::pair<int, int> >& restricts)
{
    restrict_pair = restricts;
    for(unsigned int i = 0; i < restrict_pair.size(); ) {
        bool tag = true;
        for(unsigned int j = 0; j < restrict_pair.size(); ++j) {
            // 找到三个不同的限制条件满足:r_i=(a,b),r_j=(a,c),r_k=(c,b)
            if(i != j && restrict_pair[j].first == restrict_pair[i].first) {
                for (unsigned int k = 0; k < restrict_pair.size(); ++k) {
                    if (i != k && j != k && restrict_pair[k].second == restrict_pair[i].second &&
                        restrict_pair[j].second == restrict_pair[k].first) {
                        tag = false;
                        break;
                    }
                }
            }
            if(tag == false) break;
        }
        if(tag == false) {
            // 移除其中一个限制条件
            restrict_pair.erase(restrict_pair.begin() + i);
        }
        else ++i;
    }

    int max_prefix_num = size * (size - 1) / 2;
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));
    total_restrict_num = 0;
    // 以链表的形式存储限制条件
    for (const auto& p : restrict_pair)
    {
        // p.first must be greater than p.second
        // QUE:why p.first>p.second, should p.first<p.second? ANS:这里的大小应该不是指的
        //模式图中限制条件结点的序号大小,而是在实际匹配数据图时,要求限制条件中前1个结点对应的数据图结点
        //序号大于后一结点
        // 记录限制条件的起始结点
        restrict_index[total_restrict_num] = p.first;
        // 限制条件终止结点构造链表
        restrict_next[total_restrict_num] = restrict_last[p.second];
        restrict_last[p.second] = total_restrict_num;
        ++total_restrict_num;
    }
/*
    printf("restrict list:\n");
    for(int i=0;i<size;++i){
        printf("%d: ",i);
        for(int id=restrict_last[i];id!=-1;id=restrict_next[id]){
            printf("%d ",restrict_index[id]);
        }
        putchar('\n');
    }*/
}

// 获取同构置换的数量
int Schedule::get_multiplicity() const{
    std::vector< std::vector<int> > isomorphism_vec = get_isomorphism_vec();
    return isomorphism_vec.size();
}

void Schedule::aggressive_optimize(std::vector< std::pair<int, int> >& ordered_pairs) const
{
    std::vector< std::vector<int> > isomorphism_vec = get_isomorphism_vec();

    std::vector< std::vector< std::vector<int> > > permutation_groups;
    permutation_groups.clear();
    for (const std::vector<int>& v : isomorphism_vec)
        permutation_groups.push_back(calc_permutation_group(v, size));

    ordered_pairs.clear();

    // delete permutation group which contains 1 permutation with 2 elements and some permutation with 1 elements,
    // and record the corresponding restriction.
    for (unsigned int i = 0; i < permutation_groups.size(); )
    {
        int two_element_number = 0;
        std::pair<int, int> found_pair;
        for (const std::vector<int>& v : permutation_groups[i])
            if (v.size() == 2)
            {
                ++two_element_number;
                found_pair = std::pair<int ,int>(v[0], v[1]);
            }
            else if (v.size() != 1)
            {
                two_element_number = -1;
                break;
            }
        if (two_element_number == 1)
        {
            permutation_groups.erase(permutation_groups.begin() + i);
            isomorphism_vec.erase(isomorphism_vec.begin() + i);
            ordered_pairs.push_back(found_pair);
            assert(found_pair.first < found_pair.second);
        }
        else
            ++i;
    }

    Pattern base_dag(size);
    for (const std::pair<int, int>& pair : ordered_pairs)
        base_dag.add_ordered_edge(pair.first, pair.second);

    bool changed = true;
    while (changed && isomorphism_vec.size() != 1)
    {
        // use restrictions to delete other isomophism
        for (unsigned int i = 0; i < isomorphism_vec.size(); )
        {
            Pattern test_dag(base_dag);
            const std::vector<int>& iso = isomorphism_vec[i];
            for (const std::pair<int, int>& pair : ordered_pairs)
                test_dag.add_ordered_edge(iso[pair.first], iso[pair.second]);
            if (test_dag.is_dag() == false) // is not dag means conflict
            {
                permutation_groups.erase(permutation_groups.begin() + i);
                isomorphism_vec.erase(isomorphism_vec.begin() + i);
            }
            else
                ++i;
        }

        changed = false;
        std::pair<int, int> found_pair;
        for (unsigned int i = 0; i < permutation_groups.size(); )
        {
            int two_element_number = 0;
            for (const std::vector<int>& v : permutation_groups[i])
                if (v.size() == 2)
                {
                    ++two_element_number;
                    found_pair = std::pair<int ,int>(v[0], v[1]);
                    break;
                }
            if (two_element_number >= 1)
            {
                permutation_groups.erase(permutation_groups.begin() + i);
                isomorphism_vec.erase(isomorphism_vec.begin() + i);
                assert(found_pair.first < found_pair.second);
                ordered_pairs.push_back(found_pair);
                base_dag.add_ordered_edge(found_pair.first, found_pair.second);
                changed = true;
                break;
            }
            else
                ++i;
        }
    }
    assert(isomorphism_vec.size() == 1);
}

// 获取所有可以让同构消除的多组限制条件
// Schedule::aggressive_optimize(...) can only get one valid restrictions
// but in this function, we try our best to find more restrictions
// WARNING: the restrictions in ordered_pairs_vector may NOT CORRECT
void Schedule::aggressive_optimize_get_all_pairs(std::vector< std::vector< std::pair<int, int> > >& ordered_pairs_vector)
{
    // 与模式同构的置换集
    std::vector< std::vector<int> > isomorphism_vec = get_isomorphism_vec();
//    std::cout << "isomorphism_vec:[\n";
//    for (const auto &rs: isomorphism_vec) {
//        std::cout << " [";
//        for (const auto &x: rs) {
//            std::cout << x << " ";
//        }
//        std::cout << "]\n";
//    }
//    std::cout << "]\n";

    // 模式同构对应的循环置换的乘积的集合
    //permutation_groups[i]表示一个同构对应的循环置换乘积
    //permutation_groups[i][x]为一个循环置换
    std::vector< std::vector< std::vector<int> > > permutation_groups;
    permutation_groups.clear();

    // 对每个同构计算循环置换乘积
    for (const std::vector<int>& v : isomorphism_vec)
        permutation_groups.push_back(calc_permutation_group(v, size));

    ordered_pairs_vector.clear();

    std::vector< std::pair<int,int> > ordered_pairs;    // 对换集合
    ordered_pairs.clear();

    // delete permutation group which contains 1 permutation with 2 elements and some permutation with 1 elements,
    // and record the corresponding restriction.
    for (unsigned int i = 0; i < permutation_groups.size(); )
    {
        int two_element_number = 0;
        std::pair<int, int> found_pair; //对换集合
        // 在当前循环置换乘积中寻找2-cycle即对换
        for (const std::vector<int>& v : permutation_groups[i])
            if (v.size() == 2)
            {
                ++two_element_number;
                found_pair = std::pair<int ,int>(v[0], v[1]);
            }
            else if (v.size() != 1)
            {
                two_element_number = -1;
                break;  // QUE:为什么遇到非1的循环置换直接退出当前同构置换的判断
            }
        // 若对换只有一个
        if (two_element_number == 1)
        {
            // 移除当前循环置换乘积和同构
            permutation_groups.erase(permutation_groups.begin() + i);
            isomorphism_vec.erase(isomorphism_vec.begin() + i);
            // 添加到作为限制条件的对换集合中
            ordered_pairs.push_back(found_pair);
            // 该大小关系通过calc_permutation_group()中计算算法可以保证
            assert(found_pair.first < found_pair.second);
        }
        else
            ++i;
    }
    // 构造一个只有结点的模式,作为限制条件的DAG图
    Pattern base_dag(size);
    // 将所有对换以边的形式添加到模式中
    for (const std::pair<int, int>& pair : ordered_pairs)
        base_dag.add_ordered_edge(pair.first, pair.second);
    // 深度优先的找寻一组限制条件(对换集)
    //结果记录在ordered_pairs_vector中,其中的每个元素是一组能将同构消除的限制条件
    aggressive_optimize_dfs(base_dag, isomorphism_vec, permutation_groups,
                            ordered_pairs, ordered_pairs_vector);

}

// 根据ordered_pairs中对换构造DAG图消除同构,并递归至同构消除,得到的一组限制条件记录到ordered_pairs_vector中
void Schedule::aggressive_optimize_dfs(Pattern base_dag, std::vector< std::vector<int> > isomorphism_vec,
                                       std::vector< std::vector< std::vector<int> > > permutation_groups,
                                       std::vector< std::pair<int,int> > ordered_pairs,
                                       std::vector< std::vector< std::pair<int,int> > >& ordered_pairs_vector) {
    for (unsigned int i = 0; i < isomorphism_vec.size(); )
    {
        Pattern test_dag(base_dag);
        const std::vector<int>& iso = isomorphism_vec[i];   // 同构置换
        // 添加对换置换后的边
        for (const std::pair<int, int>& pair : ordered_pairs)
            test_dag.add_ordered_edge(iso[pair.first], iso[pair.second]);
        // 若添加后成环则移除当前同构及其置换乘积
        if (test_dag.is_dag() == false) // is not dag means conflict
        {
            permutation_groups.erase(permutation_groups.begin() + i);
            isomorphism_vec.erase(isomorphism_vec.begin() + i);
        }
        else
            ++i;
    }
    // 若只剩一个同构,则ordered_pairs里的对换正好可以消除其他同构,
    //则直接添加到结果数组中
    if(isomorphism_vec.size() == 1) {
        ordered_pairs_vector.push_back(ordered_pairs);
        return;
    }

    // 若剩余多个同构
    std::pair<int, int> found_pair;
    // 遍历剩余的置换乘积
    for (unsigned int i = 0; i < permutation_groups.size(); )
    {
        int two_element_number = 0;
        // 遍历每个循环置换
        for (const std::vector<int>& v : permutation_groups[i])
        {
            // 找2-cycle
            if (v.size() == 2)
            {
                ++two_element_number;
                found_pair = std::pair<int ,int>(v[0], v[1]);
                // 记录先前的置换乘积,同构集,对换集和DAG图
                std::vector< std::vector< std::vector<int> > > next_permutation_groups = permutation_groups;
                std::vector< std::vector<int> > next_isomorphism_vec = isomorphism_vec;
                std::vector< std::pair<int,int> > next_ordered_pairs = ordered_pairs;
                Pattern next_base_dag = base_dag;
                // 移除当前置换乘积和同构
                next_permutation_groups.erase(next_permutation_groups.begin() + i);
                next_isomorphism_vec.erase(next_isomorphism_vec.begin() + i);
                assert(found_pair.first < found_pair.second);
                // 添加当前对换到对换集
                next_ordered_pairs.push_back(found_pair);
                // 添加边到DAG图
                next_base_dag.add_ordered_edge(found_pair.first, found_pair.second);
                // 递归找寻做限制条件的对换
                aggressive_optimize_dfs(next_base_dag, next_isomorphism_vec,
                                        next_permutation_groups, next_ordered_pairs,
                                        ordered_pairs_vector);
            }
        }
        // QUE:为什么若有多于一个对换则直接退出遍历置换乘积
        if( two_element_number >= 1) {
            break;
        }
        else {
           ++i;
        }
    }

}

void Schedule::GraphZero_aggressive_optimize(std::vector< std::pair<int, int> >& ordered_pairs) const {
    std::vector< std::vector<int> > Aut;
    GraphZero_get_automorphisms(Aut);

    std::vector< std::pair<int,int> > L;
    L.clear();

    for(int v = 0; v < size; ++v) { // iterate all elements in schedule
        std::vector< std::vector<int> > stabilized_aut;
        stabilized_aut.clear();

        for(int i = 0; i < Aut.size(); ++i) {
            std::vector<int>& x = Aut[i];
            if( x[v] == v) {
                stabilized_aut.push_back(x);
            }
            else {
                int x1  = v, x2 = x[v];
                if( x1 > x2) {
                    int tmp = x1;
                    x1 = x2;
                    x2 = tmp;
                }
                bool tag = true;
                std::pair<int,int> cur = std::make_pair(x1, x2);
                for(int j = 0; j < L.size(); ++j)
                    if( L[j]  == cur) {
                        tag = false;
                        break;
                    }
                if(tag) {
                    L.push_back(cur);
                }
            }
        }
        Aut = stabilized_aut;
    }

    ordered_pairs.clear(); // In GraphZero paper, this vector's name is 'L'

    for(int i = 0; i < L.size(); ++i) {
        bool tag = true;
        for(int j = 0; j < ordered_pairs.size(); ++j)
            if( L[i].second == ordered_pairs[j].second) {
                tag = false;
                if( L[i].first > ordered_pairs[j].first) ordered_pairs[j].first = L[i].first;
                break;
            }
        if(tag) ordered_pairs.push_back(L[i]);
    }
}

void Schedule::GraphZero_get_automorphisms(std::vector< std::vector<int> > &Aut) const {
    int p[size];
    Aut.clear();
    for(int i = 0; i < size; ++i) p[i] = i;
    do{
        bool tag = true;
        for(int i = 0; i < size; ++i) {
            for(int j = 0; j < size; ++j)
                if( adj_mat[INDEX(i, j, size)] != adj_mat[INDEX(p[i], p[j], size)]) {
                    tag = false;
                    break;
                }
            if( !tag ) break;
        }
        if(tag) {
            std::vector<int> tmp;
            tmp.clear();
            for(int i = 0; i < size; ++i) tmp.push_back(p[i]);
            Aut.push_back(tmp);
        }
    } while( std::next_permutation(p, p + size) );

}

// 获取与模式同构的置换集
std::vector< std::vector<int> > Schedule::get_isomorphism_vec() const
{
    unsigned int pow = 1;
    for (int i = 2; i <= size; ++i)
        pow *= i;
    std::vector< std::vector<int> > vec;
    vec.clear();
    bool use[size];
    for (int i = 0; i < size; ++i)
        use[i] = false;
    std::vector<int> tmp_vec;
    // 获取节点的全部排列
    get_full_permutation(vec, use, tmp_vec, 0);
    assert(vec.size() == pow);
    std::vector< std::vector<int> > isomorphism_vec;
    isomorphism_vec.clear();
    for (const std::vector<int>& v : vec)
    {
        bool flag = true;
        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                // 要求结点间的连接性置换前后保持一致
                if (adj_mat[INDEX(i, j, size)] != 0) {
                    if (adj_mat[INDEX(v[i], v[j], size)] == 0) // not isomorphism
                    {
                        flag = false;
                        break;
                    }
                }
            }
        }
        if (flag)
            isomorphism_vec.push_back(v);
    }
    return isomorphism_vec;
}

// 获取所有结点[0,size)的全排列到vec二维数组,用于做候选schedules
void Schedule::get_full_permutation(std::vector< std::vector<int> >& vec,
                                    bool use[], std::vector<int> tmp_vec,
                                    int depth) const
{
    if (depth == size)
    {
        vec.push_back(tmp_vec);
        return;
    }
    for (int i = 0; i < size; ++i)
        if (use[i] == false)
        {
            use[i] = true;
            tmp_vec.push_back(i);
            get_full_permutation(vec, use, tmp_vec, depth + 1);
            tmp_vec.pop_back();
            use[i] = false;
        }
}

// 获取模式同构对应的循环置换的乘积
std::vector< std::vector<int> > Schedule::calc_permutation_group(const std::vector<int> vec, int size)
{
    // vec是模式的同构
    bool use[size];
    for (int i = 0; i < size; ++i)
        use[i] = false;
    std::vector< std::vector<int> > res;
    res.clear();
    for (unsigned int i = 0; i < vec.size(); ++i)
        if (use[i] == false)
        {
            std::vector<int> tmp_vec;
            tmp_vec.clear();
            tmp_vec.push_back(i);
            use[i] = true;
            int x = vec[i];
            while (use[x] == false)
            {
                use[x] = true;
                tmp_vec.push_back(x);
                x = vec[x];
            }
            res.push_back(tmp_vec);
        }
    return res;
}

void Schedule::performance_modeling(int* best_order, std::vector< std::vector<int> > &candidates, int v_cnt, unsigned int e_cnt) {
    int* order;
    int* rank;

    double* p_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];

    double p = e_cnt * 1.0 / v_cnt / v_cnt;

    p_size[0] = v_cnt;
    for(int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p;
    }

    order = new int[size];
    rank = new int[size];

    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];
    for(const std::vector<int>& vec : candidates) {
        for(int i = 0; i < size; ++i)
            order[i] = vec[i];
        // check whether it is valid schedule
        bool is_valid = true;
        for(int i = 1; i < size; ++i) {
            bool have_edge = false;
            for(int j = 0; j < i; ++j)
                if( adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if( have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if( is_valid == false ) continue;

        for(int i = 0; i < size; ++i) rank[order[i]] = i;
        int* cur_adj_mat;
        cur_adj_mat = new int[size*size];
        for(int i = 0; i < size; ++i)
            for(int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector< std::pair<int,int> > restricts;
        //TODO BUG!!!!!
        GraphZero_aggressive_optimize(restricts);
        int restricts_size = restricts.size();
        std::sort(restricts.begin(), restricts.end());
        double* sum;
        sum = new double[restricts_size];
        for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
        int* tmp;
        tmp = new int[size];
        for(int i = 0; i < size; ++i) tmp[i] = i;
        do {
            for(int i = 0; i < restricts_size; ++i)
                if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                    sum[i] += 1;
                }
                else break;
        } while( std::next_permutation(tmp, tmp + size));
        double total = 1;
        for(int i = 2; i <= size; ++i) total *= i;
        for(int i = 0; i < restricts_size; ++i)
            sum[i] = sum[i] /total;
        for(int i = restricts_size - 1; i > 0; --i)
            sum[i] /= sum[i - 1];

        double val = 1;
        for(int i = 0; i < size; ++i) invariant_size[i].clear();
        for(int i = size - 1; i >= 0; --i) {
            int cnt_forward = 0;
            int cnt_backward = 0;
            for(int j = 0; j < i; ++j)
                if(cur_adj_mat[INDEX(j, i, size)])
                    ++cnt_forward;
            for(int j = i + 1; j < size; ++j)
                if(cur_adj_mat[INDEX(j, i, size)])
                    ++cnt_backward;

            int c = cnt_forward;
            for(int j = i - 1; j >= 0; --j)
                if(cur_adj_mat[INDEX(j, i, size)])
                    invariant_size[j].push_back(c--);

            for(int j = 0; j < invariant_size[i].size(); ++j)
                if(invariant_size[i][j] > 1)
                    val += p_size[invariant_size[i][j] - 1] + p_size[1];
            for(int j = 0; j < restricts_size; ++j)
                if(restricts[j].second == i)
                    val *=  sum[j];
            val *= p_size[cnt_forward];

        }
        if( have_best == false || val < min_val) {
            have_best = true;
            for(int i = 0; i < size; ++i)
                best_order[i] = order[i];
            min_val = val;
        }
        delete[] sum;
        delete[] tmp;
        delete[] cur_adj_mat;

    }

    delete[] order;
    delete[] rank;
    delete[] p_size;
}

void Schedule::bug_performance_modeling(int* best_order, std::vector< std::vector<int> > &candidates, int v_cnt, unsigned int e_cnt) {
    int* order;
    int* rank;

    double* p_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];

    double p = e_cnt * 1.0 / v_cnt / v_cnt;

    p_size[0] = v_cnt;
    for(int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p;
    }

    order = new int[size];
    rank = new int[size];

    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];
    for(const std::vector<int>& vec : candidates) {
        for(int i = 0; i < size; ++i)
            order[i] = vec[i];
        // check whether it is valid schedule
        bool is_valid = true;
        for(int i = 1; i < size; ++i) {
            bool have_edge = false;
            for(int j = 0; j < i; ++j)
                if( adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if( have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if( is_valid == false ) continue;

        for(int i = 0; i < size; ++i) rank[order[i]] = i;
        int* cur_adj_mat;
        cur_adj_mat = new int[size*size];
        for(int i = 0; i < size; ++i)
            for(int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector< std::vector< std::pair<int,int> > > restricts_vector;
        restricts_generate(cur_adj_mat, restricts_vector);
        for(int restricts_rank = 0; restricts_rank < restricts_vector.size(); ++restricts_rank) {
            std::vector< std::pair<int,int> >& restricts = restricts_vector[restricts_rank];
            int restricts_size = restricts.size();
            std::sort(restricts.begin(), restricts.end());
            double* sum;
            sum = new double[restricts_size];
            for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
            int* tmp;
            tmp = new int[size];
            for(int i = 0; i < size; ++i) tmp[i] = i;
            do {
                for(int i = 0; i < restricts_size; ++i)
                    if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                        sum[i] += 1;
                    }
                    else break;
            } while( std::next_permutation(tmp, tmp + size));
            double total = 1;
            for(int i = 2; i <= size; ++i) total *= i;
            for(int i = 0; i < restricts_size; ++i)
                sum[i] = sum[i] /total;
            for(int i = restricts_size - 1; i > 0; --i)
                sum[i] /= sum[i - 1];

            double val = 1;
            for(int i = 0; i < size; ++i) invariant_size[i].clear();
            for(int i = size - 1; i >= 0; --i) {
                int cnt_forward = 0;
                int cnt_backward = 0;
                for(int j = 0; j < i; ++j)
                    if(cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_forward;
                for(int j = i + 1; j < size; ++j)
                    if(cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_backward;

                int c = cnt_forward;
                for(int j = i - 1; j >= 0; --j)
                    if(cur_adj_mat[INDEX(j, i, size)])
                        invariant_size[j].push_back(c--);

                for(int j = 0; j < invariant_size[i].size(); ++j)
                    if(invariant_size[i][j] > 1)
                        val += p_size[invariant_size[i][j] - 1] + p_size[1];
                for(int j = 0; j < restricts_size; ++j)
                    if(restricts[j].second == i)
                        val *=  sum[j];
                val *= p_size[cnt_forward];

            }
            if( have_best == false || val < min_val) {
                have_best = true;
                for(int i = 0; i < size; ++i)
                    best_order[i] = order[i];
                min_val = val;
            }
            delete[] sum;
            delete[] tmp;
        }
        delete[] cur_adj_mat;

    }

    delete[] order;
    delete[] rank;
    delete[] p_size;
}

void Schedule::new_performance_modeling(int* best_order, std::vector< std::vector<int> > &candidates, int v_cnt, unsigned int e_cnt, long long tri_cnt) {
    int* order;
    int* rank;

    double* p_size;
    double* pp_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];
    pp_size = new double[max_degree];

    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt;

    p_size[0] = v_cnt;
    for(int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p0;
    }
    pp_size[0] = 1;
    for(int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i-1] * p1;
    }

    order = new int[size];
    rank = new int[size];

    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];
    for(const std::vector<int>& vec : candidates) {
        for(int i = 0; i < size; ++i)
            order[i] = vec[i];
        // check whether it is valid schedule
        bool is_valid = true;
        for(int i = 1; i < size; ++i) {
            bool have_edge = false;
            for(int j = 0; j < i; ++j)
                if( adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if( have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if( is_valid == false ) continue;

        for(int i = 0; i < size; ++i) rank[order[i]] = i;
        int* cur_adj_mat;
        cur_adj_mat = new int[size*size];
        for(int i = 0; i < size; ++i)
            for(int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector< std::vector< std::pair<int,int> > > restricts_vector;
        restricts_generate(cur_adj_mat, restricts_vector);
        for(int restricts_rank = 0; restricts_rank < restricts_vector.size(); ++restricts_rank) {
            std::vector< std::pair<int,int> >& restricts = restricts_vector[restricts_rank];
            int restricts_size = restricts.size();
            std::sort(restricts.begin(), restricts.end());
            double* sum;
            sum = new double[restricts_size];
            for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
            int* tmp;
            tmp = new int[size];
            for(int i = 0; i < size; ++i) tmp[i] = i;
            do {
                for(int i = 0; i < restricts_size; ++i)
                    if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                        sum[i] += 1;
                    }
                    else break;
            } while( std::next_permutation(tmp, tmp + size));
            double total = 1;
            for(int i = 2; i <= size; ++i) total *= i;
            for(int i = 0; i < restricts_size; ++i)
                sum[i] = sum[i] /total;
            for(int i = restricts_size - 1; i > 0; --i)
                sum[i] /= sum[i - 1];

            double val = 1;
            for(int i = 0; i < size; ++i) invariant_size[i].clear();
            for(int i = size - 1; i >= 0; --i) {
                int cnt_forward = 0;
                int cnt_backward = 0;
                for(int j = 0; j < i; ++j)
                    if(cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_forward;
                for(int j = i + 1; j < size; ++j)
                    if(cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_backward;

                int c = cnt_forward;
                for(int j = i - 1; j >= 0; --j)
                    if(cur_adj_mat[INDEX(j, i, size)])
                        invariant_size[j].push_back(c--);

                for(int j = 0; j < invariant_size[i].size(); ++j)
                    if(invariant_size[i][j] > 1)
                        val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
                val += 1;
                for(int j = 0; j < restricts_size; ++j)
                    if(restricts[j].second == i)
                        val *=  sum[j];
                val *= p_size[1] * pp_size[ cnt_forward - 1 ];

            }
            if( have_best == false || val < min_val) {
                have_best = true;
                for(int i = 0; i < size; ++i)
                    best_order[i] = order[i];
                min_val = val;
            }
            delete[] sum;
            delete[] tmp;
        }
        delete[] cur_adj_mat;

    }

    delete[] order;
    delete[] rank;
    delete[] p_size;
    delete[] pp_size;
}

// 容斥定理优化初始化函数
void Schedule::init_in_exclusion_optimize() {
    // 当前排列的Phase2中定义的k值,要求排列最后k个结点不直连
    int optimize_num = in_exclusion_optimize_num;

    assert( in_exclusion_optimize_num > 1);

    int* id;
    id = new int[ optimize_num ];

    // in_exclusion_val[2*n-1]和in_exclusion_val[2*n-2]分别表示
    //在n个结点所构成的所有图的可能中总边数为奇数和偶数的连通图的数目
    int* in_exclusion_val;
    in_exclusion_val = new int[ optimize_num * 2];

    // 遍历k值
    for(int n = 1; n <= optimize_num; ++n) {
        DisjointSetUnion dsu(n);
        int m = n * (n - 1) / 2;    // 完全图的边数

        in_exclusion_val[ 2 * n - 2 ] = 0;
        in_exclusion_val[ 2 * n - 1 ] = 0;

        if( n == 1) {
            ++in_exclusion_val[0];
            continue;
        }

        std::pair<int,int> edge[m];     // 存所有的无向边
        int e_cnt = 0;
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < i; ++j)
                edge[e_cnt++] = std::make_pair(i,j);

        // 一共m条边,每个边用1位表示是否存在
        // s则表示n个结点的一种图的情况, 此处遍历n个结点存在的2^m种情况
        // BUG:s为int类型,sizeof(int)=32位,因此m最多为30,否则会出现溢出错误
        //由于m=n*(n-1)/2,因此n最多为8,超过后下述代码将出错,一般不会出现该情况
        for(int s = 0; s < (1<<m); ++s) {
            dsu.init();
            int bit_cnt = 0;    // 记录该种图的情况下边的数目
            // 遍历所有无向边
            for(int i = 0; i < m; ++i) {
                // 判断某条边是否存在
                if (s & (1 << i)) {
                    ++bit_cnt;
                    dsu.merge(edge[i].first, edge[i].second);   // 连接结点
                }
            }
            // 若最后图是连通的,即只有1个连通分量
            if( dsu.get_set_size() == 1) {
                // 按奇偶性区分
                if( bit_cnt & 1) ++in_exclusion_val[2 * n -1];
                else ++in_exclusion_val[ 2 * n - 2];
            }
        }
    }

    in_exclusion_optimize_group.clear();
    in_exclusion_optimize_val.clear();

    get_in_exclusion_optimize_group(0, id, 0, in_exclusion_val);

    delete[] id;
    delete[] in_exclusion_val;
}

// 获取容斥定理的分组
void Schedule::get_in_exclusion_optimize_group(int depth, int* id, int id_cnt, int* in_exclusion_val) {
    // depth from 0 to k
    // id为大小为in_exclusion_optimize_num(k)的数组
    //id[i]=0~i, i from 0 to id_cnt-1, id_cnt from 0 to in_exclusion_optimize_num

    if( depth == in_exclusion_optimize_num) {
        int* size = new int[id_cnt];
        for(int i = 0; i < id_cnt; ++i)
            size[i] = 0;
        for(int i = 0; i < in_exclusion_optimize_num; ++i)
            size[ id[i] ] ++;
        int val[2];
        val[0] = in_exclusion_val[ size[0] * 2 - 2 ];
        val[1] = in_exclusion_val[ size[0] * 2 - 1 ];
        for(int i = 1; i < id_cnt; ++i) {
            int tmp0 = val[0];
            int tmp1 = val[1];
            // QUE:这么计算的原理
            val[0] = tmp0 * in_exclusion_val[ size[i] * 2 - 2] + tmp1 * in_exclusion_val[ size[i] * 2 - 1];
            val[1] = tmp0 * in_exclusion_val[ size[i] * 2 - 1] + tmp1 * in_exclusion_val[ size[i] * 2 - 2];
        }

        std::vector< std::vector<int> > group;
        group.clear();
        for(int i = 0; i < id_cnt; ++i) {
            std::vector<int> cur;
            cur.clear();
            for(int j = 0; j < in_exclusion_optimize_num; ++j)
                if( id[j] == i) cur.push_back(j);
            group.push_back(cur);
        }

        in_exclusion_optimize_group.push_back(group);
        in_exclusion_optimize_val.push_back( val[0] - val[1] );

        delete[] size;
        return;
    }

    id[depth] = id_cnt;

    get_in_exclusion_optimize_group(depth + 1, id, id_cnt + 1, in_exclusion_val);

    for(int i = 0; i < id_cnt; ++i) {
        id[depth] = i;
        get_in_exclusion_optimize_group(depth + 1, id, id_cnt, in_exclusion_val);
    }
}

void Schedule::print_schedule() const{
    printf("Schedule:\n");
    for(int i = 0; i < size; ++i) {
        for(int j = 0; j < size; ++j)
            printf("%d", adj_mat[INDEX(i,j,size)]);
        puts("");
    }
}

void Schedule::GraphZero_performance_modeling(int* best_order, int v_cnt, unsigned int e_cnt) {
    int* order;
    int* rank;

    double* p_size;
    double* anti_p;
    p_size = new double[size];
    anti_p = new double[size];

    double p = e_cnt * 2.0 / v_cnt / v_cnt;

    printf("fuck p %.6lf\n",p);
    p_size[0] = v_cnt;
    for(int i = 1; i < size; ++i) {
        p_size[i] = p_size[i-1] * p;
        printf("p %d %.6lf\n", i, p_size[i]);
    }
    anti_p[0] = 1;
    for(int i = 1; i < size; ++i) {
        anti_p[i] = anti_p[i-1] * (1-p);
        printf("anti %d %.6lf\n", i, anti_p[i]);
    }

    order = new int[size];
    rank = new int[size];

    for(int i = 0; i < size; ++i) order[i] = i;
    double min_val;
    bool have_best = false;
    do {
        // check whether it is valid schedule
        bool is_valid = true;
        for(int i = 1; i < size; ++i) {
            bool have_edge = false;
            for(int j = 0; j < i; ++j)
                if( adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if( have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if( is_valid == false ) continue;

        for(int i = 0; i < size; ++i) rank[order[i]] = i;
        int* cur_adj_mat;
        cur_adj_mat = new int[size*size];
        for(int i = 0; i < size; ++i)
            for(int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector< std::pair<int,int> > restricts;
        GraphZero_aggressive_optimize(restricts);
        int restricts_size = restricts.size();
        std::sort(restricts.begin(), restricts.end());
        double* sum;
        sum = new double[restricts_size];
        for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
        int* tmp;
        tmp = new int[size];
        for(int i = 0; i < size; ++i) tmp[i] = i;
        do {
            for(int i = 0; i < restricts_size; ++i)
                if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                    sum[i] += 1;
                }
                else break;
        } while( std::next_permutation(tmp, tmp + size));
        double total = 1;
        for(int i = 2; i <= size; ++i) total *= i;
        for(int i = 0; i < restricts_size; ++i)
            sum[i] = sum[i] /total;
        for(int i = restricts_size - 1; i > 0; --i)
            sum[i] /= sum[i - 1];

        double val = 1;
        for(int i = size - 1; i >= 0; --i) {
            int cnt_forward = 0;
            for(int j = 0; j < i; ++j)
                if(cur_adj_mat[INDEX(j, i, size)])
                    ++cnt_forward;

            for(int j = 0; j < restricts_size; ++j)
                if(restricts[j].second == i)
                    val *=  sum[j];
      //      val *= p_size[cnt_forward + 1] * anti_p[i - cnt_forward];
            val *= v_cnt;
            for(int j = 0; j < i - cnt_forward; ++j)
                val *= (1-p);
            for(int j = 0; j < cnt_forward; ++j)
                val *= p;
        }
        if( have_best == false || val <= min_val) {
            have_best = true;
            for(int i = 0; i < size; ++i)
                best_order[i] = order[i];
            min_val = val;
            printf("gz upd %.10lf\n", val);
        }

        delete[] cur_adj_mat;
        delete[] sum;
        delete[] tmp;

    } while( std::next_permutation(order, order + size) );

    delete[] order;
    delete[] rank;
    delete[] p_size;
    delete[] anti_p;
}

void Schedule::restrict_selection(int v_cnt, unsigned int e_cnt, long long tri_cnt, std::vector< std::vector< std::pair<int,int> > > ordered_pairs_vector, std::vector< std::pair<int,int> > &best_restricts) const{
    
    double* p_size;
    double* pp_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];
    pp_size = new double[max_degree];

    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt;
    
    p_size[0] = v_cnt;
    for(int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p0;
    }
    
    pp_size[0] = 1;
    for(int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i-1] * p1;
    }

    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];

    for(int cur_restricts = 0; cur_restricts < ordered_pairs_vector.size(); ++cur_restricts) {
        std::vector< std::pair<int,int> > &restricts = ordered_pairs_vector[cur_restricts];
        int restricts_size = restricts.size();
        std::sort(restricts.begin(), restricts.end());

        double* sum;
        sum = new double[restricts_size];
        for(int i = 0; i < restricts_size; ++i) sum[i] = 0;

        int* tmp;
        tmp = new int[size];
        for(int i = 0; i < size; ++i) tmp[i] = i;

        do {
            for(int i = 0; i < restricts_size; ++i)
                if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                    sum[i] += 1;
                }
                else break;
        } while( std::next_permutation(tmp, tmp + size));

        double total = 1;
        for(int i = 2; i <= size; ++i) total *= i;

        for(int i = 0; i < restricts_size; ++i)
            sum[i] = sum[i] /total;
        for(int i = restricts_size - 1; i > 0; --i)
            sum[i] /= sum[i - 1];

        double val = 1;
        for(int i = 0; i < size; ++i) invariant_size[i].clear();

        for(int i = size - 1; i >= 0; --i) {
            int cnt_forward = 0;
            int cnt_backward = 0;
            for(int j = 0; j < i; ++j)
                if(adj_mat[INDEX(j, i, size)])
                    ++cnt_forward;
            for(int j = i + 1; j < size; ++j)
                if(adj_mat[INDEX(j, i, size)])
                    ++cnt_backward;

            int c = cnt_forward;
            for(int j = i - 1; j >= 0; --j)
                if(adj_mat[INDEX(j, i, size)])
                    invariant_size[j].push_back(c--);

            for(int j = 0; j < invariant_size[i].size(); ++j)
                if(invariant_size[i][j] > 1)
                    val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
            val += 1;
            for(int j = 0; j < restricts_size; ++j)
                if(restricts[j].second == i)
                    val *=  sum[j];
            val *= p_size[1] * pp_size[cnt_forward-1];

        }
        if( have_best == false || val < min_val) {
            have_best = true;
            best_restricts = restricts;
            min_val = val;
        }

        delete[] sum;
        delete[] tmp;

    }

    delete[] p_size;
    delete[] pp_size;
    assert(have_best);
}

// 生成限制条件并进行验证
void Schedule::restricts_generate(const int* cur_adj_mat, std::vector< std::vector< std::pair<int,int> > > &restricts) {
    Schedule schedule(cur_adj_mat, get_size());
    // restricts中为多组限制条件
    schedule.aggressive_optimize_get_all_pairs(restricts);
//    std::cout << "all restricts_vec:[\n";
//    for (const auto &rs: restricts) {
//        std::cout << " [";
//        for (const auto &p: rs) {
//            std::cout << "(" << p.first << ", " << p.second << ") ";
//        }
//        std::cout << "]\n";
//    }
//    std::cout << "]\n";

    int size = schedule.get_size();
    Graph* complete;
    DataLoader* D = new DataLoader();
    // 构造一个结点数为size+1的完全图
    // QUE:此处图结点数为什么是size+1,按照论文应该是size? ANS:结点数不重要,在有无限制条件下匹配的的比值都是同构数
    assert(D->load_complete(complete, size));
    // ans=ans_without/automorphisms_count, ans_without即不使用限制条件进行的完全图模式匹配
    long long all = complete->pattern_matching( schedule, 1);
    long long ans = all / schedule.get_multiplicity();
//    std::cout << "ans: " << ans << std::endl;

    int thread_num = 1;
    for(int i = 0; i < restricts.size(); ) {
        Schedule cur_schedule(schedule.get_adj_mat_ptr(), schedule.get_size());
        // 设置限制条件
        cur_schedule.add_restrict(restricts[i]);
        // ans_with即带限定条件后的完全图模式匹配
        long long cur_ans = complete->pattern_matching( cur_schedule, thread_num);
        // 不相等则不为正确的限制条件则移除
        if( cur_ans != ans) {
            restricts.erase(restricts.begin() + i);
        }
        else {
            ++i;
        }
    }

    delete complete;
    delete D;
}

// 获取当前排列的Phase2中定义的k值,要求排列最后k个结点不直连
int Schedule::get_vec_optimize_num(const std::vector<int> &vec) {
    bool is_valid = true;
    // 判断当前候选排列是否已经移除了不满足Phase1规则的排列
    for(int i = 1; i < size; ++i) {
        bool have_edge = false;
        for(int j = 0; j < i; ++j)
            if(adj_mat[INDEX(vec[i], vec[j], size)]) {
                have_edge = true;
                break;
            }
        if(!have_edge) {
            is_valid = false;
            break;
        }
    }
    if(!is_valid) return -1;
    // 根据Phase2,要求最后k个结点任意两个不直连
    for(int k = 2; k <= size; ++k) {
        bool flag = true;
        for(int i = size - k + 1; i < size; ++i)
            // 置换后结点间有边退出
            if(adj_mat[INDEX(vec[size - k], vec[i], size)]) {
                flag = false;
                break;
            }
        if(!flag) return k - 1;
    }
    assert(0);
    return -1;
}

double Schedule::our_estimate_schedule_restrict(const std::vector<int> &order, const std::vector< std::pair<int,int> > &pairs,
                                                int v_cnt, unsigned int e_cnt, long long tri_cnt) {
    // order:一种schedule
    // pairs:限制条件对,是基于经schedule重排模式图得到的

    int max_degree = get_max_degree();

    double p_size[max_degree];
    double pp_size[max_degree];

    // 在加载图数据时边数e_cnt进行了乘2
    // 论文中的p_1=2*|E_G|/|V_G|^2
    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;
    // 论文中的p_2=tri_cnt*|V_G|/(2*|E_G|)^2
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt;

    // p_size[i]=p_size[i-1]*p0=v_cnt*(p0^i)=|V_G|*(p_1)^i
    //p_size[0]=|V_G|,p_size[1]=|V_G|*p_1
    p_size[0] = v_cnt;
    for(int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p0;
    }
    //pp_size[i]=pp_size[i-1]*p1=p1^i=p_2^i
    pp_size[0] = 1;
    for(int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i-1] * p1;
    }

    int rank[size];
    for(int i = 0; i < size; ++i) rank[order[i]] = i;

    int* cur_adj_mat;   // 根据schedule顺序转换后的模式图
    cur_adj_mat = new int[size*size];
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector< std::pair<int,int> > restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());

    double sum[restricts_size];
    for(int i = 0; i < restricts_size; ++i) sum[i] = 0;

    int tmp[size];
    for(int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for(int i = 0; i < restricts_size; ++i)
            if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                //sum[i]此时记录的是满足第i个限制条件的排列个数
                sum[i] += 1;
            }
            else break;
    } while( std::next_permutation(tmp, tmp + size));   // 遍历所有排列
    // next_permutation()用于获取下一个更大的排列

    // total=1*2*3*...size=(size)!
    double total = 1;
    for(int i = 2; i <= size; ++i) total *= i;
    // sum[i]=(s[i]/total)/(s[i-1]/total)=s[i]/s[i-1]
    // sum[i]此时记录的为在满足前i-1个限制条件的前提下满足第i个限制条件的概率,
    //相当于论文中的(1-f_i)
    for(int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] /total;
    for(int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    std::vector<int> invariant_size[size];
    for(int i = 0; i < size; ++i) invariant_size[i].clear();

    // 性能开销参数,代表当前configuration的性能开销,越小越好
    double val = 1;
    // 嵌套计算cost_i
    //此处i:n-1->0,论文中为n->1 (n==size即结点数)
    for(int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;    // 该结点与比自己序号小的结点的边数
        int cnt_backward = 0;   // 该结点与比自己序号大的阶段的边数
        for(int j = 0; j < i; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;
        for(int j = i + 1; j < size; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_backward;

        int c = cnt_forward;
        // invariant_size[j]数组的大小表示与结点j相连且比j大的结点数量
        // invariant_size[j][k]表示与结点j相连的一个结点
        for(int j = i - 1; j >= 0; --j)
            if(cur_adj_mat[INDEX(j, i, size)])
                invariant_size[j].push_back(c--);

        // +c_i
        //求解的过程为父前缀结点集个结点和当前结点邻域求交集
        //假设以该结点结尾的前缀有m个结点,m=invariant_size[i][j]
        //则父前缀结点集有m-1个结点,根据论文,其结点集是求交集的结果,基数为|V_G|*p_1*p_2^(m-1-1)
        //当前结点的邻域与结点的度数有关,可以用平均度数表示,即2|E_G|/|V_G|=|V_G|*p_1
        //因此以下循环相当于是求所有以该结点为结尾的前缀的结点集的开销
        for(int j = 0; j < invariant_size[i].size(); ++j){
            if(invariant_size[i][j] > 1) {
                //val+=|V_G|*p_1 * p_2^(invariant_size[i][j]-1-1) + |V_G|*p_1
                val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
            }
        }
        val += 1;
        // ×(1-f_i),经限制条件过滤的概率
        for(int j = 0; j < restricts_size; ++j)
            if(restricts[j].second == i)
                val *=  sum[j];
        // ×l_i,即本层循环的基数,
        //也就是该结点的前缀的结点集的大小,假设前缀有m个结点,m=cnt_forward
        //在i!=0时,前缀的结点集是个m结点邻域交集,根据论文,结点集的大小为|V_G|*p_1*p_2^(m-1)
        if( i ) {
            // val*=|V_G|*p_1 * p_2^(cnt_forward-1)
            val *= p_size[1] * pp_size[ cnt_forward - 1 ];
        }
        //在i=0时是最外层循环,即遍历所有数据图结点的循环,循环的基数是结点个数|V_G|
        else {
            // val *= |V_G|
            val *= p_size[0];
        }
    }
    delete[] cur_adj_mat;

    return val;
}

double Schedule::GraphZero_estimate_schedule_restrict(const std::vector<int> &order, const std::vector< std::pair<int,int> > &pairs, int v_cnt, unsigned int e_cnt) {
    int max_degree = get_max_degree();

    double p_size[max_degree];
    double p = e_cnt * 1.0 / v_cnt / v_cnt;

    p_size[0] = v_cnt;
    for(int i = 1;i < max_degree; ++i) {
        p_size[i] = p_size[i-1] * p;
    }

    int rank[size];
    for(int i = 0; i < size; ++i) rank[order[i]] = i;

    int* cur_adj_mat;
    cur_adj_mat = new int[size*size];
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector< std::pair<int,int> > restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());

    double sum[restricts_size];
    for(int i = 0; i < restricts_size; ++i) sum[i] = 0;

    int tmp[size];
    for(int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for(int i = 0; i < restricts_size; ++i)
            if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            }
            else break;
    } while( std::next_permutation(tmp, tmp + size));

    double total = 1;
    for(int i = 2; i <= size; ++i) total *= i;

    for(int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] /total;
    for(int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    std::vector<int> invariant_size[size];
    for(int i = 0; i < size; ++i) invariant_size[i].clear();

    double val = 1;
    for(int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        int cnt_backward = 0;
        for(int j = 0; j < i; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;
        for(int j = i + 1; j < size; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_backward;

        int c = cnt_forward;
        for(int j = i - 1; j >= 0; --j)
            if(cur_adj_mat[INDEX(j, i, size)])
                invariant_size[j].push_back(c--);

        for(int j = 0; j < invariant_size[i].size(); ++j)
            if(invariant_size[i][j] > 1)
                val += p_size[invariant_size[i][j] - 1] + p_size[1];
        for(int j = 0; j < restricts_size; ++j)
            if(restricts[j].second == i)
                val *=  sum[j];
        val *= p_size[cnt_forward];

    }

    delete[] cur_adj_mat;

    return val;
}

double Schedule::Naive_estimate_schedule_restrict(const std::vector<int> &order, const std::vector< std::pair<int,int> > &pairs, int v_cnt, unsigned int e_cnt) {

    double p = e_cnt * 2.0 / v_cnt / v_cnt;

    int rank[size];
    for(int i = 0; i < size; ++i) rank[order[i]] = i;

    int* cur_adj_mat;
    cur_adj_mat = new int[size*size];
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector< std::pair<int,int> > restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());

    double sum[restricts_size];
    for(int i = 0; i < restricts_size; ++i) sum[i] = 0;
    int tmp[size];

    for(int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for(int i = 0; i < restricts_size; ++i)
            if(tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            }
            else break;
    } while( std::next_permutation(tmp, tmp + size));

    double total = 1;
    for(int i = 2; i <= size; ++i) total *= i;

    for(int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] /total;
    for(int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    double val = 1;
    for(int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        for(int j = 0; j < i; ++j)
            if(cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;

        for(int j = 0; j < restricts_size; ++j)
            if(restricts[j].second == i)
                val *=  sum[j];
        val *= v_cnt;
        for(int j = 0; j < i - cnt_forward; ++j)
            val *= (1-p);
        for(int j = 0; j < cnt_forward; ++j)
            val *= p;
    }

    delete[] cur_adj_mat;

    return val;
}

// 按照Phase1规则移除无效排列
void Schedule::remove_invalid_permutation(
        std::vector< std::vector<int> > &candidate_permutations
        ) {
    for(unsigned int i = 0; i < candidate_permutations.size(); ) {
        const auto& vec = candidate_permutations[i];    // 当前置换
        bool tag = true;
        for(int x = 1; x < size; ++x) {
            bool have_edge = false;
            for(int y = 0; y < x; ++y) {
                // 置换后结点间有边就保留
                // Phase1:要求搜索的结点必须与至少一个先前已搜索的结点直连
                if (adj_mat[INDEX(vec[x], vec[y], size)]) {
                    have_edge = true;
                    break;
                }
            }
            if(!have_edge) {
                tag = false;
                break;
            }
        }
        if(tag) {
            ++i;
        }
        else {
            // 若置换后的该结点与前面的都不直连则移除
            candidate_permutations.erase(candidate_permutations.begin() + i);
        }
    }
}

int Schedule::get_in_exclusion_optimize_num_when_not_optimize() {
    std::vector<int> I;
    for(int i = 0; i < size; ++i) I.push_back(i);
    return get_vec_optimize_num(I);
}

// 计算容斥定理优化时的冗余值
//即论文中的x,由于对最后k个结点使用容斥定理时未考虑限制条件产生的同构,通过除以x去除
void Schedule::set_in_exclusion_optimize_redundancy() {
    int tmp = get_in_exclusion_optimize_num();
    if(tmp <= 1) {
        in_exclusion_optimize_redundancy = 1;
    }
    else {
        Graph* complete;
        DataLoader* D = new DataLoader();
        assert(D->load_complete(complete, get_size()));
        delete D;
        in_exclusion_optimize_redundancy = 1;
        // 在有k个结点的情况下使用容斥定理,得到在容斥定理的情况下满足的匹配数量
        long long ans = complete->pattern_matching( *this, 1);
        // 将k个结点值0,即不使用容斥定理时,得到的匹配数量
        set_in_exclusion_optimize_num(0);
        long long true_ans = complete->pattern_matching( *this, 1);
        set_in_exclusion_optimize_num(tmp);
        delete complete;
        // 计算出容易的同构数
        //根据论文考虑限制条件的实际数量应为ans=ans_IEP/x => x = ans_IEP/ans
        in_exclusion_optimize_redundancy = ans / true_ans;
    }
}
