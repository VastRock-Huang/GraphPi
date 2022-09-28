#pragma once

#include "pattern.h"
#include "prefix.h"
#include "disjoint_set_union.h"

#include <vector>

class Schedule
{
public:
    //TODO : more kinds of constructors to construct different Schedules from
    // one Pattern
    Schedule(const Pattern& pattern, bool& is_pattern_valid,
             int performance_modeling_type, int restricts_type,
             bool use_in_exclusion_optimize, int v_cnt, unsigned int e_cnt,
             long long tri_cnt = 0);
    // performance_modeling type = 0 : not use modeling
    //                      type = 1 : use our modeling
    //                      type = 2 : use GraphZero's modeling
    //                      type = 3 : use naive modeling
    // restricts_type = 0 : not use restricts
    //                = 1 : use our restricts
    //                = 2 : use GraphZero's restricts
    Schedule(const int* _adj_mat, int _size);
    ~Schedule();
    //! 获取前缀数组prefix的数量
    inline int get_total_prefix_num() const { return total_prefix_num;}
    //! 获取第prefix_id个前缀的父前缀在prefix数组中的索引
    inline int get_father_prefix_id(int prefix_id) const { return father_prefix_id[prefix_id];}
    //! 获取结点loop的前缀索引
    inline int get_loop_set_prefix_id(int loop) const { return loop_set_prefix_id[loop];}
    //! 获取模式的结点数
    inline int get_size() const { return size;}
    //! 返回以i结点结尾的前缀在prefix数组的索引
    inline int get_last(int i) const { return last[i];}
    //! 返回和索引为i的前缀结尾结点相同的下个前缀索引
    inline int get_next(int i) const { return next[i];}
    //! 返回当前排列的Phase2中定义的k值,要求排列最后k个结点不直连
    inline int get_in_exclusion_optimize_num() const { return in_exclusion_optimize_num;}
    inline void set_in_exclusion_optimize_num(int num) { in_exclusion_optimize_num = num; }
    int get_in_exclusion_optimize_num_when_not_optimize();
    void add_restrict(const std::vector< std::pair<int, int> >& restricts);
    //! 获取限制条件数量
    inline int get_total_restrict_num() const { return total_restrict_num;}
    //! 获取以结点v为终点的限制条件的链表头
    inline int get_restrict_last(int i) const { return restrict_last[i];}
    inline int get_restrict_next(int i) const { return restrict_next[i];}
    //! 返回索引为i的限制条件的起始结点
    inline int get_restrict_index(int i) const { return restrict_index[i];}
    inline long long get_in_exclusion_optimize_redundancy() const { return in_exclusion_optimize_redundancy; }
    int get_max_degree() const;
    int get_multiplicity() const;
    void aggressive_optimize(std::vector< std::pair<int,int> >& ordered_pairs) const;
    void aggressive_optimize_get_all_pairs(std::vector< std::vector< std::pair<int,int> > >& ordered_pairs_vector);
    void aggressive_optimize_dfs(Pattern base_dag, std::vector< std::vector<int> > isomorphism_vec, std::vector< std::vector< std::vector<int> > > permutation_groups, std::vector< std::pair<int,int> > ordered_pairs, std::vector< std::vector< std::pair<int,int> > >& ordered_pairs_vector);
    void restrict_selection(int v_cnt, unsigned int e_cnt, long long tri_cnt, std::vector< std::vector< std::pair<int,int> > > ordered_pairs_vector, std::vector< std::pair<int,int> >& best_restricts) const;
    void restricts_generate(const int* cur_adj_mat, std::vector< std::vector< std::pair<int,int> > > &restricts);

    void GraphZero_aggressive_optimize(std::vector< std::pair<int,int> >& ordered_pairs) const;
    void GraphZero_get_automorphisms(std::vector< std::vector<int> > &Aut) const;

    std::vector< std::vector<int> > get_isomorphism_vec() const;
    static std::vector< std::vector<int> > calc_permutation_group(const std::vector<int> vec, int size);
    inline const int* get_adj_mat_ptr() const {return adj_mat;}
        

    void print_schedule() const;
    //! 容斥定理分组,每个二维数组是一种情况,将k个结点分成几组
    std::vector< std::vector< std::vector<int> > >in_exclusion_optimize_group;
    std::vector< int > in_exclusion_optimize_val;
    //! 限制条件数组,每个限制条件是起始结点和终结点的结点对
    std::vector< std::pair<int,int> > restrict_pair;
private:
    //! 模式的一维化邻接矩阵
    int* adj_mat;
    //! 记录前缀的父前缀在prefix数组的索引
    int* father_prefix_id;
    //! last[v]表示以v结点结尾的前缀在prefix数组的索引,为链式结构的链表头
    int* last;
    //! 辅助last数组,last[p]表示和索引为p的前缀数组prefix[p]结尾结点相同的下一个前缀索引
    int* next;
    //! 结点的前缀对应在prefix数组的索引,loop_set_prefix_id[v]表示结点v的前缀为prefix[loop_set_prefix_id[v]]
    int* loop_set_prefix_id;
    //! 存储所有结点的前缀数组,prefix[p]表示一种前缀数组
    Prefix* prefix;
    //! restrict_last[v]表示以结点v为终点的限制条件的链表头
    int* restrict_last;
    //! 辅助restrict_last数组,restrict_last[r]表示和索引为r的限制条件的终点相同的下一个限制条件索引
    int* restrict_next;
    //! restrict_index[r]表示索引为r的限制条件的起始节点
    int* restrict_index;
    //! 结点数
    int size;
    //! 前缀数组prefix的数量
    int total_prefix_num;
    //! 限制条件数量
    int total_restrict_num;
    //! 当前排列的Phase2中定义的k值,要求排列最后k个结点不直连
    int in_exclusion_optimize_num;
    //! 使用容斥定理未考虑后k个阶段的限制条件而引入的冗余倍数
    long long in_exclusion_optimize_redundancy;

    void build_loop_invariant();
    int find_father_prefix(int data_size, const int* data);
    void get_full_permutation(std::vector< std::vector<int> >& vec, bool use[], std::vector<int> tmp_vec, int depth) const;
    void performance_modeling(int* best_order, std::vector< std::vector<int> > &candidates, int v_cnt, unsigned int e_cnt);
    void bug_performance_modeling(int* best_order, std::vector< std::vector<int> > &candidates, int v_cnt, unsigned int e_cnt);
    void new_performance_modeling(int* best_order, std::vector< std::vector<int> > &candidates, int v_cnt, unsigned int e_cnt, long long tri_cnt);
    void GraphZero_performance_modeling(int* best_order, int v_cnt, unsigned int e_cnt);

    double our_estimate_schedule_restrict(const std::vector<int> &order, const std::vector< std::pair<int,int> > &pairs, int v_cnt, unsigned int e_cnt, long long tri_cnt);
    double GraphZero_estimate_schedule_restrict(const std::vector<int> &order, const std::vector< std::pair<int,int> > &pairs, int v_cnt, unsigned int e_cnt);
    double Naive_estimate_schedule_restrict(const std::vector<int> &order, const std::vector< std::pair<int,int> > &paris, int v_cnt, unsigned int e_cnt);

    void get_in_exclusion_optimize_group(int depth, int* id, int id_cnt, int* in_exclusion_val); 
    //use principle of inclusion-exclusion to optimize
    void init_in_exclusion_optimize();
    
    int get_vec_optimize_num(const std::vector<int> &vec);

    void remove_invalid_permutation(std::vector< std::vector<int> > &candidate_permutations);

    void set_in_exclusion_optimize_redundancy();
};
