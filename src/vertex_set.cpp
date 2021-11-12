#include "../include/vertex_set.h"
#include <algorithm>

int VertexSet::max_intersection_size = -1;

VertexSet::VertexSet()
:data(nullptr), size(0), allocate(false)
{}

// 初始化分配data动态数组
void VertexSet::init()
{
    if (allocate && data != nullptr)
        size = 0; // do not reallocate
    else
    {
        size = 0;
        allocate = true;
        data = new int[max_intersection_size];
    }
}

//this function is only used for tmp_set in graph.cpp (i.e., init(Graph.max_degree))
void VertexSet::init(int input_size)
{
    if (allocate == true && data != nullptr)
        size = 0;
    else
    {
        size = 0;
        allocate = true;
        data = new int[input_size];
    }
}

// 将内部数组及大小浅拷贝为输入的参数
void VertexSet::init(int input_size, int* input_data)
{
    if (allocate && data != nullptr)
        delete[] data;
    size = input_size;
    data = input_data;
    allocate = false;
}

// 深拷贝输入数组到内部数组
void VertexSet::copy(int input_size, const int* input_data)
{
    size = input_size;
    for(int i = 0; i < input_size; ++i) data[i] = input_data[i];
}

VertexSet::~VertexSet()
{
    if (allocate== true && data != nullptr)
        delete[] data;
}

// 求两结点交集,结果存到内部数组
void VertexSet::intersection(const VertexSet& set0, const VertexSet& set1, int min_vertex, bool clique)
{
    //set0和set1中的结点应该是升序排列的

    // 若两结点集VertexSet相同则拷贝结点集到内部数组
    if (&set0 == &set1)
    {
        copy(set0.get_size(), set0.get_data_ptr());
        return;
    }
    int i = 0;
    int j = 0;
    int size0 = set0.get_size();
    int size1 = set1.get_size();

    // TODO : Try more kinds of calculation.
    // Like
    // while (true)
    //     ..., if (++i == size0) break;
    //     ..., if (++j == size1) break;
    //     ......
    // Maybe we can also use binary search if one set is very small and another is large.
    if (clique) {
        if (set0.get_data(0) >= min_vertex || set1.get_data(0) >= min_vertex)
            return;
    }
    int data0 = set0.get_data(0);
    int data1 = set1.get_data(0);
    if (clique) {   // 对于完全子图结点序号大于min_vertex就退出
        // TODO : Try more kinds of calculation.
        // For example, we can use binary search find the last element which is smaller than min_vertex, and set its index as loop_size.
        while (i < size0 && j < size1) {
            if (data0 < data1) {
                if ((data0 = set0.get_data(++i)) >= min_vertex)
                    break;
            } else if (data0 > data1) {
                if ((data1 = set1.get_data(++j)) >= min_vertex)
                    break;
            } else {
                push_back(data0);
                if ((data0 = set0.get_data(++i)) >= min_vertex)
                    break;
                if ((data1 = set1.get_data(++j)) >= min_vertex)
                    break;
            }
        }
    }
    else {  // 非完全子图 clique==false
        // 将set0和set1中所有相同的结点存入VertexSet内部数组
        while (i < size0 && j < size1) {
            data0 = set0.get_data(i);
            data1 = set1.get_data(j);
            if (data0 < data1)
                ++i;
            else if (data0 > data1)
                ++j;
            else {
                push_back(data0);
                ++i;
                ++j;
            }
        }
    }
}

// 求当前结点集和set1结点集的交集
void VertexSet::intersection_with(const VertexSet& set1) {
    if (this == &set1)
        return;
    const VertexSet& set0 = *this;
    int i = 0;
    int j = 0;
    int size0 = set0.get_size();
    int size1 = set1.get_size();

    // TODO : Try more kinds of calculation.
    // Like
    // while (true)
    //     ..., if (++i == size0) break;
    //     ..., if (++j == size1) break;
    //     ......
    // Maybe we can also use binary search if one set is very small and another is large.
    int data0;
    int data1;
    size = 0;   // 将当前集合清空
    while (i < size0 && j < size1)
    {
        data0 = set0.get_data(i);
        data1 = set1.get_data(j);
        if (data0 < data1)
            ++i;
        else if (data0 > data1)
            ++j;
        else
        {
            push_back(data0);
            ++i;
            ++j;
        }
    }
}

// 求给定prefix_id对应的父前缀的结点集vertex_set[father_id]与某个结点的领域(input_data)的交集
void VertexSet::build_vertex_set(const Schedule& schedule, const VertexSet* vertex_set, int* input_data,
                                 int input_size, int prefix_id, int min_vertex, bool clique)
{
    // input_data: &edge[l],是某个结点在edge数组中的起点. 整个数组表示该结点的所有边的终结点数组
    // input_size: input_data数组的大小,也是该结点的边数
    // min_vertex = -1
    // clique = false

    // 父前缀索引
    int father_id = schedule.get_father_prefix_id(prefix_id);
    if (father_id == -1) {  // 若无父前缀,则初始化为该结点的边数组
        init(input_size, input_data);
    }
    else  {
        init();
        VertexSet tmp_vset;
        tmp_vset.init(input_size, input_data);
        // 求vertex[father_id]和tmp_vset的结点交集
        intersection(vertex_set[father_id], tmp_vset, min_vertex, clique);
    }
}

void VertexSet::insert_ans_sort(int val)
{
    int i;
    for (i = size - 1; i >= 0; --i)
        if (data[i] >= val)
            data[i + 1] = data[i];
        else
        {
            data[i + 1] = val;
            break;
        }
    if (i == -1)
        data[0] = val;
    ++size;
}

// 判断val是否在该结点集中
bool VertexSet::has_data(int val)
{
    for (int i = 0; i < size; ++i)
        if (data[i] == val)
            return true;
    return false;
}

// set0的前size_after_restrict个结点与set1结点集作差,返回作差后集合的大小
int VertexSet::unorderd_subtraction_size(const VertexSet& set0, const VertexSet& set1, int size_after_restrict)
{
    int size0 = set0.get_size();
    int size1 = set1.get_size();
    if (size_after_restrict != -1)
        size0 = size_after_restrict;

    int ret = size0;
    const int* set0_ptr = set0.get_data_ptr();
    for (int j = 0; j < size1; ++j)
        if (std::binary_search(set0_ptr, set0_ptr + size0, set1.get_data(j)))
            --ret;
    return ret;
}
