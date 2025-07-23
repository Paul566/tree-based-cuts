//
// Created by alevran on 7/9/2025.
//

#ifndef RANGEQUERY_H
#define RANGEQUERY_H
#include <cstdint>
#include <vector>


class RangeQuery {
public:
    explicit RangeQuery(const std::vector<int64_t>& nums);
    void update(int l, int r, int64_t addend);
    int64_t query_min(int l, int r);
    int get_arg(int l, int r, int64_t v);

private:
    std::vector<int64_t> t;
    std::vector<int64_t> lazy;
    int n;
    void build(const std::vector<int64_t>& nums, int v, int tl, int tr);
    void push(int v);
    void recursive_update(int v, int tl, int tr, int l, int r, int64_t addend);
    int64_t recursive_query(int v, int tl, int tr, int l, int r);

};



#endif //RANGEQUERY_H
