//
// Created by alevran on 7/9/2025.
//

#ifndef RANGEQUERY_H
#define RANGEQUERY_H
#include <map>
#include <vector>


class RangeQuery {
public:
    explicit RangeQuery(const std::vector<int>& nums);
    explicit RangeQuery(const std::map<int, int>& nums);
    void update(int l, int r, int addend);
    int query_min(int l, int r);

    int query_min_full();

    int query_min_without(int i);

private:
    std::vector<int> t;
    std::vector<int> lazy;
    int n;
    void build(const std::vector<int>& nums, int v, int tl, int tr);
    void push(int v);
    void recursive_update(int v, int tl, int tr, int l, int r, int addend);
    int recursive_query(int v, int tl, int tr, int l, int r);

};



#endif //RANGEQUERY_H
