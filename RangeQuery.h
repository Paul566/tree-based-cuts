//
// Created by alevran on 7/9/2025.
//

#ifndef RANGEQUERY_H
#define RANGEQUERY_H
#include <vector>


class RangeQuery {
public:
    explicit RangeQuery(const std::vector<long>& nums);
    void update(int l, int r, long addend);
    long query_min(int l, int r);

private:
    std::vector<long> t;
    std::vector<long> lazy;
    int n;
    void build(const std::vector<long>& nums, int v, int tl, int tr);
    void push(int v);
    void recursive_update(int v, int tl, int tr, int l, int r, long addend);
    long recursive_query(int v, int tl, int tr, int l, int r);

};



#endif //RANGEQUERY_H
