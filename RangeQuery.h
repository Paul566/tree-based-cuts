//
// Created by alevran on 7/9/2025.
//

#ifndef RANGEQUERY_H
#define RANGEQUERY_H
#include <vector>


class RangeQuery {
public:
    explicit RangeQuery(const std::vector<float>& nums);
    void update(int l, int r, float addend);
    float query_min(int l, int r);
    int get_arg(int l, int r, float value, float tolerance);

private:
    std::vector<float> t;
    std::vector<float> lazy;
    int n;
    void build(const std::vector<float>& nums, int v, int tl, int tr);
    void push(int v);
    void recursive_update(int v, int tl, int tr, int l, int r, float addend);
    float recursive_query(int v, int tl, int tr, int l, int r);

};



#endif //RANGEQUERY_H
