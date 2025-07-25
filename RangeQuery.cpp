//
// Created by alevran on 7/9/2025.
//

#include "RangeQuery.h"

#include <cassert>
#include <climits>

RangeQuery::RangeQuery(const std::vector<int64_t> &nums) {
    n = nums.size();
    t.resize(4*n);
    lazy.resize(4*n);
    build(nums, 1, 0, n-1);

}

void RangeQuery::update(int l, int r, int64_t addend) {
    assert((0 <= l && l <= r && r < n));
    recursive_update(1, 0, n-1, l, r, addend);
}

int64_t RangeQuery::query_min(int l, int r) {
    assert((0 <= l && l <= r && r < n));
    return recursive_query(1, 0, n-1, l, r);
}

void RangeQuery::build(const std::vector<int64_t>& nums, int v, int tl, int tr) {
    if (tl == tr) {
        t[v] = nums[tl];
    } else {
        int tm = (tl + tr) / 2;
        build(nums, v*2, tl, tm);
        build(nums, v*2+1, tm+1, tr);
        t[v] = std::min(t[v*2], t[v*2 + 1]);
    }
}

void RangeQuery::push(int v) {
    t[v*2] += lazy[v];
    lazy[v*2] += lazy[v];
    t[v*2+1] += lazy[v];
    lazy[v*2+1] += lazy[v];
    lazy[v] = 0;
}


void RangeQuery::recursive_update(int v, int tl, int tr, int l, int r, int64_t addend) {
    if (l > r)
        return;
    if (l == tl && tr == r) {
        t[v] += addend;
        lazy[v] += addend;
    } else {
        push(v);
        int tm = (tl + tr) / 2;
        recursive_update(v*2, tl, tm, l, std::min(r, tm), addend);
        recursive_update(v*2+1, tm+1, tr, std::max(l, tm+1), r, addend);
        t[v] = std::min(t[v*2], t[v*2+1]);
    }
}

int64_t RangeQuery::recursive_query(int v, int tl, int tr, int l, int r) {
    if (l > r)
        return INT_MAX;
    if (l == tl && tr == r)
        return t[v];
    push(v);
    int tm = (tl + tr) / 2;
    return std::min(recursive_query(v*2, tl, tm, l, std::min(r, tm)),
               recursive_query(v*2+1, tm+1, tr, std::max(l, tm+1), r));
}

int RangeQuery::get_arg(int l, int r, int64_t v) {
    while (l < r) {
        int m = (l + r) / 2;
        if (query_min(l, m) == v) {
            r = m;
        } else {
            l = m+1;
        }
    }
    assert((query_min(l, l) == v));
    return l;
}
