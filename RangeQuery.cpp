//
// Created by alevran on 7/9/2025.
//

#include "RangeQuery.h"

#include <cassert>
#include <climits>

RangeQuery::RangeQuery(const std::vector<int> &nums) {
    n = nums.size();
    t.resize(4*n);
    lazy.resize(4*n);
    build(nums, 1, 0, n-1);

}

RangeQuery::RangeQuery(const std::map<int, int> &nums) {
    int max_element = nums.rbegin()->first;
    n = max_element + 1;
    std::vector<int> nums_arr(n, INT_MAX);
    for (auto [key, value]: nums) {
        nums_arr[key] = value;
    }
    t.resize(4*n);
    lazy.resize(4*n);
    build(nums_arr, 1, 0, n-1);
}


void RangeQuery::update(int l, int r, int addend) {
    assert((0 <= l && l <= r && r < n));
    recursive_update(1, 0, n-1, l, r, addend);
}

int RangeQuery::query_min(int l, int r) {
    assert((0 <= l && l <= r && r < n));
    return recursive_query(1, 0, n-1, l, r);
}

int RangeQuery::query_min_full() {
    return query_min(0, n-1);
}

int RangeQuery::query_min_without(int i) {
    assert((0 <= i && i < n));
    int res = INT_MAX;
    if (i > 0) {
        res = std::min(res, query_min(0, i-1));
    }
    if (i < n - 1) {
        res = std::min(res, query_min(i+1, n-1));
    }
    return res;
}

void RangeQuery::build(const std::vector<int>& nums, int v, int tl, int tr) {
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


void RangeQuery::recursive_update(int v, int tl, int tr, int l, int r, int addend) {
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

int RangeQuery::recursive_query(int v, int tl, int tr, int l, int r) {
    if (l > r)
        return INT_MAX;
    if (l == tl && tr == r)
        return t[v];
    push(v);
    int tm = (tl + tr) / 2;
    return std::min(recursive_query(v*2, tl, tm, l, std::min(r, tm)),
               recursive_query(v*2+1, tm+1, tr, std::max(l, tm+1), r));
}
