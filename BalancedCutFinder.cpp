//
// Created by alevran on 7/18/2025.
//

#include "BalancedCutFinder.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <set>

BalancedCutFinder::BalancedCutFinder(const Graph& graph, const Tree& tree, float factor) :
graph(graph), tree(tree), factor(factor) {
    n = graph.adj_list.size();
    InitializeWeights();
    global_weight = 0;
    InitializeSegments();
    for (int edge = 0; edge < n-1; edge++) {
        edge_to_subtree_size.emplace_back(tree.subtree_sizes[tree.ordered_edges[edge]]);
    }
}


void BalancedCutFinder::InitializeWeights() {
    std::set<int> seen;

    const int bound = static_cast<int>(n * factor / 2);

    for (int path : tree.path_indices) {
        if (path < 0 || seen.contains(path)) {
            continue;
        }
        seen.insert(path);

        int max_tree_subtree = tree.subtree_sizes[tree.path_upper_end[path]];
        if (max_tree_subtree < bound) {
            // no need to maintain weights array for small enough weights
            continue;
        }

        for (int t_edge_start = tree.parenting_edge_index[tree.path_lower_end[path]];
        t_edge_start <= tree.parenting_edge_index[tree.path_upper_end[path]]; t_edge_start++) {
            if (tree.subtree_sizes[tree.ordered_edges[t_edge_start]] >= bound) {
                std::vector<int64_t> path_weights;
                for (int t_edge = t_edge_start; t_edge <= tree.parenting_edge_index[tree.path_upper_end[path]]; t_edge++ ) {
                    int node1 = tree.ordered_edges[t_edge];
                    int node2 = tree.parent[node1];
                    int64_t t_weight = graph.weights.at({node1, node2});
                    path_weights.emplace_back(t_weight);
                }
                weights.insert({path, RangeQuery(path_weights)});
                paths_gaps.insert({path, t_edge_start});
                path_array_length.insert({path, path_weights.size()});
                break;
            }
        }


    }
}

void BalancedCutFinder::InitializeSegments() {
    tree_edge_to_path_in_edges.resize(tree.adj_list.size() - 1);
    tree_edge_to_path_out_edges.resize(tree.adj_list.size() - 1);
    for (int i = 0; i < static_cast<int>(graph.adj_list.size()); ++i) {
        for (auto [j, weight] : graph.adj_list[i]) {
            if (i < j) {
                HandleGraphEdge({i, j, weight});
            }
        }
    }
}

void BalancedCutFinder::HandleGraphEdge(const Edge& edge) {
    auto [node1, node2, weight] = edge;
    if ((tree.parent[node1] == node2) || (tree.parent[node2] == node1)) {
        return;
    }

    int lca = tree.LCA(node1, node2);

    std::vector<std::pair<int, int>> segments;
    auto segments1 = GetSegmentsFromHalfPath(node1, lca);
    auto segments2 = GetSegmentsFromHalfPath(node2, lca);
    segments.reserve(segments1.size() + segments2.size());
    segments.insert(segments.end(), segments1.begin(), segments1.end());
    segments.insert(segments.end(), segments2.begin(), segments2.end());
    // TODO: check if there is more efficient way to do it

    std::sort(segments.begin(), segments.end());

    //sanity check, to be removed
    for (int i = 0; i < segments.size(); ++i) {
        assert((segments[i].first <= segments[i].second));
        if (i < segments.size() - 1) {
            assert((segments[i].second < segments[i+1].first));
        }
    }

    for (int big_segment_start = 0; big_segment_start < segments.size(); ++big_segment_start) {

        int big_segment_end = big_segment_start;
        while (big_segment_end < segments.size() - 1 &&
            segments[big_segment_end].second == segments[big_segment_start].first - 1) {
            ++big_segment_end;
            }


        tree_edge_to_path_in_edges[segments[big_segment_start].first].push_back(edge);

        tree_edge_to_path_out_edges[segments[big_segment_end].second].push_back(edge);

    }

    AddPath(edge, false);
}

std::vector<std::pair<int, int>> BalancedCutFinder::GetSegmentsFromHalfPath(int node, int lca) const {
    std::vector<std::pair<int, int>> result;
    if (node == lca) {
        return result;
    }

    while (true) {

        int upperend_node = tree.path_upper_end[tree.path_indices[node]];

        // case 1: upperend_node is higher than lca, so we cut segment at lca-1
        // (-1 b/c the last edge's node is before lca)
        if (tree.depth[upperend_node] <= tree.depth[lca]) {
            result.emplace_back(tree.parenting_edge_index[node],
                                tree.parenting_edge_index[lca] - 1);
            break;
        }

        // case 2: upperend_node is lca's son, still last segment
        // (cannot be united with previous as lca may be root)
        if (tree.depth[upperend_node] == tree.depth[lca] + 1) {
            result.emplace_back(tree.parenting_edge_index[node],
                                tree.parenting_edge_index[upperend_node]);
            break;
        }

        // case 3: there are more segments, add this and continue
        result.emplace_back(tree.parenting_edge_index[node],
                            tree.parenting_edge_index[upperend_node]);

        node = tree.parent[upperend_node];
    }
    return result;
}

void BalancedCutFinder::AddPath(const Edge &edge, bool negative) {
    auto [node1, node2, weight] = edge;
    if (negative) {
        weight = -weight;
    }
    int lca = tree.LCA(node1, node2);
    for (int node: {node1, node2}) {
        if (node == lca) {
            continue;
        }

        while (true) {

            int path = tree.path_indices[node];
            int upper_end_node = tree.path_upper_end[path];

            // case 1: upper_end_node is higher than lca, so we update at lca-1
            // (-1 b/c the last edge's node is before lca)
            if (tree.depth[upper_end_node] <= tree.depth[lca]) {
                Update(path, tree.parenting_edge_index[node], tree.parenting_edge_index[lca] - 1, weight);
                break;
            }

            // case 2: upper_end_node is lca's son, update until upper_end_node
            // (cannot be united with previous as lca may be root)
            if (tree.depth[upper_end_node] == tree.depth[lca] + 1) {
                Update(path, tree.parenting_edge_index[node], tree.parenting_edge_index[upper_end_node], weight);
                break;
            }

            // case 3: there are more segments, update this and continue
            Update(path, tree.parenting_edge_index[node], tree.parenting_edge_index[upper_end_node], weight);

            node = tree.parent[upper_end_node];
        }
    }
}

void BalancedCutFinder::AddOutPath(const Edge &edge, bool negative) {
    AddPath(edge, !negative);
    auto [_, __, weight] = edge;
    if (negative) {
        weight = -weight;
    }
    global_weight += weight;
}

void BalancedCutFinder::Update(int path, int l, int r, int weight) {
    if (!weights.contains(path)) {
        return;
    }
    int gap = paths_gaps[path];
    if (r < gap) {
        return;
    }
    weights.at(path).update(std::max(l - gap, 0), r - gap, weight);

}

std::pair<int64_t, PathChunck> BalancedCutFinder::MinCutWithEdge(int edge) {
    std::pair<int64_t, PathChunck> res = {INT64_MAX, {-1, -1, -1}};
    int node = tree.ordered_edges[edge];
    int edge_subtree = tree.subtree_sizes[tree.ordered_edges[edge]];
    int lca_edge;
    PathChunck path_chuck;


    for (auto [path, query] : weights) {
        int min_edge = tree.parenting_edge_index[tree.path_lower_end[path]];
        int max_edge = tree.parenting_edge_index[tree.path_upper_end[path]];
        int lca_node_lower = tree.LCA(node, tree.path_lower_end[path]);
        int lca_node_upper = tree.LCA(node, tree.path_upper_end[path]);

        if (lca_node_lower == node && lca_node_upper == tree.path_upper_end[path]) {
            // case 1: node is inside the path
            path_chuck = {path, min_edge, edge-1};
            res = std::min(res, MinOfPath(path_chuck, n-edge_subtree, true, query));
            path_chuck = {path, edge + 1, max_edge};
            res = std::min(res, MinOfPath(path_chuck, edge_subtree, false, query));
        } else if (lca_node_upper == node) {
            //case 2: path is fully contained in subtree
            path_chuck = {path, min_edge, max_edge};
            res = std::min(res, MinOfPath(path_chuck, n-edge_subtree, true, query));
        } else if (lca_node_lower == tree.path_lower_end[path]) {
            // case 3: node is in the subtree of path
            path_chuck = {path, min_edge, max_edge};
            res = std::min(res, MinOfPath(path_chuck, edge_subtree, false, query));
        } else if ((lca_edge = tree.parenting_edge_index[lca_node_lower]) <= max_edge &&
            lca_edge >= min_edge) {
            // case 4: node is descant of something inside the path
            path_chuck = {path, min_edge, lca_edge - 1};
            res = std::min(res, MinOfPath(path_chuck, edge_subtree, true, query));
            path_chuck = {path, lca_edge, max_edge};
            res = std::min(res, MinOfPath(path_chuck, edge_subtree, false, query));
        } else {
            // case 5: none of the above
            path_chuck = {path, min_edge, max_edge};
            res = std::min(res, MinOfPath(path_chuck, edge_subtree, true, query));
        }


    }
    if (res.first == INT64_MAX) {
        return res;
    }
    auto [best, best_path] = res;
    return {best + global_weight, best_path};
}

std::pair<int64_t, PathChunck> BalancedCutFinder::MinOfPath(PathChunck path_chunck, int other_subtree,
    bool is_in_subtree, RangeQuery& query) {

    auto [path, min_edge, max_edge] = path_chunck;
    // n*factor <= other_subtree+x <= n-n*factor, x>=other_subtree
    int min_subtree = std::max(other_subtree, static_cast<int>(ceil(n*factor - other_subtree)));
    int max_subtree = n-n*factor - other_subtree;
    if (max_subtree < min_subtree) {
        return {INT64_MAX, {-1, -1, -1}};
    }
    GetMinMaxEdge(min_edge, max_edge, min_subtree, max_subtree, is_in_subtree);

    if (min_edge <= max_edge) {
        return {query.query_min(min_edge - paths_gaps[path], max_edge - paths_gaps[path]),
            {path, min_edge, max_edge}};
    }
    return {INT64_MAX, {-1, -1, -1}};
}

int BinarySearch(const std::vector<int>& nums, int left, int right, int bound, bool upper) {
    // given that nums (in range [left,right] is strickly increasing
    // if upper: return minimum i s.t. left <= i <= right and nums[i]>=bound (or right+1 if non exist)
    // otherwise return maximum i s.t. left <= i <= right and nums[i]<=bound (or left-1 if non exist)
    assert(0 <= left && left <= right && right < nums.size());
    if (upper && nums[right] < bound) {
        return right + 1;
    }
    if (!upper && nums[left] > bound) {
        return left - 1;
    }
    while (right - left > 1) {
        int mid = left + (right - left) / 2;
        if (nums[mid] > bound) {
            right = mid;
        } else {
            left = mid;
        }
    }
    if (upper) {
        return nums[left] >= bound ? left : left + 1;
    }
    return nums[right] <= bound ? right : right - 1;
}


void BalancedCutFinder::GetMinMaxEdge(int &min_edge, int &max_edge, int min_subtree,
    int max_subtree, bool is_in_subtree) {
    if (min_edge > max_edge) {
        return;
    }
    if (is_in_subtree) {
        min_edge = BinarySearch(edge_to_subtree_size, min_edge, max_edge, min_subtree, true);
        if (min_edge > max_edge) {
            return;
        }
        max_edge = BinarySearch(edge_to_subtree_size, min_edge, max_edge, max_subtree, false);
    } else {
        min_edge = BinarySearch(edge_to_subtree_size, min_edge, max_edge, n-max_subtree, true);
        if (min_edge > max_edge) {
            return;
        }
        max_edge = BinarySearch(edge_to_subtree_size, min_edge, max_edge, n-min_subtree, false);
    }

}

std::pair<std::vector<int>, int64_t> BalancedCutFinder::TwoRespectedBalancedCut() {
    std::pair<int64_t, std::pair<int, int>> res = {INT64_MAX, {-1, -1}};

    for (int t_edge = 0; t_edge < static_cast<int>(tree.ordered_edges.size()); ++t_edge) {
        for (auto g_edge: tree_edge_to_path_in_edges[t_edge]) {
            AddPath(g_edge, true);
            AddOutPath(g_edge, false);
        }

        auto [min_curr_cut, path_chunk] = MinCutWithEdge(t_edge);
        if (min_curr_cut != INT64_MAX) {
            int64_t array_val = min_curr_cut - global_weight;
            int node1 = tree.ordered_edges[t_edge];
            int node2 = tree.parent[node1];
            min_curr_cut += graph.weights.at({node1, node2}); // weight of t_edge
            if (min_curr_cut < res.first) {
                res.first = min_curr_cut;
                auto [path, min_edge, max_edge] = path_chunk;
                int t_edge2 = weights.at(path).get_arg(min_edge - paths_gaps[path],
                    max_edge - paths_gaps[path], array_val) + paths_gaps[path];
                res.second = {t_edge, t_edge2};
            }

        }

        for (auto g_edge: tree_edge_to_path_out_edges[t_edge]) {
            AddPath(g_edge, false);
            AddOutPath(g_edge, true);

        }
    }

    auto [cut_size, edges] = res;
    auto [edge1, edge2] = edges;
    if (edge1 == -1) {
        assert(cut_size == INT64_MAX);
        return {{}, INT64_MAX};
    }
    int node1 = tree.ordered_edges[edge1];
    int node2 = tree.ordered_edges[edge2];
    auto sub_trees = tree.SubtreesNodes(node1, node2);
    auto cut = sub_trees.first;
    cut.insert(cut.end(), sub_trees.second.begin(), sub_trees.second.end());
    return  {cut, cut_size};
}
