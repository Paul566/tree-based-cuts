#include "Tree.h"

#include <algorithm>
#include <iostream>
#include <queue>
#include <stack>
#include <stdexcept>
#include <stdint.h>
#include <set>
#include <cassert>

Tree::Tree(std::vector<std::pair<int, int> > _edge_list) {
    if (_edge_list.empty()) {
        throw std::runtime_error("The edge list is empty");
    }

    adj_list = std::vector(_edge_list.size() + 1, std::vector<int>());
    for (auto edge : _edge_list) {
        adj_list[edge.first].push_back(edge.second);
        adj_list[edge.second].push_back(edge.first);
    }

    InitializeTreeStructure();
    InitializePostorderNodes();
    InitializeSubtreeSizes();
    InitializeHeavyParents();
    InitializePathData();
    tree_edge_to_path_in_edges.resize(adj_list.size() - 1);
    tree_edge_to_path_out_edges.resize(adj_list.size() - 1);
    InitializeWeights();
    global_weight = 0;
}

Tree::Tree() {
}

void Tree::UpdateDeltaCut(int node1, int node2) {
    if ((parent[node1] == node2) || (parent[node2] == node1)) {
        return;
    }

    if ((parenting_edge_index[node1] == 0) || (parenting_edge_index[node2] == 0)) {
        ++cut_size_e0;
    }

    int lca = LCA(node1, node2);
    UpdateDeltaCutHalfPath(node1, lca);
    UpdateDeltaCutHalfPath(node2, lca);
}

void Tree::HandleGraphEdge(const Edge& edge) {
    auto [node1, node2] = edge;
    if ((parent[node1] == node2) || (parent[node2] == node1)) {
        return;
    }

    int lca = LCA(node1, node2);

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

    AddPath(edge, 1);
}

std::pair<std::vector<int>, int> Tree::OneRespectedMincut() {
    int min_cut_size = INT32_MAX;
    int best_edge = 0;
    for (int i = 0; i < static_cast<int>(cut_values.size()); ++i) {
        if (cut_values[i] < min_cut_size) {
            min_cut_size = cut_values[i];
            best_edge = ordered_edges[i];
        }
    }

    return {SubtreeNodes(best_edge), min_cut_size};
}

std::pair<std::vector<int>, float> Tree::OneRespectedSparsestCut() {
    float min_value = INT32_MAX;
    int best_edge = 0;
    for (int i = 0; i < static_cast<int>(cut_values.size()); ++i) {
        const float this_sparsity = static_cast<float>(cut_values[i]) / static_cast<float>(
            std::min(subtree_sizes[ordered_edges[i]],
            static_cast<int>(adj_list.size()) - subtree_sizes[ordered_edges[i]]));
        if (this_sparsity < min_value) {
            min_value = this_sparsity;
            best_edge = ordered_edges[i];
        }
    }

    return {SubtreeNodes(best_edge), min_value};
}

std::pair<std::vector<int>, int> Tree::OneRespectedBalancedCut(float ratio) {
    int min_cut_size = INT32_MAX;
    int best_edge = -1;
    for (int i = 0; i < static_cast<int>(cut_values.size()); ++i) {
        if ((static_cast<float>(subtree_sizes[ordered_edges[i]]) < static_cast<float>(adj_list.size()) * ratio) ||
            (static_cast<float>(adj_list.size()) - static_cast<float>(subtree_sizes[ordered_edges[i]]) <
                static_cast<float>(adj_list.size()) * ratio)) {
            continue;
        }

        if (cut_values[i] < min_cut_size) {
            min_cut_size = cut_values[i];
            best_edge = ordered_edges[i];
        }
    }

    if (best_edge == -1) {
        // no balanced cut
        return {std::vector<int>(), INT32_MAX};
    }
    return {SubtreeNodes(best_edge), min_cut_size};
}

int Tree::TwoRespectedMinCut() {
    int min_cut_size = INT32_MAX;
    //std::pair<int, int> best_edges = {-1, -1};
    for (int t_edge = 0; t_edge < static_cast<int>(ordered_edges.size()); ++t_edge) {
        for (auto g_edge: tree_edge_to_path_in_edges[t_edge]) {
            AddPath(g_edge, -1);
            AddOutPath(g_edge, 1);

        }
        int min_curr_cut = 1 + MinWithout(t_edge);
        min_cut_size = std::min(min_cut_size, min_curr_cut);

        for (auto g_edge: tree_edge_to_path_out_edges[t_edge]) {
            AddPath(g_edge, 1);
            AddOutPath(g_edge, -1);

        }
    }
    return  min_cut_size;
}

void Tree::InitializeTreeStructure() {
    // initialize the children, parent and depth vectors
    // parent of the root is -1

    root = 0;
    parent = std::vector<int>(adj_list.size());
    children = std::vector<std::vector<int> >(adj_list.size(), std::vector<int>());
    depth = std::vector<int>(adj_list.size(), 0);
    delta_cut = std::vector<int>(adj_list.size() - 2, 0);
    cut_size_e0 = 1;
    cut_values = std::vector<int>(adj_list.size() - 1, 0);

    std::queue<int> queue;
    std::vector<bool> visited(adj_list.size(), false);
    int num_visited = 0;
    queue.push(root);
    visited[root] = true;
    parent[root] = -1;
    while (!queue.empty()) {
        int node = queue.front();
        queue.pop();
        visited[node] = true;
        ++num_visited;
        for (int neighbor : adj_list[node]) {
            if (!visited[neighbor]) {
                queue.push(neighbor);
                children[node].push_back(neighbor);
                parent[neighbor] = node;
                depth[neighbor] = depth[node] + 1;
            }
        }
    }

    if (num_visited != static_cast<int>(adj_list.size())) {
        throw std::runtime_error("The tree must be connected");
    }
}

void Tree::InitializePostorderNodes() {
    postorder_nodes.reserve(adj_list.size());

    std::stack<int> stack;
    stack.push(root);
    while (!stack.empty()) {
        int node = stack.top();
        stack.pop();
        postorder_nodes.push_back(node);
        for (int child : children[node]) {
            stack.push(child);
        }
    }
}

void Tree::InitializeSubtreeSizes() {
    subtree_sizes = std::vector<int>(adj_list.size(), 1);

    for (int i = static_cast<int>(postorder_nodes.size()) - 1; i >= 0; --i) {
        for (int child : children[postorder_nodes[i]]) {
            subtree_sizes[postorder_nodes[i]] += subtree_sizes[child];
        }
    }
}

void Tree::InitializeHeavyParents() {
    heavy_parent = std::vector<bool>(adj_list.size(), false);

    std::queue<int> queue;
    queue.push(root);
    while (!queue.empty()) {
        int node = queue.front();
        queue.pop();
        if (children[node].empty()) {
            continue;
        }
        int max_subtree_size = 0;
        int heaviest_child = 0;
        for (int child : children[node]) {
            if (subtree_sizes[child] > max_subtree_size) {
                max_subtree_size = subtree_sizes[child];
                heaviest_child = child;
            }
        }
        heavy_parent[heaviest_child] = true;

        for (int child : children[node]) {
            queue.push(child);
        }
    }
}

void Tree::InitializePathData() {
    ordered_edges.reserve(adj_list.size() - 1);
    path_indices = std::vector(adj_list.size(), -1);

    std::vector<bool> visited(adj_list.size(), false);
    for (int i = static_cast<int>(postorder_nodes.size()) - 1; i >= 0; --i) {
        if ((!visited[postorder_nodes[i]]) && (postorder_nodes[i] != root)) {
            int node = postorder_nodes[i];

            path_lower_end.push_back(node);
            ordered_edges.push_back(node);
            path_indices[node] = static_cast<int>(path_lower_end.size() - 1);
            visited[node] = true;
            while ((heavy_parent[node]) && (depth[node] >= 2)) {
                node = parent[node];

                ordered_edges.push_back(node);
                path_indices[node] = static_cast<int>(path_lower_end.size() - 1);
                visited[node] = true;
            }
            path_upper_end.push_back(node);
        }
    }

    parenting_edge_index = std::vector<int>(adj_list.size(), -1);
    for (int i = 0; i < static_cast<int>(ordered_edges.size()); ++i) {
        parenting_edge_index[ordered_edges[i]] = i;
    }
}

void Tree::InitializeWeights() {

    for (int path : path_indices) {
        if (path < 0 || weights.contains(path)) {
            continue;
        }

        int length = depth[path_lower_end[path]] - depth[path_upper_end[path]] + 1;
        std::vector<int> path_weights(length, 1);
        weights.insert({path, RangeQuery(path_weights)});
    }
}

int Tree::LCA(int node1, int node2) const {
    if ((node1 == root) || (node2 == root)) {
        return root;
    }

    while (path_indices[node1] != path_indices[node2]) {
        if ((node1 == root) || (node2 == root)) {
            return root;
        }
        if (depth[path_upper_end[path_indices[node1]]] < depth[path_upper_end[path_indices[node2]]]) {
            std::swap(node1, node2);
        }
        node1 = parent[path_upper_end[path_indices[node1]]];
    }

    if (depth[node1] < depth[node2]) {
        return node1;
    }
    return node2;
}

void Tree::UpdateDeltaCutHalfPath(int node, int lca) {
    while (depth[node] > depth[lca]) {
        if (parenting_edge_index[node] - 1 >= 0) {
            ++delta_cut[parenting_edge_index[node] - 1];
        }

        int upperend_node = path_upper_end[path_indices[node]];

        if (depth[upperend_node] <= depth[lca]) {
            if (parenting_edge_index[lca] - 1 >= 0) {
                --delta_cut[parenting_edge_index[lca] - 1];
            }
            break;
        }

        if (parenting_edge_index[upperend_node] <= static_cast<int>(delta_cut.size()) - 1) {
            --delta_cut[parenting_edge_index[upperend_node]];
        }
        node = parent[upperend_node];
    }
}

std::vector<std::pair<int, int>> Tree::GetSegmentsFromHalfPath(int node, int lca) const {
    std::vector<std::pair<int, int>> result;
    if (node == lca) {
        return result;
    }

    while (true) {

        int upperend_node = path_upper_end[path_indices[node]];

        // case 1: upperend_node is higher than lca, so we cut segment at lca-1
        // (-1 b/c the last edge's node is before lca)
        if (depth[upperend_node] <= depth[lca]) {
            result.emplace_back(parenting_edge_index[node],
                                parenting_edge_index[lca] - 1);
            break;
        }

        // case 2: upperend_node is lca's son, still last segment
        // (cannot be united with previous as lca may be root)
        if (depth[upperend_node] == depth[lca] + 1) {
            result.emplace_back(parenting_edge_index[node],
                                parenting_edge_index[upperend_node]);
            break;
        }

        // case 3: there are more segments, add this and continue
        result.emplace_back(parenting_edge_index[node],
                            parenting_edge_index[upperend_node]);

        node = parent[upperend_node];
    }
    return result;
}

void Tree::AddPath(const Edge &edge, int weight) {
    auto [node1, node2] = edge;
    int lca = LCA(node1, node2);
    for (int node: {node1, node2}) {
        if (node == lca) {
            continue;
        }

        while (true) {

            int path = path_indices[node];
            int upper_end_node = path_upper_end[path];
            int low_edge = parenting_edge_index[path_lower_end[path]];

            // case 1: upper_end_node is higher than lca, so we update at lca-1
            // (-1 b/c the last edge's node is before lca)
            if (depth[upper_end_node] <= depth[lca]) {
                weights.at(path).update(parenting_edge_index[node] - low_edge,
                    parenting_edge_index[lca] - 1 - low_edge, weight);
                break;
            }

            // case 2: upper_end_node is lca's son, update until upper_end_node
            // (cannot be united with previous as lca may be root)
            if (depth[upper_end_node] == depth[lca] + 1) {
                weights.at(path).update(parenting_edge_index[node] - low_edge,
                    parenting_edge_index[upper_end_node] - low_edge, weight);
                break;
            }

            // case 3: there are more segments, update this and continue
            weights.at(path).update(parenting_edge_index[node] - low_edge,
                                parenting_edge_index[upper_end_node] - low_edge,
                                weight);

            node = parent[upper_end_node];
        }
    }
}

void Tree::AddOutPath(const Edge &edge, int weight) {
    AddPath(edge, -weight);
    global_weight += weight;
}

int Tree::MinWithout(int edge) {
    int res = INT32_MAX;
    int edge_path = path_indices[ordered_edges[edge]];
    for (auto [path, query] : weights) {
        if (path == edge_path) {
            int low_edge = parenting_edge_index[path_lower_end[path]];
            res = std::min(res, query.query_min_without(edge - low_edge));
        } else {
            res = std::min(res, query.query_min_full());
        }
    }
    return res + global_weight;
}


void Tree::UpdateCutValues() {
    cut_values[0] = cut_size_e0;
    for (int i = 0; i < static_cast<int>(delta_cut.size()); ++i) {
        cut_values[i + 1] = cut_values[i] + delta_cut[i];
    }
}

std::vector<int> Tree::SubtreeNodes(const int vertex) const {
    // returns a list of vertices in the subtree, including this vertex

    std::vector<int> subtree;
    subtree.reserve(subtree_sizes[vertex]);

    std::queue<int> queue;
    queue.push(vertex);
    while (!queue.empty()) {
        int node = queue.front();
        queue.pop();
        subtree.push_back(node);
        for (int child : children[node]) {
            queue.push(child);
        }
    }

    return subtree;
}

std::vector<int> Tree::OutSubtreeNodes(int vertex) const {
    // return the vertices outside the subtree
    std::vector<int> result;

    auto subtree = SubtreeNodes(vertex);
    std::set<int> subtree_set(subtree.begin(), subtree.end());
    for (int vertex2 = 0; vertex2 < adj_list.size(); ++vertex2) {
        if (!(subtree_set.count(vertex2))) {
            result.push_back(vertex2);
        }
    }
    return result;
}

bool Tree::IsAncestor(int node1, int node2) const {
    return LCA(node1, node2) == node1;
}

std::pair<std::vector<int>,std::vector<int>> Tree::SubtreesNodes(int vertex1, int vertex2) const {
    // returns a list of vertices in the subtree, including this vertex
    if (vertex1 == vertex2) {
        throw std::runtime_error("vertex1 == vertex2");
    }

    if (IsAncestor(vertex1, vertex2)) {
        return std::make_pair(OutSubtreeNodes(vertex1), SubtreeNodes(vertex2));
    }

    if (IsAncestor(vertex2, vertex1)) {
        return std::make_pair(SubtreeNodes(vertex1), OutSubtreeNodes(vertex2));
    }

    return std::make_pair(SubtreeNodes(vertex1), SubtreeNodes(vertex2));

}


int Tree::GetRoot() const {
    return root;
}

std::vector<int> Tree::GetDepths() const {
    return depth;
}
