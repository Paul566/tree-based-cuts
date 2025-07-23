#include "Tree.h"

#include <algorithm>
#include <iostream>
#include <queue>
#include <stack>
#include <stdexcept>
#include <stdint.h>
#include <set>
#include <cassert>

Tree::Tree(int num_vertices, const std::vector<std::pair<int, int> > &_edge_list) {
    if (num_vertices < 2) {
        throw std::runtime_error("In Tree: need at least 2 vertices");
    }

    adj_list = std::vector(num_vertices, std::vector<int>());
    for (auto edge : _edge_list) {
        adj_list[edge.first].push_back(edge.second);
        adj_list[edge.second].push_back(edge.first);
    }

    InitializeTreeStructure();
    InitializePostorderNodes();
    InitializeSubtreeSizes();
    InitializeHeavyParents();
    InitializePathData();
}

Tree::Tree() {
}

void Tree::UpdateDeltaCut(int node1, int node2, int weight) {
    if ((parenting_edge_index[node1] == 0) || (parenting_edge_index[node2] == 0)) {
        cut_size_e0 += weight;
    }

    int lca = LCA(node1, node2);
    UpdateDeltaCutHalfPath(node1, lca, weight);
    UpdateDeltaCutHalfPath(node2, lca, weight);
}

std::pair<std::vector<int>, int64_t> Tree::OneRespectedMincut() const {
    int64_t min_cut_size = INT32_MAX;
    int best_edge = 0;
    for (int i = 0; i < static_cast<int>(cut_values.size()); ++i) {
        if (cut_values[i] < min_cut_size) {
            min_cut_size = cut_values[i];
            best_edge = ordered_edges[i];
        }
    }

    return {SubtreeNodes(best_edge), min_cut_size};
}

std::tuple<std::vector<int>, int64_t, int> Tree::OneRespectedSparsestCut() const {
    int64_t best_cut_size = INT64_MAX;
    int best_denominator = 1;
    int best_edge = 0;
    for (int i = 0; i < static_cast<int>(cut_values.size()); ++i) {
        int64_t this_cut_size = cut_values[i];
        int this_denominator = std::min(subtree_sizes[ordered_edges[i]],
                                        static_cast<int>(adj_list.size()) - subtree_sizes[ordered_edges[i]]);

        if (static_cast<double>(this_cut_size) * static_cast<double>(best_denominator) < static_cast<double>(
            best_cut_size) * static_cast<double>(this_denominator)) {   // compared fractions
            best_cut_size = this_cut_size;
            best_denominator = this_denominator;
            best_edge = ordered_edges[i];
        }
    }

    return {SubtreeNodes(best_edge), best_cut_size, best_denominator};
}

std::pair<std::vector<int>, int64_t> Tree::OneRespectedBalancedCut(float ratio) const {
    int64_t min_cut_size = INT64_MAX;
    int best_edge = -1;
    for (int i = 0; i < static_cast<int>(cut_values.size()); ++i) {
        // TODO make a cleaner comparison
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
        return {std::vector<int>(), INT64_MAX};
    }
    return {SubtreeNodes(best_edge), min_cut_size};
}

void Tree::InitializeTreeStructure() {
    // initialize the children, parent and depth vectors
    // parent of the root is -1
    // if the adj_list gives a disconnected graph, connects the components to the root

    root = 0;
    parent = std::vector<int>(adj_list.size());
    children = std::vector<std::vector<int> >(adj_list.size(), std::vector<int>());
    depth = std::vector<int>(adj_list.size(), 0);
    delta_cut = std::vector<int64_t>(adj_list.size() - 2, 0.);
    cut_size_e0 = 0;
    cut_values = std::vector<int64_t>(adj_list.size() - 1, 0.);

    std::vector<bool> visited(adj_list.size(), false);
    int num_visited = 0;

    for (int component_root = 0; component_root < static_cast<int>(adj_list.size()); ++component_root) {
        if (visited[component_root]) {
            continue;
        }
        std::queue<int> queue;
        queue.push(component_root);
        visited[component_root] = true;
        if (component_root == root) {
            parent[component_root] = -1;
        } else {
            parent[component_root] = root;
            children[root].push_back(component_root);
        }

        while (!queue.empty()) {
            const int node = queue.front();
            queue.pop();
            visited[node] = true;
            ++num_visited;
            for (int neighbor : adj_list[node]) {
                if (!visited[neighbor]) {
                    queue.push(neighbor);
                    children[node].push_back(neighbor);
                    parent[neighbor] = node;
                    depth[neighbor] = depth[node] + 1;
                } else {
                    if (neighbor != parent[node]) {
                        throw std::runtime_error("In Tree: the provided list of edges is not acyclic");
                    }
                }
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

void Tree::UpdateDeltaCutHalfPath(int node, int lca, int weight) {
    while (depth[node] > depth[lca]) {
        if (parenting_edge_index[node] - 1 >= 0) {
            delta_cut[parenting_edge_index[node] - 1] += weight;
        }

        int upperend_node = path_upper_end[path_indices[node]];

        if (depth[upperend_node] <= depth[lca]) {
            if (parenting_edge_index[lca] - 1 >= 0) {
                delta_cut[parenting_edge_index[lca] - 1] -= weight;
            }
            break;
        }

        if (parenting_edge_index[upperend_node] <= static_cast<int>(delta_cut.size()) - 1) {
            delta_cut[parenting_edge_index[upperend_node]] -= weight;
        }
        node = parent[upperend_node];
    }
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

std::pair<std::vector<int>, std::vector<int> > Tree::SubtreesNodes(int vertex1, int vertex2) const {
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
