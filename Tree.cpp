#include "Tree.h"

#include <queue>
#include <stack>
#include <stdexcept>
#include <stdint.h>

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

int Tree::GetRoot() const {
    return root;
}
