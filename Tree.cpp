#include "Tree.h"

#include <queue>
#include <stack>


Tree::Tree(std::vector<std::pair<int, int>> _edge_list) {
    const int n = static_cast<int>(_edge_list.size()) + 1;
    adj_list = std::vector(n, std::vector<int>());
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

void Tree::AddEdge(int node1, int node2) {
    if ((parenting_edge_index[node1] == 0) || (parenting_edge_index[node2] == 0)) {
        ++cut_size_e0;
    }

    while (node1 != node2) {
        if (depth[node1] < depth[node2]) {
            std::swap(node1, node2);
        }

        if (path_indices[node1] == path_indices[node2]) {
            if (parenting_edge_index[node1] - 1 >= 0) {
                ++delta_cut[parenting_edge_index[node1] - 1];
            }
            if (parenting_edge_index[node2] - 1 >= 0) {
                --delta_cut[parenting_edge_index[node2] - 1];
            }
            break;
        }

        if (parenting_edge_index[node1] - 1 >= 0) {
            ++delta_cut[parenting_edge_index[node1] - 1];
        }
        int upperend_node = path_upper_end[path_indices[node1]];
        -- delta_cut[upperend_node];
        node1 = parent[upperend_node];
    }
}

void Tree::InitializeTreeStructure() {
    // initialize the children, parent and depth vectors
    // parent of the root is -1

    root = 0;
    parent = std::vector<int>(adj_list.size());
    children = std::vector<std::vector<int>>(adj_list.size(), std::vector<int>());
    depth = std::vector<int>(adj_list.size(), 0);
    delta_cut = std::vector<int>(adj_list.size() - 1, 0);
    cut_size_e0 = 1;

    std::queue<int> queue;
    std::vector<bool> visited(adj_list.size(), false);
    queue.push(root);
    visited[root] = true;
    parent[root] = -1;
    while (!queue.empty()) {
        int node = queue.front();
        queue.pop();
        visited[node] = true;
        for (int neighbor : adj_list[node]) {
            if (!visited[neighbor]) {
                queue.push(neighbor);
                children[node].push_back(neighbor);
                parent[neighbor] = node;
                depth[neighbor] = depth[node] + 1;
            }
        }
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
            while (heavy_parent[node]) {
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
