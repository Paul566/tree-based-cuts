#include "Tree.h"

#include <queue>


Tree::Tree(const std::vector<std::vector<int> > &_adj_list) {
    adj_list = _adj_list;

    InitializeTree();
}

Tree::Tree(std::vector<std::pair<int, int>> _edge_list) {
    const int n = static_cast<int>(_edge_list.size()) + 1;
    adj_list = std::vector(n, std::vector<int>());
    for (auto edge : _edge_list) {
        adj_list[edge.first].push_back(edge.second);
        adj_list[edge.second].push_back(edge.first);
    }

    InitializeTree();
}

Tree::Tree() {

}

void Tree::InitializeTree() {
    // initialize the children and parent vectors
    // parent of the root is -1

    parent = std::vector<int>(adj_list.size());
    children = std::vector<std::vector<int>>(adj_list.size(), std::vector<int>());

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
            }
        }
    }
}
