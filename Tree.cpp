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
    InitializeSubtreeSizes();
}

Tree::Tree() {

}

void Tree::InitializeTreeStructure() {
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

void Tree::InitializeSubtreeSizes() {
    subtree_sizes = std::vector<int>(adj_list.size(), 1);

    std::stack<int> postorder_stack;

    std::stack<int> stack;
    stack.push(root);
    while (!stack.empty()) {
        int node = stack.top();
        stack.pop();
        postorder_stack.push(node);
        for (int child : children[node]) {
            stack.push(child);
        }
    }

    while (!postorder_stack.empty()) {
        int node = postorder_stack.top();
        postorder_stack.pop();
        for (int child : children[node]) {
            subtree_sizes[node] += subtree_sizes[child];
        }
    }
}

void Tree::InitializeHeavyParents() {
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
        for (int i = 0; i < static_cast<int>(children[node].size()); ++i) {
            if (subtree_sizes[children[node][i]] > max_subtree_size) {
                max_subtree_size = subtree_sizes[children[node][i]];
                heaviest_child = children[node][i];
            }
        }
        heavy_parent[heaviest_child] = node;
    }
}
