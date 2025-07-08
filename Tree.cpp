#include "Tree.h"


Tree::Tree(const std::vector<std::vector<int> > &_adj_list) {
    adj_list = _adj_list;
}

Tree::Tree(std::vector<std::pair<int, int>> _edge_list) {
    const int n = static_cast<int>(_edge_list.size()) + 1;
    adj_list = std::vector(n, std::vector<int>());
    for (auto edge : _edge_list) {
        adj_list[edge.first].push_back(edge.second);
        adj_list[edge.second].push_back(edge.first);
    }
}

Tree::Tree() {
    adj_list = std::vector<std::vector<int>>();
}
