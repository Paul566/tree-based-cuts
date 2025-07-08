#include "Graph.h"

#include <iostream>

Graph::Graph(const std::vector<std::vector<int>> &_adj_list) {
    adj_list = _adj_list;
}

void Graph::PrintGraph() const {
    for (int i = 0; i < adj_list.size(); ++i) {
        std::cout << i << ":";
        for (int j : adj_list[i]) {
            std::cout << " " << j;
        }
        std::cout << std::endl;
    }
}
