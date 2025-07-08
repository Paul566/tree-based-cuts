#include <iostream>
#include <vector>
#include "Graph.h"
#include "Tree.h"


int main() {

    std::vector<std::vector<int>> adj_list;
    adj_list.emplace_back(std::vector<int>{1, 2, 3});
    adj_list.emplace_back(std::vector<int>{0, 3});
    adj_list.emplace_back(std::vector<int>{0, 3});
    adj_list.emplace_back(std::vector<int>{0, 1, 2});

    Graph graph(adj_list);
    graph.PrintGraph();

    return 0;
}