#include <iostream>
#include <vector>
#include "Graph.h"
#include "RangeQuery.h"
#include "Tree.h"


int main() {

    std::vector<std::vector<int>> adj_list;
    adj_list.emplace_back(std::vector<int>{1, 2, 3});
    adj_list.emplace_back(std::vector<int>{0, 3});
    adj_list.emplace_back(std::vector<int>{0, 3});
    adj_list.emplace_back(std::vector<int>{0, 1, 2});

    Graph graph(adj_list);
    graph.PrintGraph();

    std::vector<int> nums({2, 3, 5, 3, 4});
    RangeQuery rangeQuery(nums);
    std::cout << "min: " << rangeQuery.query_min(0, 4) << std::endl;
    rangeQuery.update(0, 2, 2);
    std::cout << "min after update: " << rangeQuery.query_min(0, 4) << '\n' << std::endl;
    return 0;
}