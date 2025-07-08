#ifndef GRAPH_H
#define GRAPH_H

#include <random>
#include <string>
#include <vector>
#include "Tree.h"


class Graph {
    public:
        std::vector<std::vector<int>> adj_list;
        Tree tree;

        explicit Graph(const std::vector<std::vector<int>> &_adj_list, const std::string &tree_init_type="random_mst", int seed=239);

        std::vector<int> SparsestOneRespectedCut();

        void PrintGraph() const;

    private:
        std::mt19937 generator;

        void InitRandomMST();

        std::vector<std::tuple<int, int, int>> RandomlyWeightedEdgeList();

        void InitRandomSpanningTree();
};



#endif //GRAPH_H
