#ifndef GRAPH_H
#define GRAPH_H

#include <random>
#include <string>
#include <vector>
#include "Tree.h"

class Graph {
    public:
        std::vector<std::vector<int> > adj_list;
        std::string tree_init_type;
        Tree tree;

        explicit Graph(const std::vector<std::vector<int> > &_adj_list,
                       const std::string &_tree_init_type = "random_mst",
                       int seed = 239);

        std::pair<std::vector<int>, int> OneRespectedSparsestCut();

        std::pair<std::vector<int>, int> OneRespectedMincut();

        std::pair<std::vector<int>, int> OneRespectedBalancedCut(float ratio);

        void PrintGraph() const;

    private:
        std::mt19937 generator;

        void PrepareTree();

        void InitRandomMST();

        std::vector<std::tuple<int, int, int> > RandomlyWeightedEdgeList();

        void InitRandomSpanningTree();
};

#endif //GRAPH_H
