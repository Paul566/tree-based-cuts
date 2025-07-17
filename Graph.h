#ifndef GRAPH_H
#define GRAPH_H

#include <random>
#include <string>
#include <vector>
#include "Tree.h"

class Graph {
    public:
        std::vector<std::vector<int> > adj_list;
        Tree tree;

        explicit Graph(const std::vector<std::vector<int> > &_adj_list,
                       const std::string &_tree_init_type = "random_mst",
                       int seed_for_tree = 239);

        std::pair<std::vector<int>, float> OneRespectedSparsestCut();
        // returns (list of vertices in a cut, sparsity of the cut)

        std::pair<std::vector<int>, int> OneRespectedMincut();
        // returns (list of vertices in a cut, value of the cut)

        std::pair<std::vector<int>, int> OneRespectedBalancedCut(float ratio);
        // returns (list of vertices in a cut, value of the cut)

        int TwoRespectedMinCut();

        // for testing purposes:
        std::pair<std::vector<int>, float> SlowOneRespectedSparsestCut() const;
        std::pair<std::vector<int>, int> SlowOneRespectedMincut() const;
        std::pair<std::vector<int>, int> SlowOneRespectedBalancedCut(float ratio) const;

        std::pair<std::vector<int>, int> SlowTwoRespectedMinCut() const;
        std::pair<std::vector<int>, int> SlowTwoRespectedBalancedCut(float ratio) const;

        void CalculateTree(std::string tree_init_type);

        void PrintGraph() const;

    private:
        std::mt19937 generator;

        void PrepareTree(std::string tree_init_type);

        void InitRandomMST();

        std::vector<std::tuple<int, int, int> > RandomlyWeightedEdgeList();
        // returns a list of edges (weight, vertex1, vertex2), weights are random ints

        void InitRandomSpanningTree();

        int SlowCutSize(const std::vector<int>& cut) const;
        // returns the size of the cut, the cut is a list of vertices
};

#endif //GRAPH_H
