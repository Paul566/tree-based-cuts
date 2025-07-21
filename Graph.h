#ifndef GRAPH_H
#define GRAPH_H

#include <map>
#include <random>
#include <string>
#include <vector>
#include "Tree.h"

class Graph {
    public:
        std::vector<std::vector<std::pair<int, float>> > adj_list;    // an edge is a pair (destination, weight)
        Tree tree;
        std::map<std::pair<int, int>, float> weights;

        explicit Graph(const std::vector<std::vector<std::pair<int, float>> > &_adj_list,
                       const std::string &_tree_init_type = "random_mst",
                       int seed_for_tree = 239);

        explicit Graph(const std::vector<std::vector<int> > &_adj_list,
                       const std::string &_tree_init_type = "random_mst",
                       int seed_for_tree = 239);

        void FillWeights();

        Graph() {};

        std::pair<std::vector<int>, float> OneRespectedSparsestCut();
        // returns (list of vertices in a cut, sparsity of the cut)

        std::pair<std::vector<int>, float> OneRespectedMincut();
        // returns (list of vertices in a cut, value of the cut)

        std::pair<std::vector<int>, float> OneRespectedBalancedCut(float ratio);
        // returns (list of vertices in a cut, value of the cut)

        float TwoRespectedMinCut();

        // for testing purposes:
        std::pair<std::vector<int>, float> SlowOneRespectedSparsestCut() const;
        std::pair<std::vector<int>, float> SlowOneRespectedMincut() const;
        std::pair<std::vector<int>, float> SlowOneRespectedBalancedCut(float ratio) const;

        std::pair<std::vector<int>, float> SlowTwoRespectedMinCut() const;
        std::pair<std::vector<int>, float> SlowTwoRespectedBalancedCut(float ratio) const;

        void CalculateTree(std::string tree_init_type);

        void PrintGraph() const;

    private:
        std::mt19937 generator;

        void InitMST(std::vector<std::tuple<float, int, int> > weighted_edges);

        std::vector<std::tuple<float, int, int> > RandomlyWeightedEdgeList();
        // returns a list of edges (weight, vertex1, vertex2), weights are random ints

        std::vector<std::tuple<float, int, int>> NegativelyWeightedEdgeList();

        void InitRandomSpanningTree();

        int RandomNeighbor(int vertex);

        float SlowCutSize(const std::vector<int>& cut) const;
        // returns the size of the cut, the cut is a list of vertices
};

#endif //GRAPH_H
