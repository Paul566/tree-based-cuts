#ifndef GRAPH_H
#define GRAPH_H

#include <map>
#include <random>
#include <string>
#include <vector>
#include "Tree.h"

class Graph {
    public:
        std::vector<std::vector<std::pair<int, int>> > adj_list;    // an edge is a pair (destination, weight)
        Tree tree;
        std::map<std::pair<int, int>, int> weights;

        explicit Graph(const std::vector<std::vector<std::pair<int, int>> > &_adj_list,
                       const std::string &_tree_init_type = "random_mst",
                       int seed_for_tree = 239);

        explicit Graph(const std::vector<std::vector<int> > &_adj_list,
                       const std::string &_tree_init_type = "random_mst",
                       int seed_for_tree = 239);

        void FillWeights();

        Graph() {};

        std::tuple<std::vector<int>, int64_t, int> OneRespectedSparsestCut() const;
        // returns (list of vertices in a cut, sparsity of the cut)

        std::pair<std::vector<int>, int64_t> OneRespectedMincut() const;
        // returns (list of vertices in a cut, value of the cut)

        std::pair<std::vector<int>, int64_t> OneRespectedBalancedCut(float ratio) const;
        // returns (list of vertices in a cut, value of the cut)

        // for testing purposes:
        std::tuple<std::vector<int>, int64_t, int> SlowOneRespectedSparsestCut() const;
        std::pair<std::vector<int>, int64_t> SlowOneRespectedMincut() const;
        std::pair<std::vector<int>, int64_t> SlowOneRespectedBalancedCut(float ratio) const;

        std::pair<std::vector<int>, int64_t> SlowTwoRespectedMinCut() const;
        std::pair<std::vector<int>, int64_t> SlowTwoRespectedBalancedCut(float ratio) const;

        void CalculateTree(std::string tree_init_type);

        void PrintGraph() const;

        std::vector<std::pair<int, int>> GetEdgeList();

    private:
        std::mt19937 generator;

        void InitMST(std::vector<std::tuple<int, int, int> > weighted_edges);

        std::vector<std::tuple<int, int, int> > RandomlyWeightedEdgeList();
        // returns a list of edges (weight, vertex1, vertex2), weights are random ints

        std::vector<std::tuple<int, int, int>> NegativelyWeightedEdgeList();

        void InitRandomSpanningTree();

        int RandomNeighbor(int vertex);

        int64_t SlowCutSize(const std::vector<int>& cut) const;
        // returns the size of the cut, the cut is a list of vertices
};

#endif //GRAPH_H
