//
// Created by paul on 7/21/25.
//

#ifndef CLUSTERER_H
#define CLUSTERER_H
#include <unordered_set>
#include <vector>
#include "Graph.h"

class Clusterer {
    public:
        std::vector<std::vector<float> > points;
        int min_samples;
        // min_samples nearest neighbors are used to compute the core densities and the nearest neighbor graph

        std::string finding_neighbors_method;

        Clusterer(const std::vector<std::vector<float> > &_points, int _min_samples, const std::string &_finding_neighbors_method);

        std::vector<int> ClusterLabels(int num_clusters);

    private:
        bool DimensionConsistencyCheck() const; // checks that all the points have the same dimension

        void InitializeNeighborGraph(int max_integer_weight = 1'000'000'000);
        // TODO add brute force nearest neighbors also, in case of high dimension

        void MakeKCuts(int num_cuts);

        static float WeightFunction(float distance);   // the weight of an edge is 1/distance

        float Distance(int point_index_1, int point_index_2) const;

        int SparsestCutEdge() const;  // the node of the lower end of the sparsest cut edge

        void CutEdge(int edge); // makes a cut at a given edge, edge = its lower endpoint

        void UpdateSubtreeSizes(int cut_edge);

        void UpdateComponentSizes(int cut_edge);

        int RootOfComponent(int vertex) const;

        std::vector<int> ComponentLabels();
        // returns a vector of labels for the components of the tree separated by the edges in cut_tree_edges set

        std::vector<float> core_distances;
        Graph neighbor_graph;
        std::unordered_set<int> cut_tree_edges;
        std::vector<int> component_sizes;   // each vertex remembers the size of its components
};

#endif //CLUSTERER_H
