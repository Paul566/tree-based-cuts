//
// Created by paul on 7/21/25.
//

#ifndef CLUSTERER_H
#define CLUSTERER_H
#include <vector>
#include "Graph.h"

class Clusterer {
    public:
        std::vector<std::vector<float> > points;
        int min_samples;
        // min_samples nearest neighbors are used to compute the core densities and the nearest neighbor graph

        std::string finding_neighbors_method;
        std::vector<float> core_distances;
        Graph neighbor_graph;

        Clusterer(const std::vector<std::vector<float> > &_points, int _min_samples, const std::string &_finding_neighbors_method);

        std::vector<int> ClusterLabels(int num_clusters);

    private:
        bool DimensionConsistencyCheck() const; // checks that all the points have the same dimension

        void InitializeNeighborGraph();
        // TODO add brute force nearest neighbors also, in case of high dimension

        static float WeightFunction(float distance);   // the weight of an edge is 1/distance

        float Distance(int point_index_1, int point_index_2) const;
};

#endif //CLUSTERER_H
