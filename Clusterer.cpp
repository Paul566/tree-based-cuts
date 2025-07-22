//
// Created by paul on 7/21/25.
//

#include "Clusterer.h"
#include "KDTree.hpp"

Clusterer::Clusterer(const std::vector<std::vector<float> > &_points,
                     int _min_samples,
                     const std::string &_finding_neighbors_method) : points(_points),
                                                                     finding_neighbors_method(
                                                                         _finding_neighbors_method),
                                                                     min_samples(_min_samples) {
    if (!DimensionConsistencyCheck()) {
        throw std::runtime_error("In Clusterer: Dimensions of all the points mush be the same");
    }
    InitializeNeighborGraph();
}

std::vector<int> Clusterer::ClusterLabels(const int num_clusters) {
    if (num_clusters != 2) {
        throw std::runtime_error("In Clusterer: not implemented num_clusters != 2");
    }

    std::vector<int> labels(points.size(), 0);
    auto [cut, sparsity] = neighbor_graph.OneRespectedSparsestCut();
    for (int index : cut) {
        labels[index] = 1;
    }
    return labels;
}

bool Clusterer::DimensionConsistencyCheck() const {
    if (points.empty()) {
        return true;
    }

    size_t dimension = points[0].size();
    for (int i = 1; i < static_cast<int>(points.size()); ++i) {
        if (points[i].size() != dimension) {
            return false;
        }
    }
    return true;
}

void Clusterer::InitializeNeighborGraph() {
    auto kdtree = KDTree(points);
    core_distances = std::vector<float>(points.size(), 0.);

    std::vector<std::vector<std::pair<int, float> > > adj_list = std::vector<std::vector<std::pair<int, float> > >(
        points.size(),
        std::vector<std::pair<int, float> >());

    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        auto neighbors = kdtree.nearest_indices(points[i], min_samples);
        float max_distance = 0.;
        for (int neighbor : neighbors) {
            if (neighbor == i) {
                continue;
            }
            float this_distance = Distance(i, neighbor);
            adj_list[i].emplace_back(neighbor, this_distance);
            // so far use distances as weights, correct wrt mutual reachability distance in the next loop

            if (this_distance > max_distance) {
                max_distance = this_distance;
            }
        }
        core_distances[i] = max_distance;
    }

    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        auto neighbors = kdtree.nearest_indices(points[i], min_samples);
        for (int j = 0; j < static_cast<int>(adj_list[i].size()); ++j) {
            adj_list[i][j].second = WeightFunction(std::max(adj_list[i][j].second, std::max(
                                                            core_distances[i],
                                                            core_distances[adj_list[i][j].first])));
        }
    }

    neighbor_graph = Graph(adj_list, "maximum_spanning_tree");
}

float Clusterer::WeightFunction(const float distance) {
    return 1. / distance;
}

float Clusterer::Distance(const int point_index_1, const int point_index_2) const {
    float norm_squared = 0.;
    for (int i = 0; i < static_cast<int>(points[point_index_1].size()); ++i) {
        norm_squared += (points[point_index_1][i] - points[point_index_2][i]) * (points[point_index_1][i] - points[
            point_index_2][i]);
    }
    return std::sqrt(norm_squared);
}
