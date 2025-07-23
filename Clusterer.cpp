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

    component_sizes = std::vector<int>(points.size(), static_cast<int>(points.size()));
}

std::vector<int> Clusterer::ClusterLabels(const int num_clusters) {
    if (num_clusters != 2) {
        throw std::runtime_error("In Clusterer: not implemented num_clusters != 2");
    }

    std::vector<int> labels(points.size(), 0);
    auto [cut, cut_size, denominator] = neighbor_graph.OneRespectedSparsestCut();
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

    // TODO do the rescaling
    std::vector<std::vector<std::pair<int, int> > > adj_list = std::vector<std::vector<std::pair<int, int> > >(
        points.size(),
        std::vector<std::pair<int, int> >());

    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        auto neighbors = kdtree.nearest_indices(points[i], min_samples);
        float max_distance = 0.;
        for (int neighbor : neighbors) {
            float this_distance = Distance(i, neighbor);

            adj_list[i].emplace_back(neighbor, this_distance);
            // so far use distances as weights, correct wrt mutual reachability distance in the next loop
            bool i_in_neigbor_list = false;
            // check if adj_list[neighbor] already contains i
            // TODO this leads to O(m * maxdeg^2) time, maybe change that
            for (auto [node, dist] : adj_list[neighbor]) {
                if (node == i) {
                    i_in_neigbor_list = true;
                    break;
                }
            }
            if (!i_in_neigbor_list) {
                adj_list[neighbor].emplace_back(i, this_distance);
            }

            if (this_distance > max_distance) {
                max_distance = this_distance;
            }
        }
        core_distances[i] = max_distance;
    }

    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        auto neighbors = kdtree.nearest_indices(points[i], min_samples);
        for (int j = 0; j < static_cast<int>(adj_list[i].size()); ++j) {
                adj_list[i][j].second = WeightFunction(std::max(static_cast<float>(adj_list[i][j].second), std::max(
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

int Clusterer::SparsestCutEdge() const {
    float min_value = INT32_MAX;
    int best_edge = 0;
    for (int i = 0; i < static_cast<int>(neighbor_graph.tree.cut_values.size()); ++i) {
        const float this_sparsity = neighbor_graph.tree.cut_values[i] / static_cast<float>(
            std::min(neighbor_graph.tree.subtree_sizes[neighbor_graph.tree.ordered_edges[i]],
            component_sizes[neighbor_graph.tree.ordered_edges[i]] - neighbor_graph.tree.subtree_sizes[neighbor_graph.tree.ordered_edges[i]]));
        if ((this_sparsity < min_value) && (!cut_tree_edges.contains(neighbor_graph.tree.ordered_edges[i]))) {
            min_value = this_sparsity;
            best_edge = neighbor_graph.tree.ordered_edges[i];
        }
    }

    return best_edge;
}

// void Clusterer::CutEdge(const int edge) {
//     cut_tree_edges.insert(edge);
//     UpdateSubtreeSizes(edge);
//     UpdateComponentSizes(edge);
// }
//
// void Clusterer::UpdateSubtreeSizes(int cut_edge) {
//     int local_root = RootOfComponent(cut_edge);
//     int subtree_size = neighbor_graph.tree.subtree_sizes[cut_edge];
//     int cur_vertex = neighbor_graph.tree.parent[cut_edge];
//     while (cur_vertex != local_root) {
//         // TODO
//     }
// }
