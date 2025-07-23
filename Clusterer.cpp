//
// Created by paul on 7/21/25.
//

#include "Clusterer.h"

#include <iostream>
#include <queue>

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

void Clusterer::InitializeNeighborGraph(int max_integer_weight) {
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
            bool i_in_neighbor_list = false;
            // check if adj_list[neighbor] already contains i
            // TODO this leads to O(m * maxdeg^2) time, maybe change that
            for (auto [node, dist] : adj_list[neighbor]) {
                if (node == i) {
                    i_in_neighbor_list = true;
                    break;
                }
            }
            if (!i_in_neighbor_list) {
                adj_list[neighbor].emplace_back(i, this_distance);
            }

            if (this_distance > max_distance) {
                max_distance = this_distance;
            }
        }
        core_distances[i] = max_distance;
    }

    float max_weight = 0.;
    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        for (int j = 0; j < static_cast<int>(adj_list[i].size()); ++j) {
            adj_list[i][j].second = WeightFunction(std::max(adj_list[i][j].second, std::max(
                                                                core_distances[i],
                                                                core_distances[adj_list[i][j].first])));
            if (adj_list[i][j].second > max_weight) {
                max_weight = adj_list[i][j].second;
            }
        }
    }

    std::vector<std::vector<std::pair<int, int>>> adj_list_for_graph(adj_list.size(), std::vector<std::pair<int, int> >());
    for (int i = 0; i < static_cast<int>(adj_list.size()); ++i) {
        adj_list_for_graph[i].reserve(adj_list[i].size());
        for (auto [neighbor, weight] : adj_list[i]) {
            if (neighbor > i) {
                int this_weight = static_cast<int>(weight * max_integer_weight / max_weight);
                adj_list_for_graph[i].emplace_back(neighbor, this_weight);
                adj_list_for_graph[neighbor].emplace_back(i, this_weight);
            }
        }
    }
    neighbor_graph = Graph(adj_list_for_graph, "maximum_spanning_tree");
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
    long best_cut_size = INT64_MAX;
    int best_denominator = 1;
    int best_edge = 0;
    for (int i = 0; i < static_cast<int>(neighbor_graph.tree.cut_values.size()); ++i) {
        long this_cut_size = neighbor_graph.tree.cut_values[i];
        int this_edge = neighbor_graph.tree.ordered_edges[i];
        int this_denominator = std::min(neighbor_graph.tree.subtree_sizes[this_edge],
                                        component_sizes[this_edge] - neighbor_graph.tree.subtree_sizes[this_edge]);

        if ((static_cast<double>(this_cut_size) * static_cast<double>(best_denominator) <
            static_cast<double>(best_cut_size) * static_cast<double>(this_denominator)) &&
            (!cut_tree_edges.contains(neighbor_graph.tree.ordered_edges[i]))) {
            best_cut_size = this_cut_size;
            best_denominator = this_denominator;
            best_edge = this_edge;
            }
    }

    return best_edge;
}

void Clusterer::CutEdge(const int edge) {
    cut_tree_edges.insert(edge);
    UpdateSubtreeSizes(edge);
    UpdateComponentSizes(edge);
}

void Clusterer::UpdateSubtreeSizes(int cut_edge) {
    int local_root = RootOfComponent(cut_edge);
    int subtree_size = neighbor_graph.tree.subtree_sizes[cut_edge];

    int cur_vertex = neighbor_graph.tree.parent[cut_edge];
    while (cur_vertex != local_root) {
        neighbor_graph.tree.subtree_sizes[cur_vertex] -= subtree_size;
        cur_vertex = neighbor_graph.tree.parent[cur_vertex];
    }
    neighbor_graph.tree.subtree_sizes[local_root] -= subtree_size;
}

void Clusterer::UpdateComponentSizes(int cut_edge) {
    int subtree_size = neighbor_graph.tree.subtree_sizes[cut_edge];

    std::queue<int> queue;
    queue.push(cut_edge);
    while (!queue.empty()) {
        int cur_vertex = queue.front();
        queue.pop();
        component_sizes[cur_vertex] = subtree_size;
        for (int child : neighbor_graph.tree.children[cur_vertex]) {
            if (!cut_tree_edges.contains(child)) {
                queue.push(child);
            }
        }
    }

    queue.push(RootOfComponent(neighbor_graph.tree.parent[cut_edge]));
    while (!queue.empty()) {
        int cur_vertex = queue.front();
        queue.pop();
        component_sizes[cur_vertex] -= subtree_size;
        for (int child : neighbor_graph.tree.children[cur_vertex]) {
            if (!cut_tree_edges.contains(child)) {
                queue.push(child);
            }
        }
    }
}

int Clusterer::RootOfComponent(const int vertex) const {
    int cur_vertex = vertex;
    while (!cut_tree_edges.contains(cur_vertex)) {
        if (cur_vertex == neighbor_graph.tree.root) {
            return neighbor_graph.tree.root;
        }
        cur_vertex = neighbor_graph.tree.parent[cur_vertex];
    }
    return cur_vertex;
}
