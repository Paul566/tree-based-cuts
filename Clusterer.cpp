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
    if (min_samples > static_cast<int>(_points.size())) {
        throw std::runtime_error("In Clusterer: min_samples is greater than the number of points");
    }

    InitializeNeighborGraph();

    component_sizes = std::vector<int>(points.size(), static_cast<int>(points.size()));
}

std::vector<int> Clusterer::ClusterLabels(const int num_clusters) {
    if (num_clusters < 1) {
        throw std::runtime_error("num_clusters must be at least 1");
    }

    if (!cut_tree_edges.empty()) {
        throw std::runtime_error("In ClusterLabels: the tree has already been partitioned");
    }

    MakeKCuts(num_clusters - 1);

    return ComponentLabels();
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
    auto neighbors_list = GetNearestNeighbors();
    core_distances = std::vector<float>(points.size(), 0.);

    std::vector<std::vector<std::pair<int, float> > > adj_list = std::vector<std::vector<std::pair<int, float> > >(
        points.size(),
        std::vector<std::pair<int, float> >());

    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        float max_distance = 0.;
        for (int neighbor : neighbors_list[i]) {
            if (neighbor == i) {
                continue;
            }

            float this_distance = Distance(i, neighbor);

            adj_list[i].emplace_back(neighbor, this_distance);
            // so far use distances as weights, correct wrt mutual reachability distance in the next loop
            bool i_in_neighbor_list = false;
            // check if adj_list[neighbor] already contains i
            // TODO this is O(m * maxdeg^2) time, change this
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

std::vector<std::vector<int>> Clusterer::GetNearestNeighbors() {
    if (finding_neighbors_method == "kdtree") {
        return GetNearestNeighborsKDT();
    }
    if (finding_neighbors_method == "brute_force") {
        return GetNearestNeighborsBrute();
    }
    throw std::runtime_error("Unsupported finding_neighbors_method: " + finding_neighbors_method);
}

std::vector<std::vector<int>> Clusterer::GetNearestNeighborsKDT() {
    auto kdtree = KDTree(points);
    std::vector<std::vector<int>> adj_list(points.size(), std::vector<int>());
    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        auto neighbors = kdtree.nearest_indices(points[i], min_samples);
        adj_list[i].reserve(min_samples);
        for (int neighbor : neighbors) {
            adj_list[i].push_back(neighbor);
        }
    }

    return adj_list;
}

std::vector<std::vector<int>> Clusterer::GetNearestNeighborsBrute() {
    std::vector<std::vector<int>> adj_list(points.size(), std::vector<int>());
    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        adj_list[i] = KNearestNeighborsBrute(i, min_samples);
    }
    return adj_list;
}

std::vector<int> Clusterer::KNearestNeighborsBrute(int point_index, int num_neighbors) {
    // TODO this is O(n log(n)), can make it O(n) with quick select
    std::vector<std::pair<float, int>> dist_to_points;
    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        dist_to_points.emplace_back(Distance(point_index, i), i);
    }
    std::sort(dist_to_points.begin(), dist_to_points.end());

    std::vector<int> neighbors;
    neighbors.reserve(num_neighbors);
    for (int i = 0; i < num_neighbors; ++i) {
        neighbors.push_back(dist_to_points[i].second);
    }

    return neighbors;
}

void Clusterer::MakeKCuts(const int num_cuts) {
    for (int i = 0; i < num_cuts; ++i) {
        const int next_edge = SparsestCutEdge();
        CutEdge(next_edge);
    }
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

        if (((static_cast<double>(this_cut_size) * static_cast<double>(best_denominator) <
            static_cast<double>(best_cut_size) * static_cast<double>(this_denominator))) ||
            ((this_cut_size * best_denominator == best_cut_size * this_denominator) && (this_denominator > best_denominator))) {
            if (!cut_tree_edges.contains(neighbor_graph.tree.ordered_edges[i])) {
                // handles the case of disconnected graphs, we need the most balanced 0-cut then
                best_cut_size = this_cut_size;
                best_denominator = this_denominator;
                best_edge = this_edge;
            }
            }
    }

    return best_edge;
}

void Clusterer::CutEdge(const int edge) {
    // TODO subtract the cut edges
    cut_tree_edges.insert(edge);
    UpdateSubtreeSizes(edge);
    UpdateComponentSizes(edge);
}

void Clusterer::UpdateSubtreeSizes(int cut_edge) {
    int subtree_size = neighbor_graph.tree.subtree_sizes[cut_edge];

    int cur_vertex = neighbor_graph.tree.parent[cut_edge];
    int local_root = RootOfComponent(cur_vertex);
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

std::vector<int> Clusterer::ComponentLabels() {
    std::vector<int> labels(points.size(), 0);
    int component_counter = 0;

    std::queue<int> queue;
    queue.push(neighbor_graph.tree.root);
    while (!queue.empty()) {
        int cur_vertex = queue.front();
        queue.pop();

        for (int child : neighbor_graph.tree.children[cur_vertex]) {
            queue.push(child);

            if (!cut_tree_edges.contains(child)) {
                labels[child] = labels[cur_vertex];
            } else {
                ++component_counter;
                labels[child] = component_counter;
            }
        }
    }

    return labels;
}
