#include <algorithm>
#include <iostream>
#include <queue>
#include <set>
#include <vector>
#include "Graph.h"
#include "RangeQuery.h"
#include "Tree.h"

bool IsConnected(std::vector<std::vector<int>> adj_list) {
    if (adj_list.empty()) {
        return true;
    }

    std::vector<bool> visited(adj_list.size(), false);
    std::queue<int> queue;

    queue.push(0);
    visited[0] = true;
    while (!queue.empty()) {
        int node = queue.front();
        queue.pop();

        for (int neighbor : adj_list[node]) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                queue.push(neighbor);
            }
        }
    }

    for (bool value : visited) {
        if (!value) {
            return false;
        }
    }
    return true;
}

std::vector<std::vector<int> > RandomGraph(const int num_vertices,
                                                     const int num_edges,
                                                     std::mt19937 &generator) {
    std::vector<std::vector<int> > adj_list(num_vertices, std::vector<int>());

    if (static_cast<long>(num_edges) > static_cast<long>(num_vertices) * static_cast<long>(num_vertices - 1) / 2) {
        throw std::runtime_error("In RandomGraph: too many edges");
    }

    std::set<std::pair<int, int> > edges;

    std::uniform_int_distribution<> dist(0, num_vertices - 1);
    while (edges.size() < static_cast<size_t>(num_edges)) {
        int first_vertex = dist(generator); // Random vertex u
        int second_vertex = dist(generator);

        if (first_vertex != second_vertex) {
            auto edge = std::minmax(first_vertex, second_vertex);
            if (!edges.contains(edge)) {
                edges.insert(edge);
                adj_list[edge.first].push_back(edge.second);
                adj_list[edge.second].push_back(edge.first);
            }
        }
    }

    return adj_list;
}

std::vector<std::vector<int> > RandomConnectedGraph(const int num_vertices,
                                                     const int num_edges,
                                                     std::mt19937 &generator) {
    // takes a random graph with num_vertices and num_edges, and adds random edges until the graph becomes connected

    std::vector<std::vector<int> > adj_list(num_vertices, std::vector<int>());
    std::set<std::pair<int, int> > edges;

    std::uniform_int_distribution<> dist(0, num_vertices - 1);
    while ((edges.size() < static_cast<size_t>(num_edges)) || (!IsConnected(adj_list))) {
        int first_vertex = dist(generator); // Random vertex u
        int second_vertex = dist(generator);

        if (first_vertex != second_vertex) {
            auto edge = std::minmax(first_vertex, second_vertex);
            if (!edges.contains(edge)) {
                edges.insert(edge);
                adj_list[edge.first].push_back(edge.second);
                adj_list[edge.second].push_back(edge.first);
            }
        }
    }

    return adj_list;
}

void RunRandomUnweightedTests(int max_vertices,
                              int num_tests,
                              std::mt19937 &generator) {
    std::uniform_int_distribution<> dist_vertices(2, max_vertices);

    for (int i = 0; i < num_tests; ++i) {
        if (i % 1 == 0) {
            std::cout << "test " << i << "\n";
        }

        int num_vertices = dist_vertices(generator);
        std::uniform_int_distribution<> dist_edges(num_vertices - 1, num_vertices * (num_vertices - 1) / 2);
        int num_edges = dist_edges(generator);
        std::uniform_real_distribution<> dist_sparsity(0, 0.5);
        float sparsity = dist_sparsity(generator);

        auto adj_list = RandomConnectedGraph(num_vertices, num_edges, generator);
        Graph graph(adj_list, "random_mst",239);

        auto [_sc, sparsest_cut_value] = graph.OneRespectedSparsestCut();
        auto [_mc, min_cut_value] = graph.OneRespectedMincut();
        auto [_bc, balanced_cut_value] = graph.OneRespectedBalancedCut(sparsity);
        auto [_ssc, slow_sparsest_cut_value] = graph.SlowOneRespectedSparsestCut();
        auto [_smc, slow_min_cut_value] = graph.SlowOneRespectedMincut();
        auto [_sbc, slow_balanced_cut_value] = graph.SlowOneRespectedBalancedCut(sparsity);

        if (sparsest_cut_value != slow_sparsest_cut_value) {
            graph.PrintGraph();
            throw std::runtime_error("sparsest cut test failed");
        }
        if (min_cut_value != slow_min_cut_value) {
            graph.PrintGraph();
            throw std::runtime_error("min cut test failed");
        }
        if (balanced_cut_value != slow_balanced_cut_value) {
            graph.PrintGraph();
            throw std::runtime_error("balanced cut test failed");
        }
    }

    std::cout << "All tests passed!\n";
}


int main() {
    std::random_device rd;
    // std::mt19937 generator(rd());
    std::mt19937 generator(239);

    RunRandomUnweightedTests(100, 100000, generator);

    /*std::vector<std::vector<int>> adj_list;
    adj_list.emplace_back(std::vector<int>{4});
    adj_list.emplace_back(std::vector<int>{4, 2});
    adj_list.emplace_back(std::vector<int>{4, 1, 3});
    adj_list.emplace_back(std::vector<int>{4, 2});
    adj_list.emplace_back(std::vector<int>{2, 0, 1, 3});

    Graph graph(adj_list, "random_mst",239);
    auto [_sc, sparsest_cut_value] = graph.OneRespectedSparsestCut();
    std::cout << sparsest_cut_value << std::endl;*/

    return 0;
}