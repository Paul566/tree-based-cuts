#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <sstream>
#include <unordered_set>
#include <vector>

#include "BalancedCutFinder.h"
#include "Clusterer.h"
#include "Graph.h"
#include "Tree.h"
#include "KDTree.hpp"

std::vector<std::vector<int> > ReadGraph(const std::string &path) {
    std::vector<std::vector<int> > adj_list;
    std::fstream input_file(path);

    if (input_file.is_open()) {
        std::string line;
        std::getline(input_file, line);

        int num_vertices, num_edges;
        std::istringstream stream(line);
        stream >> num_vertices >> num_edges;

        for (int i = 0; i < num_vertices; ++i) {
            std::getline(input_file, line);
            std::vector<int> neighbors;

            stream = std::istringstream(line);
            int next_vertex;
            while (stream >> next_vertex) {
                neighbors.push_back(next_vertex - 1);
            }

            adj_list.push_back(neighbors);
        }
    } else {
        throw std::runtime_error("input file " + path + " not found");
    }
    input_file.close();

    return adj_list;
}

std::vector<std::vector<float>> ReadPointsFromCSV(const std::string& filename) {
    std::vector<std::vector<float>> points;
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("input file " + filename + " not found");
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<float> row;
        std::stringstream line_stream(line);
        std::string cell;

        while (std::getline(line_stream, cell, ',')) {
            row.push_back(std::stof(cell));
        }

        if (!row.empty()) {
            points.push_back(row);
        }
    }

    file.close();

    return points;
}

template<typename T>
void ExportVector(const std::string &path, const std::vector<T> &vector) {
    std::ofstream output(path);
    std::cout << std::fixed;
    std::cout << std::setprecision(4);

    for (T element : vector) {
        output << std::setprecision(4) << element << "\n";
    }

    output.close();
}

bool IsConnected(std::vector<std::vector<int> > adj_list) {
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

    if (static_cast<int64_t>(num_edges) > static_cast<int64_t>(num_vertices) * static_cast<int64_t>(num_vertices - 1) / 2) {
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

std::vector<std::vector<std::pair<int, int>>> AddRandomWeights(std::vector<std::vector<int>> adj_list, int min_weight, int max_weight, std::mt19937 &generator) {
    std::vector<std::vector<std::pair<int, int>>> weighted_adj_list(adj_list.size(), std::vector<std::pair<int, int>>());

    std::uniform_int_distribution<> dist_weights(min_weight, max_weight);
    for (int i = 0; i < static_cast<int>(adj_list.size()); ++i) {
        weighted_adj_list[i].reserve(adj_list[i].size());
        for (int neighbor : adj_list[i]) {
            if (neighbor > i) {
                float weight = dist_weights(generator);
                weighted_adj_list[i].emplace_back(neighbor, weight);
                weighted_adj_list[neighbor].emplace_back(i, weight);
            }
        }
    }

    return weighted_adj_list;
}

std::vector<std::vector<int> > SquareGridGraph(int size) {
    // returns an adjacency list for a size*size square grid
    std::vector<std::vector<int> > adj_list(size * size, std::vector<int>());
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i > 0) {
                adj_list[i * size + j].push_back((i - 1) * size + j);
            }
            if (j > 0) {
                adj_list[i * size + j].push_back(i * size + j - 1);
            }
            if (i < size - 1) {
                adj_list[i * size + j].push_back((i + 1) * size + j);
            }
            if (j < size - 1) {
                adj_list[i * size + j].push_back(i * size + j + 1);
            }
        }
    }

    return adj_list;
}

void RunRandomTests(int max_vertices,
                              int num_tests,
                              std::mt19937 &generator,
                              std::string tree_init_type,
                              int min_weight,
                              int max_weight) {
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

        auto unweighted_adj_list = RandomConnectedGraph(num_vertices, num_edges, generator);
        auto adj_list = AddRandomWeights(unweighted_adj_list, min_weight, max_weight, generator);
        Graph graph(adj_list, tree_init_type, 239);

        auto [_sc, sparsest_cut_size, denominator] = graph.OneRespectedSparsestCut();
        auto [_mc, min_cut_value] = graph.OneRespectedMincut();
        auto [_bc, balanced_cut_value] = graph.OneRespectedBalancedCut(sparsity);
        auto [_ssc, slow_sparsest_cut_size, slow_denominator] = graph.SlowOneRespectedSparsestCut();
        auto [_smc, slow_min_cut_value] = graph.SlowOneRespectedMincut();
        auto [_sbc, slow_balanced_cut_value] = graph.SlowOneRespectedBalancedCut(sparsity);

        auto [_strc, slow_two_respected_cut] = graph.SlowTwoRespectedBalancedCut(sparsity);
        BalancedCutFinder b(graph, graph.tree, sparsity);
        auto [_trc, two_respected_cut] = b.TwoRespectedBalancedCut();

        if (sparsest_cut_size * slow_denominator != slow_sparsest_cut_size * denominator) {
            graph.PrintGraph();
            throw std::runtime_error("sparsest cut test failed");
        }
        if (min_cut_value != slow_min_cut_value) {
            graph.PrintGraph();
            std::cout << min_cut_value << " vs " << slow_min_cut_value << std::endl;
            throw std::runtime_error("min cut test failed");
        }
        if (balanced_cut_value != slow_balanced_cut_value) {
            graph.PrintGraph();
            std::cout << balanced_cut_value << " vs " << slow_balanced_cut_value << std::endl;
            throw std::runtime_error("balanced cut test failed");
        }
        if (two_respected_cut != slow_two_respected_cut) {
            graph.PrintGraph();
            std::cout << two_respected_cut << " vs " << slow_two_respected_cut << std::endl;
            throw std::runtime_error("two balanced cut test failed");
        }
    }

    std::cout << "All tests passed!\n";
}

int BalancedCutSize(float ratio,
                    int num_tries,
                    std::vector<std::vector<int> > adj_list,
                    std::string tree_init_type = "random_mst",
                    int seed_for_tree = 239) {
    Graph graph(adj_list, tree_init_type, seed_for_tree);
    auto [best_cut, best_value] = graph.OneRespectedBalancedCut(ratio);

    for (int i = 0; i < num_tries - 1; ++i) {
        graph.CalculateTree(tree_init_type);
        auto [cut, value] = graph.OneRespectedBalancedCut(ratio);
        if (value < best_value) {
            best_value = value;
        }
    }

    return best_value;
}

std::vector<std::string> AllFiles(const std::string &directory_path) {
    std::vector<std::string> files;

    try {
        for (const auto &entry : std::filesystem::directory_iterator(directory_path)) {
            if (std::filesystem::is_regular_file(entry.status())) {
                files.push_back(entry.path().string());
            }
        }
    } catch (const std::filesystem::filesystem_error &e) {
        std::cerr << "Cannot access directory: " << e.what() << std::endl;
    }

    std::sort(files.begin(), files.end());

    return files;
}

void RunBalancedCutBenchmark(float ratio,
                             int num_tries_per_mln_vertices,
                             std::string tree_init_type = "random_mst",
                             int seed_for_tree = 239) {
    std::string prefix = std::filesystem::current_path().string() + "/../graphs/";

    auto files = AllFiles(prefix);
    for (const auto &file : files) {
        size_t pos = file.find_last_of('/');
        if ((file.substr(pos + 1) == "bcsstk29.graph") || (file.substr(pos + 1) == "bcsstk31.graph") || (file.
            substr(pos + 1) == "fe_body.graph") || (file.substr(pos + 1) == "fe_pwt.graph")) {
            // these graphs are disconnected
            continue;
        }

        auto adj_list = ReadGraph(file);
        int num_tries = num_tries_per_mln_vertices * 1000000 / adj_list.size();
        std::cout << file.substr(pos + 1) << "\tnum_tries: " << num_tries << std::endl;
        std::cout << BalancedCutSize(ratio, num_tries, adj_list, tree_init_type, seed_for_tree) << '\n' << std::endl;
    }
}

void ExportGridCut(int size,
                   std::string filename_cut,
                   std::string filename_depths,
                   int seed,
                   std::string init_tree_type) {
    auto adj_list = SquareGridGraph(size);
    Graph graph(adj_list, init_tree_type, seed);
    auto [cut, cut_size, denominator] = graph.OneRespectedSparsestCut();
    std::unordered_set<int> cut_set;
    for (int element : cut) {
        cut_set.insert(element);
    }
    std::vector<int> cut_mask(size * size, 0);
    for (int i = 0; i < size * size; ++i) {
        if (cut_set.contains(i)) {
            cut_mask[i] = 1;
        }
    }

    ExportVector(filename_cut, cut_mask);
    ExportVector(filename_depths, graph.tree.GetDepths());
}

int main() {
    std::random_device rd;
    // std::mt19937 generator(rd());
    std::mt19937 generator(239);

    auto points = ReadPointsFromCSV("../clustering/points.csv");

    Clusterer clusterer(points, 5, "kdtree");
    auto labels = clusterer.ClusterLabels(2);
    ExportVector("../clustering/labels", labels);
    for (auto label : labels) {
        std::cout << label << std::endl;
    }

    // RunRandomTests(100, 100000, generator, "random_mst", 0, 10);

    return 0;
}
