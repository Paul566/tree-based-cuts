#include "Graph.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include "DisjointSets.h"

Graph::Graph(const std::vector<std::vector<int>> &_adj_list, const std::string &tree_init_type, int seed) {
    adj_list = _adj_list;
    // generator = std::mt19937(seed);
    std::random_device rd;
    generator = std::mt19937(rd());

    if (tree_init_type == "random_mst") {
        InitRandomMST();
    } else {
        if (tree_init_type == "random_spanning_tree") {
            InitRandomSpanningTree();
        } else {
            throw std::invalid_argument("Unsupported tree initialization type");
        }
    }
}

std::vector<int> Graph::SparsestOneRespectedCut() {
    // TODO
    return std::vector<int>();
}

void Graph::PrintGraph() const {
    std::cout << "Adjacency list:" << std::endl;
    for (int i = 0; i < adj_list.size(); ++i) {
        std::cout << i << ":";
        for (int j : adj_list[i]) {
            std::cout << " " << j;
        }
        std::cout << std::endl;
    }

    std::cout << "Tree adjacency list:" << std::endl;
    for (int i = 0; i < adj_list.size(); ++i) {
        std::cout << i << ":";
        for (int j : tree.adj_list[i]) {
            std::cout << " " << j;
        }
        std::cout << std::endl;
    }
}

void Graph::InitRandomMST() {
    auto weighted_edges = RandomlyWeightedEdgeList();
    std::sort(weighted_edges.begin(), weighted_edges.end());

    DisjointSets disjoint_sets(static_cast<int>(adj_list.size()));
    std::vector<std::pair<int, int>> mst_edges;
    mst_edges.reserve(static_cast<int>(weighted_edges.size()) - 1);

    for (const auto edge : weighted_edges) {
        if (!disjoint_sets.InSameSet(std::get<1>(edge), std::get<2>(edge))) {
            disjoint_sets.Unite(std::get<1>(edge), std::get<2>(edge));
            mst_edges.emplace_back(std::get<1>(edge), std::get<2>(edge));
        }
    }

    tree = Tree(mst_edges);
}

std::vector<std::tuple<int, int, int>> Graph::RandomlyWeightedEdgeList() {
    // returns a list of edges (weight, vertex1, vertex2)
    // weights are random ints

    std::vector<std::tuple<int, int, int>> weighted_edges;

    for (int i = 0; i < adj_list.size(); ++i) {
        for (int j = 0; j < adj_list[i].size(); ++j) {
            if (i < adj_list[i][j]) {
                weighted_edges.emplace_back(generator(), i, adj_list[i][j]);
            }
        }
    }

    return weighted_edges;
}

void Graph::InitRandomSpanningTree() {
    // TODO
}
