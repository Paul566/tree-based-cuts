#include "Graph.h"

#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "DisjointSets.h"

Graph::Graph(const std::vector<std::vector<int> > &_adj_list, const std::string &_tree_init_type, int seed_for_tree) {
    if (_adj_list.size() <= 1) {
        throw std::runtime_error("Graph::Graph(): need at least two vertices");
    }

    adj_list = _adj_list;
    generator = std::mt19937(seed_for_tree);
    // std::random_device rd;
    // generator = std::mt19937(rd());

    CalculateTree(_tree_init_type);
}

std::pair<std::vector<int>, float> Graph::OneRespectedSparsestCut() {
    return tree.OneRespectedSparsestCut();
}

std::pair<std::vector<int>, int> Graph::OneRespectedMincut() {
    return tree.OneRespectedMincut();
}

std::pair<std::vector<int>, int> Graph::OneRespectedBalancedCut(float ratio) {
    return tree.OneRespectedBalancedCut(ratio);
}

int Graph::TwoRespectedMinCut() {
    return tree.TwoRespectedMinCut();
}


std::pair<std::vector<int>, float> Graph::SlowOneRespectedSparsestCut() const {
    float min_sparsity = INT32_MAX;
    std::vector<int> best_cut;

    for (int vertex = 0; vertex < static_cast<int>(adj_list.size()); ++vertex) {
        if (vertex == tree.GetRoot()) {
            continue;
        }

        int cut_size = SlowCutSize(tree.SubtreeNodes(vertex));
        float part_size = static_cast<float>(tree.SubtreeNodes(vertex).size());
        float this_sparsity = static_cast<float>(cut_size) / static_cast<float>(
            std::min(part_size,static_cast<float>(adj_list.size()) - part_size));
        if (this_sparsity < min_sparsity) {
            min_sparsity = this_sparsity;
            best_cut = tree.SubtreeNodes(vertex);
        }
    }

    return {best_cut, min_sparsity};
}

std::pair<std::vector<int>, int> Graph::SlowOneRespectedMincut() const {
    int min_cut_size = INT32_MAX;
    std::vector<int> best_cut;

    for (int vertex = 0; vertex < static_cast<int>(adj_list.size()); ++vertex) {
        if (vertex == tree.GetRoot()) {
            continue;
        }

        int cut_size = SlowCutSize(tree.SubtreeNodes(vertex));
        if (cut_size < min_cut_size) {
            min_cut_size = cut_size;
            best_cut = tree.SubtreeNodes(vertex);
        }
    }

    return {best_cut, min_cut_size};
}

std::pair<std::vector<int>, int> Graph::SlowOneRespectedBalancedCut(float ratio) const {
    int min_cut_size = INT32_MAX;
    std::vector<int> best_cut;

    for (int vertex = 0; vertex < static_cast<int>(adj_list.size()); ++vertex) {
        if (vertex == tree.GetRoot()) {
            continue;
        }

        int cut_size = SlowCutSize(tree.SubtreeNodes(vertex));
        float part_size = static_cast<float>(tree.SubtreeNodes(vertex).size());

        if ((part_size < static_cast<float>(adj_list.size()) * ratio) ||
            (static_cast<float>(adj_list.size()) - part_size < static_cast<float>(adj_list.size()) * ratio)) {
            continue;
            }

        if (cut_size < min_cut_size) {
            min_cut_size = cut_size;
            best_cut = tree.SubtreeNodes(vertex);
        }
    }

    if (min_cut_size == INT32_MAX) {
        // no balanced cut
        return {std::vector<int>(), INT32_MAX};
    }
    return {best_cut, min_cut_size};
}

std::pair<std::vector<int>, int> Graph::SlowTwoRespectedMinCut() const{
    int min_cut_size = INT32_MAX;
    std::vector<int> best_cut;
    int graph_size = static_cast<int>(adj_list.size());

    for (int vertex1 = 0; vertex1 < graph_size; ++vertex1) {
        if (vertex1 == tree.GetRoot()) {
            continue;
        }

        for (int vertex2 = vertex1 + 1; vertex2 < graph_size; ++vertex2) {
            if (vertex2 == tree.GetRoot()) {
                continue;
            }

            auto sub_trees = tree.SubtreesNodes(vertex1, vertex2);
            auto cut = sub_trees.first;
            cut.insert(cut.end(), sub_trees.second.begin(), sub_trees.second.end());

            int cut_size = SlowCutSize(cut);

            if (cut_size < min_cut_size) {
                min_cut_size = cut_size;
                best_cut = cut;
            }
        }
    }

    if (min_cut_size == INT32_MAX) {
        // no balanced cut
        return {std::vector<int>(), INT32_MAX};
    }
    return {best_cut, min_cut_size};
}

std::pair<std::vector<int>, int> Graph::SlowTwoRespectedBalancedCut(float ratio) const{
    int min_cut_size = INT32_MAX;
    std::vector<int> best_cut;
    int graph_size = static_cast<int>(adj_list.size());

    for (int vertex1 = 0; vertex1 < graph_size; ++vertex1) {
        if (vertex1 == tree.GetRoot()) {
            continue;
        }

        for (int vertex2 = vertex1 + 1; vertex2 < graph_size; ++vertex2) {
            if (vertex2 == tree.GetRoot()) {
                continue;
            }

            auto sub_trees = tree.SubtreesNodes(vertex1, vertex2);
            auto cut = sub_trees.first;
            cut.insert(cut.end(), sub_trees.second.begin(), sub_trees.second.end());

            int cut_size = SlowCutSize(cut);
            float part_size = cut.size();

            if ((part_size < static_cast<float>(adj_list.size()) * ratio) ||
                (static_cast<float>(adj_list.size()) - part_size < static_cast<float>(adj_list.size()) * ratio)) {
                continue;
                }

            if (cut_size < min_cut_size) {
                min_cut_size = cut_size;
                best_cut = cut;
            }
        }
    }

    if (min_cut_size == INT32_MAX) {
        // no balanced cut
        return {std::vector<int>(), INT32_MAX};
    }
    return {best_cut, min_cut_size};
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
}

void Graph::CalculateTree(std::string tree_init_type) {
    if (tree_init_type == "random_mst") {
        InitRandomMST();
    } else {
        if (tree_init_type == "random_spanning_tree") {
            InitRandomSpanningTree();
        } else {
            throw std::invalid_argument("Unsupported tree initialization type");
        }
    }

    for (int i = 0; i < static_cast<int>(adj_list.size()); ++i) {
        for (int j : adj_list[i]) {
            if (i < j) {
                tree.UpdateDeltaCut(i, j);
                tree.HandleGraphEdge(std::make_pair(i, j));
            }
        }
    }
    tree.UpdateCutValues();
}

void Graph::InitRandomMST() {
    auto weighted_edges = RandomlyWeightedEdgeList();
    std::sort(weighted_edges.begin(), weighted_edges.end());

    DisjointSets disjoint_sets(static_cast<int>(adj_list.size()));
    std::vector<std::pair<int, int> > mst_edges;
    mst_edges.reserve(static_cast<int>(weighted_edges.size()));

    for (const auto edge : weighted_edges) {
        if (!disjoint_sets.InSameSet(std::get<1>(edge), std::get<2>(edge))) {
            disjoint_sets.Unite(std::get<1>(edge), std::get<2>(edge));
            mst_edges.emplace_back(std::get<1>(edge), std::get<2>(edge));
        }
    }

    for (int i = 1; i < static_cast<int>(adj_list.size()); ++i) {
        if (!disjoint_sets.InSameSet(0, i)) {
            throw std::runtime_error("The graph is disconnected");
        }
    }

    tree = Tree(mst_edges);
}

std::vector<std::tuple<int, int, int> > Graph::RandomlyWeightedEdgeList() {
    std::vector<std::tuple<int, int, int> > weighted_edges;

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
    std::unordered_set<int> in_tree;
    std::vector<std::pair<int, int>> tree_edges;
    std::uniform_int_distribution<int> vertex_dist(0, static_cast<int>(adj_list.size()) - 1);

    // Start from a random node
    int start = vertex_dist(generator);
    in_tree.insert(start);

    while (in_tree.size() < static_cast<int>(adj_list.size())) {
        int next_vertex;
        // Pick a random node not in the tree yet
        do {
            next_vertex = vertex_dist(generator);
        } while (in_tree.contains(next_vertex));

        std::unordered_map<int, int> parent;
        std::unordered_set<int> visited;
        int cur_vertex = next_vertex;

        // Loop-erased random walk
        while (!in_tree.contains(cur_vertex)) {
            visited.insert(cur_vertex);
            int next = RandomNeighbor(cur_vertex);
            parent[cur_vertex] = next;

            // Loop erasure: if a cycle is formed, overwrite earlier path
            if (visited.contains(next)) {
                while (visited.contains(cur_vertex)) {
                    visited.erase(cur_vertex);
                    cur_vertex = parent[cur_vertex];
                }
                cur_vertex = next;
            } else {
                cur_vertex = next;
            }
        }

        // Add the path from u to the tree
        cur_vertex = next_vertex;
        while (!in_tree.count(cur_vertex)) {
            int cur_parent = parent[cur_vertex];
            tree_edges.emplace_back(cur_vertex, cur_parent);
            in_tree.insert(cur_vertex);
            cur_vertex = cur_parent;
        }
    }

    tree = Tree(tree_edges);
}

int Graph::RandomNeighbor(const int vertex) {
    const auto& neighbors = adj_list[vertex];
    std::uniform_int_distribution<int> dist(0, neighbors.size() - 1);
    return neighbors[dist(generator)];
}

int Graph::SlowCutSize(const std::vector<int>& cut) const {
    std::unordered_set<int> cut_set;
    for (int vertex : cut) {
        cut_set.insert(vertex);
    }

    int cut_size = 0;
    for (int vertex : cut) {
        for (int neighbor : adj_list[vertex]) {
            if (!cut_set.contains(neighbor)) {
                ++cut_size;
            }
        }
    }

    return cut_size;
}
