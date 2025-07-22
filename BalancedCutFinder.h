//
// Created by alevran on 7/18/2025.
//

#ifndef BALANCEDCUTFINDER_H
#define BALANCEDCUTFINDER_H
#include <map>

#include "Graph.h"


class BalancedCutFinder {
    public:
    BalancedCutFinder(const Graph& graph, const Tree& tree, float factor);
    std::pair<std::vector<int>, float> TwoRespectedBalancedCut();

    private:
    const Graph& graph;
    const Tree& tree;
    float factor;
    int n;
    // dict from i to list of edges (u,v) s.t. e[i-1] isn't in u->v, but e[i] is (or e[i] is and i==0)
    std::vector<std::vector<Edge>> tree_edge_to_path_in_edges;
    // dict from i to list of edges (u,v) s.t. e[i] is in u->v, but e[i+1] isn't (or e[i] is and i==n-1)
    std::vector<std::vector<Edge>> tree_edge_to_path_out_edges;
    std::map<int, RangeQuery> weights;
    float global_weight;

    // path p-> how much gap we have between e_i and weights[p]
    std::map<int, int> paths_gaps;
    std::map<int, int> path_array_length;
    std::vector<int> edge_to_subtree_size;

    void InitializeWeights();

    void InitializeSegments();

    void HandleGraphEdge(const Edge &edge);

    std::vector<std::pair<int, int>> GetSegmentsFromHalfPath(int node, int lca) const;

    // add `weight` among all tree edges in the path connecting nodes of `edge`
    void AddPath(const Edge& edge, bool negative);

    void Update(int path, int l, int r, float weight);

    // add `weight` among all tree edges out the path connecting nodes of `edge`
    void AddOutPath(const Edge& edge, bool negative);

    // find min cut (within [n*factor, n-n*factor]) through `edge` and another tree edge
    std::pair<float, int> MinCutWithEdge(int edge);

    // find the min edge inside the path with index `path` in the segment [min_edge, max_edge].
    // together with another subtree of size `other_subtree`, the subtrees together need to have size
    // within [n*factor, n-n*factor]. If is_in_subtree, calculate the subtree rooted in it; otherwise, the outside subtree
    std::pair<float, int> MinOfPath(int path, int min_edge, int max_edge, int other_subtree, bool is_in_subtree, RangeQuery &query);

    // adjust min_tree, max_tree to be in the window size we want to check
    void GetMinMaxEdge(int &min_edge, int &max_edge, int min_subtree, int max_subtree, bool is_in_subtree);

};



#endif //BALANCEDCUTFINDER_H
