#ifndef TREE_H
#define TREE_H
#include <vector>

#include "RangeQuery.h"

typedef std::tuple<int, int, int> Edge;   // (vertex1, vertex2, weight)

class Tree {
    public:
        explicit Tree(int num_vertices, const std::vector<std::pair<int, int>> &_edge_list);

        Tree();

        void UpdateDeltaCut(int node1, int node2, int weight);

        // updates the delta_cut vector, when handling an edge of the graph (node1, node2)

        void UpdateCutValues();
        // computes the sizes of 1-respected cuts based on delta_cut

        std::pair<std::vector<int>, int64_t> OneRespectedMincut() const;

        std::tuple<std::vector<int>, int64_t, int> OneRespectedSparsestCut() const;
        // returns (cut, size of the cut, size of the smallest part)

        std::pair<std::vector<int>, int64_t> OneRespectedBalancedCut(float ratio) const;
        // minimum ratio-balanced cut, if no such cut, returns ({}, inf)

        std::vector<int> SubtreeNodes(int vertex) const;
        std::pair<std::vector<int>,std::vector<int>> SubtreesNodes(int vertex1, int vertex2) const;

        int GetRoot() const;

        std::vector<int> GetDepths() const;

        std::vector<std::pair<int, int>> GetEdgeList() const;

    private:
        int root;
        std::vector<std::vector<int>> adj_list;
        std::vector<std::vector<int>> children;
        std::vector<int> parent;
        std::vector<int> subtree_sizes; // subtree size of a leaf vertex is 1
        std::vector<bool> heavy_parent; // for a vertex v, heavy_parent[v] is true if (v, parent[v]) is a heavy edge
        std::vector<int> ordered_edges;
        // edges_ordering[i] is e_i'th lowest end vertex (so e_i is an edge (edges_ordering[i], parent[edges_ordering[i]]))
        std::vector<int> parenting_edge_index; // index of the edge (i, parent[i]) in the ordering
        std::vector<int> postorder_nodes;
        std::vector<int> depth; // depth of the root is 0
        std::vector<int> path_indices;
        // path_indices[i] is the index of a path that the edge (i, parent[i]) beint64_ts to, -1 if i is root
        std::vector<int> path_upper_end; // path_upper_end[path_index] is the vertex v such that (v, parent[v]) is the uppermost edge of the path
        std::vector<int> path_lower_end; // path_lower_end[path_index] is the vertex v such that (v, parent[v]) is the lowest edge of the path
        std::vector<int64_t> delta_cut; // delta_cut[i] = cut_size(e_{i+1}) - cut_size(e_i)
        int64_t cut_size_e0; // the size of the cut if we cut e[0]
        std::vector<int64_t> cut_values; // cut_values[i] is the size of the 1-respected cut defined by e_i
        std::vector<int64_t> weights; // weights[i] - weight of e[i]

        void InitializeTreeStructure();

        void InitializePostorderNodes();

        void InitializeSubtreeSizes();

        void InitializeHeavyParents();

        void InitializePathData();

        int LCA(int node1, int node2) const;
        bool IsAncestor(int node1, int node2) const;

        std::vector<int> OutSubtreeNodes(int vertex) const;

        void UpdateDeltaCutHalfPath(int node, int lca, int weight);



        friend class BalancedCutFinder;
        friend class Clusterer;
};



#endif //TREE_H
