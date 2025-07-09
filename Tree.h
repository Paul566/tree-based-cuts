#ifndef TREE_H
#define TREE_H
#include <vector>

class Tree {
    public:
        int root;
        std::vector<std::vector<int>> adj_list;
        std::vector<std::vector<int>> children;
        std::vector<int> parent;

        explicit Tree(std::vector<std::pair<int, int>> _edge_list);

        Tree();

    private:
        std::vector<int> subtree_sizes; // subtree size of a leaf vertex is 1
        std::vector<bool> heavy_parent; // for a vertex v, heavy_parent[v] is true if (v, parent[v]) is a heavy edge

        void InitializeTreeStructure();

        void InitializeSubtreeSizes();

        void InitializeHeavyParents();
};



#endif //TREE_H
