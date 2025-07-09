#ifndef TREE_H
#define TREE_H
#include <vector>

class Tree {
    public:
        int root;
        std::vector<std::vector<int>> adj_list;
        std::vector<std::vector<int>> children;
        std::vector<int> parent;

        explicit Tree(const std::vector<std::vector<int>> &_adj_list);

        explicit Tree(std::vector<std::pair<int, int>> _edge_list);

        Tree();

    private:
        void InitializeTree();
};



#endif //TREE_H
