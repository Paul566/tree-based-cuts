#ifndef TREE_H
#define TREE_H
#include <vector>

class Tree {
    public:
        std::vector<std::vector<int>> adj_list;

        explicit Tree(const std::vector<std::vector<int>> &_adj_list);

        explicit Tree(std::vector<std::pair<int, int>> _edge_list);

        Tree();
};



#endif //TREE_H
