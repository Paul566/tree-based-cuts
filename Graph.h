#ifndef GRAPH_H
#define GRAPH_H

#include <vector>


class Graph {
    public:
        std::vector<std::vector<int>> adj_list;

        explicit Graph(const std::vector<std::vector<int>> &_adj_list);

        void PrintGraph() const;


};



#endif //GRAPH_H
