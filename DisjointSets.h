#ifndef K_FORESTS_DISJOINTSETS_H
#define K_FORESTS_DISJOINTSETS_H


#include <vector>

class DisjointSets {
    public:
        explicit DisjointSets(int size);

        bool InSameSet(int first, int second);

        void Unite(int first, int second);

    private:
        std::vector<int> ranks;
        std::vector<int> parents;

        int Representative(int element);

        void UniteRepresentatives(int first, int second);
};


#endif //K_FORESTS_DISJOINTSETS_H
