cmake_minimum_required(VERSION 3.31)
project(tree_based_cuts)

set(CMAKE_CXX_STANDARD 20)

add_executable(tree_based_cuts main.cpp Graph.cpp Graph.h Tree.cpp Tree.h DisjointSets.cpp DisjointSets.h
        RangeQuery.cpp
        RangeQuery.h)
