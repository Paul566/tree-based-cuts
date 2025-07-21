// This file is taken from the repository https://github.com/crvs/KDTree/tree/master
// The repository contains the following disclaimer:

// BSD 3-Clause License
//
// Copyright (c) 2018, J. Frederico Carvalho
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#pragma once

/// @file KDTree.hpp
/// @author J. Frederico Carvalho
///
/// This is an adaptation of the KD-tree implementation in rosetta code
///  https://rosettacode.org/wiki/K-d_tree
/// It is a reimplementation of the C code using C++.
/// It also includes a few more queries than the original

#include <algorithm>
#include <functional>
#include <list>
#include <memory>
#include <vector>

/// The point type (vector of float precision floats)
using point_t = std::vector<float>;

/// Array of indices
using indexArr = std::vector<size_t>;

/// Pair of point and Index
using pointIndex = typename std::pair<std::vector<float>, size_t>;

class KDNode {
  public:
    using KDNodePtr = std::shared_ptr<KDNode>;
    size_t index;
    point_t x;
    KDNodePtr left;
    KDNodePtr right;

    // initializer
    KDNode();
    KDNode(point_t const&, size_t const&, KDNodePtr const&, KDNodePtr const&);
    KDNode(pointIndex const&, KDNodePtr const&, KDNodePtr const&);
    ~KDNode();

    // getter
    float coord(size_t const&);

    // conversions
    explicit operator bool();
    explicit operator point_t();
    explicit operator size_t();
    explicit operator pointIndex();
};

using KDNodePtr = std::shared_ptr<KDNode>;

KDNodePtr NewKDNodePtr();

// square euclidean distance
inline float dist2(point_t const&, point_t const&);
inline float dist2(KDNodePtr const&, KDNodePtr const&);

// Need for sorting
class comparer {
  public:
    size_t idx;
    explicit comparer(size_t idx_);
    inline bool compare_idx(std::pair<std::vector<float>, size_t> const&, //
                            std::pair<std::vector<float>, size_t> const&  //
    );
};

using pointIndexArr = typename std::vector<pointIndex>;

inline void sort_on_idx(pointIndexArr::iterator const&, //
                        pointIndexArr::iterator const&, //
                        size_t idx);

using pointVec = std::vector<point_t>;

class KDTree {

  public:
    KDTree() = default;

    /// Build a KDtree
    explicit KDTree(pointVec point_array);

    /// Get the point which lies closest to the input point.
    /// @param pt input point.
    point_t nearest_point(point_t const& pt);

    /// Get the index of the point which lies closest to the input point.
    ///
    /// @param pt input point.
    size_t nearest_index(point_t const& pt);

    /// Get the point and its index which lies closest to the input point.
    ///
    /// @param pt input point.
    pointIndex nearest_pointIndex(point_t const& pt);

    /// Get both the point and the index of the points closest to the input
    /// point.
    ///
    /// @param pt input point.
    /// @param num_nearest Number of nearest points to return.
    ///
    /// @returns a vector containing the points and their respective indices
    /// which are at a distance smaller than rad to the input point.
    pointIndexArr nearest_pointIndices(point_t const& pt,
                                       size_t const& num_nearest);

    /// Get the nearest set of points to the given input point.
    ///
    /// @param pt input point.
    /// @param num_nearest Number of nearest points to return.
    ///
    /// @returns a vector containing the points which are at a distance smaller
    /// than rad to the input point.
    pointVec nearest_points(point_t const& pt, size_t const& num_nearest);

    /// Get the indices of points closest to the input point.
    ///
    /// @param pt input point.
    /// @param num_nearest Number of nearest points to return.
    ///
    /// @returns a vector containing the indices of the points which are at a
    /// distance smaller than rad to the input point.
    indexArr nearest_indices(point_t const& pt, size_t const& num_nearest);

    /// Get both the point and the index of the points which are at a distance
    /// smaller than the input radius to the input point.
    ///
    /// @param pt input point.
    /// @param rad input radius.
    ///
    /// @returns a vector containing the points and their respective indices
    /// which are at a distance smaller than rad to the input point.
    pointIndexArr neighborhood(point_t const& pt, float const& rad);

    /// Get the points that are at a distance to the input point which is
    /// smaller than the input radius.
    ///
    /// @param pt input point.
    /// @param rad input radius.
    ///
    /// @returns a vector containing the points which are at a distance smaller
    /// than rad to the input point.
    pointVec neighborhood_points(point_t const& pt, float const& rad);

    /// Get the indices of points that are at a distance to the input point
    /// which is smaller than the input radius.
    ///
    /// @param pt input point.
    /// @param rad input radius.
    ///
    /// @returns a vector containing the indices of the points which are at a
    /// distance smaller than rad to the input point.
    indexArr neighborhood_indices(point_t const& pt, float const& rad);

  private:
    KDNodePtr make_tree(pointIndexArr::iterator const& begin,
                        pointIndexArr::iterator const& end,
                        size_t const& level);

    void knearest_(KDNodePtr const& branch, point_t const& pt,
                   size_t const& level, size_t const& num_nearest,
                   std::list<std::pair<KDNodePtr, float>>& k_nearest_buffer);

    void node_query_(KDNodePtr const& branch, point_t const& pt,
                     size_t const& level, size_t const& num_nearest,
                     std::list<std::pair<KDNodePtr, float>>& k_nearest_buffer);

    // default caller
    KDNodePtr nearest_(point_t const& pt);

    void neighborhood_(KDNodePtr const& branch, point_t const& pt,
                       float const& rad2, size_t const& level,
                       pointIndexArr& nbh);

    KDNodePtr root_;
    KDNodePtr leaf_;
};
