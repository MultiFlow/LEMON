/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2003-2009
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#ifndef LEMON_MULTICOMMODITY_FLOW_H
#define LEMON_MULTICOMMODITY_FLOW_H

/// \ingroup approx
///
/// \file
/// \brief Approximation algorithms for solving multicommodity flow problems.

#include <lemon/adaptors.h>
#include <lemon/bfs.h>
#include <lemon/dijkstra.h>
#include <lemon/math.h>
#include <lemon/path.h>
#include <lemon/preflow.h>
#include <lemon/tolerance.h>
#include <vector>

//
// TODO:
//
//   - The algorithms should provide primal and dual solution as well:
//       primal objective value, primal solution (flow), primal solution (flow)
//       for each commodity dual objective value, dual solution
//     Note that these algorithms are approximation methods, so the found primal
//     and dual obj. values are not necessarily optimal, ie. the found primal
//     obj. value could be smaller than the found dual objective value and their
//     ratio provides a proof for the approximation factor of both the primal
//     and the dual solution.
//
//   - compute an upper bound for the optimal primal obj. value, e.g. using
//       x_i = min ( out_cap(s_i), in_cap(t_i) )  for each commodity i
//     the sum of x_i or the minimum of x_i/d_i ratios can be utilized.
//
//   - give back as good approximation guarantee r as possible for which:
//     dual = r*primal
//
//   - revise the implementations and their API
//

namespace lemon {

/// \addtogroup approx
/// @{

/// \brief Implementation of an approximation algorithm for finding
/// a maximum multicommodity flow.
///
/// \ref MaxMulticommodityFlow implements Lisa Fleischer's algorithm
/// for finding a maximum multicommodity flow.
///
/// \tparam Digraph The type of the digraph the algorithm runs on.
/// \tparam CAP The type of the capacity map.
/// The default map type is \ref concepts::Digraph::ArcMap
/// "GR::ArcMap<double>".
/// \tparam V The value type used in the algorithm.
/// The default value type is \c CAP::Value.
#ifdef DOXYGEN
template <typename GR, typename CAP, typename V>
#else
template <typename GR, typename CAP = typename GR::template ArcMap<double>,
          typename V = typename CAP::Value>
#endif
class MaxMulticommodityFlow {
  TEMPLATE_DIGRAPH_TYPEDEFS(GR);

  typedef SimplePath<GR> SPath;
  typedef typename SPath::ArcIt PathArcIt;

public:
  ///  The type of the digraph the algorithm runs on.
  typedef GR Digraph;
  /// The type of the capacity map.
  typedef CAP CapacityMap;
  /// The value type of the algorithm.
  typedef V Value;

#ifdef DOXYGEN
  /// The type of the flow map (primal solution).
  typedef Digraph::ArcMap<Value> FlowMap;
  /// The type of the length map (dual soluiton).
  typedef Digraph::ArcMap<Value> LengthMap;
#else
  /// The type of the flow map (primal solution).
  typedef typename Digraph::template ArcMap<Value> FlowMap;
  /// The type of the length map (dual soluiton).
  typedef typename Digraph::template ArcMap<Value> LengthMap;
#endif

private:
  // The digraph the algorithm runs on
  const Digraph &_g;
  // The capacity map
  const CapacityMap &_cap;

  // The number of commdities
  int _k;
  // The source nodes for each commodity
  std::vector<Node> _s;
  // The sink nodes for each commodity
  std::vector<Node> _t;
  // The weight of each commodity
  // std::vector<Weight> _w;

  // The flow map (the primal solution).
  FlowMap _flow;
  // The length map (the dual solution).
  LengthMap _len;

  struct PathData {
    SPath path;
    Value value;
    Value min_cap;
    int hash;

    PathData() {}
    PathData(const SPath &p, Value v, Value mc, int h = 0)
        : path(p), value(v), min_cap(mc), hash(h) {}
  };

public:
  /// \brief Constructor.
  ///
  /// Constructor.
  ///
  /// \param digraph The digraph the algorithm runs on.
  /// \param capacity The capacities of the edges.
  MaxMulticommodityFlow(const Digraph &digraph, const CapacityMap &capacity)
      : _g(digraph), _cap(capacity), _k(0), _flow(_g), _len(_g) {}

  /// \brief Constructor.
  ///
  /// Constructor.
  ///
  /// \param gr The digraph the algorithm runs on.
  /// \param capacity The capacities of the edges.
  /// \param k The number of commodities.
  /// \param sources A vector containing the source nodes
  /// (its size must be at least \c k).
  /// \param sinks A vector containing the sink nodes
  /// (its size must be at least \c k).
  MaxMulticommodityFlow(const Digraph &gr, const CapacityMap &capacity, int k,
                        const std::vector<Node> &sources,
                        const std::vector<Node> &sinks)
      : _g(gr), _cap(capacity), _k(k), _s(sources), _t(sinks), _flow(_g),
        _len(_g) {
    _s.resize(k);
    _t.resize(k);
  }

  /// Destructor.
  ~MaxMulticommodityFlow() {}

  /// \brief Add a new commodity.
  ///
  /// This function adds a new commodity with the given source and
  /// sink nodes.
  void addCommodity(const Node &s, const Node &t) {
    ++_k;
    _s.push_back(s);
    _t.push_back(t);
  }

  /// \name Execution control

  /// @{

  /// \brief Run the algorithm for finding a maximum multicommodity flow.
  ///
  /// This function runs the algorithm for finding a maximum multicommodity
  /// flow.
  /// \param r The approximation factor. It must be greater than 1.
  /// \return The effective approximation factor (at least 1), i.e. the
  /// ratio of the values of the found dual and primal solutions.
  /// Therefore both the primal and dual solutions are r-approximate.
  ///
  /// \note The smaller parameter \c r is, the slower the algorithm runs,
  /// but it might find better solution. According to our tests, using
  /// the default value the algorithm performs well and it usually manage
  /// to find the exact optimum due to heuristic improvements.
  Value run(Value r = 2.0, Value opt = 0) {

#define MMF_FLEISCHER_IMPR
#define MMF_HEURISTIC_IMPR
//#define MMF_EXP_UPDATE
#define MMF_DEBUG_OUTPUT

    // Initialize epsilon and delta parameters
    Tolerance<Value> ftol, ctol;
    if (!ftol.less(1.0, r))
      return 0;
#ifdef MMF_FLEISCHER_IMPR
#ifdef MMF_EXP_UPDATE
    Value epsilon = (r + 1 - 2 * pow(r, 0.5)) / (r - 1);
#else
    Value epsilon = (2 * r + 1 - pow(8 * r + 1, 0.5)) / (2 * r);
#endif
#else
#ifdef MMF_EXP_UPDATE
    Value epsilon = (2 * r + 1 - pow(8 * r + 1, 0.5)) / (2 * r);
#else
    Value epsilon = 1 - pow(r, -0.5);
#endif
#endif

#ifdef MMF_DEBUG_OUTPUT
    std::cout << "r = " << r << "; epsilon = " << epsilon << "\n";
#endif

    int node_num = countNodes(_g);
    Value delta =
        pow(1 + epsilon, 1 - 1 / epsilon) / pow(node_num, 1 / epsilon);
    ctol.epsilon(delta * 0.0001);

    // Initialize flow values and lengths
    for (ArcIt a(_g); a != INVALID; ++a) {
      _flow[a] = 0;
      _len[a] = delta;
    }

    // Check commodities
    for (int i = 0; i < _k; ++i) {
      if (_s[i] == _t[i] || !bfs(_g).run(_s[i], _t[i])) {
        _s[i] = _s[_k - 1];
        _t[i] = _t[_k - 1];
        --_k;
      }
    }
    if (_k <= 0)
      return 0;

    // Successive improvement along shortest paths
    std::vector<std::vector<PathData>> paths(_k);
#ifdef MMF_FLEISCHER_IMPR
    Value limit = delta * (1 + epsilon);
    int i = 0;
#else
    std::vector<SPath> curr_paths(_k);
#endif
    SPath curr_path;
    Value curr_dist;
    int iter = 0;
    while (true) {
#ifdef MMF_FLEISCHER_IMPR
      // Find the shortest path for the i-th commodity
      // (w.r.t the current length function)
      dijkstra(_g, _len).path(curr_path).dist(curr_dist).run(_s[i], _t[i]);

      // Check if the path can be used and increase i or the
      // length limit if necessary
      if (!ctol.less(curr_dist, limit)) {
        if (++i == _k) {
          if (!ctol.less(limit, 1.0))
            break;
          i = 0;
          limit *= (1 + epsilon);
          if (ctol.less(1.0, limit))
            limit = 1.0;
        }
        continue;
      }
#else
      // Find the shortest path
      Value min_dist = 1.0;
      int i = -1;
      for (int j = 0; j < _k; ++j) {
        if (!dijkstra(_g, _len)
                 .path(curr_paths[j])
                 .dist(curr_dist)
                 .run(_s[j], _t[j])) {
          std::cerr << "\n\nERROR!!! " << j << ":" << _k << "\n\n";
        }
        if (curr_dist < min_dist) {
          min_dist = curr_dist;
          i = j;
        }
      }
      if (i < 0)
        break;
      curr_path = curr_paths[i];
#endif

      ++iter;

      // Check if the path is used before
      int h = hash(curr_path);
      int j;
      for (j = 0; j < int(paths[i].size()); ++j) {
        bool b = (paths[i][j].hash == h) &&
                 (paths[i][j].path.length() == curr_path.length());
        if (b) {
          for (PathArcIt a1(paths[i][j].path), a2(curr_path);
               b && a1 != INVALID && a2 != INVALID; ++a1, ++a2)
            b = _g.id(a1) == _g.id(a2);
          if (b)
            break;
        }
      }

      // Set/modify the flow value for the found path
      Value gamma;
      if (j == int(paths[i].size())) {
        gamma = std::numeric_limits<Value>::infinity();
        for (PathArcIt a(curr_path); a != INVALID; ++a) {
          if (_cap[a] < gamma)
            gamma = _cap[a];
        }
        paths[i].push_back(PathData(curr_path, gamma, gamma, h));
      } else {
        gamma = paths[i][j].min_cap;
        paths[i][j].value += gamma;
      }

      // Modify the arc lengths along the current path
#ifdef MMF_EXP_UPDATE
      for (PathArcIt a(curr_path); a != INVALID; ++a)
        _len[a] *= exp(epsilon * gamma / _cap[a]);
#else
      for (PathArcIt a(curr_path); a != INVALID; ++a)
        _len[a] *= 1 + epsilon * gamma / _cap[a];
#endif
    } // end while

    // Setting flow values
    for (int i = 0; i < _k; ++i) {
      for (int j = 0; j < int(paths[i].size()); ++j) {
        Value val = paths[i][j].value;
        for (PathArcIt a(paths[i][j].path); a != INVALID; ++a)
          _flow[a] += val;
      }
    }

    // Scaling the solution
    Value lambda = 0;
    for (ArcIt a(_g); a != INVALID; ++a) {
      if (_flow[a] / _cap[a] > lambda)
        lambda = _flow[a] / _cap[a];
    }
    if (ftol.positive(lambda)) {
      lambda = 1 / lambda;
      for (ArcIt a(_g); a != INVALID; ++a)
        _flow[a] *= lambda;
      for (int i = 0; i < _k; ++i) {
        for (int j = 0; j < int(paths[i].size()); ++j)
          paths[i][j].value *= lambda;
      }
    }

#ifdef MMF_DEBUG_OUTPUT
    Value value1 = 0;
    for (int i = 0; i < _k; ++i) {
      for (InArcIt a(_g, _t[i]); a != INVALID; ++a)
        value1 += _flow[a];
      for (OutArcIt a(_g, _t[i]); a != INVALID; ++a)
        value1 -= _flow[a];
    }
    Value max_opt = value1 / (1 - epsilon);

    std::cout << "\n";
    std::cout << "Computation finished. Iterations: " << iter << "\n";
    std::cout << "Solution value:   " << value1;
    if (opt > 0)
      std::cout << "  \t" << value1 / opt * 100 << "%";
    std::cout << std::endl;
#else
    opt = opt; // avoid warning
#endif

#ifdef MMF_DEBUG_OUTPUT
    Value tmp = 0;
    for (ArcIt a(_g); a != INVALID; ++a)
      tmp += _cap[a] * _len[a];
    std::cout << "Final dual: " << tmp << "\n";
#endif

#ifdef MMF_HEURISTIC_IMPR
    // Heuristic improvement 1: improve the solution along one of the
    // found paths if it is possible
    int opi = 0;
    for (int i = 0; i < _k; ++i) {
      for (int j = 0; j < int(paths[i].size()); ++j) {
        Value d = std::numeric_limits<Value>::infinity();
        for (PathArcIt a(paths[i][j].path); a != INVALID; ++a)
          if (_cap[a] - _flow[a] < d)
            d = _cap[a] - _flow[a];
        if (ftol.positive(d)) {
          paths[i][j].value += d;
          for (PathArcIt a(paths[i][j].path); a != INVALID; ++a)
            _flow[a] += d;
          ++opi;
        }
      }
    }

    // Heuristic improvement 2: improve the solution along new paths
    // if it is possible
    int npi = 0;
    BoolArcMap filter(_g);
    for (ArcIt a(_g); a != INVALID; ++a)
      filter[a] = ftol.positive(_cap[a] - _flow[a]);
    FilterArcs<const Digraph, BoolArcMap> sub_graph(_g, filter);
    SPath p;
    for (int i = 0; i < _k; ++i) {
      while (bfs(sub_graph).path(p).run(_s[i], _t[i])) {
        Value d = std::numeric_limits<Value>::infinity();
        for (PathArcIt a(p); a != INVALID; ++a)
          if (_cap[a] - _flow[a] < d)
            d = _cap[a] - _flow[a];
        paths[i].push_back(PathData(p, d, d));
        for (PathArcIt a(p); a != INVALID; ++a) {
          _flow[a] += d;
          filter[a] = ftol.positive(_cap[a] - _flow[a]);
        }
        ++npi;
      }
    }
#endif

#ifdef MMF_DEBUG_OUTPUT
    Value value2 = 0;
    for (int i = 0; i < _k; ++i) {
      for (InArcIt a(_g, _t[i]); a != INVALID; ++a)
        value2 += _flow[a];
      for (OutArcIt a(_g, _t[i]); a != INVALID; ++a)
        value2 -= _flow[a];
    }
    epsilon = 1 - value2 / max_opt;

#ifdef MMF_HEURISTIC_IMPR
    std::cout << "\nHeuristics applied. Old impr: " << opi
              << "  new impr: " << npi << "\n";
#endif
    std::cout << "Solution value:   " << value2;
    if (opt > 0)
      std::cout << "  \t" << value2 / opt * 100 << "%";
    std::cout << std::endl << std::endl;
#endif

    // TODO: better return value
    return r;
  }

  /// @}

private:
  int hash(const SPath &path) {
    int h = path.length();
    for (PathArcIt a(path); a != INVALID; ++a)
      h = (h << 4) ^ (h >> 28) ^ _g.id(a);
    return h;
  }

public:
  /// \name Query Functions
  /// The result of the algorithm can be obtained using these
  /// functions.\n
  /// \ref run() must be called before using them.

  /// @{

  /// \brief Return the flow value on the given arc.
  ///
  /// This function returns the flow value on the given arc.
  ///
  /// \pre \ref run() must be called before using this function.
  Value flow(const Arc &arc) const { return _flow[arc]; }

  /// \brief Return a const reference to an arc map storing the
  /// found flow.
  ///
  /// This function returns a const reference to an arc map storing
  /// the found flow.
  ///
  /// \pre \ref run() must be called before using this function.
  const FlowMap &flowMap() const { return _flow; }

  /// \brief Return the length of the given arc (the dual value).
  ///
  /// This function returns the length of the given arc (the dual value).
  ///
  /// \pre \ref run() must be called before using this function.
  Value length(const Arc &arc) const { return _len[arc]; }

  /// \brief Return a const reference to an arc map storing the
  /// found lengths (the dual solution).
  ///
  /// This function returns a const reference to an arc map storing
  /// the found lengths (the dual solution).
  ///
  /// \pre \ref run() must be called before using this function.
  const LengthMap &lengthMap() const { return _len; }

  /// @}

}; // class MaxMulticommodityFlow

///@}

} // namespace lemon

#endif // LEMON_MULTICOMMODITY_FLOW_H
