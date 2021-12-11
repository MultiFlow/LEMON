/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2003-2021
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
 * Additional Copyright file multicommodity_flow.h: Lancaster University (2021)
 */

#ifndef LEMON_MULTICOMMODITY_FLOW_H_
#define LEMON_MULTICOMMODITY_FLOW_H_

/// \ingroup approx
///
/// \file
/// \brief Fleischer's algorithms for finding approximate multicommodity flows.

#include <algorithm>  // min
#include <limits>     // numeric_limits
#include <random>     // random_device, mt19937
#include <set>
#include <vector>

#include "lemon/adaptors.h"
#include "lemon/bfs.h"
#include "lemon/concepts/maps.h"
#include "lemon/dijkstra.h"
#include "lemon/maps.h"
#include "lemon/math.h"
#include "lemon/path.h"
#include "lemon/preflow.h"
#include "lemon/static_graph.h"
#include "lemon/tolerance.h"

namespace lemon {

/// \brief Default traits class of \ref ApproxMCF class.
///
/// Default traits class of \ref ApproxMCF class.
/// \tparam GR The type of the digraph.
/// \tparam CAPType The type of the capacity values.
/// \tparam CAPMap The type of the capacity map.
/// \tparam V The type for the flows and costs.
/// It must conform to the \ref concepts::ReadMap "ReadMap" concept.
#ifdef DOXYGEN
template<typename GR, typename CAPType, typename CAPMap, typename V>
#else
template<
    typename GR,
    typename CAPMap = typename GR::template ArcMap<double>,
    typename V      = double>
#endif
struct ApproxMCNFDefaultTraits {
  /// The type of the digraph
  typedef GR Digraph;
  /// The type of the capacity map
  typedef CAPMap CapacityMap;
  /// The type of the capacity map values. By default is double.
  typedef typename CapacityMap::Value CapacityValue;

  /// \brief The type used for flow values and costs.
  ///
  /// \c CapacityValue must be convertible to \c Value.
  /// By default is \c double.
  typedef V                                        Value;
  typedef typename Digraph::template ArcMap<Value> ValueMap;

  /// The tolerance type used for internal computations
  typedef lemon::Tolerance<Value> Tolerance;

  /// \brief The path type of the decomposed flows
  ///
  /// It must conform to the \ref lemon::concepts::Path "Path" concept
  typedef lemon::Path<Digraph>                 Path;
  typedef typename lemon::Path<Digraph>::ArcIt PathArcIt;

  /// \brief Instantiates a ValueMap.
  ///
  /// This function instantiates a \ref ValueMap.
  /// \param digraph The digraph for which we would like to define
  /// the large cost map. This can be used for flows or otherwise.
  static ValueMap* createValueMapPtr(
      const Digraph& digraph,
      const Value&   value = 0) {
    return new ValueMap(digraph, value);
  }
  // Same as above just const
  static const ValueMap* createValueMapConstPtr(
      const Digraph& digraph,
      const Value&   value = 0) {
    return new ValueMap(digraph, value);
  }
};

/// \addtogroup approx
/// @{
/// \brief Fleischer's algorithms for the multicommodity flow problem.
///
/// Default traits class of \ref ApproxMCF class.
/// \tparam GR The type of the digraph.
/// \tparam TR The type of the traits class.
#ifdef DOXYGEN
template<typename GR, typename TR>
#else
template<
    typename GR,
    typename CM = typename GR::template ArcMap<double>,
    typename TR = ApproxMCNFDefaultTraits<GR, CM>>
#endif
class ApproxMCF {
 public:
  /// The type of the digraph
  typedef typename TR::Digraph Digraph;
  /// The type of the capacity map
  typedef typename TR::CapacityMap CapacityMap;
  /// The type of the capacity values
  typedef typename TR::CapacityValue CapacityValue;

  /// \brief The large cost type
  ///
  /// The large cost type used for internal computations.
  /// By default is \c double.
  typedef typename TR::Value Value;
  /// Arc Map with values of type \ref Value
  typedef typename TR::ValueMap ValueMap;

  /// The tolerance type
  typedef typename TR::Tolerance Tolerance;

  /// \brief The path type of the found cycles
  ///
  /// The path type of the found cycles.
  /// Using the \ref lemon::ApproxMCNFDefaultTraits "default traits class",
  /// it is \ref lemon::Path "Path<Digraph>".
  typedef typename TR::Path Path;
  /// The arc iterator for paths.
  typedef typename TR::PathArcIt PathArcIt;

  /// The \ref lemon::ApproxMCNFDefaultTraits "traits class" of the
  /// algorithm
  typedef TR Traits;

  /// Types of approximation type.
  enum ApproxType {
    /// Maximum flow using Fleischer's algorithm.
    FLEISCHER_MAX,
    /// Maximum concurrent flow using Fleischer's algorithm.
    FLEISCHER_MAX_CONCURRENT,
    /// Minimum cost concurrent flow using binary search with cost-bounded
    /// \ref FLEISCHER_MAX_CONCURRENT
    BINARY_SEARCH_MIN_COST
  };

 private:
  struct PathData {
    Path path;
    // The amount of flow through the path
    Value value;
    // The minimum capacity along all edges in the path
    CapacityValue min_cap;
    // The actual cost of the path
    Value cost;
    int   hash;

    PathData() {}
    PathData(const Path& p, Value v, CapacityValue m, Value c, int h = 0)
        : path(p), value(v), min_cap(m), cost(c), hash(h) {}
  };

  int hash(const Path& path) const {
    int h = path.length();
    for (PathArcIt a(path); a != INVALID; ++a)
      h = (h << 4) ^ (h >> 28) ^ _graph.id(a);
    return h;
  }

  // Matrix indexed by commodity with all the paths found. i.e. [i][p] where i
  // is the commodity index and p is the path index.
  typedef std::vector<std::vector<PathData>> PathMatrix;
  // Vector indexed by commodity with \ref ValueMap pointers.
  typedef std::vector<ValueMap*>       ValueMapVector;
  typedef std::vector<const ValueMap*> ConstValueMapVector;

  TEMPLATE_DIGRAPH_TYPEDEFS(Digraph);

  const Digraph& _graph;
  const int      _num_nodes;
  const int      _num_arcs;
  // Input capacity map
  const CapacityMap& _cap;

  // The number of commodities
  int _num_commodities;
  // The source nodes for each commodity
  std::vector<Node> _source;
  // The target nodes for each commodity
  std::vector<Node> _target;
  // The demands for each commodity (may be 0, or undefined).
  std::vector<CapacityValue> _demand;
  // The cost for each commodity (may be 0, or undefined)
  // See \ref ValueMapVector.
  ConstValueMapVector _cost;

  /**
   * Algorithm parameters
   */
  // type of approximation algorithm. See \ref ApproxType.
  ApproxType _approx_type;
  // Stopping criterion
  Value _epsilon;
  // Whether the internal heuristics will be applied at the end of the
  // algorithms
  bool _heuristic;
  // Whether there is feasible unsplitable flow
  bool _feasible_unsplit;
  // Objective function of the primal problem
  Value _obj;
  // Objective function of the dual problem
  Value _dual_obj;
  bool  _local_cost;
  //
  Value _delta;
  // Expected maximum number of iterations
  Value _max_iter;
  //
  Value _limit;
  // Number of iterations
  int _iter;
  // Tolerance function. See \ref lemon::Tolerance.
  Tolerance _tol;

  // For the BINARY_SEARCH_MIN_COST case
  // Cost of the final solution.
  Value _total_cost;
  // Cost bounding factor.
  Value _phi;
  // Cost bound.
  Value _cost_bound;

  // The aggregated flow per arc over all commodities
  ValueMap _total_flow;
  // The value of the length per arc (the dual solution).
  ValueMap _len;
  // The flow for each commodity. See \ref ValueMapVector.
  ValueMapVector _flow;
  // The paths found for different commodities. See \ref PathMatrix.
  PathMatrix _paths;

  // Solution Attributes
  std::vector<Path>  _final_paths;
  std::vector<Value> _final_flows;
  std::vector<Value> _final_costs;

 public:
  /// \brief Default constructor
  ApproxMCF(const Digraph& gr, const CapacityMap& cp)
      : _graph(gr),
        _num_nodes(countNodes(_graph)),
        _num_arcs(countArcs(_graph)),
        _cap(cp),
        _num_commodities(0),
        _feasible_unsplit(true),
        _obj(0),
        _dual_obj(0),
        _local_cost(true),
        _iter(0),
        _cost_bound(0),
        _total_flow(_graph),
        _len(_graph) {}

  ~ApproxMCF() {
    for (int i = 0; i < _num_commodities; ++i) {
      delete _flow[i];
      _flow[i] = nullptr;
      if (_local_cost) {
        delete _cost[i];
        _cost[i] = nullptr;
      }
    }
    _flow.clear();
    _cost.clear();
  }

  /// \name Parameters
  /// The parameters of the algorithm can be specified using these
  /// functions.

  /// @{

  /// \brief Set single commodity parameters.
  ///
  /// \param s The source node.
  /// \param t The target node.
  /// \param d The required amount of flow from node \c s to node \c t
  /// (i.e. the supply of \c s and the demand of \c t).
  /// \param cost_map The cost arc map.
  ///
  /// \return <tt>(*this)</tt>
  ApproxMCF& addCommodity(
      const Node&          s,
      const Node&          t,
      const CapacityValue& d,
      const ValueMap&      cost_map) {
    _source.push_back(s);
    _target.push_back(t);
    if (d > 0)
      _demand.push_back(d);
    _local_cost = false;
    _cost.push_back(&cost_map);
    ++_num_commodities;
    return *this;
  }

  /// @overload no cost map, optional demand
  ApproxMCF&
  addCommodity(const Node& s, const Node& t, const CapacityValue& d = 0) {
    const ValueMap* cost_map_ptr = Traits::createValueMapConstPtr(_graph);
    addCommodity(s, t, d, *cost_map_ptr);
    _local_cost = true;
    return *this;
  }

  /// @}

  /// \name Execution control

  /// @{

  /// \brief Run the algorithm for finding a multicommodity flow.
  ///
  /// This function runs the algorithm for finding a maximum/minimum cost
  /// multicommodity flow.
  ///
  /// \param at  The type of approximation algorithm to use. @see \ref
  /// ApproxMCF::ApproxType.
  /// \param heuristic Whether the internal heuristics will be applied at the
  // end of Fleischer's algorithm.
  /// \param epsilon The stopping criterion. Must be positive and larger than
  /// 0. The smaller the value, the more accurate the result, but the longer the
  /// execution time. By default is 0.1.
  void run(
      const ApproxType& at        = FLEISCHER_MAX,
      const bool&       heuristic = true,
      const Value&      epsilon   = 0.1) {
    _approx_type = at;
    _epsilon     = epsilon;
    _heuristic   = heuristic;
    check();
    switch (_approx_type) {
      case FLEISCHER_MAX:
        runFleischerMax();
        break;
      case BINARY_SEARCH_MIN_COST:
        runBinarySearchMinCost();
        break;
      case FLEISCHER_MAX_CONCURRENT:
        runFleischerMaxConcurrent();
        break;
    }
  }

  /// @}

 private:
  // Initialize flow values and lengths
  void init() {
    _iter     = 0;
    _obj      = 0;
    _dual_obj = 0;
    // Allocate paths
    _paths.resize(_num_commodities);
    // Solution paths
    _final_paths.resize(_num_commodities);
    _final_flows.resize(_num_commodities);
    _final_costs.resize(_num_commodities);

    for (int i = 0; i < _num_commodities; ++i)
      _paths[i].clear();
    _delta =
        pow(1 + _epsilon, 1 - 1 / _epsilon) / pow(_num_nodes, 1 / _epsilon);
    _tol.epsilon(_delta * 1e-5);
    _limit = _delta * (1 + _epsilon);

    // Allocate/reset flow per commodity and arc
    _flow.resize(_num_commodities);
    for (int i = 0; i < _num_commodities; ++i)
      if (!_flow[i]) {
        _flow[i] = Traits::createValueMapPtr(_graph);
      } else {
        for (ArcIt a(_graph); a != INVALID; ++a) {
          (*_flow[i])[a] = 0;
        }
      }

    // Initialise length and flow map
    for (ArcIt a(_graph); a != INVALID; ++a) {
      _len[a]        = _delta;
      _total_flow[a] = 0;
    }

    // Initialise phi for cost bounding
    if (_approx_type == BINARY_SEARCH_MIN_COST && _cost_bound > 0)
      _phi = _delta / _cost_bound;
    else
      _phi = 0;
  }

  // Preprocessing checks to ensure that:
  //
  // 1. Commodities are well-defined
  // 2. Input sizes are consistent
  // 3. Concurrent case: Feasible split flow exists for each commodity that
  // satisfies demand. In the case this doesn't pass, will throw an exception.
  //
  // TODO(): Group commodities by source, so we only have to solve the shortest
  // path once for each group (and length values).
  void check() {
    LEMON_ASSERT(_tol.positive(_epsilon), "Epsilon must be greater than 0");

    // Remove commodities with source = sink or unreachable
    for (int i = 0; i < _num_commodities; ++i) {
      if (_source[i] == _target[i] ||
          !bfs(_graph).run(_source[i], _target[i])) {
        _source[i] = _source[_num_commodities - 1];
        _target[i] = _target[_num_commodities - 1];
        --_num_commodities;
      }
    }

    const int& s_size = _source.size();
    const int& t_size = _target.size();
    const int& d_size = _demand.size();

    switch (_approx_type) {
      case FLEISCHER_MAX: {
        LEMON_ASSERT(_tol.positive(_num_commodities), "No commodities");
        // Check consistency of vector sizes
        LEMON_ASSERT(
            (_num_commodities == s_size && s_size == t_size &&
             !_tol.positive(d_size)),
            "All vectors should be of size equal to the number of commodities, "
            "and have no demand.");
        break;
      }
      case FLEISCHER_MAX_CONCURRENT:
      case BINARY_SEARCH_MIN_COST: {
        // Check consistency of vector sizes
        LEMON_ASSERT(
            (_num_commodities == s_size && s_size == t_size &&
             s_size == d_size),
            "All vectors should be of size equal to the number of commodities");
        // Remove commodities with 0 demand.
        for (int i = 0; i < _num_commodities; ++i) {
          if (!_tol.positive(_demand[i])) {
            _source[i] = _source[_num_commodities - 1];
            _target[i] = _target[_num_commodities - 1];
            --_num_commodities;
          }
        }
        LEMON_ASSERT(_tol.positive(_num_commodities), "No commodities");

        // Check if feasible routing of demand exists for each commodity
        for (int i = 0; i < _num_commodities; ++i) {
          // Only check commodity if the demand in non-zero
          // Check if feasible split flow exists
          // Total capacity outgoing at the source
          CapacityValue total_cap_source = 0;
          CapacityValue max_cap_source   = 0;
          for (OutArcIt a(_graph, _source[i]); a != INVALID; ++a) {
            total_cap_source += _cap[a];
            if (max_cap_source < _cap[a])
              max_cap_source = _cap[a];
          }

          // Total capacity incoming at the target
          CapacityValue total_cap_target = 0;
          CapacityValue max_cap_target   = 0;
          for (InArcIt a(_graph, _target[i]); a != INVALID; ++a) {
            total_cap_target += _cap[a];
            if (max_cap_target < _cap[a])
              max_cap_target = _cap[a];
          }

          // Check if enough capacity to satisfy the demand
          LEMON_ASSERT(
              (total_cap_source >= _demand[i] &&
               total_cap_target >= _demand[i]),
              "No feasible (split) flow to satisfy demand for commodity" +
                  std::to_string(i));

          // Check if feasible unsplittable flow exists that satisfies demand
          if (max_cap_source < _demand[i] || max_cap_target < _demand[i])
            _feasible_unsplit = false;
        }
      }
    }
  }

  // Run standard max concurrent algorithm iteratively updating the budget, or
  // cost bound, a la binary search to obtain an approximation to the
  // multicommodity flow minimum cost
  void runBinarySearchMinCost() {
    // Run standard max concurrent
    runFleischerMaxConcurrent();

    Value lower_bound = 0, upper_bound = getTotalCost();

    // Best solution attributes
    Value          curr_best_cost = std::numeric_limits<Value>::infinity();
    ValueMapVector curr_best_flow;

    while ((upper_bound - lower_bound) / upper_bound > _epsilon &&
           upper_bound >= lower_bound) {
      // Run Max Concurrent algorithm with cost bound
      _cost_bound = (upper_bound + lower_bound) / 2;
      runFleischerMaxConcurrent();
      // Update upper/lower bounds
      if (_total_cost < _cost_bound)
        upper_bound = _cost_bound;
      else
        lower_bound = _cost_bound;
      // Save current solution if appropriate
      if (_total_cost < curr_best_cost) {
        curr_best_cost = _total_cost;
        curr_best_flow = _flow;
      }
    }
    // Update flow
    _total_cost = curr_best_cost;
    _flow       = curr_best_flow;
  }

  void runFleischerMaxConcurrent() {
    // Successive improvement along shortest paths
    init();
    Path  curr_path;
    Value curr_dist;
    while (true) {
      // Find the shortest path for the i-th commodity
      // (w.r.t the current length function)
      for (int i = 0; i < _num_commodities; ++i) {
        // Route demand
        CapacityValue d_i = _demand[i];
        while (true) {  // while (phase)
          dijkstra(_graph, _len)
              .path(curr_path)
              .dist(curr_dist)
              .run(_source[i], _target[i]);
          // Get
          const CapacityValue& routed_flow = updatePaths(curr_path, i, d_i);
          // Update demand
          d_i -= routed_flow;
          updateLen(curr_path, routed_flow, i);
          updateDualObjective();
          if (d_i <= 0 || !_tol.less(_dual_obj, 1))
            break;
        }  // end while (phase)
      }
      ++_iter;
      updateDualObjective();
      if (!_tol.less(_dual_obj, 1))
        break;
    }  // end while (main)

    // Clean up
    scaleFinalFlowsMaxConcurrent();

    if (_heuristic)
      randomisedRounding();

    setFlowPerCommodity();
    setFinalCost();
    saveFinalPaths();
  }

  void runFleischerMax() {
    // Successive improvement along shortest paths
    int   i = 0;
    Path  curr_path;
    Value curr_dist;
    init();
    while (true) {  // start while
      dijkstra(_graph, _len)
          .path(curr_path)
          .dist(curr_dist)
          .run(_source[i], _target[i]);

      // Check if the path can be used and increase i or the length limit if
      // necessary
      if (!_tol.less(curr_dist, _limit)) {
        if (++i == _num_commodities) {
          if (!_tol.less(_limit, 1.0))
            break;
          i = 0;
          _limit *= (1 + _epsilon);
          if (_tol.less(1.0, _limit))
            _limit = 1.0;
        }
        continue;
      }
      const CapacityValue& routed_flow = updatePaths(curr_path, i);
      // Update dual values
      updateLen(curr_path, routed_flow, i);
      ++_iter;
    }  // end while

    scaleFinalFlowsMaxFlow();

    if (_heuristic) {
      runHeurMax1();
      runHeurMax2();
    }

    setFlowPerCommodity();

    // Get dual objective
    _dual_obj = getDualObjective();
    // Get objective value
    _obj = 0;
    for (int i = 0; i < _num_commodities; ++i) {
      for (InArcIt a(_graph, _target[i]); a != INVALID; ++a)
        _obj += (*_flow[i])[a];
      for (OutArcIt a(_graph, _target[i]); a != INVALID; ++a)
        _obj -= (*_flow[i])[a];
    }
  }

  // Randomised rounding - if flow is split over multiple paths, randomised
  // rounding randomly rounds up the largest proportion of flow/cost.
  //
  // Only runs if flow can be routed without splitting.
  // If costs are given, to make expensive paths less likely to occur, we use,
  // for each path j a value[j] = flow / cost for the discrete distribution.
  // This distribution samples path j with a probability value[j] / sum values
  void randomisedRounding() {
    // Only works if flow can be routed without splitting.
    if (!_feasible_unsplit)
      return;
    for (int i = 0; i < _num_commodities; ++i) {
      const int&         p_size       = static_cast<int>(_paths[i].size());
      Value              total_flow_i = 0;
      std::vector<Value> values(p_size, 0);

      for (int j = 0; j < p_size; ++j) {
        const Value& val  = _paths[i][j].value;
        const Value& cost = _paths[i][j].cost;
        total_flow_i += val;
        values[j] = val;
        if (_tol.positive(cost))
          values[j] /= cost;
      }
      // Get random index wrt values[j] / sum of values
      int                          chosen_p_idx;
      std::random_device           rd;
      std::mt19937                 gen(rd());
      std::discrete_distribution<> d(values.begin(), values.end());
      chosen_p_idx = d(gen);
      // Check routing all demand through chosen_p_idx remains capacity feasible
      bool capacity_feasible = true;
      // cost is all demand routed through here
      for (PathArcIt a(_paths[i][chosen_p_idx].path); a != INVALID; ++a) {
        if (_total_flow[a] - _paths[i][chosen_p_idx].value + _demand[i] >
            _cap[a]) {
          capacity_feasible = false;
          break;
        }
      }
      if (capacity_feasible) {
        for (int j = 0; j < p_size; ++j) {
          // Update total flows and path flow
          if (j == chosen_p_idx) {  // path j routes demand
            for (PathArcIt a(_paths[i][j].path); a != INVALID; ++a) {
              _total_flow[a] = _total_flow[a] - _paths[i][j].value + _demand[i];
            }
            _paths[i][j].value = _demand[i];
          } else {  // path j flow is 0
            for (PathArcIt a(_paths[i][j].path); a != INVALID; ++a) {
              _total_flow[a] = _total_flow[a] - _paths[i][j].value;
            }
            _paths[i][j].value = 0;
          }
        }
      }
    }
  }

  // Heuristic improvement 1 (for FLEISCHER_MAX only).
  //
  // Attempts to improve the solution along one of the found paths if it is
  // possible
  void runHeurMax1() {
    for (int i = 0; i < _num_commodities; ++i) {
      for (int j = 0; j < static_cast<int>(_paths[i].size()); ++j) {
        CapacityValue d = std::numeric_limits<CapacityValue>::infinity();
        for (PathArcIt a(_paths[i][j].path); a != INVALID; ++a)
          if (_cap[a] - _total_flow[a] < d)
            d = _cap[a] - _total_flow[a];
        if (_tol.positive(d)) {
          _paths[i][j].value += d;
          for (PathArcIt a(_paths[i][j].path); a != INVALID; ++a)
            _total_flow[a] += d;
        }
      }
    }
  }

  // Heuristic improvement 2 (for FLEISCHER_MAX only).
  //
  // Attempts to improve the solution along new paths
  void runHeurMax2() {
    BoolArcMap filter(_graph);
    for (ArcIt a(_graph); a != INVALID; ++a)
      filter[a] = _tol.positive(_cap[a] - _total_flow[a]);
    FilterArcs<const Digraph, BoolArcMap> sub_graph(_graph, filter);
    Path                                  p;
    for (int i = 0; i < _num_commodities; ++i) {
      while (bfs(sub_graph).path(p).run(_source[i], _target[i])) {
        CapacityValue d = std::numeric_limits<CapacityValue>::infinity();
        for (PathArcIt a(p); a != INVALID; ++a)
          if (_cap[a] - _total_flow[a] < d)
            d = _cap[a] - _total_flow[a];
        _paths[i].push_back(PathData(p, d, d, 0, 0));
        for (PathArcIt a(p); a != INVALID; ++a) {
          _total_flow[a] += d;
          filter[a] = _tol.positive(_cap[a] - _total_flow[a]);
        }
      }
    }
  }

  // Update paths for current commodity and return flow routed
  Value updatePaths(
      const Path&          curr_path,
      const int&           i,
      const CapacityValue& d_i =
          std::numeric_limits<CapacityValue>::infinity()) {
    int h = hash(curr_path);
    int j;  // path index
    for (j = 0; j < static_cast<int>(_paths[i].size()); ++j) {
      bool b = (_paths[i][j].hash == h) &&
               (_paths[i][j].path.length() == curr_path.length());
      if (b) {
        for (PathArcIt a1(_paths[i][j].path), a2(curr_path);
             b && a1 != INVALID && a2 != INVALID;
             ++a1, ++a2)
          b = _graph.id(a1) == _graph.id(a2);
        if (b)
          break;
      }
    }
    const bool path_seen = !(j == static_cast<int>(_paths[i].size()));

    // Set/modify the flow value for the found path
    Value         routed_flow = 0;
    CapacityValue min_cap;
    if (!path_seen) {
      // Add new path
      Value path_cost = 0;
      min_cap         = std::numeric_limits<CapacityValue>::infinity();
      for (PathArcIt a(curr_path); a != INVALID; ++a) {
        path_cost += (*_cost[i])[a];
        if (_cap[a] < min_cap)
          min_cap = _cap[a];
      }
      routed_flow = std::min(min_cap, d_i);
      // Save path
      _paths[i].push_back(
          PathData(curr_path, routed_flow, min_cap, path_cost, h));
    } else {
      // Update path flow
      min_cap     = _paths[i][j].min_cap;
      routed_flow = std::min(min_cap, d_i);
      _paths[i][j].value += routed_flow;
    }
    return routed_flow;
  }

  // Modify the arc lengths along the current path
  void updateLen(const Path& p, const Value& f, const int& i) {
    Value path_cost = 0;
    for (PathArcIt a(p); a != INVALID; ++a) {
      const Value& c_a = (*_cost[i])[a];
      _len[a] *= 1 + _epsilon * f / _cap[a] + c_a * _phi;
      path_cost += c_a;
    }
    // Update phi if needed
    if (_tol.positive(_phi))
      _phi *= (1 + (_epsilon * path_cost * f) / _cost_bound);
  }

  void updateDualObjective() {
    _dual_obj = 0;
    for (ArcIt a(_graph); a != INVALID; ++a)
      _dual_obj += _cap[a] * _len[a];
  }

  void setFinalCost() {
    _total_cost = 0;
    for (int i = 0; i < _num_commodities; ++i) {
      for (int j = 0; j < static_cast<int>((_paths[i].size())); ++j) {
        for (PathArcIt a(_paths[i][j].path); a != INVALID; ++a)
          _total_cost += (*_cost[i])[a] * (*_flow[i])[a];
      }
    }
  }

  // Scale flows in the \ref FLEISCHER_MAX_CONCURRENT.
  // Flows are scaled by the number of iterations.
  void scaleFinalFlowsMaxConcurrent() {
    int scaler = _iter;
    for (int i = 0; i < _num_commodities; ++i) {
      for (int j = 0; j < static_cast<int>(_paths[i].size()); ++j) {
        _paths[i][j].value /= scaler;
        for (PathArcIt a(_paths[i][j].path); a != INVALID; ++a) {
          _total_flow[a] += _paths[i][j].value;
        }
      }
    }
    _obj = (_iter) / scaler;
  }

  // Update total flow map
  void updateTotalFlow() {
    for (int i = 0; i < _num_commodities; ++i) {
      for (int j = 0; j < static_cast<int>(_paths[i].size()); ++j) {
        CapacityValue val = _paths[i][j].value;
        for (PathArcIt a(_paths[i][j].path); a != INVALID; ++a)
          _total_flow[a] += val;
      }
    }
  }

  // Scale flows per commodity for the FLEISCHER_MAX case.
  // scaler = max {total_flow_a / cap_a} such that all flows are capacity
  // feasible
  void scaleFinalFlowsMaxFlow() {
    updateTotalFlow();
    Value scaler = 0;
    for (ArcIt a(_graph); a != INVALID; ++a) {
      if (_total_flow[a] / _cap[a] > scaler)
        scaler = _total_flow[a] / _cap[a];
    }
    if (_tol.positive(scaler)) {
      for (ArcIt a(_graph); a != INVALID; ++a) {
        _total_flow[a] /= scaler;
      }
      for (int i = 0; i < _num_commodities; ++i) {
        for (int j = 0; j < static_cast<int>(_paths[i].size()); ++j)
          _paths[i][j].value /= scaler;
      }
    }
  }

  // Set flow decomposed solution for exposure
  void setFlowPerCommodity() {
    // Setting flow values using saved paths
    for (int i = 0; i < _num_commodities; ++i) {
      for (int j = 0; j < static_cast<int>(_paths[i].size()); ++j) {
        Value val = _paths[i][j].value;
        for (PathArcIt a(_paths[i][j].path); a != INVALID; ++a)
          (*_flow[i])[a] += val;
      }
    }
  }

  void saveFinalPaths() {
    for (int i = 0; i < _num_commodities; ++i) {
      const int& p_size     = static_cast<int>(_paths[i].size());
      bool       saved_path = false;
      Value      path_flow  = 0;
      Value      path_cost  = std::numeric_limits<Value>::infinity();
      Path       path;

      for (int j = 0; j < p_size; ++j) {
        const Value& v = _paths[i][j].value;
        const Value& c = _paths[i][j].cost;
        if (_tol.positive(v)) {
          if (saved_path)
            std::cout << "[WARNING] multiple paths with positive flow"
                      << std::endl;
          // Keep path with largest flow
          if (path_flow < v && c < path_cost) {
            path       = _paths[i][j].path;
            path_flow  = v;
            path_cost  = c;
            saved_path = true;
          }
        }
      }
      _final_paths[i] = path;
      _final_flows[i] = path_flow;
      _final_costs[i] = path_cost;
    }
  }

 public:
  const Value& flow(const int& k, const Arc& arc) const {
    return (*_flow[k])[arc];
  }

  const Value& flow(const int& k, const int& arc_id) const {
    return (*_flow[k])[_graph.arcFromId(arc_id)];
  }

  const Value& length(const Arc& arc) const { return _len[arc]; }

  const ValueMap& lengthMap() const { return _len; }

  const Value& getObjective() const { return _obj; }

  const Value& getDualObjective() const { return _dual_obj; }

  const Value& getTotalCost() const { return _total_cost; }

  const ValueMap& getFlowMap(const int& k) const { return *_flow[k]; }

  const Path& getPath(const int& k) const { return _final_paths[k]; }

  const Value& getPathCost(const int& k) const { return _final_costs[k]; }

  const Value& getPathFlow(const int& k) const { return _final_flows[k]; }
};  // namespace lemon
/// @}

}  // namespace lemon

#endif  // LEMON_MULTICOMMODITY_FLOW_H_
