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
 */

// Author: David Torres Sanchez <d.torressanchez@lancaster.ac.uk>
// Lancaster University (2021)

#ifndef LEMON_MULTICOMMODITY_FLOW_H_
#define LEMON_MULTICOMMODITY_FLOW_H_

/// \ingroup approx
///
/// \file
/// \brief Fleischer's algorithms for finding approximate multicommodity flows.

#include <lemon/adaptors.h>
#include <lemon/bfs.h>
#include <lemon/concepts/maps.h>
#include <lemon/dijkstra.h>
#include <lemon/maps.h>
#include <lemon/math.h>
#include <lemon/path.h>
#include <lemon/preflow.h>
#include <lemon/static_graph.h>
#include <lemon/tolerance.h>

#include <algorithm>  // min
#include <limits>     // numeric_limits
#include <random>     // random_device, mt19937
#include <set>
#include <vector>

namespace lemon {

/// \brief Default traits class of \ref ApproxMCF class.
///
/// Default traits class of \ref ApproxMCF class.
/// \tparam GR The type of the digraph.
/// \tparam CM The type of the capacity map.
/// It must conform to the \ref concepts::ReadMap "ReadMap" concept.
#ifdef DOXYGEN
template<typename GR, typename CAP, typename LC>
#else
template<
    typename GR,
    typename CAP = typename GR::template ArcMap<double>,
    typename LC  = double>
#endif
struct ApproxMCNFDefaultTraits {
  /// The type of the digraph
  typedef GR Digraph;
  /// The type of the capacity map
  typedef CAP CapacityMap;
  /// The type of the capacity map values. By default is double.
  typedef typename CapacityMap::Value Value;

  /// \brief The large cost type used for internal computations
  ///
  /// The large cost type used for internal computations.
  /// \c Value must be convertible to \c LargeCost.
  /// By default is \c double.
  typedef LC                                           LargeCost;
  typedef typename Digraph::template ArcMap<LargeCost> LargeCostMap;

  /// The tolerance type used for internal computations
  typedef lemon::Tolerance<LargeCost> Tolerance;

  /// \brief The path type of the decomposed flows
  ///
  /// It must conform to the \ref lemon::concepts::Path "Path" concept
  typedef lemon::Path<Digraph>                 Path;
  typedef typename lemon::Path<Digraph>::ArcIt PathArcIt;

  /// \brief Instantiates a LargeCostMap.
  ///
  /// This function instantiates a \ref LargeCostMap.
  /// \param digraph The digraph for which we would like to define
  /// the large cost map. This can be used for flows or otherwise.
  static LargeCostMap* createLargeCostMapPtr(
      const Digraph&   digraph,
      const LargeCost& value = 0) {
    return new LargeCostMap(digraph, value);
  }
  // Same as above just const
  static const LargeCostMap* createLargeCostMapConstPtr(
      const Digraph&   digraph,
      const LargeCost& value = 0) {
    return new LargeCostMap(digraph, value);
  }
};

/// \addtogroup approx
/// @{
/// \brief .
///
/// Default traits class of \ref ApproxMCF class.
/// \tparam GR The type of the digraph.
/// \tparam CM The type of the capacity map.
/// It must conform to the \ref concepts::ReadMap "ReadMap" concept.
#ifdef DOXYGEN
template<typename GR, typename CM, typename TR>
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
  typedef typename TR::Value Value;

  /// \brief The large cost type
  ///
  /// The large cost type used for internal computations.
  /// By default is \c double.
  typedef typename TR::LargeCost LargeCost;
  /// Arc Map with values of type \ref LargeCost
  typedef typename TR::LargeCostMap LargeCostMap;

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
    Path      path;
    LargeCost value;
    Value     min_cap;
    LargeCost cost;
    int       hash;

    PathData() {}
    PathData(const Path& p, LargeCost v, Value m, LargeCost c, int h = 0)
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
  // Vector indexed by commodity with \ref LargeCostMap pointers.
  typedef std::vector<LargeCostMap*>       LargeCostMapVector;
  typedef std::vector<const LargeCostMap*> ConstLargeCostMapVector;

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
  std::vector<Value> _demand;
  // The cost for each commodity (may be 0, or undefined)
  // See \ref LargeCostMapVector.
  std::vector<const LargeCostMap*> _cost;

  /**
   * Algorithm parameters
   */
  // type of approximation algorithm. See \ref ApproxType.
  ApproxType _approx_type;
  // Stopping criterion
  LargeCost _epsilon;
  // Whether the internal heuristics will be applied at the end of the
  // algorithms
  bool _heuristic;
  // Whether there is feasible unsplitable flow
  bool _feasible_unsplit;
  // Objective function of the primal problem
  LargeCost _obj;
  // Objective function of the dual problem
  LargeCost _dual_obj;
  bool      _local_cost;
  //
  LargeCost _delta;
  // Expected maximum number of iterations
  LargeCost _max_iter;
  //
  LargeCost _limit;
  // Number of iterations
  int _iter;
  // Tolerance function. See \ref lemon::Tolerance.
  Tolerance _tol;

  // For the BINARY_SEARCH_MIN_COST case
  // Cost of the final solution.
  LargeCost _total_cost;
  // Cost bounding factor.
  LargeCost _phi;
  // Cost bound.
  LargeCost _cost_bound;

  // The aggregated flow per arc over all commodities
  LargeCostMap _total_flow;
  // The value of the length per arc (the dual solution).
  LargeCostMap _len;
  // The flow for each commodity. See \ref LargeCostMapVector.
  LargeCostMapVector _flow;
  // The paths found for different commodities. See \ref PathMatrix.
  PathMatrix _paths;

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

  /// \brief Add a new commodity.
  ///
  /// This function adds a new commodity with the given source and
  /// sink nodes.
  ApproxMCF& addCommodity(
      const Node&         s,
      const Node&         t,
      const Value&        d,
      const LargeCostMap& cost_map) {
    _source.push_back(s);
    _target.push_back(t);
    if (d > 0)
      _demand.push_back(d);
    _local_cost = false;
    _cost.push_back(&cost_map);
    ++_num_commodities;
    return *this;
  }

  // @overload no cost
  ApproxMCF& addCommodity(const Node& s, const Node& t, const Value& d = 0) {
    const LargeCost* cost_map_ptr = Traits::createLargeCostMapConstPtr(_graph);
    addCommodity(s, t, d, *cost_map_ptr);
    _local_cost = true;
    return *this;
  }

  ApproxMCF& addCommodities(
      const int&                        k,
      const std::vector<Node>&          s,
      const std::vector<Node>&          t,
      const std::vector<Value>&         d = {},
      const std::vector<LargeCostMap*>& c = {}) {
    _num_commodities = k;
    _source          = s;
    _target          = t;

    if (!d.empty())
      _demand = d;

    _cost.resize(_num_commodities);
    for (int i = 0; i < _num_commodities; ++i) {
      if (c.empty()) {
        _local_cost = true;
        _cost[i]    = Traits::createLargeCostMapConstPtr(_graph);
      } else {
        _cost[i] = c[i];
      }
    }
    return *this;
  }

  /// \name Execution control

  /// @{

  /// \brief Run the algorithm for finding a maximum multicommodity flow.
  ///
  /// This function runs the algorithm for finding a maximum multicommodity
  /// flow.
  /// \param epsilon The approximation factor. It must be greater than 1.
  /// \return The effective approximation factor (at least 1), i.e. the
  /// ratio of the values of the found dual and primal solutions.
  /// Therefore both the primal and dual solutions are r-approximate.
  ///
  /// \note The smaller parameter \c r is, the slower the algorithm runs,
  /// but it might find better solution. According to our tests, using
  /// the default value the algorithm performs well and it usually manage
  /// to find the exact optimum due to heuristic improvements.
  void run(
      const ApproxType& at        = FLEISCHER_MAX,
      const bool&       heuristic = true,
      const LargeCost&  epsilon   = 0.1) {
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

 private:
  // Initialize flow values and lengths
  void init() {
    // Restart
    _iter     = 0;
    _obj      = 0;
    _dual_obj = 0;
    // Allocate paths
    _paths.resize(_num_commodities);
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
        _flow[i] = Traits::createLargeCostMapPtr(_graph);
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

  // Checks whether:
  // 1. Commodities are well-defined
  // 2. Input sizes are consistent
  // 3. Concurrent case: Feasible split flow exists for each commodity that
  // satisfies demand. In the case this doesn't pass, will throw an exception.
  //
  // TODO(torressa): Group commodities by source!
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
          Value total_cap_source = 0;
          Value max_cap_source   = 0;
          for (OutArcIt a(_graph, _source[i]); a != INVALID; ++a) {
            total_cap_source += _cap[a];
            if (max_cap_source < _cap[a])
              max_cap_source = _cap[a];
          }

          // Total capacity incoming at the target
          Value total_cap_target = 0;
          Value max_cap_target   = 0;
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

    LargeCost lower_bound = 0, upper_bound = getTotalCost();

    // Best solution attributes
    LargeCost curr_best_cost = std::numeric_limits<LargeCost>::infinity();
    LargeCostMapVector curr_best_flow;

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
    Path      curr_path;
    LargeCost curr_dist;
    while (true) {
      // Find the shortest path for the i-th commodity
      // (w.r.t the current length function)
      for (int i = 0; i < _num_commodities; ++i) {
        // Route demand
        Value d_i = _demand[i];
        while (true) {  // while (phase)
          dijkstra(_graph, _len)
              .path(curr_path)
              .dist(curr_dist)
              .run(_source[i], _target[i]);
          // Get
          const Value& routed_flow = updatePaths(curr_path, i, d_i);
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
  }

  void runFleischerMax() {
    // Successive improvement along shortest paths
    int       i = 0;
    Path      curr_path;
    LargeCost curr_dist;
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
      const Value& routed_flow = updatePaths(curr_path, i);
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
      const int&             p_size       = static_cast<int>(_paths[i].size());
      LargeCost              total_flow_i = 0;
      std::vector<LargeCost> values(p_size, 0);

      for (int j = 0; j < p_size; ++j) {
        const LargeCost& val  = _paths[i][j].value;
        const LargeCost& cost = _paths[i][j].cost;
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
        Value d = std::numeric_limits<Value>::infinity();
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
  // Attempts to improve the solution along one of the found paths if it is
  // possible
  void runHeurMax2() {
    BoolArcMap filter(_graph);
    for (ArcIt a(_graph); a != INVALID; ++a)
      filter[a] = _tol.positive(_cap[a] - _total_flow[a]);
    FilterArcs<const Digraph, BoolArcMap> sub_graph(_graph, filter);
    Path                                  p;
    for (int i = 0; i < _num_commodities; ++i) {
      while (bfs(sub_graph).path(p).run(_source[i], _target[i])) {
        Value d = std::numeric_limits<Value>::infinity();
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
      const Path&  curr_path,
      const int&   i,
      const Value& d_i = std::numeric_limits<Value>::infinity()) {
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
    Value routed_flow = 0, min_cap;
    if (path_seen) {
      // Reuse info
      min_cap     = _paths[i][j].min_cap;
      routed_flow = std::min(min_cap, d_i);
      _paths[i][j].value += routed_flow;
    } else {
      // Extract relevant info
      Value path_cost = 0;
      min_cap         = std::numeric_limits<Value>::infinity();
      // Gamma = min capacity of all arcs in the path
      for (PathArcIt a(curr_path); a != INVALID; ++a) {
        path_cost += (*_cost[i])[a];
        if (_cap[a] < min_cap)
          min_cap = _cap[a];
      }
      routed_flow = std::min(min_cap, d_i);
      // Save path
      _paths[i].push_back(
          PathData(curr_path, routed_flow, min_cap, path_cost * min_cap, h));
    }
    return routed_flow;
  }

  // Modify the arc lengths along the current path
  void updateLen(const Path& p, const Value& f, const int& i) {
    LargeCost path_cost = 0;
    for (PathArcIt a(p); a != INVALID; ++a) {
      _len[a] *= 1 + _epsilon * f / _cap[a] + (*_cost[i])[a] * _phi;
      path_cost += (*_cost[i])[a];
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
        Value val = _paths[i][j].value;
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
    LargeCost scaler = 0;
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
        LargeCost val = _paths[i][j].value;
        for (PathArcIt a(_paths[i][j].path); a != INVALID; ++a)
          (*_flow[i])[a] += val;
      }
    }
  }

 public:
  const LargeCost& flow(const int& k, const Arc& arc) const {
    return (*_flow[k])[arc];
  }

  const LargeCost& flow(const int& k, const int& arc_id) const {
    return (*_flow[k])[_graph.arcFromId(arc_id)];
  }

  const LargeCost& length(const Arc& arc) const { return _len[arc]; }

  const LargeCostMap& lengthMap() const { return _len; }

  const LargeCost& getObjective() const { return _obj; }

  const LargeCost& getDualObjective() const { return _dual_obj; }

  const LargeCost& getTotalCost() const { return _total_cost; }
};  // namespace lemon
/// @}

}  // namespace lemon

#endif  // LEMON_MULTICOMMODITY_FLOW_H_
