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

#include <lemon/concepts/digraph.h>
#include <lemon/concepts/maps.h>
#include <lemon/lgf_reader.h>
#include <lemon/multicommodity_flow.h>
#include <lemon/smart_graph.h>

#include <iostream>
#include <sstream>

#include "test/test_tools.h"

char test_lgf[] =
    "@nodes\n"
    "label\n"
    "0\n"
    "1\n"
    "2\n"
    "3\n"
    "@arcs\n"
    "    label  cap   cost  solution_1  solution_2 solution_1_max  "
    "solution_2_max\n"
    "0 1   0    10.0  1.0      10.0        0.0         10.0    0.0\n"
    "0 2   1    10.0  1.0       0.0        10.0        8.3     1.7\n"
    "1 3   2    10.0  1.0      10.0        0.0         10.0    0.0\n"
    "2 3   3    10.0  1000.0   0.0         0.0         8.3     0.0\n"
    "@attributes\n"
    "source_1 0\n"
    "target_1 3\n"
    "source_2 0\n"
    "target_2 2\n";

using namespace lemon;

// Returns string for test-naming
template<typename GR, typename A = typename GR::Arc>
std::string getEdgeName(const GR& graph, const A& edge) {
  return std::to_string(graph.id(graph.source(edge))) + "->" +
         std::to_string(graph.id(graph.target(edge)));
}

template<class Digraph>
void checkApproxMCF() {
  TEMPLATE_DIGRAPH_TYPEDEFS(Digraph);

  Digraph                                           graph;
  Node                                              s1, t1, s2, t2;
  typedef typename Digraph::template ArcMap<double> DoubleMap;
  DoubleMap cap(graph), cost(graph), sol1(graph), sol2(graph), sol1_max(graph),
      sol2_max(graph);
  Tolerance<double> tol;
  tol.epsilon(1e-1);  // 1 decimal place

  std::istringstream input(test_lgf);
  DigraphReader<Digraph>(graph, input)
      .arcMap("cap", cap)
      .arcMap("cost", cost)
      .arcMap("solution_1", sol1)  // flow in the solution for commodity 1
      .arcMap("solution_2", sol2)  // flow in the solution for commodity 1
      .arcMap(
          "solution_1_max",
          sol1_max)  // flow in the solution for commodity (max-flow) 1
      .arcMap(
          "solution_2_max",
          sol2_max)          // flow in the solution for commodity (max-flow) 1
      .node("source_1", s1)  // source node commodity 1
      .node("target_1", t1)  // target node commodity 1
      .node("source_2", s2)  // source node commodity 2
      .node("target_2", t2)  // target node commodity 2
      .run();

  // Max
  {
    ApproxMCF<Digraph> max(graph, cap);
    max.addCommodity(s1, t1).addCommodity(s2, t2).run();
    for (ArcIt e(graph); e != INVALID; ++e) {
      const std::string test_name =
          "[FLEISCHER_MAX] flow on " + getEdgeName<Digraph>(graph, e);
      // Check flow for commodity 1 (index 0) to 1 decimal place
      check(
          !tol.different(max.flow(0, e), sol1_max[e]),
          test_name + " commodity 1");
      // Check flow for commodity 2 (index 1) to 1 decimal place
      check(
          !tol.different(max.flow(1, e), sol2_max[e]),
          test_name + " commodity 2");
    }
  }

  // Max Concurrent
  {
    ApproxMCF<Digraph> max(graph, cap);
    max.addCommodity(s1, t1, 10.0)
        .addCommodity(s2, t2, 10.0)
        .run(ApproxMCF<Digraph>::FLEISCHER_MAX_CONCURRENT);
    for (ArcIt e(graph); e != INVALID; ++e) {
      const std::string test_name = "[FLEISCHER_MAX_CONCURRENT] flow on " +
                                    getEdgeName<Digraph>(graph, e);
      // Check commodity 1 (index 0)
      check(max.flow(0, e) == sol1[e], test_name + " commodity 1");
      // Check commodity 2 (index 1)
      check(max.flow(1, e) == sol2[e], test_name + " commodity 2");
    }
  }

  // Min cost
  {
    ApproxMCF<Digraph> max(graph, cap);
    max.addCommodity(s1, t1, 10.0, cost)
        .addCommodity(s2, t2, 10.0, cost)
        .run(ApproxMCF<Digraph>::BINARY_SEARCH_MIN_COST);
    for (ArcIt e(graph); e != INVALID; ++e) {
      const std::string test_name =
          "[BINARY_SEARCH_MIN_COST] flow on " + getEdgeName<Digraph>(graph, e);
      // Check commodity 1 (index 0)
      check(max.flow(0, e) == sol1[e], test_name + " commodity 1");
      // Check commodity 2 (index 1)
      check(max.flow(1, e) == sol2[e], test_name + " commodity 2");
    }
    // Check total cost
    check(max.getTotalCost() == 30.0, "min-cost failed");
  }

  // Max Concurrent - 3 commodities
  {
    ApproxMCF<Digraph> max(graph, cap);
    max.addCommodity(s1, t1, 5.0, cost)
        .addCommodity(s1, t2, 10.0, cost)
        .addCommodity(s1, t1, 5.0, cost)
        .run(ApproxMCF<Digraph>::FLEISCHER_MAX_CONCURRENT);

    const DoubleMap& flow = max.getFlowMap(0);
    check(
        max.getTotalCost() == 30.0,
        "3 commodities FLEISCHER_MAX_CONCURRENT failed");
  }
}

template<class Digraph>
void checkApproxMCFDisconnected() {
  TEMPLATE_DIGRAPH_TYPEDEFS(Digraph);

  Digraph                                           graph;
  typedef typename Digraph::template ArcMap<double> DoubleMap;
  DoubleMap cap(graph), cost(graph), sol1(graph), sol2(graph);

  Node n0 = graph.addNode();
  Node n1 = graph.addNode();
  Node n2 = graph.addNode();
  Node n3 = graph.addNode();

  Arc a;
  a       = graph.addArc(n0, n1);
  cap[a]  = 10.0;
  cost[a] = 1.0;
  sol1[a] = 10.0;
  sol2[a] = 0.0;

  a       = graph.addArc(n2, n3);
  cap[a]  = 10.0;
  cost[a] = 1.0;
  sol1[a] = 0.0;
  sol2[a] = 10.0;

  // Max
  {
    ApproxMCF<Digraph> max(graph, cap);
    max.addCommodity(n0, n1).addCommodity(n2, n3).run();
    for (ArcIt e(graph); e != INVALID; ++e) {
      const std::string test_name =
          "[FLEISCHER_MAX] flow on " + getEdgeName<Digraph>(graph, e);
      // Check commodity 1
      check(max.flow(0, e) == sol1[e], test_name + " commodity 1");
      // Check commodity 2
      check(max.flow(1, e) == sol2[e], test_name + " commodity 2");
    }
  }
  // Max Concurrent
  {
    ApproxMCF<Digraph> max(graph, cap);
    max.addCommodity(n0, n1, 10.0, cost)
        .addCommodity(n2, n3, 10.0, cost)
        .run(ApproxMCF<Digraph>::FLEISCHER_MAX_CONCURRENT);
    for (ArcIt e(graph); e != INVALID; ++e) {
      const std::string test_name = "[FLEISCHER_MAX_CONCURRENT] flow on " +
                                    getEdgeName<Digraph>(graph, e);
      // Check commodity 1
      check(max.flow(0, e) == sol1[e], test_name + " commodity 1");
      // Check commodity 2
      check(max.flow(1, e) == sol2[e], test_name + " commodity 2");
    }
  }
}

int main() {
  checkApproxMCF<SmartDigraph>();
  checkApproxMCF<ListDigraph>();
  checkApproxMCFDisconnected<SmartDigraph>();
  checkApproxMCFDisconnected<ListDigraph>();
}
