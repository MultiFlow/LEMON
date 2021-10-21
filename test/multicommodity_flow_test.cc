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

#include "test_tools.h"

const std::string lgf =
    "@nodes\n"
    "label\n"
    "0\n"
    "1\n"
    "2\n"
    "3\n"
    "@arcs\n"
    "    label  cap   cost \n"
    "0 1   0    10.0  1.0\n"
    "0 2   1    10.0  1.0\n"
    "1 3   2    10.0  1.0\n"
    "2 3   3    10.0  1000.0\n"
    "@attributes\n"
    "source_1 0\n"
    "target_1 3\n"
    "source_2 0\n"
    "target_2 2\n";

using namespace lemon;

void testLgf() {
  DIGRAPH_TYPEDEFS(SmartDigraph);

  Tolerance<double> tol;
  tol.epsilon(1e-5);

  SmartDigraph                 graph;
  Node                         s1, t1, s2, t2;
  SmartDigraph::ArcMap<double> cap(graph), cost(graph);

  std::istringstream input(lgf);
  DigraphReader<SmartDigraph>(graph, input)
      .arcMap("cap", cap)
      .arcMap("cost", cost)
      .node("source_1", s1)
      .node("target_1", t1)
      .node("source_2", s2)
      .node("target_2", t2)
      .run();

  int number_commodities = 2;
  // Max
  {
    ApproxMCF<SmartDigraph> max(
        graph, cap, cost, number_commodities, {s1, s2}, {t1, t2});
    max.run();
  }
  // Max Concurrent
  {
    ApproxMCF<SmartDigraph> max(
        graph, cap, cost, number_commodities, {s1, s2}, {t1, t2}, {10.0, 10.0});
    max.run(0.1, ApproxMCF<SmartDigraph>::FLEISCHER_MAX_CONCURRENT);
    check(max.getTotalCost() == 30.0, "max concurrent failed");
  }
  // Min cost
  {
    ApproxMCF<SmartDigraph> max(
        graph, cap, cost, number_commodities, {s1, s2}, {t1, t2}, {10.0, 10.0});
    max.run(0.1, ApproxMCF<SmartDigraph>::BINARY_SEARCH_MIN_COST);
    check(max.getTotalCost() == 30.0, "min-cost failed");
  }
  // Min cost (using constructor with cost vector)
  {
    std::vector<ApproxMCF<SmartDigraph>::Traits::LargeCostMap*> cost_vector;
    cost_vector.resize(number_commodities);
    for (int i = 0; i < number_commodities; ++i) {
      cost_vector[i] =
          ApproxMCF<SmartDigraph>::Traits::createLargeCostMap(graph);
      for (ArcIt a(graph); a != INVALID; ++a) {
        (*cost_vector[i])[a] = cost[a];
      }
    }
    ApproxMCF<SmartDigraph> max(
        graph,
        cap,
        cost_vector,
        number_commodities,
        {s1, s2},
        {t1, t2},
        {10.0, 10.0});
    max.run(0.1, ApproxMCF<SmartDigraph>::BINARY_SEARCH_MIN_COST);
    check(max.getTotalCost() == 30.0, "min-cost failed");
  }
  // Min cost (using default constructor and addCommodity)
  {
    ApproxMCF<SmartDigraph>                       max(graph, cap, cost);
    ApproxMCF<SmartDigraph>::Traits::LargeCostMap cost1(graph, 0.0);
    ApproxMCF<SmartDigraph>::Traits::LargeCostMap cost2(graph, 0.0);
    mapCopy(graph, cost, cost1);
    mapCopy(graph, cost, cost2);
    max.addCommodity(s1, t1, 10.0, cost1).addCommodity(s2, t2, 10.0, cost2);
    max.run(0.1, ApproxMCF<SmartDigraph>::BINARY_SEARCH_MIN_COST);
    check(max.getTotalCost() == 30.0, "min-cost failed");
  }
  // Max Concurrent - 3 commodities
  {
    ApproxMCF<SmartDigraph> max(
        graph, cap, cost, 3, {s1, s1, s1}, {t1, t2, t1}, {5.0, 10.0, 5.0});
    max.run(0.1, ApproxMCF<SmartDigraph>::FLEISCHER_MAX_CONCURRENT);
    std::cout << "max.getTotalCost() " << max.getTotalCost() << "\n";
    check(
        !tol.different(max.getTotalCost(), 30.0),
        "max concurrent with 3 commodities failed");
  }
};

void testDisconnected() {
  DIGRAPH_TYPEDEFS(SmartDigraph);

  SmartDigraph                 g;
  SmartDigraph::ArcMap<double> cap(g);
  SmartDigraph::ArcMap<double> cost(g);

  Node n0 = g.addNode();
  Node n1 = g.addNode();
  Node n2 = g.addNode();
  Node n3 = g.addNode();

  Arc a;
  a       = g.addArc(n0, n1);
  cap[a]  = 10.0;
  cost[a] = 1.0;

  a       = g.addArc(n2, n3);
  cap[a]  = 10.0;
  cost[a] = 1.0;

  // Max
  {
    ApproxMCF<SmartDigraph> max(g, cap, cost, 2, {n0, n2}, {n1, n3});
    max.run();
    for (int k = 0; k < 2; ++k) {
      std::cout << "Commodity " << k << "\n";
      for (ArcIt e(g); e != INVALID; ++e)
        std::cout << "flow on " << g.id(g.source(e)) << "->"
                  << g.id(g.target(e)) << " is = " << max.flow(k, e) << "\n";
    }
  }
  // Max Concurrent
  {
    ApproxMCF<SmartDigraph> max(
        g, cap, cost, 2, {n0, n2}, {n1, n3}, {10.0, 10.0});
    max.run(0.1, ApproxMCF<SmartDigraph>::FLEISCHER_MAX_CONCURRENT);
    for (int k = 0; k < 2; ++k) {
      std::cout << "Commodity " << k << "\n";
      for (ArcIt e(g); e != INVALID; ++e)
        std::cout << "flow on " << g.id(g.source(e)) << "->"
                  << g.id(g.target(e)) << " is = " << max.flow(k, e) << "\n";
    }
  }
};

int main() {
  testLgf();
  testDisconnected();
};
