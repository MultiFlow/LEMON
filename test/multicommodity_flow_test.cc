/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2003-2013
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

#include <iostream>
#include <lemon/concepts/digraph.h>
#include <lemon/concepts/maps.h>
#include <lemon/multicommodity_flow.h>
#include <lemon/smart_graph.h>

using namespace lemon;

void maxFlowTest() {
  DIGRAPH_TYPEDEFS(SmartDigraph);

  SmartDigraph g;
  SmartDigraph::ArcMap<double> cap(g);
  Node s = g.addNode();
  Node t = g.addNode();
  Node n1 = g.addNode();
  Node n2 = g.addNode();
  Node n3 = g.addNode();
  Arc a;
  a = g.addArc(s, n1);
  cap[a] = 20.0;
  a = g.addArc(s, n2);
  cap[a] = 10.0;
  a = g.addArc(s, n3);
  cap[a] = 10.0;
  a = g.addArc(n1, n2);
  cap[a] = 20.0;
  a = g.addArc(n2, n3);
  cap[a] = 20.0;
  a = g.addArc(n1, t);
  cap[a] = 20.0;
  a = g.addArc(n2, t);
  cap[a] = 20.0;
  a = g.addArc(n3, t);
  cap[a] = 20.0;

  MaxMulticommodityFlow<SmartDigraph> max(g, cap, 1, {s}, {t});
  max.run();
};

int main() {
  maxFlowTest();
  return 0;
}
