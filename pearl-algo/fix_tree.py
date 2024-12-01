import sys
sys.path.append('/mnt/d/Egyetem/Routing Cikk/fast-failover')
sys.path.append('/mnt/d/Egyetem/Routing Cikk/fast-failover/pearl_algo')

import networkx as nx
import argparse

from graph_analise import find_pearls, neighbor_of_subgraph
from graph_analise import cut_cutting_edges, is_3_connected, read_in_graph
from graph_analise import *

from util import GreedyArborescenceDecomposition
from util import contract_paths_keep_root, leaf_pearls, totally_contract_pearl
from util import fixing_pearl_parent_pointers, separate_by_cutting_nodes
from util import create_big_2_conn_from_Ion, create_simple_2_pearl_multigraph
from util import reset_arb_attribute, get_arborescence_list, edge_connectivity_maxflow

from heapq import heappush, heappop
import random

import networkx as nx
import pandas as pd
from pearl_routing import *


#'''
# compute the k^th arborescence of g greedily
def FindTree(g, k):
    print(f'Ezt használom a {k} fenyő keresésénél')
    T = nx.MultiDiGraph()
    T.add_node(g.graph['root'])
    R = {g.graph['root']}
    dist = dict()
    dist[g.graph['root']] = 0

    # heap of all border edges in form [(edge metric, (e[0], e[1]),...]
    h = []
    # heap of all border edges in form [(edge metric, (e[0], e[1], e[k])),...]
    in_edges = sorted(g.in_edges(
        g.graph['root'], keys=True), key=lambda k: random.random())
    for e in in_edges:
        heappush(h, (0, e))
        if k > 1:
            continue

    while len(h) > 0:
        (d, e) = heappop(h)           # d valami távolság, e az él, itt ugye e[0] az új csúcs e[1] már kész
        g.remove_edge(*e)             # g-ből kivesszük az élet
        if (e[0] not in R) and (k == 1 or edge_connectivity_maxflow(g, e[0], g.graph['root']) >= k-1):   #ha ez az él jó nekünk, berakjuk a fába
            dist[e[0]] = d+1               # ez a csúcs eggyel messzebb van, mint az e[1]
            R.add(e[0])                    # a kész csúcsokhoz hozzáadjuk, amit hozzá kell adni
            in_edges = sorted(g.in_edges(e[0], keys=True), key=lambda k: random.random())
            for e in in_edges:
                heappush(h, (d+1, e))
                if k > 1:
                    continue
            T.add_edge(*e)
        else:
            g.add_edge(*e)               # ha nem jó nekünk, visszarakjuk a gráfba

    if len(R) < len(g.nodes()):
        print(
            f"Couldn't find next edge for tree with {g.graph['root']}, ", k, len(R))
        sys.stdout.flush()

    print(T)
    return T
'''

# compute the k^th arborescence of g greedily
def FindTree(g, k):
    T = nx.MultiDiGraph()
    T.add_node(g.graph['root'])
    R = {g.graph['root']}
    dist = dict()
    dist[g.graph['root']] = 0
    # heap of all border edges in form [(edge metric, (e[0], e[1])),...]
    h = []
    preds = sorted(g.predecessors(
        g.graph['root']), key=lambda k: random.random())
    for x in preds:
        heappush(h, (0, (x, g.graph['root'])))
        if k > 1:
            continue
    while len(h) > 0:
        (d, e) = heappop(h)
        g.remove_edge(*e)
        if e[0] not in R and (k == 1 or nx.edge_connectivity(g, e[0], g.graph['root']) >= k-1):
            dist[e[0]] = d+1
            R.add(e[0])
            preds = sorted(g.predecessors(e[0]), key=lambda k: random.random())
            for x in preds:
                if x not in R:
                    heappush(h, (d+1, (x, e[0])))
            T.add_edge(*e)
        else:
            g.add_edge(*e)
    if len(R) < len(g.nodes()):
        print(
            "Couldn't find next edge for tree with g.graph['root'], ", k, len(R))
        sys.stdout.flush()
    return T
'''


# associate a greedy arborescence decomposition with g
def GreedyArborescenceDecomposition(g):
    reset_arb_attribute(g)
    gg = g.to_directed()
    K = g.graph['k']
    k = K
    while k > 0:
        T = FindTree(gg, k)
        if T is None:
            return None
        for (u, v, key) in T.edges(keys=True):
            g[u][v][key]['arb'] = K-k
        gg.remove_edges_from(T.edges(keys=True))
        k = k-1
    return get_arborescence_list(g)


def main():
    graph = nx.read_graphml('/mnt/d/Egyetem/Routing Cikk/fast-failover/pearl-algo/graph_sets/example_graphs/Colt_pearl.graphml')
    graph.graph['k'] = 3
    print(is_3_connected(graph), graph.graph['root'])
    multidigraph = graph.to_directed()
    print(multidigraph)

    #find the arborescences
    unordered_arborescence_list = GreedyArborescenceDecomposition(multidigraph)
    arborescence_list = []

    for arb in unordered_arborescence_list:
        #print(arb)
        pass



if __name__ == '__main__':
    main()