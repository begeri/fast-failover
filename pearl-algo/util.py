import sys
sys.path.append('/mnt/d/Egyetem/Routing Cikk/fast-failover')

import networkx as nx
import argparse
from graph_analise import find_pearls, neighbor_of_subgraph
from graph_analise import cut_cutting_edges, is_3_connected, read_in_graph
from graph_analise import *
from graph_gen import add_ear

from heapq import heappush, heappop
import random

import networkx as nx
import pandas as pd

random.seed(42)


### ARBORESCENCE UTILS

# reset the arb attribute for all edges to -1, i.e., no arborescence assigned yet
def reset_arb_attribute(g):
    for (u, v, key) in g.edges(keys= True):
        g[u][v][key]['arb'] = -1

# given a graph return the arborescences in a dictionary with indices as keys
def get_arborescence_dict(g):
    arbs = {}
    for (u, v, key) in g.edges(keys=True):
        index = g[u][v][key]['arb']
        if index not in arbs:
            arbs[index] = nx.MultiDiGraph()
            arbs[index].graph['root'] = g.graph['root']
            arbs[index].graph['index'] = index
        arbs[index].add_edge(u, v, key)
    return arbs

# given a graph return a list of its arborescences
def get_arborescence_list(g):
    arbs = get_arborescence_dict(g)
    sorted_indices = sorted([i for i in arbs.keys() if i >= 0])
    return [arbs[i] for i in sorted_indices]

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
        #print(
        #    "Couldn't find next edge for tree with g.graph['root'], ", k, len(R))
        sys.stdout.flush()
    return T

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



### PEARL UITLS

#contract each path into an edge, if 'root' is not on the path
def contract_paths_keep_root(G):
    '''
    Nomen est omen, no 2 degree nodes
    Works only with MultiGraph

    If graph has a root property, then keeps it, if it is degree 2, adds an edge
    '''
    # list of all degree 2 nodes
    nodes_to_remove = [v for v in G.nodes if G.degree(v) == 2]

    #if one of the degree 2 nodes is the root, don't remove it
    if 'root' in G.graph.keys(): 
        if G.graph['root'] in nodes_to_remove:
            nodes_to_remove.remove(G.graph['root'])
            
            #neighbors = list(G.neighbors(G.graph['root']))
            #for neighbor in neighbors:
            #    G.add_edge(neighbor, G.graph['root'])


            
    # remove the ones we want to remove
    for v in nodes_to_remove:
        neighbors = list(G.neighbors(v))
        G.remove_node(v)
        if len(neighbors)==2:
            G.add_edge(neighbors[0], neighbors[1])
        # preserve connectivity 
        else: G.add_edge(neighbors[0], neighbors[0])

    return G

#modify the output of find pearls
def leaf_pearls(multi_graph, level='unknown'):
    '''
    Creates a list of pearls, that are the current leafs of a pearl arborescence.
    Returns a list of dictionaries, that contains the information of the pearls.
    PEARL[i]['nodes'] - list of nodes for pearl i
    PEARL[i]['cut_edges'] - the two edges defining pearl i
    PEARL[i]['boundary] - the two nodes being the boundary of pearl i
    PEARL[i]['neighbor_nodes'] = neighbors
    '''

    pearls = find_pearls(multi_graph)
    #this is the output list of dictionaries
    decorated_pearls = []
    for pearl in pearls:
        #setting up the dict
        dict = {}
        #the nodes in a pearl
        dict['nodes'] = pearl

        #boundary
        comp = list(pearl)
        neighbors = list(neighbor_of_subgraph(multi_graph, comp))
        side_of_pearl = []
        for node in comp:
            if node in neighbor_of_subgraph(multi_graph, node_set=neighbors):
                side_of_pearl.append(node)
        dict['boundary'] = side_of_pearl

        #list the edges
        cut_edges = []
        for x in side_of_pearl:
            for v in neighbors:
                for edge in multi_graph.subgraph([x,v]).edges:
                    cut_edges.append(edge)
        dict['cut_edges'] = cut_edges
        dict['neighbor_nodes'] = neighbors
        dict['level_of_pearl'] = level
        dict['parent_pearl'] = 'unknown'
        dict['id'] = 'unknown'
        
        decorated_pearls.append(dict)

    return decorated_pearls

# given a graph and a pearl in it, delete the nodes of the pearl and 
# connect the neighbor nodes, if there are multiple
def totally_contract_pearl(multi_graph, pearl):
    '''
    nomen est omen
    multi_graph - the graph we are working on
    comp - the nodes of the pearl to contract

    We replace the pearl with one edge between the neighbors. 
    If there is only one neighbor, we just delete the pearl.
    '''
    comp = pearl['nodes']
    #contract pearls
    neighbors = list(neighbor_of_subgraph(multi_graph, comp))
    
    multi_graph.remove_nodes_from(comp)
    if len(neighbors)==2:
        key = multi_graph.add_edge(neighbors[0],neighbors[1])
        pearl['substitute_edge'] = (neighbors[0],neighbors[1], key)
 


    return multi_graph

# sets pearl parint pointers
def fixing_pearl_parent_pointers(simple_pearl_list, graph):
    '''
    Given a 2-connected graph and its pearls, where graph has already pearls and pearls have already 
    levels on them find the parents of each non-root pearl
    INPUT:
    graph - nx.MultiGraph
    simple_pearl_list - list of dictionaries, dictionaries being the pearls

    OUTPUT:
    simple_pearl_list - updated pearl['parent_pearl'] datas
    '''

    
    #setting the parents well
    for pearl in reversed(simple_pearl_list):
        if pearl['parent_pearl']!='this is the root':
            pearl['parent_pearl'] = 'unknown'
        
    #setting the parents well
    for pearl in reversed(simple_pearl_list):
        if pearl['neighbor_nodes'] != 'this is the root':
            neighbor_pearl_ids = []
            for neighbor in pearl['neighbor_nodes']:
                neighbor_pearl_ids.append(graph.nodes[neighbor]['pearl_id'])
            # f there is only one neighbour, that is gonna be the parent pearl
            if len(neighbor_pearl_ids)==1:
                pearl['parent_pearl'] = neighbor_pearl_ids[0]
            elif neighbor_pearl_ids[0]==neighbor_pearl_ids[1]:    
                #if the 2 neighbors are the same, that is the parent
                pearl['parent_pearl'] = neighbor_pearl_ids[0]
            # if the two neighbors are not the same
            elif neighbor_pearl_ids[0]!=neighbor_pearl_ids[1]:
                # checking if one of them is the parent
                # this can happen only, if the yet bad level (pearl height) is bigger
                if simple_pearl_list[neighbor_pearl_ids[0]]['level_of_pearl']>pearl['level_of_pearl']:
                    pearl['parent_pearl'] = neighbor_pearl_ids[0]
                elif simple_pearl_list[neighbor_pearl_ids[1]]['level_of_pearl']>pearl['level_of_pearl']:
                    pearl['parent_pearl'] = neighbor_pearl_ids[1]
                else:   
                    # finding the parent:
                    neighbor_1_parent_id = simple_pearl_list[neighbor_pearl_ids[0]]['parent_pearl']
                    neighbor_2_parent_id = simple_pearl_list[neighbor_pearl_ids[1]]['parent_pearl']
                    # if one of the neighbors is the pearl, that is the parent
                    if   neighbor_1_parent_id == 'this is the root': 
                        pearl['parent_pearl'] = neighbor_pearl_ids[0]
                    elif neighbor_2_parent_id == 'this is the root':
                        pearl['parent_pearl'] = neighbor_pearl_ids[1]
                    # parents can be 'unknown', if we haven't set them yet.
                    else:
                        if neighbor_1_parent_id != 'unknown':
                            if neighbor_2_parent_id != 'unknown':
                                # in theory this 'if' is unnecessary, the pearls are in such order, that one of the neighbors
                                if neighbor_1_parent_id == neighbor_2_parent_id:
                                    print(simple_pearl_list[neighbor_1_parent_id]['level_of_pearl'])
                                    # in theory this if is unnecessary, if the two neighbors are different, either one of them is 
                                    # the parent of all 3 of them have the same parent
                                    if simple_pearl_list[neighbor_1_parent_id]['level_of_pearl']>pearl['level_of_pearl']:
                                        pearl['parent_pearl'] = neighbor_1_parent_id
                                else: 
                                    print('Unexpected behaviour_1')

                            
                            else:   #neighbor_2_parent is unknown
                                if simple_pearl_list[neighbor_1_parent_id]['level_of_pearl']>pearl['level_of_pearl']:                                        
                                    pearl['parent_pearl'] = neighbor_1_parent_id
                        else: 
                            print('Unexpected behaviour_2')

    return simple_pearl_list

# find the cutting nodes
def cutting_nodes(G):
    if len(list(nx.connected_components(G)))!=1:
        print('not connected graph')
        return []
    
    cuts = []
    for node in G.nodes:
        test_graph = G.copy()
        test_graph.remove_node(node)
        if len(list(nx.connected_components(test_graph)))>1:
            cuts.append(node)
    return cuts

# create the 2-node-connected subgraphs of a graph 
def separate_by_cutting_nodes(G):
    '''
    A recursive algorithm to compute the 2-node-connected subgraphs of a graph

    INPUT:
    G is a nx.MultiGraph with a G.graph['root'] property

    OUTPUT:
    a list of 2-node-connected graphs, that are a partition of E(G), each graph is a nx.MultiGraph
    '''
    root = G.graph['root']
    cuts = list(cutting_nodes(G))
    

    list_of_2_node_conn_comps = []

    # if we are 2-node-connected, return the whole graph 
    if cuts==[]:
        return [G]

    #if we are not, recursion by one of the cutting nodes
    curr_cutting_node = cuts[0]
    graph_pieces = G.copy()
    graph_pieces.remove_node(curr_cutting_node)
    comps = nx.connected_components(graph_pieces)

    list(G.neighbors(curr_cutting_node))
    for comp in comps:
        smaller_graph = G.subgraph([curr_cutting_node, *comp]).copy()
        #if root is not in the smaller comp
        if G.graph['root'] not in comp:        # not comp but smaller graph
            lista = list(comp)
            if type(G.graph['root']) == 'str':
                lista.sort(key = lambda node: len(nx.shortest_path(G, str(source=G.graph['root']), target=str(node))))
            else: 
                lista.sort(key = lambda node: len(nx.shortest_path(G, source=G.graph['root'], target=node)))
            closest_node = lista[0]
            smaller_graph.graph['root'] = closest_node
        
        # recursively compute all the 2-node-connected components
        list_of_2_conn_comps_of_this_comp = separate_by_cutting_nodes(smaller_graph)
        for result_comp in list_of_2_conn_comps_of_this_comp:
            list_of_2_node_conn_comps.append(result_comp)
    

    return list_of_2_node_conn_comps





### TEST GRAPHS    

# 2 opened K4-s 
def create_simple_2_pearl_multigraph():
    # Create a MultiGraph
    G = nx.MultiGraph()

    # Add the first complete graph (nodes 0 to 3)
    G1_nodes = range(4)
    G.add_nodes_from(G1_nodes)
    G.add_edges_from((i, j) for i in G1_nodes for j in G1_nodes if (i < j and not (i==1 and j==2)))

    # Add the second complete graph (nodes 4 to 7)
    G2_nodes = range(4, 8)
    G.add_nodes_from(G2_nodes)
    G.add_edges_from((i, j) for i in G2_nodes for j in G2_nodes if (i < j and not (i==4 and j==5)))

    # Connect the two graphs with two edges
    G.add_edge(1, 4)  # First connecting edge
    G.add_edge(2, 5)  # Second connecting edge

    # Add problem specific meta info to graph
    G.graph['k'] = 3           
    G.graph['root'] = 0

    return G

# Ion has one big 2 connected component, good for 
# testing multi-leveled pearl stuff
def create_big_2_conn_from_Ion():
    '''
    Reads in Ion.graphml, takes the big component of it (that has pearl depth of 5 and 7 pearls),
    discards the degree 2 nodes and returns the graph
    '''
    #read in graph 
    file_path='/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/Ion.graphml'
    G, is_list = read_in_graph(file_path)
    Graph = nx.MultiGraph(G)
 
    #discard degree 2 nodes, so found pearls are not littered by them
    trimmed_graph = Graph.copy()
    trimmed_graph = contract_paths_keep_root(trimmed_graph)

    #for consistency lets assume, that the root is a degree 3 node:
    if 'root' not in trimmed_graph.graph.keys():                       #if not specified earlier, choose a random root
        trimmed_graph.graph['root'] = random.choice(list(trimmed_graph.nodes))



    #dissect into 2-connected components
    #comps = cut_cutting_edges(trimmed_graph)
    comps = separate_by_cutting_nodes(trimmed_graph)

    #discard single nodes and single edges
    proper_comps = [comp for comp in comps if len(comp.nodes)>2]
    if proper_comps == []:             #if every 2 connected component is a single node, then it is a tree
        return 1, 0, [], False
    depth_list = []

    # now it is only one comp
    # measure each component
    for comp in proper_comps:
        #get the graph structure we need
        graph = trimmed_graph.subgraph(comp).copy()

        # finding the new root
        # if the original root is not in the component, then it is the closest one to the root
        if trimmed_graph.graph['root'] not in comp:
            #we need to find the unique closest node to the root
            lista = list(comp)
            lista.sort(key = lambda node: len(nx.shortest_path(Graph, source=str(Graph.graph['root']), target=str(node))))
            closest_node = lista[0]
            graph.graph['root'] = closest_node

        #discard degree 2 nodes, so found pearls are not littered by them
        #if this is before finding the new root, that can make problems
        graph = contract_paths_keep_root(graph)
    return graph
    




def create_2_conn_example_graph():
    G, is_list = read_in_graph('/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/HiberniaGlobal.graphml')

    trimmed_graph = G.copy()   
    trimmed_graph = contract_paths_keep_root(trimmed_graph)

    if 'root' not in trimmed_graph.graph.keys():                       #if not specified earlier, choose a random root
        trimmed_graph.graph['root'] = random.choice(list(trimmed_graph.nodes)) 

    #dissect into 2-connected components
    comps = cut_cutting_edges(trimmed_graph)
    #discard single nodes
    proper_comps = [comp for comp in comps if len(comp)>1]
    if proper_comps == []:             #if every 2 connected component is a single node, then it is a tree
        return 1, 0    
    #get the graph structure we need
    graph = trimmed_graph.subgraph(proper_comps[0]).copy()
    #discard degree 2 nodes, so found pearls are not littered by them
    graph = contract_paths_keep_root(graph)

    return graph



##### GARBAGE/NOT USED INTERESTING STUFF:

def cactus_from_pearls(simple_pearl_list, graph):
    cactus = nx.MultiGraph()
    for pearl in reversed(simple_pearl_list):
        cactus.add_node(pearl['id'])
        if pearl['neighbor_nodes'] != 'this is the root':
            neighbor_pearl_ids = []
            for neighbor in pearl['neighbor_nodes']:
                neighbor_pearl_ids.append(graph.nodes[neighbor]['pearl_id'])
                
            if neighbor_pearl_ids[0]!=neighbor_pearl_ids[1]:
            # creating the cactus
                if not cactus.has_edge(neighbor_pearl_ids[0], pearl['id']):
                    cactus.add_edge(neighbor_pearl_ids[0], pearl['id'])
                if not cactus.has_edge(neighbor_pearl_ids[1], pearl['id']):
                    cactus.add_edge(neighbor_pearl_ids[1], pearl['id'])

    return cactus






