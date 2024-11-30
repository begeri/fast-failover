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

### UTILS
#(at least what is not imported from graph_analise)

#
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


#TODO
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
    graph_pieces = G.copy()

    list_of_2_node_conn_comps = []

    # if we are 2-node-connected, return the whole graph 
    if cuts==[]:
        return [G]

    #if we are not, recursion by one of the cutting nodes
    curr_cutting_node = cuts[0]
    graph_pieces.remove_node(curr_cutting_node)
    comps = nx.connected_components(graph_pieces)

    list(G.neighbors(curr_cutting_node))
    for comp in comps:
        smaller_graph = G.subgraph([curr_cutting_node, *comp]).copy()
        #if root is not in the smaller comp
        if G.graph['root'] not in comp:
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



#ARBORESCENCES

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

##PEARLS

#contract each path into an edge, if 'root' is not on the path
def contract_paths_keep_root(G):
    '''
    Nomen est omen, no 2 degree nodes
    Works only with MultiGraph
    '''
    nodes_to_remove = [v for v in G.nodes if G.degree(v) == 2]

    #if the root is of degree 2, doubling the edges means, that we get a routing, that first tries to reach d
    #in one of the edges and then on the other.
   
    if G.graph['root'] in nodes_to_remove:
        nodes_to_remove.remove(G.graph['root'])
        neighbors = list(G.neighbors(G.graph['root']))
        for neighbor in neighbors:
            G.add_edge(G.graph['root'], neighbor )

    for v in nodes_to_remove:
        neighbors = list(G.neighbors(v))
        G.remove_node(v)
        if len(neighbors)==2:
            G.add_edge(neighbors[0], neighbors[1])
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

def totally_contract_pearl(multi_graph, comp):
    '''
    nomen est omen
    multi_graph - the graph we are working on
    comp - the nodes of the pearl to contract

    We replace the pearl with one edge between the neighbors. 
    If there is only one neighbor, we just delete the pearl.
    '''
    #contract pearls
    neighbors = list(neighbor_of_subgraph(multi_graph, comp))
    
    multi_graph.remove_nodes_from(comp)
    if len(neighbors)==2:
        multi_graph.add_edge(neighbors[0],neighbors[1])
 


    return multi_graph



#create the arborescence of the pearls in a 2-edge-connected graph
def partition_2_conn_into_pearls(multi_graph):
    '''
    Creates the partition of the pearls, indexes the pearls, levels the pearls. Stores a pointer to the corresponding pearl dict in each
    node. Computes the parent-child relations of the pearls.

    INPUT: multi_graph - nx.MultiGraph with multi_graph.graph['root'] as the root of the graph 

    RETURNS: the level of the root pearl, which is the deepest
    '''
    root = multi_graph.graph['root']
    list_of_all_pearls = []
    remaining_graph = multi_graph.copy()
    current_pearl_level = 1
    newest_pearl_index = -1
    while is_3_connected(remaining_graph) == False:
        #there are still leaves!
        #print('Ezt nézed')
        #print(remaining_graph.edges)
        #print(remaining_graph.graph['root'])
        current_level_pearls = leaf_pearls(remaining_graph, level=current_pearl_level)
        root_pearl = None

        for pearl in current_level_pearls:
            #if root is in it: its not in this level, remove from the list, then pass
            if root in pearl['nodes']:
                root_pearl = pearl
            else:
                #for the rest of the pearls:
                newest_pearl_index+=1
                pearl['id'] = newest_pearl_index
                for node in pearl['nodes']:
                    multi_graph.nodes[node]['pearl_id'] = pearl['id']
                #remove the pearl
                remaining_graph = totally_contract_pearl(remaining_graph, pearl['nodes'])
        #if is_3_connected(remaining_graph): 
                

            
        if root_pearl != None: current_level_pearls.remove(root_pearl)


            
        #end of loop accounting
        current_pearl_level+=1
        list_of_all_pearls.append(current_level_pearls)
        #if current_pearl_level>=3: break


    #the last 3 connected graph is the highest level pearl - the root pearl
    comp = nx.connected_components(remaining_graph)
    #setting up the root_dict
    root_dict = {}
    #the nodes in a pearl
    root_dict['nodes'] = list(*comp)
    newest_pearl_index+=1
    root_dict['id'] = newest_pearl_index
    for node in root_dict['nodes']:
        multi_graph.nodes[node]['pearl_id'] = root_dict['id']
    root_dict['level_of_pearl'] = current_pearl_level

    #irrelevant stuff at the root
    root_dict['boundary'] = 'this is the root'
    root_dict['cut_edges'] = 'this is the root'
    root_dict['neighbor_nodes'] = 'this is the root'
    root_dict['parent_pearl'] = 'this is the root'

    list_of_all_pearls.append([root_dict])

    for level_list in list_of_all_pearls:
        for pearl in level_list:
            if pearl['neighbor_nodes'] != 'this is the root':
                for neighbor in pearl['neighbor_nodes']:
                    id_s = [multi_graph.nodes[neighbor]['pearl_id'] for neighbor in pearl['neighbor_nodes'] ]
                max(id_s)
                pearl['parent_pearl'] = max(id_s)

    num_pearls = newest_pearl_index+1
    return current_pearl_level, num_pearls, multi_graph, list_of_all_pearls

#measure the pearl depth of a graph        
def pearl_depth_of_graph(simple_graph, verbose = False):
    #read in graph 
    Graph = nx.MultiGraph(simple_graph)
    if 'root' not in Graph.graph.keys():                       #if not specified earlier, choose a random root
        Graph.graph['root'] = random.choice(list(Graph.nodes))

    #print('Original root', Graph.graph['root'])

    if verbose: 
        print('and the edges are')
        print(Graph.edges)


    #discard degree 2 nodes, so found pearls are not littered by them
    trimmed_graph = Graph.copy()
    trimmed_graph = contract_paths_keep_root(trimmed_graph)
    #print('Trimmed graph nodes', trimmed_graph.nodes)

    #dissect into 2-connected components
    #comps = cut_cutting_edges(trimmed_graph)
    comps = separate_by_cutting_nodes(trimmed_graph)
    #discard single nodes
    proper_comps = [comp for comp in comps if len(comp.nodes)>2]
    if proper_comps == []:             #if every 2 connected component is a single node, then it is a tree
        return 1, 0, []
    depth_list = []

    if verbose:
        print('We have proper comps')
    #measure each component
    for comp in proper_comps:
        #if verbose:
        #print('Examining', comp.nodes)
        #get the graph structure we need
        graph = trimmed_graph.subgraph(comp).copy()
        #print(graph.graph['root'])
        #if the original root is not in the component, then it is the closest one to the root
        if Graph.graph['root'] not in comp:
            #we need to find the unique closest node to the root

            #print('Graph nodes: ', Graph.nodes)
            lista = list(comp)
            lista.sort(key = lambda node: len(nx.shortest_path(Graph, source=str(Graph.graph['root']), target=str(node))))
            closest_node = lista[0]
            graph.graph['root'] = closest_node
            #print('Original one', Graph.graph['root'], 'closest node', closest_node )
        
        #discard degree 2 nodes, so found pearls are not littered by them
        #if this is before finding the new root, that can make problems
        graph = contract_paths_keep_root(graph)

        #Check those pearls
        #print('3as jel')
        #print(graph.nodes)
        depth, num_pearls, graph_with_pearl_data, list_of_all_pearls = partition_2_conn_into_pearls(graph)
        pearl_depth = nx.diameter(extract_pearl_tree(graph))+1
        #print(depth, pearl_depth)
        
        depth_list.append(pearl_depth)

    return max(depth_list), num_pearls, list_of_all_pearls


### MEASUREMENTS


def check_a_graph(file_path, verbose = False):
    '''
    prints and returns different graph properties
    '''
    #pearl relatd infos
    G, is_list = read_in_graph(file_path)
    pearl_depth, num_pearls, pearl_list = pearl_depth_of_graph(G, verbose=False)

    # compute the average pearl size, min size, max size
    pearl_size_list = []
    min_pearl_size = -200
    max_pearl_size = -200
    if pearl_list != []:
        for level in pearl_list:
            for pearl in level:
                curr_size = len(list(pearl['nodes']))
                pearl_size_list.append(curr_size)
                #updating the min, max
                if min_pearl_size < 0:         min_pearl_size = curr_size
                if min_pearl_size > curr_size: min_pearl_size = curr_size
                if max_pearl_size < curr_size: max_pearl_size = curr_size
        avg_pearl_size = sum(pearl_size_list)/len(pearl_size_list)
    else:
        min_pearl_size = 0
        max_pearl_size = 0
        avg_pearl_size = 0
    
    #get basic infos
    Graph = nx.MultiGraph(G)
    node_num = len(Graph.nodes)
    edge_num = len(Graph.edges)

    #check if subdivision graph
    is_subdivision = False
    if pearl_depth == 1:
        for node in G.nodes:
            if G.degree(node) == 2:
                is_subdivision = True


    if verbose:
        print('node num = ', node_num, ', edge num =', edge_num, ', pearl depth = ', pearl_depth,' num pearls = ', num_pearls,', at ', file_path)
        print('average pearl size = ', avg_pearl_size, 'min = ', min_pearl_size, 'max = ', max_pearl_size)
    return pearl_depth, node_num, edge_num, num_pearls, avg_pearl_size, min_pearl_size, max_pearl_size, is_subdivision

def check_everything():

    fnames = [ '/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/'+fname for fname in os.listdir('/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/') if fname.endswith('.graphml')] 
    fnames = filter(lambda f: nx.edge_connectivity(nx.read_graphml(f))>0,fnames)

    graph_data = pd.DataFrame(columns=['NUM NODES', 'NUM EDGES', 'PEARL DEPTH', 'GRAPH NAME', 'NUM PEARLS'])
    subdivison_graphs_data = pd.DataFrame(columns=['NUM NODES', 'NUM EDGES', 'PEARL DEPTH', 'GRAPH NAME', 'NUM PEARLS'])
    print(graph_data)
    index_of_graph = 0
    deepest = ''
    depth = 0
    for name in fnames:
        print(name) 
        end_of_name = name.split(sep = '/')[-1]
        #getting the data
        this_depth, node_num, edge_num, num_pearls, avg_pearl_size, min_pearl_size, max_pearl_size, is_subdivision = check_a_graph(name)
        print(this_depth)
        #converting the data to dataframe
        this_data = pd.DataFrame( index=[index_of_graph], data =  {'NUM NODES': node_num, 'NUM EDGES': edge_num, 'PEARL DEPTH': this_depth, 'GRAPH NAME':end_of_name, 'NUM PEARLS': num_pearls, 'AVG PEARL SIZE': avg_pearl_size, 'MIN PEARL SIZE': min_pearl_size, 'MAX PEARL SIZE': max_pearl_size})
        graph_data = pd.concat([graph_data, this_data], axis=0)
        if this_depth>=depth:
            deepest = name
            depth = this_depth 

        #listing the subdivision graphs
        if is_subdivision:
            subdivison_graphs_data = pd.concat([subdivison_graphs_data, this_data], axis=0)

        #next graph
        index_of_graph+=1
    print( deepest, depth)

    return graph_data, subdivison_graphs_data

#reads in a graph and for each possible choice of root computes the pearl depth 
def analyze_a_graph(file_path):
    '''
    Reads in the input graph and checks it num_experiences times. It chooses a new random root every time. 
    '''
    G, is_list = read_in_graph(file_path)

    for node in G.nodes:
        G.graph['root'] = node
        pearl_depth, num_pearls, pearl_list = pearl_depth_of_graph(G, verbose=False)
        print(pearl_depth)

    
    Graph = nx.MultiGraph(G)
    node_num = len(Graph.nodes)
    edge_num = len(Graph.edges)

    #print('node num = ', node_num, ', edge num =', edge_num, ', pearl depth = ', pearl_depth,' num pearls = ', num_pearls,', at ', file_path)
    
    return pearl_depth, node_num, edge_num, num_pearls

#given a 2 connected multigraph partitioned into pearls, creates the tree of pearls
def extract_pearl_tree(multi_graph):
    '''
    Expects a 2 connected multi_graph, returns the arborescence of pearls
    '''
    current_pearl_level, newest_pearl_index, multi_graph, list_of_all_pearls = partition_2_conn_into_pearls(multi_graph)
    pearl_tree = nx.Graph()
    for level in reversed( list_of_all_pearls):
        for pearl in level:
            pearl_tree.add_node(pearl['id'])
            pearl_tree.nodes[pearl['id']]['pearl_data'] = pearl
            if pearl['parent_pearl'] != 'this is the root':
                pearl_tree.add_edge(pearl['id'], pearl['parent_pearl'])

    #print('List of all pearls')
    #print(list_of_all_pearls)
    print("Tree edges, info")
    print(pearl_tree.edges)
    print(list(nx.connected_components(pearl_tree)))
    return pearl_tree





##ROUTING

#TODO
# routes a general graph
def route_a_multigraph(multigraph):
    '''
    todo
    '''
    #

    # partition into 2-connected components
    # list the connections
    
    # route each component

    # route the connections of the connected components
    return 0

#TODO
# routes a 2-connected graph
def route_a_2_conn_multigraph(multigraph):
    '''
    Partition the edges of a 2 connected multigraph into pearls, then route each pearl
    '''
    # for testing the code working:
    # list the edges to check if the are routed. 

    # partition the edges of multigraph into pearls
    current_pearl_level, num_pearls, multi_graph, list_of_all_pearls = partition_2_conn_into_pearls(multi_graph)

    # starting from backwards, route the pearls backwards
      # route the pearl for basic getting out
      # route the pearl for traversing
      # insert it into the previous routing
        # find the edge, that was replacing it
        # modify the neighbor nodes, the previous replacing edge now is the cutting edges 

    return 0

#TODO
#routing of closed_pearls   #rout_3_conn_graph() + route_a_pearl_for_traversing
def route_a_pearl(P, multigraph, multidigraph):
    '''
    Given a pearl, with boundaries, etc, it generates a arborescense based routing with rout_3_conn_graph(), and 
    route_a_pearl_for_traversing
    INPUT: 
    P - pearl dictionary
    multigraph - nx.MultiGraph
    multidigraph - nx.MultiDiGraph, to write the routing on, multigraph.to_directed()
    '''
    print('Here we start to route the pearl')
    #this is what we are working with
    pearl_graph = multigraph.subgraph(P['nodes']).copy()

    # add boundary_1, boundary_2 as an edge
    boundary_1 = P['boundary'][0]
    boundary_2 = P['boundary'][1]
    pearl_graph.add_edge(boundary_1, boundary_2)
    # If we use a specific key, the GreedyArborescences crashes
    edge_data = pearl_graph[boundary_1][boundary_2]  
    keys = list(edge_data.keys())  
    boundary_edge_key = keys[-1] 
    P['boundary_edge_key'] = boundary_edge_key 

    #we need to get to one of the edges #here we have a freedom of choice
    pearl_graph.graph['root'] = boundary_1
    pearl_digraph = pearl_graph.to_directed() 

    # route the pearl with rout_3_conn_graph
    rout_3_conn_graph(P, pearl_graph, pearl_digraph, (boundary_1, boundary_2, boundary_edge_key) )

    # route the pearl for traversing 
    # for cut edges to be in the graph, pearl_graph is not enough
    route_a_pearl_for_traversing(P, multidigraph, multigraph)

    '''
    # check the result
    print('root', pearl_digraph.graph['root']) 
    for node in pearl_digraph.nodes:
        print('node ', node)
        if node is not pearl_digraph.graph['root']:
            print(pearl_digraph.nodes[node]['R'])
    '''
    

#TODO: bouncing bit!
#routing of 3-connected graphs via arborescences
def rout_3_conn_graph(P, multigraph, multidigraph, edge):
    '''
    given a 3 connected graph, that has an extra attribute, g.graph['root'] and a specific 'edge', make a 2 resilient routing based on
    circular arborescences, where the first arborescence contains 'edge'
    '''
    # recreate the pearl graph
    pearl_graph = multigraph.subgraph(P['nodes']).copy()
    boundary_1 = P['boundary'][0]
    pearl_graph.graph['root'] = boundary_1
    pearl_digraph = pearl_graph.to_directed()



    #find the arborescences
    unordered_arborescence_list = GreedyArborescenceDecomposition(multidigraph)
    arborescence_list = []
    # pick the arb containing the boundary edge
    last_arb = None
    for arb in unordered_arborescence_list:
        if edge in arb.edges(keys=True):              #find the arb, that contains the edge that 
            last_arb = arb
            unordered_arborescence_list.remove(arb)
    for arb in unordered_arborescence_list:
        arborescence_list.append(arb)
    arborescence_list.append(last_arb)
    

    for node in multidigraph.nodes:
        if node is not multidigraph.graph['root']:
            R = {}
            # if you started under this pearl, you need to get out of here
            R_getting_out = {}
            R_bouncing    = {}
            for edge in multidigraph.in_edges(nbunch= [node], keys=True):
                # reindexing the arborescences
                if edge in arborescence_list[0].edges:
                    tree_index = 0
                elif edge in arborescence_list[1].edges:
                    tree_index = 1
                elif edge in arborescence_list[2].edges:
                    tree_index = 2
                else: 
                    print('edge not in any tree', edge)
                    tree_index = None

                if tree_index is not None:
                    curr_tree = arborescence_list[tree_index]
                    next_tree = arborescence_list[(tree_index+1)%3]
                    prev_tree = arborescence_list[(tree_index+2)%3]
                    if len(list(curr_tree.out_edges([node], keys = True))) != 1:
                        print('Tree ', tree_index, ' has more than one out edge at node ', node, '!')

                    out_edge_1 = list(curr_tree.out_edges([node], keys = True))[0]
                    out_edge_2 = list(next_tree.out_edges([node], keys = True))[0]
                    out_edge_3 = list(prev_tree.out_edges([node], keys = True))[0]
                    R_getting_out[edge] = [out_edge_1, out_edge_2, out_edge_3]
                    #R_bouncing = #TODO
            #if the package starts at 'node'
            arb_0_edges_list = list(arborescence_list[0].edges)
            arb_1_edges_list = list(arborescence_list[1].edges)
            arb_2_edges_list = list(arborescence_list[2].edges)
            out_edge_in_arb_0 = [edge for edge in arb_0_edges_list if (edge[0]==node)][0]
            out_edge_in_arb_1 = [edge for edge in arb_1_edges_list if (edge[0]==node)][0]
            out_edge_in_arb_2 = [edge for edge in arb_2_edges_list if (edge[0]==node)][0]
            R_getting_out['start'] = [out_edge_in_arb_0, out_edge_in_arb_1, out_edge_in_arb_2]

            #saving the local dict to the big dict
            R['getting_out'] = R_getting_out
            multidigraph.nodes[node]['R'] = R

    return None 

#TODO: bouncing bit! (here we actually need to think about it)
#routing of parallel routes in a graph
def route_a_pearl_for_traversing(P, multidigraph, pearl_graph):
    '''
    find 2 (or more) edge-disjoint paths between the boundaries of the pearl, then do the circular routing for 
    parallel paths
    INPUT:
    P - dictionary for pearl data
    multigraph     - multigraph, that has been partitioned into pearls, and P is one of them
    multidigraph   - multidigraph, thats edges will be used for partitioning
    '''
    # find the boundaries
    if len(P['boundary'])==1:                          # if we are only using 2-edge-connectivity, this is needed
        return 'Nope, we do not traverse this'
    else:
        boundary_1 = P['boundary'][0]
        boundary_2 = P['boundary'][1]

    #create the pearl subgraph
    pearl_id = P['id']
    if (boundary_1, boundary_2, f'boundary_edge_for_pearl_{pearl_id}') in pearl_graph.edges:
        print('Ki kellett venni')
        boundary_edge_key = P['boundary_edge_key']
        pearl_graph.remove_edge(boundary_1, boundary_2, key = boundary_edge_key)


    # find 2 arc-disjointed pathes
    paths = list(nx.edge_disjoint_paths(pearl_graph, boundary_1, boundary_2, cutoff=2))
    path_1 = paths[0]
    path_2 = paths[1]
    # create the 2-path-graph 
    path_graph_1 = pearl_graph.subgraph(path_1).copy()
    path_graph_2 = pearl_graph.subgraph(path_2).copy()
    parallel_path_graph = path_graph_1.copy()
    for node in path_graph_2:
        if node not in parallel_path_graph.nodes: parallel_path_graph.add_node(node)
    for edge in path_graph_2.edges(keys=True): 
        if edge not in parallel_path_graph.edges(keys=True): parallel_path_graph.add_edge(*edge, keys=True)
    parallel_path_graph = parallel_path_graph.to_directed()

    
    # route_degree_2_node_for_traversing for the inner ones
    for node in path_1[1:-1]:
        route_degree_2_node_for_traversing(parallel_path_graph, node)
    for node in path_2[1:-1]:
        route_degree_2_node_for_traversing(parallel_path_graph, node)
    
   
    

    # route the boundaries
    all_edges_list = list(multidigraph.edges)

    # path_1 goes primarily from path_1[0] to path_1[-1] (boundary_1 -> boundary_2)
    # path_2 goes primarily from path_1[-1] to path_1[0] (boundary_2 -> boundary_1)
    # need the 3 edges from boundaries, 3rd is in P['cut_edges']

    #choosing the appropriate in and out edges
    #boundary_1
    out_cut_edge_1    =  [edge for edge in all_edges_list if (edge[0]==boundary_1 and (edge[1] in P['neighbor_nodes']))][0]
    in_cut_edge_1     =  [edge for edge in all_edges_list if (edge[1]==boundary_1 and (edge[0] in P['neighbor_nodes']))][0]
    p1_out_edge_1     =  [edge for edge in all_edges_list if (edge[0]==boundary_1 and (edge[1] in path_1))][0]
    p1_in_edge_1      =  [edge for edge in all_edges_list if (edge[1]==boundary_1 and (edge[0] in path_1))][0]
    p2_out_edge_1     =  [edge for edge in all_edges_list if (edge[0]==boundary_1 and (edge[1] in path_2))][0]
    p2_in_edge_1      =  [edge for edge in all_edges_list if (edge[1]==boundary_1 and (edge[0] in path_2))][0]
    #boundary_2
    out_cut_edge_2    =  [edge for edge in all_edges_list if (edge[0]==boundary_2 and (edge[1] in P['neighbor_nodes']))][0]
    in_cut_edge_2     =  [edge for edge in all_edges_list if (edge[1]==boundary_2 and (edge[0] in P['neighbor_nodes']))][0]
    p1_out_edge_2     =  [edge for edge in all_edges_list if (edge[0]==boundary_2 and (edge[1] in path_1))][0]
    p1_in_edge_2      =  [edge for edge in all_edges_list if (edge[1]==boundary_2 and (edge[0] in path_1))][0]
    p2_out_edge_2     =  [edge for edge in all_edges_list if (edge[0]==boundary_2 and (edge[1] in path_2))][0]
    p2_in_edge_2      =  [edge for edge in all_edges_list if (edge[1]==boundary_2 and (edge[0] in path_2))][0]

   
    #creating dicts for boundary_1
    R_1 = {}                            # routing table for boundary_1 
    R1_traversing_dict = {} 

    R1_traversing_dict[in_cut_edge_1]    = [p1_out_edge_1, p2_out_edge_1, out_cut_edge_1]
    R1_traversing_dict[p1_in_edge_1] = [p2_out_edge_1, out_cut_edge_1, p1_out_edge_1]
    R1_traversing_dict[p2_in_edge_1] = [out_cut_edge_1, p1_out_edge_1, p2_out_edge_1]
    R_1['traversing'] = R1_traversing_dict

    #creating dicts for boundary_2
    R_2 = {}
    R2_traversing_dict = {} 
    R2_traversing_dict[in_cut_edge_2]    = [p2_out_edge_2, p1_out_edge_2, out_cut_edge_2]
    R2_traversing_dict[p2_in_edge_2] = [p1_out_edge_2, out_cut_edge_2, p2_out_edge_2]
    R2_traversing_dict[p1_in_edge_2] = [out_cut_edge_2, p2_out_edge_2, p1_out_edge_2]
    R_2['traversing'] = R2_traversing_dict

    #TODO is one of the cut_edges has failed, need to switch to bouncing mode here.

    # saving the dicts
    #if we already stored data previously, just add some new keys
    if 'R' not in parallel_path_graph.nodes[boundary_1].keys():
        parallel_path_graph.nodes[boundary_1]['R'] = R_1
    else: 
        parallel_path_graph.nodes[boundary_1]['R']['traversing'] = R1_traversing_dict
    
    if 'R' not in parallel_path_graph.nodes[boundary_2].keys():
        parallel_path_graph.nodes[boundary_2]['R'] = R_2
    else: 
        parallel_path_graph.nodes[boundary_2]['R']['traversing'] = R2_traversing_dict
    

    
    print('parallel path node data')
    for node in parallel_path_graph:
        print(node)
        multidigraph.nodes[node]['R'] = parallel_path_graph.nodes[node]['R']
        print(multidigraph.nodes[node])
    

#routing of degree 2 nodes, if we are traversing
def route_degree_2_node_for_traversing(multidigraph, node):
    '''
    It just forwards, if possible.
    '''
    R = {}
    dict_for_traversing = {}
    for in_edge in multidigraph.in_edges(node, keys=True):
        out_edge_list = list(multidigraph.out_edges(node, keys=True) )
        if out_edge_list[0][1] == in_edge[0]:                                    # just lists the nodes, decides which out_edge
            dict_for_traversing[in_edge] = [out_edge_list[1],out_edge_list[0]]   # goes forward, and puts that to first, the other one second
        else:                                                                   
            dict_for_traversing[in_edge] = [out_edge_list[0],out_edge_list[1]]
    
    R['traversing'] = dict_for_traversing 

    if 'R' not in multidigraph.nodes[node].keys():
        multidigraph.nodes[node]['R'] = R
    else: 
        multidigraph.nodes[node]['R']['traversing'] = dict_for_traversing

    multidigraph.nodes[node]['R'] = R
    
    return None

#TODO: bouncing bit! (here we need to think about it)
# routing of degree 2 nodes if we are not traversing
def route_degree_2_node_outside_of_pearl(multidigraph, node):
    pass       
        

# The header info of a package now can be represented as a dict.
# 
# package_info['bouncing'],
# package_info['pearl_level'] 
# package_info['num_steps']
#
# A routing table in given node v looks like a multi levelled dict, keys are possible header configs, 
# items are the routing permutation on the out-edges
#
# R['incoming_edge']['in_pearl'] = 
# R['incoming_edge']['traversing'] = 
#
# Now we write it onto the nodes of the graph, so we don't need the R['v']['incoming_edge']['traversing'] notation.



#routing a whole graph
def pearl_based_routing(g):
    '''
    INPUT: 
    g - nx.MultiGraph with a designated root, denoted at g.graph['root']
    '''
    C = cactus(g)

    pearl_tree = arborescence_of_pearls(C)

    #create a list of 3 connected graphs for each node of the pearl_tree

    #the remember the parent-edges of deeper pearls, and boundary nodes/edges for the pearls

    #create routing for the closed pearls.

    #match them together -> hard part
    return g


### TESTS, SIMULATIONS

#TODO
def compute_primary_path():
    pass

#TODO
def compute_path_under_1_fail():
    pass

#TODO
def compute_path_under_2_fails():
    pass

#TODO
# simple failure experiments, path lengths, etc.
# running time experiments
# 3 failure experiments
# are there subdivision graphs, if yes, how many

### RUNNING STUFF ########################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

def main(args):
    ex_multigraph = create_simple_2_pearl_multigraph()
    current_pearl_level, num_pearls, multi_graph, list_of_all_pearls = partition_2_conn_into_pearls(ex_multigraph)

    ex_multidigraph = ex_multigraph.to_directed()
    example_pearl = list_of_all_pearls[0][0]
    print(ex_multigraph.subgraph(example_pearl['nodes']).edges)
    print('Start of the pearl routing')
    route_a_pearl(example_pearl, ex_multigraph, ex_multidigraph)



    '''
    print('Result')
    for node in multidigraph.nodes:
        print(node, multidigraph.nodes[node])
    '''
        
#iterate through topology zoo and write node num, edge num, pearl num and pearl depth
def pearl_depth_experiment():
#def main(args):
    graph_data, subdivision_data = check_everything()
    graph_data.to_csv(path_or_buf='/mnt/d/Egyetem/Routing Cikk/fast-failover/pearl-algo/topology_zoo_statistics_with_pearl_sizes.csv')
    subdivision_data.to_csv(path_or_buf='/mnt/d/Egyetem/Routing Cikk/fast-failover/pearl-algo/topology_zoo_statistics_subdivision_stats.csv')

def main_1(args):
    pearl_depth, node_num, edge_num, num_pearls, avg_pearl_size, min_pearl_size, max_pearl_size, is_subdivision = check_a_graph('/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/York.graphml')  
    print(pearl_depth, node_num, edge_num, num_pearls)


def separation_tests(args):
    file_path='/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/Aarnet.graphml'
    G, is_list = read_in_graph(file_path)

    Graph = nx.MultiGraph(G)
    Graph.graph['root'] = random.choice(list(Graph.nodes))
    Graph = contract_paths_keep_root(Graph)

    print(separate_by_cutting_nodes(Graph))


    


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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--file_path', type=str, default="/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/HiberniaGlobal.graphml",
                        help='Where is the graph')
    
    args = parser.parse_args()
    main(args)


def create_2_conn_example_graph():
    G, is_list = read_in_graph('/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/HiberniaGlobal.graphml')

    trimmed_graph = G.copy()
    if 'root' not in trimmed_graph.graph.keys():                       #if not specified earlier, choose a random root
        trimmed_graph.graph['root'] = random.choice(list(trimmed_graph.nodes))    
    trimmed_graph = contract_paths_keep_root(trimmed_graph)

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






    
    '''
    reset_arb_attribute(g)
    print(nx.edge_connectivity(g, 0,1))

    arbs = GreedyArborescenceDecomposition(g)
    for tree in arbs:
        print(tree.edges)

    #you can put a 'routing' dictionary onto the nodes.
    for node in g.nodes:
        g.nodes[node]['routing'] = {}
        #g.nodes[node]['routing']['edge_1'] = ['edge_2', 'edge_3']

    for item in g.nodes.items():
        print(item)

    #iterating through the neihbors of a node:
    for neighbor in g.neighbors(1):
        pass
        #print(neighbor)
        #print(g[1][neighbor][0])
    

    #arb_list = arbs.GreedyArborescenceDecomposition(g)
    #print(arb_list)
    
    #for arb in arb_list:
        #print(arb.edges)
    '''