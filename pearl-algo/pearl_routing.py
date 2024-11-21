import sys
sys.path.append('/mnt/d/Egyetem/Routing Cikk/fast-failover')

import networkx as nx
import argparse
from graph_analise import find_pearls, neighbor_of_subgraph
from graph_analise import *
from graph_gen import add_ear

from heapq import heappush, heappop
import random

import networkx as nx
import pandas as pd

random.seed(42)


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



#TODO
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
    current_pearl_level = 0
    newest_pearl_index = -1
    while is_3_connected(remaining_graph) == False:
        #there are still leaves!
        #print('pearl level is ', current_pearl_level)
        #print(remaining_graph.edges)
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
            for neighbor in pearl['neighbor_nodes']:
                if pearl['neighbor_nodes'] != 'this is the root':
                    pearl['parent_pearl'] = multi_graph.nodes[neighbor]['pearl_id']

    return current_pearl_level, newest_pearl_index

#measure the pearl depth of a graph        
def pearl_depth_of_graph(file_path, verbose = False):
    #read in graph 
    G, is_list = read_in_graph(file_path)
    Graph = nx.MultiGraph(G)
    Graph.graph['root'] = random.choice(list(Graph.nodes))
    #for debugging purposes:
    #Graph.graph['root'] = 7
    if verbose:
        print(Graph.nodes)
        print('root = ', Graph.graph['root'], ' is degree ', nx.degree(Graph, '7') )
        print('and the edges are')
        print(Graph.edges)


    #discard degree 2 nodes, so found pearls are not littered by them
    trimmed_graph = Graph.copy()
    trimmed_graph = contract_paths_keep_root(trimmed_graph)

    #dissect into 2-connected components
    comps = cut_cutting_edges(trimmed_graph)
    #discard single nodes
    proper_comps = [comp for comp in comps if len(comp)>1]
    if proper_comps == []:             #if every 2 connected component is a single node, then it is a tree
        return 0, 0
    depth_list = []

    #measure each component
    for comp in proper_comps:
        #get the graph structure we need
        graph = trimmed_graph.subgraph(comp).copy()
        #discard degree 2 nodes, so found pearls are not littered by them
        graph = contract_paths_keep_root(graph)

        #if the original root is not in the component, then it is the closest one to the root
        if Graph.graph['root'] not in comp:
            #we need to find the unique closest node to the root
            
            #print('Graph nodes: ', Graph.nodes)
            #print(len(nx.shortest_path(Graph, source='7', target='0')))
            lista = list(comp)
            lista.sort(key= lambda node: len(nx.shortest_path(Graph, source=str(Graph.graph['root']), target=str(node))))
            closest_node = lista[0]
            graph.graph['root'] = closest_node

        #Check those pearls
        depth, num_pearls = partition_2_conn_into_pearls(graph)
        depth_list.append(depth)

    return max(depth_list), num_pearls

def check_a_graph(file_path, verbose = False):
    pearl_depth, num_pearls = pearl_depth_of_graph(file_path, verbose=False)

    G, is_list = read_in_graph(file_path)
    Graph = nx.MultiGraph(G)
    node_num = len(Graph.nodes)
    edge_num = len(Graph.edges)
    print('node num = ', node_num, ', edge num =', edge_num, ', pearl depth = ', pearl_depth,' num pearls = ', num_pearls,', at ', file_path)
    return pearl_depth, node_num, edge_num, num_pearls

def check_everything():

    fnames = [ '/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/'+fname for fname in os.listdir('/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/') if fname.endswith('.graphml')] 
    fnames = filter(lambda f: nx.edge_connectivity(nx.read_graphml(f))>0,fnames)

    graph_data = pd.DataFrame(columns=['NUM NODES', 'NUM EDGES', 'PEARL DEPTH', 'GRAPH NAME', 'NUM PEARLS'])
    print(graph_data)
    index_of_graph = 0
    deepest = ''
    depth = 0
    for name in fnames:
        #try:
        print(name) 
        end_of_name = name.split(sep = '/')[-1]
        this_depth, node_num, edge_num, num_pearls = check_a_graph(name)
        print(this_depth)
        this_data = pd.DataFrame( index=[index_of_graph], data =  {'NUM NODES': node_num, 'NUM EDGES': edge_num, 'PEARL DEPTH': this_depth, 'GRAPH NAME':end_of_name, 'NUM PEARLS': num_pearls})
        graph_data = pd.concat([graph_data, this_data], axis=0)
        if this_depth>=depth:
            deepest = name
            depth = this_depth        
        #except:
        #    print('baj ', name)
        index_of_graph+=1
    print( deepest, depth)

    return graph_data

#TODO
def analyze_a_graph(file_path, num_experiences):
    '''
    Reads in the input graph and checks it num_experiences times. It chooses a new random root every time. 
    '''
    for _ in range(num_experiences):
        check_a_graph(file_path)




##ROUTING

#TODO
#routing of 3-connected graphs via arborescences
def rout_3_conn_graph(g):
    '''
    given a 3 connected graph, that has an extra attribute, g.graph['root'], make a 2 resilient routing based on
    circular arborescences
    '''

#TODO
#routing of closed_pearls   #might be the same as rout_3_conn_graph()
def route_a_pearl(P):
    '''
    todo
    '''

#TODO
#routing of parallel routes in a graph
def route_a_pearl_for_traversing(P):
    '''
    find 2 (or more) edge-disjoint paths between the boundaries of the pearl, then do the circular routing for 
    parallel paths
    '''

#TODO
#routing of degree 2 nodes
#might not be necessary to be a distinct function
def route_a_path(G):
    '''
    
    '''


#routing a whole graph
def pearl_based_routing(g):
    '''
    INPUT: 
    g - nx.MultiGraph with a designated root, denoted at g.graph['root']
    '''
    C = cactus(g)

    pearl_arb = arborescence_of_pearls(C)

    #create a list of 3 connected graphs for each node of the pearl_arb

    #the remember the parent-edges of deeper pearls, and boundary nodes/edges for the pearls

    #create routing for the closed pearls.

    #match them together -> hard part
    return g


### RUNNING STUFF

def main(args):
    #graph_data = check_everything()
    #print(graph_data)
    #graph_data.to_csv(path_or_buf='/mnt/d/Egyetem/Routing Cikk/fast-failover/pearl-algo/topology_zoo_statistics.csv')

    #depth, node_num, edge_num, num_pearls = check_a_graph('/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/Bellcanada.graphml', verbose=True)
    #print(depth, node_num, edge_num, num_pearls)
    
    analyze_a_graph(file_path='/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/HiberniaGlobal.graphml', num_experiences=10)
    


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--file_path', type=str, default="/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/HiberniaGlobal.graphml",
                        help='Where is the graph')
    
    args = parser.parse_args()
    main(args)






    
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