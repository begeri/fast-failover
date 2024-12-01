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

from heapq import heappush, heappop
import random

import networkx as nx
import pandas as pd

random.seed(42)



##PEARLS

#partition the nodes into pearls in a 2-edge-connected graph
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
    current_pearl_level =  1
    newest_pearl_index  = -1
    while is_3_connected(remaining_graph) == False:

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
                remaining_graph = totally_contract_pearl(remaining_graph, pearl)
                    
        if root_pearl != None: current_level_pearls.remove(root_pearl)
 
        #end of loop accounting
        current_pearl_level+=1
        list_of_all_pearls.append(current_level_pearls)


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

    
    #creating the simple list, where they are not in levels, so we can iterate through them
    simple_pearl_list = []
    for level in list_of_all_pearls:
        for pearl in level:
            simple_pearl_list.append(pearl)
    
    # from Gomory-Hu type parents to cactus type parents
    fixing_pearl_parent_pointers(simple_pearl_list, multi_graph)


    
    # at this point we have a pearl height, and it is only the same as pearl depth for the root pearl and the max length
    # pathes to leafs
    #if the parents are good, this gonna solve the problem of pearl levels
    for pearl in reversed(simple_pearl_list):
        if pearl['parent_pearl'] != 'this is the root':
            parent_pearl_id = pearl['parent_pearl']
            pearl['level_of_pearl'] = simple_pearl_list[parent_pearl_id]['level_of_pearl']-1


    num_pearls = newest_pearl_index+1
    return current_pearl_level, num_pearls, multi_graph, list_of_all_pearls

#measure the pearl depth of a graph        
def pearl_depth_of_graph(simple_graph, verbose = False):
    #convert to multigraph 
    OriginalGraph = nx.MultiGraph(simple_graph)


    is_subdivision = False
    list_of_all_pearls = []

    # For debugging/understanding print stuff
    if verbose: 
        print('At the beginning the edges are')
        print(OriginalGraph.edges)


    #discard degree 2 nodes, so found pearls are not littered by them
    trimmed_graph = OriginalGraph.copy()
    print(trimmed_graph.edges)
    trimmed_graph = contract_paths_keep_root(trimmed_graph)
    if len(list(trimmed_graph.nodes)) <1:    # if every node was less then degree 2, then we are a circle or a path or sg else not interesting
        return 1, 0, [], False

    #for consistency lets assume, that the root is a degree 3 node:
    if 'root' not in trimmed_graph.graph.keys():                       #if not specified earlier, choose a random root
        trimmed_graph.graph['root'] = random.choice(list(trimmed_graph.nodes))

    
    
    
    if verbose:
        print('The graph before dissecting')
        print(trimmed_graph.edges)

    #dissect into 2-connected components
    #comps = cut_cutting_edges(trimmed_graph)
    comps = separate_by_cutting_nodes(trimmed_graph)
    #discard single nodes and single edges and circles
    proper_comps = []
    for comp in comps:
        if len(list(comp.nodes))>2 and len(list(comp.edges))>len(list(comp.nodes)):
            proper_comps.append(comp)


    #proper_comps = [comp for comp in comps if len(comp.nodes)>2] 
    
    if proper_comps == []:             #if every 2 connected component is a single node, edge or circle, then it is a cactus
        return 1, 0, [], False
    depth_list = []
    big_graph_num_pearls = 0

    if verbose:
        print('We have proper comps')
    #measure each component
    for comp in proper_comps:
        #get the graph structure we need
        graph = trimmed_graph.subgraph(comp).copy()

        # finding the new root
        # if the original root is not in the component, then it is the closest one to the root
        if trimmed_graph.graph['root'] not in comp:
            #we need to find the unique closest node to the root
            lista = list(comp)
            lista.sort(key = lambda node: len(nx.shortest_path(OriginalGraph, source=str(trimmed_graph.graph['root']), target=str(node))))
            closest_node = lista[0]
            graph.graph['root'] = closest_node
        

        # in order for it to be a subdivision graph, there need to be a degree 2 node in a real pearl
        there_is_degree_2_node = False
        for node in graph.nodes:
            if graph.degree(node) == 2:
                there_is_degree_2_node = True


        #discard degree 2 nodes, so found pearls are not littered by them
        #if this is before finding the new root, that can make problems
        graph = contract_paths_keep_root(graph)

        # just to double check for subdivison graphs
        is_this_comp_3_connected = is_3_connected(graph)

        #Check those pearls
        depth, num_pearls, graph_with_pearl_data, list_of_pearls_in_this_comp = partition_2_conn_into_pearls(graph)
        
        # Bookkeeping of global properties, lists
        pearl_depth = nx.diameter(extract_pearl_tree(graph))+1  # depth of this comp
        depth_list.append(pearl_depth)                       
        list_of_all_pearls.append(list_of_pearls_in_this_comp)  # list of all pearls 
        big_graph_num_pearls += num_pearls                      # number of all pearls worth keeping tab

        #for the subdivision case, we need to have a component, that is exactly 1 deep and has a degree 2 node
        if there_is_degree_2_node and (pearl_depth==1):
            is_subdivision = True

        if not is_this_comp_3_connected and (pearl_depth==1):
            print('Gebasz')

    # subdivision:
    # if we found a component, that is a subdivision graph, but the max depth >1, then it is not a subdivision
    if max(depth_list)>1:
        is_subdivision = False

    return max(depth_list), big_graph_num_pearls, list_of_all_pearls, is_subdivision


### MEASUREMENTS


def check_a_graph(file_path, verbose = False):
    '''
    prints and returns different graph properties
    '''
    #pearl relatd infos
    G, is_list = read_in_graph(file_path)
    pearl_depth, num_pearls, pearl_list, is_subdivision = pearl_depth_of_graph(G, verbose=False)

    # compute the average pearl size, min size, max size
    pearl_size_list = []
    min_pearl_size = -200
    max_pearl_size = -200
    if pearl_list != []:
        for comp_list in pearl_list:
            for level in comp_list:
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
    Graph    = nx.MultiGraph(G)
    node_num = len(Graph.nodes)
    edge_num = len(Graph.edges)


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
        print(this_depth, node_num, edge_num, num_pearls, avg_pearl_size, min_pearl_size, max_pearl_size, is_subdivision)
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


    return pearl_tree




### FORWARDING

def choose_first_available(out_edge_list, failure_set):
    '''
    From a list of out going edges choose the first, that is not failed
    '''
    
    for (u,v,k) in reversed(out_edge_list):
        if ((u,v,k) not in failure_set) and ((v,u,k) not in failure_set):
            out_edge = (u,v,k)
    return out_edge
        


def send_a_package(multidigraph, starting_node, failure_set, max_step):
    '''
    Given a multidigraph properly labeled with routing tables and a starting node from that graph, 
    go after the routing tables
    It cannot use the edge set 'failure set'.
    INPUT:
    multidigraph - a nx.MultiDiGraph containing the routing infos and the pearl infos
    starting_node - a node from multidigraph
    failure_set - a subset of edges from multidigraph, which cannot be used
    max_step - stop after max_step step
    '''

    # data to keep tab
    bouncing = False
    pearl_level = None #later: pearl level of starting node
    steps_taken = 0
    # the goal
    d = multidigraph.graph['root']

    # start the walk:
    out_edge_list = multidigraph.nodes[starting_node]['R']['start']
    first_step = choose_first_available(out_edge_list, failure_set)
    #get ready for the general step
    prev_edge = first_step
    curr_node = first_step[1]
    steps_taken = 1
    #go
    while curr_node != d:
        out_edge_list = multidigraph.nodes[curr_node]['R']['basic forwarding'][prev_edge]
        print(out_edge_list)
        step_edge = choose_first_available(out_edge_list, failure_set)
        #move there
        curr_node = step_edge[1]
        prev_edge = step_edge

        steps_taken+=1
        if steps_taken>max_step: return 'Fail, maximum number of steps reached'
    return f'Success, {d} reached in {steps_taken} steps.'
        





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

    # for 2-node-connected components:
        # handle size 2 components
        # handle pearl roots - in edges in their root should be given permutations in the one bigger pearl

    return 0

#TODO
#routing of closed_pearls   #route_3_conn_graph() + route_a_pearl_for_traversing
def route_a_pearl(P, multigraph, multidigraph):
    '''
    Given a pearl, with boundaries, etc, it generates a arborescense based routing with route_3_conn_graph(), and 
    route_a_pearl_for_traversing
    INPUT: 
    P - pearl dictionary
    multigraph - nx.MultiGraph of the whole 2-conn-component
    multidigraph - nx.MultiDiGraph, to write the routing on, multigraph.to_directed()
    '''
    print('Here we start to route the pearl')
    #this is what we are working with
    pearl_graph = multigraph.subgraph(P['nodes']).copy()

    # add boundary_1, boundary_2 as an edge
    boundary_1 = P['boundary'][0]
    boundary_2 = P['boundary'][1]
    boundary_edge_key = pearl_graph.add_edge(boundary_1, boundary_2)
    P['boundary_edge_key'] = boundary_edge_key 

    # we need to get to one of the edges #here we have a freedom of choice
    pearl_graph.graph['root'] = boundary_1
    pearl_digraph = pearl_graph.to_directed() 

    # route the pearl with route_3_conn_graph
    pearl_digraph  = route_3_conn_graph(pearl_digraph, (boundary_1, boundary_2, boundary_edge_key))

    # route the pearl for traversing 
    # for cut edges to be in the graph, pearl_graph is not enough
    traversing_route_digraph = route_a_pearl_for_traversing(P, multidigraph, multigraph)

    # merging the two routing into one and finishing it into a complete one

    # if node is boundary_1: traversing + sending everything out on the cut edge, if that has failed
    # traverse back with bouncing, forwarding stuff to outside.
    # if node is boundary_2: change 'boundary_edge' to the appropriate cut edge
    # if node is else: getting out + traversing, if applicable 

    # for this, we gonna need to find edges sometimes
    all_edges_list = list(multidigraph.edges)

    for node in P['nodes']:
        if node == boundary_1:
            # traversing + sending everything out on the cut edge, if that has failed
            # traverse back with bouncing, forwarding stuff to outside.
            cut_edge  =  [edge for edge in all_edges_list if (edge[0]==boundary_1 and (edge[1] in P['neighbor_nodes']))][0]
            print(node)
            print(print(traversing_route_digraph.nodes[node]['R']))
            print(cut_edge)
            

            
        elif node == boundary_2:
            pass

         
        else:                                #if it is not a boundary getting out + traversing, if applicable 
            pass
        #    R = {}
        #    multidigraph.nodes[node]['R']

    
    '''
    for node in P['nodes']:
        print('node is ', node)
        if node in traversing_route_digraph.nodes:
            print(traversing_route_digraph.nodes[node]['R'])
        if node != pearl_graph.graph['root']:
            print(pearl_digraph.nodes[node]['R'])
    '''

#TODO: bouncing bit!
#routing of 3-connected graphs via arborescences
def route_3_conn_graph(multidigraph, edge):
    '''
    given a 3 connected graph, that has an extra attribute, g.graph['root'] and a specific 'edge', make a 2 resilient routing based on
    circular arborescences, where the first arborescence contains 'edge'
    '''

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
            R_basic_forwarding = {}
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
                    R_basic_forwarding[edge] = [out_edge_1, out_edge_2, out_edge_3]
                    #R_bouncing = #TODO
            #if the package starts at 'node'
            arb_0_edges_list = list(arborescence_list[0].edges)
            arb_1_edges_list = list(arborescence_list[1].edges)
            arb_2_edges_list = list(arborescence_list[2].edges)
            if node == '15':
                print(multidigraph.graph['root'])
                print('Debugging')
                print([edge for edge in arb_2_edges_list if (edge[0]=='15')])

            out_edge_in_arb_0 = [edge for edge in arb_0_edges_list if (edge[0]==node)][0]
            out_edge_in_arb_1 = [edge for edge in arb_1_edges_list if (edge[0]==node)][0]
            out_edge_in_arb_2 = [edge for edge in arb_2_edges_list if (edge[0]==node)][0]

            #saving the local dicts to the big dict
            R['basic forwarding'] = R_basic_forwarding
            R['start'] = [out_edge_in_arb_0, out_edge_in_arb_1, out_edge_in_arb_2]
            multidigraph.nodes[node]['R'] = R


    return multidigraph 

#TODO: bouncing bit! (here we actually need to think about it)
#routing of parallel routes in a graph
def route_a_pearl_for_traversing(P, multidigraph, multigraph):
    '''
    find 2 (or more) edge-disjoint paths between the boundaries of the pearl, then do the circular routing for 
    parallel paths
    INPUT:
    P - dictionary for pearl data
    multigraph     - multigraph, that has been partitioned into pearls, and P is one of them
    multidigraph   - multidigraph, thats edges will be used for partitioning

    OUTPUT: 
    - parallel_path_graph -  nx.MultiDiGraph, which is 2 parallel paths, a subraph of multidigraph and
    contains the routing info at 
    parallel_path_graph.nodes[node]['R']['traversing'] 
    '''
    #list the boundaries
    boundary_1 = P['boundary'][0]
    boundary_2 = P['boundary'][1]

    #create the pearl subgraph
    pearl_graph = multigraph.subgraph(P['nodes']).copy()


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

    #TODO is one of the cut_edges has failed and we are on which level?,  need to change the bouncing bit

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
    
    
    return parallel_path_graph

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
        


### RUNNING STUFF ########################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

# running tests of package sending
def package_send():
#def main(args):
    graph = nx.read_graphml('/mnt/d/Egyetem/Routing Cikk/fast-failover/pearl-algo/graph_sets/example_graphs/biggest_pearl.graphml')
    graph = contract_paths_keep_root(graph)
    for node in graph.nodes:
        if graph.degree(node)<3:
            print('Small degree node')
            print(node, graph.degree(node))
    
    graph.graph['k'] = 3
    digraph = graph.to_directed()

    # route it
    route_3_conn_graph(digraph, (1,0,0))

    #start the walk
    result = send_a_package(digraph, 3, [], 10)
    print(result)


#K4 package sending
def test_3_conn_package_sending():

    # take the one of the easiest possible graphs to try
    graph = nx.MultiGraph(nx.complete_graph(n=4))
    graph.graph['root'] = 0
    graph.graph['k'] = 3
    digraph = graph.to_directed()

    # route it
    route_3_conn_graph(digraph, (1,0,0))

    #start the walk
    result = send_a_package(digraph, 3, [], 10)
    print(result)

    
    
    




#def main(args):
def routing_tests():
    ex_multigraph = create_simple_2_pearl_multigraph()
    current_pearl_level, num_pearls, multi_graph, list_of_all_pearls = partition_2_conn_into_pearls(ex_multigraph)

    ex_multidigraph = ex_multigraph.to_directed()
    example_pearl = list_of_all_pearls[0][0]
    print(ex_multigraph.subgraph(example_pearl['nodes']).edges)
    print('Start of the pearl routing')
    route_a_pearl(example_pearl, ex_multigraph, ex_multidigraph)
        
#iterate through topology zoo and data to .csv files
def main(args):
#def pearl_depth_experiment():

    graph_data, subdivision_data = check_everything()
    graph_data.to_csv(path_or_buf='/mnt/d/Egyetem/Routing Cikk/fast-failover/pearl-algo/topology_zoo_statistics_with_pearl_sizes.csv')
    subdivision_data.to_csv(path_or_buf='/mnt/d/Egyetem/Routing Cikk/fast-failover/pearl-algo/topology_zoo_statistics_subdivision_stats.csv')

#def main(args):
def check_a_single_graph():
    pearl_depth, node_num, edge_num, num_pearls, avg_pearl_size, min_pearl_size, max_pearl_size, is_subdivision = check_a_graph('/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/Colt.graphml')  
    print(pearl_depth, node_num, edge_num, num_pearls)

#def main(args):
def separation_tests():
    file_path='/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/Aarnet.graphml'
    G, is_list = read_in_graph(file_path)

    Graph = nx.MultiGraph(G)
    Graph.graph['root'] = random.choice(list(Graph.nodes))
    Graph = contract_paths_keep_root(Graph)

    print(separate_by_cutting_nodes(Graph))

#def main(args):
def pearl_decomposition_and_cactus():
    graph = create_big_2_conn_from_Ion()

    tree = extract_pearl_tree(graph)
    print(tree)

    #Check those pearls
    depth, num_pearls, graph_with_pearl_data, list_of_all_pearls = partition_2_conn_into_pearls(graph)
    # this should be done with proper pearl depth
    pearl_depth = nx.diameter(extract_pearl_tree(graph))+1

    #creating the simple list, where they are not in levels, so we can iterate through them
    print(list_of_all_pearls[0][0])
    simple_pearl_list = []
    for level in list_of_all_pearls:
        for pearl in level:
            simple_pearl_list.append(pearl)

    

    #creating the simple list, where they are not in levels, so we can iterate through them
    simple_pearl_list = []
    for level in list_of_all_pearls:
        for pearl in level:
            simple_pearl_list.append(pearl)
        
    fixing_pearl_parent_pointers(simple_pearl_list, graph)


### UTILS, that are using functions from this file 
def max_sized_pearl_of_a_graph(filepath):
    '''
    Returns the maximum sized pearl as a 3 connected graph of 
    '''
    G, is_list = read_in_graph(filepath)
    OriginalGraph = nx.MultiGraph(G)


    list_of_all_pearls = []
    list_of_all_pearl_subgraphs = []

    #discard degree 2 nodes, so found pearls are not littered by them
    trimmed_graph = OriginalGraph.copy()
    trimmed_graph = contract_paths_keep_root(trimmed_graph)

    #for consistency lets assume, that the root is a degree 3 node:
    if 'root' not in trimmed_graph.graph.keys():                       #if not specified earlier, choose a random root
        trimmed_graph.graph['root'] = random.choice(list(trimmed_graph.nodes))

    #dissect into 2-connected components
    #comps = cut_cutting_edges(trimmed_graph)
    comps = separate_by_cutting_nodes(trimmed_graph)
    #discard single nodes and single edges and circles
    proper_comps = []
    for comp in comps:
        if len(list(comp.nodes))>2 and len(list(comp.edges))>len(list(comp.nodes)):
            proper_comps.append(comp)
    
    if proper_comps == []:             #if every 2 connected component is a single node, edge or circle, then it is a cactus
        return 1, 0, [], False


    #measure each component
    for comp in proper_comps:
        #get the graph structure we need
        graph = trimmed_graph.subgraph(comp).copy()

        # finding the new root
        # if the original root is not in the component, then it is the closest one to the root
        if trimmed_graph.graph['root'] not in comp:
            #we need to find the unique closest node to the root
            lista = list(comp)
            lista.sort(key = lambda node: len(nx.shortest_path(OriginalGraph, source=str(trimmed_graph.graph['root']), target=str(node))))
            closest_node = lista[0]
            graph.graph['root'] = closest_node

        #discard degree 2 nodes, so found pearls are not littered by them
        #if this is before finding the new root, that can make problems
        graph = contract_paths_keep_root(graph)

        #Check those pearls
        depth, num_pearls, graph_with_pearl_data, list_of_pearls_in_this_comp = partition_2_conn_into_pearls(graph)
        for level in list_of_pearls_in_this_comp:
            for pearl in level:
                subgraph = graph.subgraph(pearl['nodes']).copy()
                subgraph.add_edge(pearl['boundary'][0], pearl['boundary'][-1])
                list_of_all_pearl_subgraphs.append(subgraph)

                    
        list_of_all_pearls.append(list_of_pearls_in_this_comp)  # list of all pearls 
    
    #finding the biggest pearl
    max_pearl_size = 0
    max_pearl = None
    for pearl_graph in list_of_all_pearl_subgraphs:
        if len(list(pearl_graph.nodes))>max_pearl_size:
            max_pearl_size = len(list(pearl_graph.nodes))
            max_pearl  = pearl_graph
    

    #creating the subgraph
    return  max_pearl





if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--file_path', type=str, default="/mnt/d/Egyetem/Routing Cikk/SyPeR/topology-zoo-original/HiberniaGlobal.graphml",
                        help='Where is the graph')
    
    args = parser.parse_args()
    main(args)





#DESCRIPTION OF DATA STRUCTURES

## PEARL
#pearls are stored in a list of lists, where 

#PEARL[i]['nodes'] - list of nodes for pearl i
#PEARL[i]['cut_edges'] - the two edges defining pearl i
#PEARL[i]['boundary'] - the two nodes being the boundary of pearl i
#PEARL[i]['neighbor_nodes'] = neighbors

## PACKAGE

# package_info['bouncing'],
# package_info['pearl_level'] 
# package_info['num_steps']

## ROUTING TABLE

# A routing table in given node v looks like a multi levelled dict, keys are possible header configs, 
# items are the routing permutation on the out-edges

# R['incoming_edge']['in_pearl'] = 
# R['incoming_edge']['traversing'] = 

# Now we write it onto the nodes of the graph, so we don't need the R['v']['incoming_edge']['traversing'] notation.


#TODO
# simple failure experiments, path lengths, etc.
# running time experiments
# 3 failure experiments