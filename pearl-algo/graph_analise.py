import networkx as nx
import argparse

import os

def read_in_graph(filepath):
    '''
    Reads in .graphml files to nx graphs
    '''
    G = nx.read_graphml(filepath)

    is_list = isinstance(G, list)

    return G, is_list

#measure graph submethods
def find_pearls(multi_graph):
    '''
    Generates the pearls of a graph using the Gomory-Hu tree.
    These are the leafs of pearls arborescence
    Input:nx.MultiGraph graph
    '''
    #print("Start of pearl finding")
    #find pearls
    simple_graph = nx.Graph(multi_graph)               #gomory-hu tree is not implemented for MultiGraphs
    for edge in simple_graph.edges:
        i = edge[0]
        j = edge[1]
        simple_graph[i][j]['capacity'] = multi_graph.number_of_edges(i,j)

    tree = nx.gomory_hu_tree(simple_graph)
    #list 2-cuts, that should be removed
    important_two_cuts = []
    for cut in tree.edges:
        i = cut[0]
        j = cut[1]
        #print(i, j, tree[i][j]['weight'])
        if tree[i][j]['weight'] == 2:
            important_two_cuts.append(cut)

    #print('important two cuts', important_two_cuts)
    #remove 2-cuts from gomory-hu tree
    tree.remove_edges_from(important_two_cuts)
    #these components are the pearls
    comps = list(nx.connected_components(tree))
    #print('comps ',comps)

    #keep the real pearls
    pearls = []
    for comp in comps:
        #print('this comp ',comp)
        comp = list(comp)
        neighbors = list(neighbor_of_subgraph(multi_graph, comp))
        side_of_pearl = []
        for node in comp:
            if node in neighbor_of_subgraph(multi_graph, node_set=neighbors):
                side_of_pearl.append(node)
        #print('side of pearls ',side_of_pearl)
        cut_size = 0
        for x in side_of_pearl:
            for v in neighbors:
                cut_size += multi_graph.number_of_edges(x,v)
        if cut_size==2:                           #between comp and the rest of the graph there is actually only 2 edges
            # copy out the pearl
            pearl = multi_graph.subgraph(comp).copy()
            # connect the edges
            if len(side_of_pearl)==2:
                pearl.add_edge(side_of_pearl[0], side_of_pearl[1])
            # check if 3-connected 
            if is_3_connected(pearl):
                pearls.append(comp)

    #print('pearls ', list(pearls))
    return pearls


def neighbor_of_subgraph(graph, node_set):
    '''
    
    '''
    # Extract the subgraph
    subgraph = graph.subgraph(node_set)

    # Initialize a set to store the neighbors of the subgraph
    subgraph_neighbors = set()

    # Loop through each node in the subgraph
    for node in subgraph.nodes():
        # Get the neighbors of the node in the original graph
        neighbors = graph.neighbors(node)
        # Add neighbors that are not part of the subgraph
        subgraph_neighbors.update(n for n in neighbors if n not in subgraph.nodes())
    return subgraph_neighbors

def contract_pearl(multi_graph, comp):
    '''
    nomen est omen
    multi_graph - the graph we are working on
    comp - the nodes of the pearl to contract

    We keep one node of the original ones.
    '''
    #contract pearls
    #print('comp to remove: ', comp)
    staying_node = comp[0]
    #print('node to stay: ', staying_node)
    neighbors = list(neighbor_of_subgraph(multi_graph, comp))
    #print(neighbors[0], staying_node)
    
    multi_graph.remove_nodes_from(comp)
    multi_graph.add_node(staying_node)
    #print('szomszédok ', neighbors)
    if len(neighbors)==2:
        multi_graph.add_edge(neighbors[0],staying_node)
        multi_graph.add_edge(neighbors[1],staying_node)
    else:
        multi_graph.add_edge(neighbors[0],staying_node)
        multi_graph.add_edge(neighbors[0],staying_node)

    return multi_graph

def contract_paths(G):
    '''
    Nomen est omen, no 2 degree nodes
    Works only with MultiGraph
    '''
    nodes_to_remove = [v for v in G.nodes if G.degree(v) == 2]

    for v in nodes_to_remove:
        neighbors = list(G.neighbors(v))
        G.remove_node(v)
        if len(neighbors)==2:
            G.add_edge(neighbors[0], neighbors[1])
        else: G.add_edge(neighbors[0], neighbors[0])

    return G

def cut_cutting_edges(multi_graph):
    '''
    dissects the graph into 2-connected components
    '''
    graph = nx.Graph(multi_graph)               #gomory-hu tree is not implemented for MultiGraphs
    for edge in graph.edges:
        i = edge[0]
        j = edge[1]
        graph[i][j]['capacity'] = multi_graph.number_of_edges(i,j)

    #give capacities to edges, so that the algo understands - technicality
    for edge in graph.edges:
        i = edge[0]
        j = edge[1]
        graph[i][j]['capacity'] = 1
    #create gomory-hu tree    
    tree = nx.gomory_hu_tree(graph)

    #list cutting edges
    cutting_edges = []
    for cut in tree.edges:
        i = cut[0]
        j = cut[1]
        if tree[i][j]['weight'] == 1:
            cutting_edges.append(cut)

    #remove cutting edges
    graph.remove_edges_from(cutting_edges)
    #return the remaining connected components
    comps = nx.connected_components(graph)
    return comps

def one_round_of_contracts(multi_graph):
    #contract the pathes, so there are pearls
    multi_graph = contract_paths(multi_graph)

    if len(multi_graph.nodes)<3:
        print('ja') 
        return multi_graph
    
    #find pearls
    comps = list(find_pearls(multi_graph))
    #contract pearls
    for comp in comps:
        comp = list(comp)
        multi_graph = contract_pearl(multi_graph, comp)

    return multi_graph

def is_3_connected(multi_graph):
    """
    We check the if the smallest cut is smaller than 3.
    """
    if len(multi_graph.nodes)<=1: 
        #print('ezt még átnézni')
        return True

    simple_graph = nx.Graph(multi_graph)               #gomory-hu tree is not implemented for MultiGraphs
    for edge in simple_graph.edges:
        i = edge[0]
        j = edge[1]
        simple_graph[i][j]['capacity'] = multi_graph.number_of_edges(i,j)

    tree = nx.gomory_hu_tree(simple_graph)
    #check through the cuts
    for cut in tree.edges:
        i = cut[0]
        j = cut[1]
        #print(i, j, tree[i][j]['weight'])
        if tree[i][j]['weight'] <= 2: return False

    return True

def repeat_contracts(multi_graph, verbose= False):
    """
    Repeats contracts until multi_graph is 3-connected
    """
    number_of_contracts = 0
    #print(multi_graph, number_of_contracts)
    while is_3_connected(multi_graph)==False:
        multi_graph = one_round_of_contracts(multi_graph)
        number_of_contracts+=1
        #print('num conn components ',nx.number_connected_components(multi_graph))
        #print('This is the ', number_of_contracts, '. contract')
        #print(multi_graph)
        if verbose: 
            print(number_of_contracts, multi_graph)
            print(multi_graph.edges)
        

    return number_of_contracts

def measure_graph(file_path, verbose=False):
    #read in graph 
    Graph, is_list = read_in_graph(file_path)
    #print(graph)

    #dissect into 2-connected components
    comps = cut_cutting_edges(Graph)
    #discard single nodes
    proper_comps = [comp for comp in comps if len(comp)>1]

    max_depth = 0
    #measure each component
    for comp in proper_comps:
        graph = Graph.subgraph(comp).copy()
        #create multi-graph, so connectedness is not lost    
        multi_graph = nx.MultiGraph(graph)
        depth_of_comp = repeat_contracts(multi_graph, verbose)
        if depth_of_comp>max_depth: max_depth=depth_of_comp

    return max_depth
  
def measure_everything():
    fnames = [ 'topology-zoo-original/'+fname for fname in os.listdir('topology-zoo-original/') if fname.endswith('.graphml')] 
    fnames = filter(lambda f: nx.edge_connectivity(nx.read_graphml(f))>0,fnames)


    deepest = ''
    depth = 0
    for name in fnames:
        try:
            print(name) 
            this_depth = measure_graph(name)
            print(this_depth)
            if this_depth>=depth:
                deepest = name
                depth = this_depth

        
        except:
            print('baj ', name)
    return deepest, depth

def main(args):
    #print(measure_graph(args.file_path, verbose=False))
    print(measure_everything())

    #graph = nx.MultiGraph(nx.cycle_graph(5))
    #print(contract_paths(graph).edges)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--file_path', type=str, default="topology-zoo-original/Bellcanada.graphml",
                        help='Where is the graph')
    
    args = parser.parse_args()
    main(args)