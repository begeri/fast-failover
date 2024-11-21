import networkx as nx
import argparse
import matplotlib.pyplot as plt



def add_ear(n, graph, node_1, node_2):
    '''
    n \in {0,1,2,...}

    'node_1' and 'node_2' are nodes of 'graph'. Adds an ear with n inner nodes to 'graph',
    starting from 'node_1' and 'node_2'. 
    '''
    if node_1 not in graph.nodes:
        raise Exception("Legyen már benne a node_1 csúcs a gráfban...")
    if node_2 not in graph.nodes:
        raise Exception("Legyen már benne a node_2 csúcs a gráfban...")

    num_nodes = len(list(graph.nodes))
    new_nodes = [num_nodes + i for i in range(n)]
    #print(new_nodes)
    graph.add_nodes_from(new_nodes)
    for i in range(n-1):
        graph.add_edge(num_nodes+i, num_nodes+i+1)
    
    graph.add_edge(node_1, new_nodes[0])
    graph.add_edge(new_nodes[-1], node_2)

def write_graph(G, destination):
    f = open(str(destination), "w")
    for line in nx.generate_graphml(G):
        f.write(line)
        f.write('\n')

def read_in_graphml_graph(filepath):
    '''
    Reads in .graphml files to nx graphs
    '''
    G = nx.read_graphml(filepath)

    is_list = isinstance(G, list)

    return G, is_list

def plot_and_save_graph(graph, filepath):
    '''
    Plots the graph with edge attribute
    '''
    # Check if the graph has 'weight' attribute on edges
    has_weight = nx.get_edge_attributes(graph, 'weight')

    # Draw the graph
    pos = nx.spring_layout(graph)  # positions for all nodes
    plt.figure(figsize=(8, 6))
    nx.draw(graph, pos, with_labels=True, node_color='skyblue', node_size=2000, edge_color='gray', font_size=15, font_color='black')

    # If the graph has 'weight' attribute, draw edge labels
    if has_weight:
        edge_labels = nx.get_edge_attributes(graph, 'weight')
        nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels, font_color='red')

    # Save the plot to a file
    plt.savefig(filepath)
    plt.close()

def main(args):
    
    G = nx.Graph()
    G.add_nodes_from([0,1,2,3,4,5,6,7,8,9,10,11])
    G.add_edge(0,1)
    G.add_edge(1,2)
    G.add_edge(2,3)
    G.add_edge(2,5)
    G.add_edge(3,4)
    G.add_edge(4,5)
    G.add_edge(4,6)
    G.add_edge(6,7)
    G.add_edge(7,8)
    G.add_edge(7,10)
    G.add_edge(8,9)
    G.add_edge(9,10)
    G.add_edge(9,11)
    G.add_edge(11,0)

    

    print(list(G.nodes))
    print(list(G.edges))
    write_graph(G, 'topology-zoo/counter_try.graphml')

    G, #is_list = read_in_graphml_graph(args.filepath)
    '''
    graph = G
    for edge in graph.edges:
        i = edge[0]
        j = edge[1]
        graph[i][j]['capacity'] = 1
    tree = nx.gomory_hu_tree(graph)

    for cut in tree.edges:
        i = cut[0]
        j = cut[1]
        print(tree[i][j]['weight'])

    for cut in tree.edges:
        tree_copy = tree.copy()
        tree_copy.remove_edge(*cut)
        comps = nx.connected_components(tree_copy)
        for c in comps:
            print(c)
    '''

    
    


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--instance_path', type=str, default="topology-zoo/example_graph.graphml",
                        help='Where is the graph')
    
    args = parser.parse_args()
    main(args)










