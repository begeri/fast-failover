import pandas as pd
import matplotlib.pyplot as plt
import tikzplotlib
import networkx as nx


from matplotlib.patches import FancyArrowPatch
from matplotlib.collections import LineCollection





#def plot_columns_by_index(dataframe, col1, col2, col3, col4, output_file="topology_zoo_pearl_data.tex"):
def plot_columns_by_index(dataframe, col1, col2, col3, output_file="topology_zoo_pearl_data.tex"):
    '''
    Given a dataframe, 
    '''
    # Sort the DataFrame by col1
    sorted_df = dataframe.sort_values(by=col1).reset_index(drop=True)

    # Plot the data
    plt.figure(figsize=(10, 6))
    plt.scatter(sorted_df.index, sorted_df[col1], label=col1, color='b', s=5)
    plt.scatter(sorted_df.index, sorted_df[col2], label=col2, color='r', s=5)
    plt.scatter(sorted_df.index, sorted_df[col3], label=col3, color='y', s=5)
    #plt.scatter(sorted_df.index, sorted_df[col4], label=col4, color='k', s=5)
    plt.xlabel('Index')
    plt.ylabel('Values')
    #plt.title(f'{col1}, {col2}, {col3} and {col4}')
    plt.title(f'{col1}, {col2} and {col3}')
    plt.legend()
    plt.grid()

    # Export the plot to a TikZ file
    tikzplotlib.save(output_file)

    plt.show()

#def plot_columns_by_an_other_column(dataframe, col_x, col1, col2, col3, output_file=None):
def plot_columns_by_an_other_column(dataframe, col_x, col1, output_file=None): 
    '''
    Given a dataframe, 
    '''
    # Sort the DataFrame by col1
    sorted_df = dataframe.sort_values(by=col1).reset_index(drop=True)

    # Plot the data
    plt.figure(figsize=(10, 6))
    plt.scatter(sorted_df[col_x], sorted_df[col1], label=col1, color='b', s=5)
    #plt.scatter(sorted_df[col_x], sorted_df[col2], label=col2, color='r', s=5)
    #plt.scatter(sorted_df[col_x], sorted_df[col3], label=col3, color='y', s=5)
    
    plt.xlabel('Index')
    plt.ylabel('Values')
    plt.title(f'{col1} vs {col_x}')
    plt.legend()
    plt.grid()

    # Export the plot to a TikZ file
    tikzplotlib.save(output_file)

    plt.show()
    


def draw_graph_with_labels(graph, output_tex_file):
    """
    Draws a graph using matplotlib, colors nodes based on their labels, and exports the figure to a .tex file.

    Parameters:
    - graph: networkx.Graph or networkx.MultiGraph
    - output_tex_file: Path to save the resulting .tex file
    """
    # Extract labels from node attributes
    labels = nx.get_node_attributes(graph, 'label')
    
    # Assign unique colors to each label
    unique_labels = list(set(labels.values()))
    color_map = {label: plt.cm.tab10(i % 10) for i, label in enumerate(unique_labels)}
    
    # Map colors to nodes
    node_colors = [color_map[labels[node]] for node in graph.nodes if node in labels]
    
    # Draw the graph
    pos = nx.spring_layout(graph)  # Positioning algorithm for node layout
    nx.draw(graph, pos, with_labels=True, node_color=node_colors, node_size=300, font_size=10, font_color='white')
    
    # Save the plot to a .tex file using tikzplotlib
    tikzplotlib.save(output_tex_file)
    print(f"Graph saved to {output_tex_file}")


def draw_multigraph_with_labels(graph, output_tex_file = None):
    """
    Draws a MultiGraph using matplotlib, including parallel edges, 
    colors nodes based on their labels, and exports the figure to a .tex file.

    Parameters:
    - graph: networkx.MultiGraph
    - output_tex_file: Path to save the resulting .tex file
    """
    # Extract labels from node attributes
    labels = nx.get_node_attributes(graph, 'pearl_id')
    
    # Assign unique colors to each label
    unique_labels = list(set(labels.values()))
    color_map = {label: plt.cm.tab10(i % 10) for i, label in enumerate(unique_labels)}
    
    # Map colors to nodes
    node_colors = [color_map[labels[node]] for node in graph.nodes if node in labels]
    
    # Node positions using spring layout
    pos = nx.spring_layout(graph)
    
    # Draw nodes and labels
    nx.draw_networkx_nodes(graph, pos, node_color=node_colors, node_size=500)
    nx.draw_networkx_labels(graph, pos, font_size=10, font_color='white')

    # Draw parallel edges manually
    ax = plt.gca()
    seen_edges = {}
    for u, v, *rest in graph.edges(data=True, keys=True):
        # Determine an offset for parallel edges
        edge = (min(u, v), max(u, v))
        if edge not in seen_edges:
            seen_edges[edge] = 0
        offset = seen_edges[edge] * 0.1
        seen_edges[edge] += 1

        # Compute edge midpoints for Bezier curves
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        ctrl_x, ctrl_y = (x1 + x2) / 2 + offset, (y1 + y2) / 2 + offset

        # Create a curved edge
        line = FancyArrowPatch((x1, y1), (x2, y2),
                               connectionstyle=f"arc3,rad={offset}",
                               arrowstyle="-", linewidth=1.5, color="gray")
        ax.add_patch(line)

    # Adjust plot limits
    ax.set_aspect('equal')
    ax.autoscale_view()

    # Save the plot to a .tex file using tikzplotlib
    if output_tex_file is not None:
        tikzplotlib.save(output_tex_file)
        print(f"MultiGraph saved to {output_tex_file}")


