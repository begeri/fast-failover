import networkx as nx
import itertools


import networkx as nx

def get_faces(embedding):
    """
    Given a planar embedding (a dict-like structure with cyclic orders), 
    return a list of faces. Each face is a list of vertices (in order along its boundary).
    """
    faces = []
    visited = set()  # keep track of visited directed edges (u, v)
    # Iterate over every directed edge in the embedding.
    for u in embedding:
        for v in embedding[u]:
            if (u, v) in visited:
                continue
            face = []
            current_edge = (u, v)
            # Walk around the face until we return to the starting edge.
            while current_edge not in visited:
                visited.add(current_edge)
                cur_u, cur_v = current_edge
                face.append(cur_u)
                # Convert the AtlasView to a list for proper indexing.
                nbrs = list(embedding[cur_v])
                idx = nbrs.index(cur_u)
                next_idx = (idx + 1) % len(nbrs)
                next_u = cur_v
                next_v = nbrs[next_idx]
                current_edge = (next_u, next_v)
            faces.append(face)
    return faces

def face_edges(face):
    """
    Given a face (list of vertices in order), return its set of (undirected) edges.
    Each edge is represented as a frozenset of its two endpoints.
    """
    edges = set()
    n = len(face)
    for i in range(n):
        a = face[i]
        b = face[(i+1) % n]
        edges.add(frozenset((a, b)))
    return edges

def share_edge(face1, face2):
    """
    Returns True if the two faces share at least one edge.
    """
    return len(face_edges(face1).intersection(face_edges(face2))) > 0


def check_component_outerplanarity(G_component):
    """
    Checks if the connected component (G_component) is planar and, if so,
    whether it has a face that contains all nodes (outerplanar under the given definition).
    Returns a tuple: (is_planar, is_outerplanar, embedding, faces)
    """
    is_planar, embedding = nx.check_planarity(G_component)
    if not is_planar:
        return (False, False, None, [])
    
    # Retrieve all faces from the planar embedding.
    faces = list(get_faces(embedding))
    # The set of all nodes in the component.
    all_nodes = set(G_component.nodes())
    
    # Check if any face contains all the nodes.
    is_outerplanar = any(all_nodes.issubset(set(face)) for face in faces)
    
    return (True, is_outerplanar, embedding, faces)


def check_for_outerplanarity(G):
    """
    Checks if every connected component in the graph G is outerplanar.
    Returns True if every component is outerplanar, otherwise False.
    """
    # Iterate through each connected component in G.
    for component in nx.connected_components(G):
        subG = G.subgraph(component)
        is_planar, is_outerplanar, _, _ = check_component_outerplanarity(subG)
        # If any component is either non-planar or not outerplanar, return False.
        if not (is_planar and is_outerplanar):
            return False
    return True


# === Example Usage ===
if __name__ == "__main__":
    # Outerplanar: A cycle graph on 4 vertices.
    G1 = nx.cycle_graph(4)
    print("Cycle graph outerplanar?", check_for_outerplanarity(G1))  # Expected: True

    # Not outerplanar: K4 (planar but not outerplanar).
    G2 = nx.complete_graph(4)
    print("K4 outerplanar?", check_for_outerplanarity(G2))  # Expected: False

    # Not outerplanar: K2,3.
    G3 = nx.complete_bipartite_graph(2, 3)
    print("K2,3 outerplanar?", check_for_outerplanarity(G3))  # Expected: False
