from collections import defaultdict
# import networkx as nx
# import plotly.graph_objs as go
# from plotly.offline import plot

def find_homology(inputt):
    d = {}
    arr = []
    with open(inputt, 'r', encoding='utf-8') as inp:
        for line in inp:
            stripped_line = line.strip().split('\t')
            first_id = stripped_line[0]
            second_id = stripped_line[1]
            align_procent = float(stripped_line[2])
            if first_id != second_id:
                key1, key2 = f'{first_id} {second_id}', f'{second_id} {first_id}'
                if key1 not in d and key2 not in d and align_procent > 70.0:
                    d[key1] = None
                    d[key2] = None
                    arr.append(line)
    return arr


def cluster_homology(inputt):
    d = defaultdict(list)
    ar = find_homology(inputt)
    for line in ar:
        stripped_line = line.strip().split('\t')
        first_id = stripped_line[0]
        second_id = stripped_line[1]
        d[first_id].append(second_id)
    return d



def similarity(cluster1, cluster2):
    set1, set2 = set(cluster1), set(cluster2)
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    return len(intersection) / len(union) * 100

def merge_clusters(clusters):
    # First pass: merge clusters with 50% or more similarity
    i = 0
    while i < len(clusters):
        j = i + 1
        while j < len(clusters):
            if similarity(clusters[i], clusters[j]) >= 50.0:
                # Merge the clusters
                clusters[i] = list(set(clusters[i] + clusters[j]))
                del clusters[j]
            else:
                j += 1
        i += 1

    # Second pass: merge newly formed clusters with the rest
    merged = True
    while merged:
        merged = False
        i = 0
        while i < len(clusters) - 1:
            j = i + 1
            while j < len(clusters):
                if similarity(clusters[i], clusters[j]) >= 50.0:
                    # Merge the clusters
                    clusters[i] = list(set(clusters[i] + clusters[j]))
                    del clusters[j]
                    merged = True
                else:
                    j += 1
            i += 1

    return clusters


# def plot_clusters(clusters):
#     G = nx.Graph()
#
#     for idx, cluster in enumerate(clusters):
#         for node in cluster:
#             if node not in G:
#                 G.add_node(node)
#         for node1 in cluster:
#             for node2 in cluster:
#                 if node1 != node2:
#                     G.add_edge(node1, node2)
#
#     pos = nx.spring_layout(G)
#
#     edge_trace = go.Scatter(
#         x=[],
#         y=[],
#         line=dict(width=0.15, color='#000000'),
#         hoverinfo='none',
#         mode='lines')
#
#     for edge in G.edges():
#         x0, y0 = pos[edge[0]]
#         x1, y1 = pos[edge[1]]
#         edge_trace['x'] += tuple([x0, x1, None])
#         edge_trace['y'] += tuple([y0, y1, None])
#     node_trace = go.Scatter(
#         x=[],
#         y=[],
#         text=[],
#         mode='markers',
#         hoverinfo='text',
#         marker=dict(
#             showscale=True,
#             colorscale='Jet',
#             reversescale=True,
#             color=[],
#             size=10,
#             colorbar=dict(
#                 thickness=15,
#                 title='Num. of Connections',
#                 xanchor='left',
#                 titleside='right'
#             ),
#             line=dict(width=2)))
#
#     for node, adjacencies in enumerate(G.adjacency()):
#         node_trace['x'] += tuple([pos[adjacencies[0]][0]])
#         node_trace['y'] += tuple([pos[adjacencies[0]][1]])
#         node_trace['text'] += tuple([f'{adjacencies[0]} (connections: {len(adjacencies[1])})'])
#
#     for node, adjacencies in enumerate(G.adjacency()):
#         node_trace['marker']['color'] += tuple([len(adjacencies[1])])
#         node_info = f'{adjacencies[0]} (connections: {len(adjacencies[1])})'
#         node_trace['text'] += tuple([node_info])
#     fig = go.Figure(data=[edge_trace, node_trace],
#                     layout=go.Layout(
#                         title='<br>Exon Network Graph',
#                         titlefont=dict(size=16),
#                         showlegend=False,
#                         hovermode='closest',
#                         margin=dict(b=20, l=5, r=5, t=40),
#                         annotations=[dict(
#                             text="",
#                             showarrow=False,
#                             xref="paper", yref="paper",
#                             x=0.005, y=-0.002)],
#                         xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
#                         yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
#
#     plot(fig, filename='network.html')



if __name__ == "__main__":
    arra = []
    input_file_path = '/Users/cankale/Downloads/exon_homology_deneme_msa_input_tblastx'
    homology_data = cluster_homology(input_file_path)
    top_10_homology = sorted(homology_data.items(), key=lambda item: len(item[1]), reverse=True)[:10]
    for key, value_list in top_10_homology:
        homology_lists = [key] + value_list
        arra.append(homology_lists)
    merged_clusters = merge_clusters(arra)
    for cluster in merged_clusters:
        print(cluster)
