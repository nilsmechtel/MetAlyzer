import networkx as nx
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import ListedColormap


def draw_network_py(nodes_df, edge_df, annotation_df, output, f_width, f_height, colbar_width, cmap, bg_color,
                    unassigned_color, node_size, font_size, font_color, edge_color, ann_font_size, ann_font_color):
                        
    cmap = ListedColormap(cmap)
    cmap.set_bad(color=unassigned_color)
    colbar_width = int(colbar_width)
    
    # max value for colormap
    max_v = np.ceil(np.nanmax(nodes_df['log2FC'].abs()))

    # build networkx
    nodes_list = [(nodes_df['Label'][row],
                   {'pos': (nodes_df['x'][row], nodes_df['y'][row])}) for row in nodes_df.index]
    edges_list = np.array([(edge_df['Node1'][row], edge_df['Node2'][row]) for row in edge_df.index])
    annotation_list = [(annotation_df['Annotation'][row],
                        {'pos': (annotation_df['x'][row], annotation_df['y'][row]),
                         }) for row in annotation_df.index]
                         

    G = nx.DiGraph()
    G.add_nodes_from(nodes_list)
    G.add_nodes_from(annotation_list)
    G.add_edges_from(edges_list)

    # create figure and define axes
    fig = plt.figure()
    fig.set_figheight(f_height)
    fig.set_figwidth(f_width)

    ax_network = plt.subplot2grid(shape=(5, 100), loc=(0, 0), rowspan=5, colspan=100 - colbar_width)
    ax_colorbar = plt.subplot2grid(shape=(5, 100), loc=(2, 100 - colbar_width), colspan=colbar_width)
    
    fig.set_facecolor(bg_color)
    ax_network.set_facecolor(bg_color)
    ax_network.axis('off')

    # draw network nodes
    for shape in nodes_df['Shape'].unique():
        shape_spec_nodes = nodes_df['Label'].values[nodes_df['Shape'] == shape]
        shape_spec_log2FC = nodes_df['log2FC'].values[nodes_df['Shape'] == shape]

        nx.draw_networkx_nodes(G, nx.get_node_attributes(G, 'pos'), ax=ax_network, nodelist=shape_spec_nodes,
                               node_size=node_size, node_shape=shape,
                               node_color=shape_spec_log2FC, cmap=cmap, vmin=-max_v, vmax=max_v)

    # draw network edges
    for rad in edge_df['Rad'].unique():
        rad_spec_edges = edges_list[edge_df['Rad'].values == rad]
        nx.draw_networkx_edges(G, nx.get_node_attributes(G, 'pos'), ax=ax_network, edgelist=rad_spec_edges,
                               edge_color=edge_color, arrowstyle='-', connectionstyle=f'arc3,rad={rad}')

    # draw network node labels
    nx.draw_networkx_labels(G, nx.get_node_attributes(G, 'pos'), ax=ax_network,
                            labels={label: re.sub('#[0-9]+$', '', label) for
                                    label, pos in nx.get_node_attributes(G, 'pos').items()
                                    if label in nodes_df['Label'].values},
                            font_size=font_size, font_color=font_color)

    # draw network annotations
    nx.draw_networkx_labels(G, nx.get_node_attributes(G, 'pos'), ax=ax_network,
                            labels={label: label for label, pos in nx.get_node_attributes(G, 'pos').items()
                                    if label in annotation_df['Annotation'].values},
                            font_size=ann_font_size, font_color=ann_font_color)

    # draw colorbar
    ColorbarBase(ax=ax_colorbar, cmap=cmap, values=np.arange(-max_v - 1e-5, max_v + 1e-5, 2 * max_v / 10000),
                 label='log2FC',
                 orientation='vertical')
    fig.tight_layout()

    # save figure
    fig.savefig(output)

