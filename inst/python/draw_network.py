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
    #ax_colorbar = plt.subplot2grid(shape=(5, 100), loc=(2, 100 - colbar_width), colspan=colbar_width)
    
    fig.set_facecolor(bg_color)
    ax_network.set_facecolor(bg_color)
    ax_network.axis('off')
    
    #defining the discrete colors according to the log2fc value
    color = []
    for val in nodes_df['log2FC']:
        if val <= -3 :    color.append('#0066CC')
        elif -3 < val < -1.5:   color.append('#3399FF')
        elif -1.5 < val < -0.5:  color.append('#00CCCC')
        elif 0.5 < val < 1.5:  color.append('#FF9999')
        elif 1.5 < val < 3:  color.append('#FF6666')
        elif 3 < val:  color.append('#FF0000')
        elif -0.5 < val < 0.5:  color.append('#A0A0A0')
    nodes_df['color'] = color
    color.append('#0066CC')
    color.append('#3399FF')
    color.append('#00CCCC')
    color.append('#FF9999')
    color.append('#FF6666')
    color.append('#FF0000')
    color.append('#A0A0A0')
        
    # draw network nodes
    for shape in nodes_df['Shape'].unique():
        shape_spec_nodes = nodes_df['Label'].values[nodes_df['Shape'] == shape]
        shape_spec_log2FC = nodes_df['log2FC'].values[nodes_df['Shape'] == shape]
        shape_spec_color = nodes_df['color'].values[nodes_df['Shape'] == shape]
        print(shape_spec_color)
        print(shape_spec_nodes)
        #print(nodes_df)
        nx.draw_networkx_nodes(G, nx.get_node_attributes(G, 'pos'), ax=ax_network, nodelist=shape_spec_nodes,
                               node_size=node_size, node_shape=shape,
                               node_color=shape_spec_color)
                               
    #defining and drawing legend nodes, attributes and labels                          
    leg_nodes = ['FC > 3', '1.5 < FC < 3', '0.5 < FC < 1.5', '-1.5 < FC < -0.5','-3 < FC < -1.5','FC < -3', 'No fold change']
    
    leg_position = {'FC > 3': (2100, 1500),'1.5 < FC < 3': (2100, 1400),'0.5 < FC < 1.5': (2100, 1300),
                '-1.5 < FC < -0.5': (2100, 1100),'-3 < FC < -1.5': (2100, 1000),'FC < -3': (2100, 900),
            'No fold change': (2100, 1200)}
            
    leg_color = ['#FF0000', '#FF6666', '#FF9999', '#00CCCC', '#3399FF', '#0066CC', '#A0A0A0']
    
    nx.draw_networkx_nodes(G, pos = leg_position, ax=ax_network, nodelist=leg_nodes,
                           node_size=node_size, node_shape='o',
                           node_color=leg_color)
    nx.draw_networkx_labels(G, pos = leg_position, labels={n: n for n in leg_nodes},font_size=15, font_color='k', font_family='sans-serif', font_weight='bold')

    
    
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
                            font_size=font_size, font_color=font_color, font_weight='bold')

    # draw network annotations
    nx.draw_networkx_labels(G, nx.get_node_attributes(G, 'pos'), ax=ax_network,
                            labels={label: label for label, pos in nx.get_node_attributes(G, 'pos').items()
                                    if label in annotation_df['Annotation'].values},
                            font_size=ann_font_size, font_color=ann_font_color)

    fig.tight_layout()

    # save figure
    fig.savefig(output)

