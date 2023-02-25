import networkx as nx
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import ListedColormap


def draw_network_py(
    nodes_df,
    edge_df,
    annotation_df,
    output,
    f_width,
    f_height,
    dpi,
    bg_color,
    node_size,
    font_size,
    font_color,
    edge_color,
    ann_font_size,
    ann_font_color,
    font_family,
):

    # build networkx
    nodes_list = [
        (nodes_df["Label"][row],
        {"pos": (nodes_df["x"][row], nodes_df["y"][row])})
        for row in nodes_df.index
    ]
    annotation_list = [
        (annotation_df["Annotation"][row],
        {"pos": (annotation_df["x"][row], annotation_df["y"][row])})
        for row in annotation_df.index
    ]
    edges_list = np.array(
        [(edge_df["Node1"][row],
        edge_df["Node2"][row])
        for row in edge_df.index]
    )
    
    G = nx.DiGraph()
    G.add_nodes_from(nodes_list)
    G.add_nodes_from(annotation_list)
    G.add_edges_from(edges_list)

    # create figure
    fig, ax = plt.subplots(figsize=(f_width, f_height))
    ax.set_facecolor(bg_color)
    ax.axis('off')

    # draw network nodes
    for shape in nodes_df["Shape"].unique():
        shape_spec_nodes = nodes_df["Label"].values[nodes_df["Shape"] == shape]
        shape_spec_color = nodes_df["Color"].values[nodes_df["Shape"] == shape]
        nx.draw_networkx_nodes(
            G,
            nx.get_node_attributes(G, "pos"),
            nodelist=shape_spec_nodes,
            node_size=node_size,
            node_shape=shape,
            node_color=shape_spec_color,
        )

    # draw network edges
    for rad in edge_df["Rad"].unique():
        rad_spec_edges = edges_list[edge_df["Rad"].values == rad]
        nx.draw_networkx_edges(
            G,
            nx.get_node_attributes(G, "pos"),
            edgelist=rad_spec_edges,
            edge_color=edge_color,
            arrowstyle="-",
            connectionstyle=f"arc3,rad={rad}",
        )

    # draw network node labels
    nx.draw_networkx_labels(
        G,
        nx.get_node_attributes(G, "pos"),
        labels={
            label: re.sub("#[0-9]+$", "", label)
            for label, pos in nx.get_node_attributes(G, "pos").items()
            if label in nodes_df["Label"].values
        },
        font_size=font_size,
        font_color=font_color,
        font_weight="bold",
        font_family=font_family,
    )
    
    # draw network annotations
    nx.draw_networkx_labels(
        G,
        nx.get_node_attributes(G, "pos"),
        labels={
            label: label
            for label, pos in nx.get_node_attributes(G, "pos").items()
            if label in annotation_df["Annotation"].values
        },
        font_size=ann_font_size,
        font_color=ann_font_color,
        font_family=font_family,
    )

    # save figure
    fig.savefig(output, dpi=dpi)
