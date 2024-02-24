# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 17:18:58 2023
@author: camlo
"""
import os.path
import sys
import hbonds_saltbridges
import Bio.PDB.PDBParser
from plotly.offline import init_notebook_mode, iplot
import plotly.graph_objects as go
import pandas as pd
import glob
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import math
import logomaker
import numpy as np
import helper_functions
from Bio import *
import plotly.offline as py
import seaborn as sns
import matplotlib.colors as mcolors
from sklearn.cluster import DBSCAN
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
aa_mapping = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLN': 'Q',
    'GLU': 'E',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'
}
amino_acid_classes = {
    'G': 'Hydrophobic',
    'A': 'Hydrophobic',
    'V': 'Hydrophobic',
    'C': 'Hydrophobic',
    'P': 'Hydrophobic',
    'L': 'Hydrophobic',
    'I': 'Hydrophobic',
    'M': 'Hydrophobic',
    'W': 'Aromatic',
    'F': 'Aromatic',
    'S': 'Hydrophilic',
    'T': 'Hydrophilic',
    'Y': 'Aromatic',
    'N': 'Hydrophilic',
    'Q': 'Hydrophilic',
    'D': 'Charged interactions',
    'E': 'Charged interactions',
    'K': 'Charged interactions',
    'R': 'Charged interactions',
    'H': 'Charged interactions'
}

def create_3d_graph(df1, df2, is_ligand, ligand_bonds={}, name_file="3d_scatter"):
    # Get XYZ positions from the DataFrame columns
    x1, y1, z1 = df1['X'], df1['Y'], df1['Z']
    x2, y2, z2 = df2['X'], df2['Y'], df2['Z']
    color_shapely = df1['shapely'].values.tolist()
    color_polar = df1['polar'].values.tolist()
    color_hb = df1['H-bond'].values.tolist()
    color_cluster = df1["cluster_color"].values.tolist()
    cluster_names = df1["cluster_index"].values.tolist()
    if is_ligand:
        names = df2['atom_name'].values.tolist()
        color_df2 = df2["color"].values.tolist()
        size2 = 15
    else:
        names = df2['residue_index'].values.tolist()
        color_df2 = "black"
        size2 = 15
        color_salt = df1['Salt Bridge'].values.tolist()
    # init_notebook_mode(connected=True)
    scatter_trace1 = go.Scatter3d(
        x=x1,
        y=y1,
        z=z1,
        mode='markers',
        marker=dict(
            size=9,
            color=color_shapely,
            opacity=1,
            line=dict(color='black', width=2)
        ),
        text=df1['AA'],
        hoverinfo='text',
        hoverlabel=dict(bgcolor='yellow', bordercolor='black'),
        name="Binding Residues"
    )
    scatter_trace2 = go.Scatter3d(
        x=x2,
        y=y2,
        z=z2,
        mode='markers',
        marker=dict(
            size=size2,
            color=color_df2,
            opacity=1,
            line=dict(color='white', width=5)
        ),
        text=names,
        hoverinfo='text',
        hoverlabel=dict(bgcolor='gray', bordercolor='white'),
        name="Target"

    )
    buttons = []
    buttons.append(dict(label='Shapely Colours', method='restyle',
                        args=[{'text': [df1['AA']], "marker.color": [color_shapely]}, [0]]))
    buttons.append(
        dict(label='Amino Colours', method='restyle', args=[{'text': [df1['AA']], 'marker.color': [color_polar]}, [0]]))
    buttons.append(
        dict(label='H-Bond', method='restyle', args=[{'text': [df1['AA']], 'marker.color': [color_hb]}, [0]]))
    if not is_ligand:
        buttons.append(dict(label='Salt-Bridges', method='restyle',
                            args=[{'text': [df1['AA']], 'marker.color': [color_salt]}, [0]]))
    buttons.append(dict(label='Hotspots', method='restyle',
                        args=[{'text': [cluster_names], 'marker.color': [color_cluster]}, [0]]))
    updatemenus = [
        dict(buttons=buttons, showactive=True),
        dict(direction='down', x=0.1, xanchor='left', y=1.1, yanchor='top'),
    ]
    # Create the layout
    layout = go.Layout(
        scene=dict(
            xaxis=dict(title='X'),
            yaxis=dict(title='Y'),
            zaxis=dict(title='Z'),
            camera=dict(
                eye=dict(x=2, y=-2, z=1.5)  # Adjust the eye position to view all eight regions
            )
        ),
        title='MAGPIE'
    )
    graphs = [scatter_trace1, scatter_trace2]
    if not is_ligand:
        line_trace_target = go.Scatter3d(
            x=x2,
            y=y2,
            z=z2,
            mode='lines',
            line=dict(
                color='black',  # Choose a color that stands out
                width=8  # Adjust line width as needed
            ),
            hoverinfo='skip',  # Optionally disable hover info for the lines
            showlegend=False  # Optionally hide the line trace from the legend
        )
        graphs.append(line_trace_target)
    else:
        for bond in ligand_bonds:
            for pair in ligand_bonds[bond]:
                atom_coords1 = df2.loc[df2['atom_serial_number'] == str(bond)]
                atom_coords2 = df2.loc[df2['atom_serial_number'] == str(pair)]
                point1 = atom_coords1[['X', 'Y', 'Z']].values
                point2 = atom_coords2[['X', 'Y', 'Z']].values
                line_trace = go.Scatter3d(
                    x=[point1[0][0], point2[0][0]],
                    y=[point1[0][1], point2[0][1]],
                    z=[point1[0][2], point2[0][2]],
                    mode='lines',
                    line=dict(color='black', width=8),
                    hoverinfo='skip',
                    showlegend=False
                )
                graphs.append(line_trace)
    # Create the figure and add the traces
    fig = go.Figure(data=graphs, layout=layout)
    fig.update_layout(updatemenus=updatemenus)
    # Show the interactive plot
    py.plot(fig, filename=f'{name_file}.html')


def find_nearest_points(target, binders, radius):
    target_points_in_contact = []
    binder_points_in_contact = []
    for i in range(0, len(target)):
        # Create KDTree for binder_points
        # Query binders DataFrame within the specified radius of target_point
        current_target = target.iloc[i]
        current_binder = binders.iloc[i]
        pd.set_option('display.max_colwidth', None)
        binder_points = []
        target_points = []
        for row in current_binder:
            if row is not None and len(row.keys()) != 0:
                binder_points.append([row['X'], row['Y'], row['Z']])
        for row in current_target:
            if row is not None and len(row.keys()) != 0:
                target_points.append([row['X'], row['Y'], row['Z']])
        tree = KDTree(np.array(binder_points))
        indices = tree.query_ball_point(np.array(target_points), r=radius)
        targets_in_contact = []
        binder_in_contact = []
        for i, index in enumerate(indices):
            if len(index) == 0:
                continue
            targets_in_contact.append(i)
            binder_in_contact += index
        target_points_in_contact.append(sorted(list(set(targets_in_contact))))
        binder_points_in_contact.append(sorted(list(set(binder_in_contact))))
    return target_points_in_contact, binder_points_in_contact


def get_coords_from_row(row):
    return [row['X'], row['Y'], row['Z']]


def calculate_arrow_position(subplot_index):
    arrow_x = (subplot_index - 1) * 0.25 + 0.1
    arrow_tail_x = arrow_x + 0.1
    return arrow_x, arrow_tail_x


def euclidean_distance(x1, y1, z1, x2, y2, z2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)


def transform_to_1_letter_code(amino_acids_3_letter):
    # Mapping dictionary for 3-letter to 1-letter code
    amino_acids_1_letter = [aa_mapping[aa] for aa in amino_acids_3_letter]
    return amino_acids_1_letter


def find_points_within_radius(binder_point, target_points, radius):
    # Convert the binder_point DataFrame to a NumPy array
    binder_np = np.array(binder_point[['X', 'Y', 'Z']])
    # Convert the target_points DataFrame to a NumPy array
    target_np = np.array(target_points[['X', 'Y', 'Z']])
    # Build a KDTree from the binder point
    tree = KDTree(binder_np)
    # Query the KDTree to find indices of points within the specified radius
    indices = tree.query_ball_point(target_np, r=radius)
    # Create a list to store the points within radius for each target point
    points_within_radius = []
    # Iterate through the indices and extract the corresponding target points
    for i, index_list in enumerate(indices):
        if len(index_list) == 0:
            continue
        points_within_radius.append(i)
    residues_points = target_points.iloc[points_within_radius]
    return residues_points


def find_closest_point(df1, df2, is_ligand):
    # Filter df2 to exclude rows where "atom_name" contains "H" if is_ligand is True
    if is_ligand:
        df2 = df2[~df2['atom_name'].str.contains("H", na=False)]
    # Calculate the average position in df1
    avg_position = df1[['X', 'Y', 'Z']].mean().values
    # Convert df2 to a NumPy array
    df2_np = np.array(df2[['X', 'Y', 'Z']])
    # Build a KDTree from df2
    tree = KDTree(df2_np)
    # Query the KDTree to find the index of the closest point in df2 to the average position in df1
    closest_index = tree.query(avg_position)[1]
    # Get the closest point from df2
    closest_point = df2.iloc[closest_index]
    return closest_point


def calculate_frequency(character, lst):
    count = lst.count(character)
    frequency = count / len(lst)
    return frequency


def calculate_bits(list_of_AA, sequence_list):
    list_of_frequencies = []
    for AA in list_of_AA:
        list_of_frequencies.append(calculate_frequency(AA, sequence_list))
    heights = []
    for f in list_of_frequencies:
        heights.append(np.abs(f * 100))
    return heights


def remove_items(test_list, item):
    # remove the item for all its occurrences
    c = test_list.count(item)
    for i in range(c):
        test_list.remove(item)
    return test_list


def create_sequence_logo_list(data, only_combined, is_ligand):
    # Titles for each type of graph
    axis_label_fontsize = 30
    title_fontsize = 36
    xtick_label_fontsize = 30
    y_lable_size = 34
    titles = ["Residues in Contact", "H-Bonds", "Salt Bridges"]
    for j, row in enumerate(data):
        # Filter out the [None, None] graphs
        valid_graphs = []
        title_row = []
        for i, graph in enumerate(row):
            if graph[1] != None:
                valid_graphs.append(graph)
                title_row.append(titles[i])
        num_graphs = len(valid_graphs)
        # Only create subplots for the valid graphs
        if num_graphs > 0:
            if j == 0:
                for i, graph in enumerate(valid_graphs):
                    fig, ax = plt.subplots(figsize=(2 * len(graph[1]), 12))  # Adjust the size as needed
                    # Create the logo with logomaker
                    logo = logomaker.Logo(graph[0], ax=ax, color_scheme='NajafabadiEtAl2017', shade_below=0.5)
                    logo.ax.set_ylabel('Frequency', fontsize=y_lable_size)
                    # Set x-ticks and labels
                    positions = [k for k in range(len(graph[1]))]
                    logo.ax.set_xticklabels(graph[1], fontsize=xtick_label_fontsize)  # Rotate labels for visibility
                    logo.ax.set_xticks(positions)
                    # Set the title for this subplot
                    ax.set_title(title_row[i], fontsize=title_fontsize)
                    plt.tight_layout()
                    plt.show()
                if only_combined:
                    break
                else:
                    continue
            fig, axs = plt.subplots(1, num_graphs, figsize=(5 * num_graphs, 6))
            fig.subplots_adjust(hspace=1.0)  # Increase the height space between rows
            if num_graphs == 1:
                axs = [axs]  # Make sure axs is iterable
            for i, graph in enumerate(valid_graphs):
                ax = axs[i]  # Select the subplot
                # Create the logo with logomaker
                logo = logomaker.Logo(graph[0], ax=ax, color_scheme='NajafabadiEtAl2017', shade_below=0.5)
                logo.ax.set_ylabel('Frequency', fontsize=axis_label_fontsize)
                # Set x-ticks and labels
                positions = [k for k in range(len(graph[1]))]
                logo.ax.set_xticklabels(graph[1], fontsize=xtick_label_fontsize)  # Rotate labels for visibility
                logo.ax.set_xticks(positions)
                # Set the title for this subplot
                ax.set_title(title_row[i], fontsize=title_fontsize)
            plt.tight_layout()
            plt.show()  # Display the plot for this row
        if only_combined:
            break


def assing_class_for_cluster(AA_list):
    cluster_info = {}
    count_hydropho = 0
    count_hydrophi = 0
    aromatic = 0
    charged = 0
    for AA in AA_list:
        if amino_acid_classes[AA] == "Hydrophobic":
            count_hydropho += 1
        elif amino_acid_classes[AA] == "Hydrophilic":
            count_hydrophi += 1
        elif amino_acid_classes[AA] == "Aromatic":
            aromatic += 1
        else:
            charged += 1
    total = len(AA_list)
    hydrophi_width = count_hydrophi / total
    hydropho_width = count_hydropho / total
    aromatic_width = aromatic / total
    charged_width = charged / total
    return [hydropho_width, hydrophi_width, aromatic_width, charged_width]


def add_mapped_column(df, column_to_map, mapping_dict, new_column_name):
    # Use the map method to create a new column based on the existing one
    df[new_column_name] = df[column_to_map].map(mapping_dict)
    return df


def map_colors(df, color_key):
    # Function to apply the color mapping with the exception
    def color_for_row(row):
        if row['cluster_index'] == -1:
            return '#FFFFFF'  # White color for the exception case
        else:
            return color_key.get(row['AA Class'], "#FFFFFF")  # You can set a default color

    # Apply the function to each row
    df['cluster_color'] = df.apply(color_for_row, axis=1)
    return df


def assign_cluster_colors(dataframe):
    # Get unique clusters excluding the noise (-1)
    unique_clusters = dataframe['cluster_index'].unique().tolist()
    widths = {}
    dataframe = add_mapped_column(dataframe, "AA", aa_mapping, "AA1")
    dataframe = add_mapped_column(dataframe, "AA1", amino_acid_classes, "AA Class")
    dataframe = map_colors(dataframe, color_key)
    for unique_cluster in unique_clusters:
        if unique_cluster == -1:
            continue
        cluster_df = dataframe[dataframe["cluster_index"] == unique_cluster]
        AA_list = transform_to_1_letter_code(cluster_df['AA'].values.tolist())
        colors = assing_class_for_cluster(AA_list)
        widths[unique_cluster] = colors
    return [dataframe, widths]


def map_column_with_dict(dataframe, mapping_dict, target_column, new_column_name):
    """
    Maps values in a specified column of a DataFrame according to a provided dictionary,
    and creates a new column with these mapped values.
    Parameters:
    - dataframe (pd.DataFrame): The DataFrame to operate on.
    - mapping_dict (dict): A dictionary where keys are the values to be mapped from the target_column
                           and values are the corresponding new values for the new_column.
    - target_column (str): The name of the column in the DataFrame whose values are to be mapped.
    - new_column_name (str): The name of the new column that will be created with mapped values.
    Returns:
    - pd.DataFrame: The original DataFrame with the new column added.
    """
    # Check if the target column exists in the DataFrame
    if target_column not in dataframe.columns:
        raise ValueError(f"Column '{target_column}' not found in DataFrame")
    # Map the target column values to the new column using the provided dictionary
    dataframe[new_column_name] = dataframe[target_column].map(mapping_dict)
    return dataframe


def cluster_3d_points(dataframe, advance_options=[2, 10]):
    # Ensure the DataFrame has the correct column names
    if not set(['X', 'Y', 'Z']).issubset(dataframe.columns):
        raise ValueError("Dataframe must contain 'X', 'Y', and 'Z' columns")
    # Extract the points from the dataframe
    points = dataframe[['X', 'Y', 'Z']].values
    # Apply DBSCAN clustering
    # The parameters `eps` and `min_samples` need to be chosen based on the dataset.
    # These are common starting values, but for optimal results they should be fine-tuned
    dbscan = DBSCAN(eps=advance_options[0], min_samples=advance_options[1], algorithm="kd_tree")
    dbscan.fit(points)
    # Add the cluster index to the dataframe
    dataframe['cluster_index'] = dbscan.labels_
    # Return the dataframe with the new 'cluster_index' column
    return dataframe


def amino_acid_statistics_per_list(list_of_lists):
    # Define all 20 amino acids (1-letter codes)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    # Initialize dictionaries to hold counts and frequencies for each sublist
    counts_per_list = {f'n{aa}': [] for aa in amino_acids}
    frequencies_per_list = {f'f{aa}': [] for aa in amino_acids}
    total_counts_per_list = []
    # Process each sublist
    for sublist in list_of_lists:
        sublist_count = {aa: 0 for aa in amino_acids}
        total_count = 0
        # Count amino acids in the sublist
        for aa in sublist:
            if aa in sublist_count:
                sublist_count[aa] += 1
                total_count += 1
        total_counts_per_list.append(total_count)
        # Calculate frequencies and update the counts for this sublist
        for aa in amino_acids:
            counts_per_list[f'n{aa}'].append(sublist_count[aa])
            frequency = sublist_count[aa] / total_count if total_count > 0 else 0
            frequencies_per_list[f'f{aa}'].append(frequency)
    # Combine the results into one dictionary
    stats_per_list = {**counts_per_list, **frequencies_per_list, 'total_counts': total_counts_per_list}
    return stats_per_list

def plot_cluster_compositions(cluster_compositions_dict):
    # Create a single figure and axis
    fig, ax = plt.subplots(figsize=(15, len(cluster_compositions_dict)))
    # Define the color for each composition entry
    composition_colors = [color_key['Hydrophobic'], color_key['Hydrophilic'],
                          color_key['Aromatic'], color_key['Charged interactions']]
    # Sort the clusters to maintain order
    sorted_cluster_indexes = sorted(cluster_compositions_dict.keys())
    # Calculate the height of each bar so that there is no space between them
    bar_height = 1.0
    # Loop through each cluster and its composition
    for cluster_index in sorted_cluster_indexes:
        composition = cluster_compositions_dict[cluster_index]
        # The starting point for each bar segment
        left = 0
        # Create each segment of the bar
        for comp, color in zip(composition, composition_colors):
            ax.barh(cluster_index, comp, height=bar_height, color=color, left=left, edgecolor='white')
            left += comp  # Increment the left position
    # Set the labels and ticks
    ax.set_xlabel('Proportion')
    ax.set_yticks(sorted_cluster_indexes)  # Set y-ticks to cluster indexes
    ax.set_yticklabels([f'Cluster {i}' for i in sorted_cluster_indexes])
    # Remove the spines for a cleaner look
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    # Create a legend, moved further down by adjusting the bbox_to_anchor parameter
    handles = [plt.Rectangle((0, 0), 1, 1, color=color_key[label]) for label in color_key]
    labels = ["Hydrophobic (G,A,V,C,P,L,I,M)", "Hydrophilic (S,T,N,Q)", "Aromatic (W,F,S)", "Charged (D,E,K,R,H)"]
    plt.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1), ncol=len(color_key))
    # Adjust the layout
    plt.tight_layout()
    # Show the plot
    plt.show()