# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 17:18:58 2023

@author: camlo
"""
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


def create_3d_graph(df1, df2,is_ligand, ligand_bonds = {}):
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
        color_df2=  df2["color"].values.tolist()
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
        hoverlabel = dict(bgcolor='yellow', bordercolor='black'),
        name = "Binding Residues"
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
        text = names,
        hoverinfo='text',
        hoverlabel=dict(bgcolor='gray', bordercolor='white'),
        name = "Target"
        
    )
    buttons = []
    buttons.append(dict(label='Shapely Colours', method='restyle',  args=[{'text' : [df1['AA']],"marker.color": [color_shapely]}, [0]]))
    buttons.append(dict(label='Amino Colours', method='restyle', args=[{'text' : [df1['AA']],'marker.color': [color_polar]}, [0]]))
    buttons.append(dict(label='H-Bond', method='restyle', args=[{'text' : [df1['AA']],'marker.color': [color_hb]}, [0]]))
    if not is_ligand:
        buttons.append(dict(label='Salt-Bridges', method='restyle', args=[{'text' : [df1['AA']],'marker.color': [color_salt]}, [0]]))
    buttons.append(dict(label='Color Cluster', method='restyle', args=[{'text' : [cluster_names],'marker.color': [color_cluster]}, [0]]))

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
                width=8 # Adjust line width as needed
            ),
            hoverinfo='skip',  # Optionally disable hover info for the lines
            showlegend=False  # Optionally hide the line trace from the legend
        )
        graphs.append(line_trace_target)

    else:
        # Prepare line data
        line_x, line_y, line_z = [], [], []
        # Check distances and prepare line data
        distance_threshold = 2.0
        n_points = len(x2)
        for i in range(n_points):
            for j in range(i + 1, n_points):
                if "H" in names[i] and "H" in names[j]:
                    continue
                if euclidean_distance(x2[i], y2[i], z2[i], x2[j], y2[j], z2[j]) <= distance_threshold:
                    line_x.extend([x2[i], x2[j], None])
                    line_y.extend([y2[i], y2[j], None])
                    line_z.extend([z2[i], z2[j], None])

        # Create line plot for connections
        lines = go.Scatter3d(x=line_x, y=line_y, z=line_z, mode='lines',
                line=dict(
                color='black',  # Choose a color that stands out
                width=8 # Adjust line width as needed
            ),
            hoverinfo='skip',  # Optionally disable hover info for the lines
            showlegend=False  # Optionally hide the line trace from the legend
        )
        graphs.append(lines)

    # Create the figure and add the traces
    fig = go.Figure(data=graphs, layout=layout)
    fig.update_layout(updatemenus=updatemenus)
    # Show the interactive plot
    iplot(fig)
    # py.plot(fig, filename='3D-scatter.html')



def find_nearest_points(target, binders, radius):
    target_points_in_contact = []
    binder_points_in_contact  = []

    for i  in range(0,  len(target)):
        # Create KDTree for binder_points
        # Query binders DataFrame within the specified radius of target_point
        current_target = target.iloc[i]
        current_binder = binders.iloc[i]
        pd.set_option('display.max_colwidth', None)


        binder_points= []
        target_points = []

        for row in current_binder:

            if row is not None and len(row.keys() )!= 0:
                binder_points.append([row['X'], row['Y'], row['Z']])
        for row in current_target:
            if row is not None and len(row.keys() )!= 0:
                target_points.append([row['X'], row['Y'], row['Z']])

        tree = KDTree(np.array(binder_points))

        indices = tree.query_ball_point(np.array(target_points), r=radius)
        targets_in_contact = []
        binder_in_contact = []

        for i,index in enumerate(indices):
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
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

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
        if len(index_list) == 0 :
            continue
        points_within_radius.append(i)

    residues_points = target_points.iloc[points_within_radius]
    return residues_points

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
def create_sequence_logo_list(data,only_combined, is_ligand):

    # Titles for each type of graph
    axis_label_fontsize = 30
    title_fontsize = 36
    xtick_label_fontsize = 30
    y_lable_size = 34
    titles = ["Residues in Contact", "H-Bonds", "Salt Bridges"]
    for j, row in enumerate(data):
        # Filter out the [None, None] graphs
        valid_graphs = []
        title_row= []
        for i,graph in enumerate(row):
            if graph[1] != None:
                valid_graphs.append(graph)
                title_row.append(titles[i])
        num_graphs = len(valid_graphs)

        # Only create subplots for the valid graphs
        if num_graphs > 0:
            if j == 0:
                for i, graph in enumerate(valid_graphs):
                    fig, ax = plt.subplots(figsize=(2*len(graph[1]), 12 ))  # Adjust the size as needed

                    # Create the logo with logomaker
                    logo = logomaker.Logo(graph[0], ax=ax, color_scheme='NajafabadiEtAl2017', shade_below=0.5)
                    logo.ax.set_ylabel('Frequency',fontsize=y_lable_size    )

                    # Set x-ticks and labels
                    positions = [k for k in range(len(graph[1]))]
                    logo.ax.set_xticklabels(graph[1], fontsize=xtick_label_fontsize)  # Rotate labels for visibility
                    logo.ax.set_xticks(positions)

                    # Set the title for this subplot
                    ax.set_title(title_row[i], fontsize = title_fontsize)
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




def assign_cluster_colors(dataframe):
    # Get unique clusters excluding the noise (-1)
    unique_clusters = dataframe['cluster_index'].unique().tolist()
    n_clusters = len(unique_clusters)
    # Generate a list of distinct colors in HEX format
    base_colors = sns.color_palette('Spectral', n_clusters)
    base_colors = base_colors.as_hex()
    mapping_dict = {}
    mapping_dict[-1] = "#FFFFFF"
    for i in range(n_clusters):
        mapping_dict[i] = base_colors[i]
    return mapping_dict

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

def cluster_3d_points(dataframe):
    # Ensure the DataFrame has the correct column names
    if not set(['X', 'Y', 'Z']).issubset(dataframe.columns):
        raise ValueError("Dataframe must contain 'X', 'Y', and 'Z' columns")

    # Extract the points from the dataframe
    points = dataframe[['X', 'Y', 'Z']].values

    # Apply DBSCAN clustering
    # The parameters `eps` and `min_samples` need to be chosen based on the dataset.
    # These are common starting values, but for optimal results they should be fine-tuned
    dbscan = DBSCAN(eps=1.7   , min_samples=5)
    dbscan.fit(points)

    # Add the cluster index to the dataframe
    dataframe['cluster_index'] = dbscan.labels_

    # Return the dataframe with the new 'cluster_index' column
    return dataframe

def plot_cluster_colors(cluster_color_map):
    """
    Plots a simple figure to show the color associated with each cluster index,
    with adjustments as requested by the user: wider bars with black lines dividing them.
    """
    # Sort clusters to ensure -1 (if present) is first
    clusters = sorted(cluster_color_map.keys(), key=lambda x: (x == -1, x))
    colors = [cluster_color_map[cluster] for cluster in clusters]

    # Plot each cluster index with its corresponding color
    plt.figure(figsize=(12, 2))  # Adjust the figure size to make bars wider
    bars = plt.bar(clusters, [1] * len(clusters), color=colors, edgecolor='black', linewidth=1, width=1.05)

    # Adding labels inside bars
    for bar, cluster in zip(bars, clusters):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2, str(cluster),
                 ha='center', va='center', color='black', fontsize=11)

    plt.title('Cluster Color Representation')
    plt.xticks([])  # Hide x-axis ticks
    plt.yticks([])  # Hide y-axis ticks
    plt.show()
def plot(list_of_paths, target_id_chain, binder_id_chain, is_ligand, distance, download_meta):

    chains_target = []
    chains_binder = []

    parser  = Bio.PDB.PDBParser(QUIET=True)
    reference_id = 0
    current_len = 0
    for i,path in enumerate(list_of_paths):
        current_structure  =  parser.get_structure(i, path)
        chains_list = list(current_structure.get_chains())  # Convert iterator to list once here
        if len(chains_list) != 2:
            print(f"{path.split('/')[0]} contains {len(chains_list)} chains instead of expected 2.")
            exit(0)
        for chain in chains_list:
            current_chain_id =  chain.get_id()
            if current_chain_id != target_id_chain and current_chain_id != binder_id_chain:
                print(f"{path.split('/')[0]} contains chain {chain.get_id()}, which is not defined, please remove or rename chain.")
                sys.exit()
            elif current_chain_id == target_id_chain:
                chains_target.append(chain)
                if len([x for x in chain]) > current_len:
                    reference_id = i
            elif current_chain_id == binder_id_chain:
                chains_binder.append(chain)
    target_chain_ca_coords = []
    binder_chain_ca_coords = []

    for i in range (0,len(chains_target)):
        if is_ligand:
            target_chain_ca_coords.append(helper_functions.extract_info_ligand(list_of_paths[i],target_id_chain))
        else:
            target_chain_ca_coords.append(helper_functions.extract_atoms_from_chain(chains_target[i],"CA", list_of_paths[i]))
        binder_chain_ca_coords.append(helper_functions.extract_atoms_from_chain(chains_binder[i],"CA",list_of_paths[i]))

    target_chain_data_frame = pd.DataFrame(target_chain_ca_coords)
    binder_chain_data_frame = pd.DataFrame(binder_chain_ca_coords)


    target_in_contact, binder_in_contact = find_nearest_points(target_chain_data_frame, binder_chain_data_frame,
                                                               distance)
    residues_to_plot = []

    for i in range (len(target_in_contact)): #Every file
        residues_in_contact_binder = []
        residues_in_contact_target = []
        for k,residue in enumerate(chains_binder[i]):
            if k in binder_in_contact[i]:
                residues_in_contact_binder.append(residue)
        if is_ligand:
            residue = chains_target[i].get_residues().__next__()
            residues_in_contact_target.append(residue)
        else:
            for k,residue in enumerate(chains_target[i]):
                if k in target_in_contact[i]:
                    residues_in_contact_target.append(residue)
        for l in range (len(residues_in_contact_binder)):
            binder_h_bond_sc = False
            binder_h_bond_bb = False
            found_salt_bridge = False
            sb_info = ["",""]
            bb_info = ["",""]
            binder_residue = residues_in_contact_binder[l]
            for k in range (len(residues_in_contact_target)):
                target_residue = residues_in_contact_target[k]
                if not binder_h_bond_sc:
                    binder_h_bond_sc, sc_info = hbonds_saltbridges.find_hydrogen_bond(target_residue, binder_residue, is_ligand, False)
                if not is_ligand:
                     if not binder_h_bond_bb:
                        binder_h_bond_bb, bb_info = hbonds_saltbridges.find_hydrogen_bond(target_residue,binder_residue, is_ligand, True)
                     if not found_salt_bridge:
                        found_salt_bridge, sb_info = hbonds_saltbridges.find_salt_bridge(target_residue,binder_residue)
                if k == len(residues_in_contact_target)-1:
                    residue_found = binder_chain_data_frame.iloc[i][binder_in_contact[i][l]]
                    residue_found["H-Bond SC Info"] = f"Acceptor: {sc_info[0]}, donor {sc_info[1]}"
                    residue_found["H-Bond BB Info"] = f"Acceptor: {bb_info[0]}, donor {bb_info[1]}"
                    if (not sb_info[0] == "" or not bb_info[0] == "" or not sb_info[0] =="" )and not is_ligand:
                        residue_found["Residue in target (Same PDB)"] = target_residue.get_full_id()[3][1]
                    if  binder_h_bond_sc:
                        residue_found["H-Bond SC"] = '#FF0000'
                    else:
                        residue_found["H-Bond SC"] = '#FFFFFF'
                    if binder_h_bond_bb:
                        residue_found["H-Bond BB"] = '#FF0000'
                    else:
                        residue_found["H-Bond BB"] = '#FFFFFF'
                    if binder_h_bond_sc or  binder_h_bond_bb:
                        residue_found["H-bond"] = "#FF0000"
                    else:
                        residue_found["H-bond"] = "#FFFFFF"
                    if not is_ligand and found_salt_bridge:
                        residue_found["Salt Bridge"] = "#FF0000"
                        residue_found["Salt Bridge Info "] = f"atom1: {sb_info[0]}, atom2  {sb_info[1]}"
                    if not is_ligand and not found_salt_bridge:
                        residue_found["Salt Bridge"] = "#FFFFFF"
                        residue_found["Salt Bridge Info"] = f"atom1: {sb_info[0]}, atom2  {sb_info[1]}"
                    residues_to_plot.append(residue_found)
                if binder_h_bond_sc and is_ligand:
                    residue_found = binder_chain_data_frame.iloc[i][binder_in_contact[i][l]]
                    residue_found["H-Bond SC"] = '#FF0000'
                    residue_found["H-Bond SC Info"] = f"Acceptor: {sc_info[0]}, donor {sc_info[1]}"
                    residue_found["H-Bond BB"] = '#FFFFFF'
                    residue_found["H-Bond BB Info"] = f"Acceptor: {sc_info[0]}, donor {sc_info[1]}"
                    residue_found["H-bond"] = "#FF0000"
                    residues_to_plot.append(residue_found)
                    break
                if binder_h_bond_sc and binder_h_bond_bb and found_salt_bridge:
                    residue_found = binder_chain_data_frame.iloc[i][binder_in_contact[i][l]]
                    residue_found["H-Bond BB"] = '#FF0000'
                    residue_found["H-Bond BB Info"] = f"Acceptor: {bb_info[0]}, donor {bb_info[1]}"
                    residue_found["H-Bond BB"] = '#FF0000'
                    residue_found["H-Bond BB Info"] = f"Acceptor: {bb_info[0]}, donor {bb_info[1]}"
                    residue_found["H-bond"] = '#FF0000'
                    residue_found["Residue in target (Same PDB)"] = target_residue.get_full_id()[3][1]
                    if not is_ligand:
                        residue_found["Salt Bridge"] = "#FF0000"
                        residue_found["Salt Bridge Info"] = f"atom1: {sb_info[0]}, atom2  {sb_info[1]}"
                    residues_to_plot.append(residue_found)
                    break

    polar_int_map = {}
    polar_int_map["#FF0000"] = "True"
    polar_int_map["#FFFFFF"] = "False"

    residue_found_df   = pd.DataFrame(residues_to_plot)
    residue_found_df = cluster_3d_points(residue_found_df)
    color_map = assign_cluster_colors(residue_found_df)
    plot_cluster_colors(color_map)
    residue_found_df = map_column_with_dict(residue_found_df, color_map, "cluster_index", "cluster_color")

    if download_meta:
        residue_found_df = map_column_with_dict(residue_found_df, polar_int_map, "H-Bond BB", "H-Bond BB Flag")
        residue_found_df = map_column_with_dict(residue_found_df, polar_int_map, "H-Bond SC", "H-Bond SC Flag")
        if not is_ligand:
            residue_found_df = map_column_with_dict(residue_found_df, polar_int_map, "Salt Bridge", "Salt Bridge Flag")
        name = list_of_paths[0].split("/")[-2]
        residue_found_df.to_csv(f"meta_data_from_{name}.csv")
    target_to_to_plot = []
    for x in target_chain_data_frame.iloc[reference_id]:
        if x is not None:
            target_to_to_plot.append(x)

    create_3d_graph(residue_found_df,pd.DataFrame(target_to_to_plot), is_ligand)

    return residue_found_df,pd.DataFrame(target_chain_ca_coords[reference_id])

def sequence_logos(residues_found, target_residues, sequence_logo_targets, is_ligand, only_combined_logo, radius ):

    warnings.filterwarnings("ignore")
    model = logomaker.get_example_matrix('ww_information_matrix',
                                         print_description=False)
    list_of_AA = model.columns.to_list()
    rows_bits_all_sq= []
    rows_bits_hb= []
    rows_bits_sc= []
    rows_bits_pc= []

    resi_combined_all_contacts = []
    resi_combined_sc_bb = []
    resi_combined_bb_hb = []
    resi_combined_pc = []

    plots = []
    plots_rows = []

    for i,target in enumerate(sequence_logo_targets):

        if is_ligand:
            point =  target_residues[target_residues['atom_name'] == target]
            all_contacts = find_points_within_radius(point, residues_found, radius)

        else:
            point = target_residues[target_residues['residue_index'] == target]
            res_nem = point["AA"].values

            target = f"{aa_mapping[res_nem[0]]}{target}"
            all_contacts = find_points_within_radius(point, residues_found, radius)
            hydrogen_bonds = all_contacts.loc[
                 (residues_found['H-bond'] == '#FF0000')]
            polar_contacts = all_contacts.loc[
                (residues_found['Salt Bridge'] == '#FF0000')]

        AA_all_contacts = transform_to_1_letter_code(all_contacts['AA'].values.tolist())
        if not is_ligand:
            AA_hb = transform_to_1_letter_code(hydrogen_bonds['AA'].values.tolist())
            AA_pc = transform_to_1_letter_code(polar_contacts['AA'].values.tolist())
        else:
            AA_pc = []
            AA_hb = []

        # Check if any transformed lists are empty and handle accordingly
        if len( AA_all_contacts) == 0:
            print(f"No AA within {radius} Å of target id: {target}")
            continue
        else:
            residue_num = f' {target} \n n = {len(AA_all_contacts)} '
            resi_combined_all_contacts.append(residue_num)
            all_bits = calculate_bits(list_of_AA, AA_all_contacts)
            rows_bits_all_sq.append(all_bits)
            df_all_contact = pd.DataFrame(columns=model.columns)
            df_all_contact = pd.concat([df_all_contact, pd.DataFrame([all_bits], columns=df_all_contact.columns)],
                                       ignore_index=True)
            plots.append([df_all_contact, [residue_num]])

        # Repeat the same check for other contact types
        if len( AA_hb) == 0:
            plots_rows.append([None,None])
        else:
            bb_residue_num = f' {target} \n n = {len(AA_hb)} '
            resi_combined_bb_hb.append(bb_residue_num)
            hb_bits = calculate_bits(list_of_AA, AA_hb)
            rows_bits_hb.append(hb_bits)
            hb_bb = pd.DataFrame(columns=model.columns)
            hb_bb = pd.concat([hb_bb, pd.DataFrame([hb_bits], columns=hb_bb.columns)], ignore_index=True)
            plots_rows.append([hb_bb, [bb_residue_num]])


        if len( AA_pc) == 0:
            plots_rows.append([None,None])
        else:
            pc_residue_num = f' {target} \n n = {len(AA_pc)} '
            resi_combined_pc.append(pc_residue_num)
            pc_bits = calculate_bits(list_of_AA, AA_pc)
            rows_bits_pc.append(pc_bits)
            df_pc = pd.DataFrame(columns=model.columns)
            df_pc = pd.concat([df_pc, pd.DataFrame([pc_bits], columns=df_pc.columns)], ignore_index=True)
            plots_rows.append([df_pc, [pc_residue_num]])

    if not len(resi_combined_all_contacts) == 1:
        df_all_contact = pd.DataFrame(columns=model.columns)
        df_all_contact = pd.concat([df_all_contact, pd.DataFrame(rows_bits_all_sq, columns=df_all_contact.columns)], ignore_index=True)
        plots.append([df_all_contact, resi_combined_all_contacts])
        if len(rows_bits_hb) == 0 :
            plots_rows.append([None,None])
        else:
            hb_bb= pd.DataFrame(columns=model.columns)
            hb_bb = pd.concat([hb_bb, pd.DataFrame(rows_bits_hb, columns=hb_bb.columns)],
                                       ignore_index=True)
            plots_rows.append([hb_bb, resi_combined_bb_hb])
        if len(rows_bits_pc) == 0:
            plots_rows.append([None,None])
        else:
            df_pc= pd.DataFrame(columns=model.columns)
            df_pc = pd.concat([df_pc, pd.DataFrame(rows_bits_pc, columns=df_pc.columns)],
                                       ignore_index=True)
            plots_rows.append([df_pc, resi_combined_pc])

    plots_by_rows = []
    for i,plot in enumerate(plots):
        plots_by_rows.append([plot, plots_rows[0+i*2], plots_rows[1+i*2]])
    plots_by_rows.insert(0, plots_by_rows.pop())
    create_sequence_logo_list(plots_by_rows,only_combined_logo, is_ligand)
