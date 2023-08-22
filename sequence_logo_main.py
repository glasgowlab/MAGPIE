# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 17:18:58 2023

@author: camlo
"""
import os.path
import plotly.subplots
from plotly.offline import init_notebook_mode, iplot
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import pandas as pd
import glob
from scipy.spatial import KDTree
from sklearn.cluster import KMeans
import math
import logomaker
import numpy as np
import helper_functions
import warnings
from matplotlib.patches import FancyArrowPatch

def create_3d_graph(df1, df2,is_ligand, ligand_bonds = {}):
    # Get XYZ positions from the DataFrame columns
    x1, y1, z1 = df1['X'], df1['Y'], df1['Z']
    x2, y2, z2 = df2['X'], df2['Y'], df2['Z']
    color_shapely = df1['shapely'].values.tolist()
    color_polar = df1['polar'].values.tolist()
    if is_ligand:
        names = df2['atom_name'].values.tolist()
        color_df2=  df2["color"].values.tolist()
        size2 = 15
    else:
        names = df2['residue_index'].values.tolist()
        color_df2 = "black"
        size2 = 15
    init_notebook_mode(connected=True)

    scatter_trace1 = go.Scatter3d(
        x=x1,
        y=y1,
        z=z1,
        mode='markers',
        marker=dict(
            size=9,
            color=color_shapely,
            opacity=0,
            line=dict(color='black', width=2)
        ),
        text=df1['AA'],  
        hoverinfo='text',
        hoverlabel = dict(bgcolor='yellow', bordercolor='black')
    )

    scatter_trace2 = go.Scatter3d(
        x=x2,
        y=y2,
        z=z2,
        mode='markers',
        marker=dict(
            size=size2,
            color=color_df2,
            opacity=0,
            line=dict(color='white', width=5)
        ),
        text = names,
        hoverinfo='text',
        hoverlabel=dict(bgcolor='gray', bordercolor='white')
    )
    buttons = []
    buttons.append(dict(label='Shapely Colours', method='restyle',  args=[{'marker.color': [color_shapely]}, [0]]))
    buttons.append(dict(label='Amino Colours', method='restyle', args=[{'marker.color': [color_polar]}, [0]]))
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
    
    
    bonds_made = []
    graphs = [scatter_trace1, scatter_trace2]
    for bond in ligand_bonds:
        for pair in ligand_bonds[bond]:
            comb = [str(bond), str(pair)]
            if comb not in bonds_made or comb not in bonds_made:
                bonds_made.append(comb)
                atom_coords = df2[df2['atom_serial_number'].isin(comb)]
                if len(atom_coords) >= 2:
                    point1 = atom_coords[['X', 'Y', 'Z']].values[0]
                    point2 = atom_coords[['X', 'Y', 'Z']].values[1]
                    line_trace = go.Scatter3d(
                        x=[point1[0], point2[0]],
                        y=[point1[1], point2[1]],
                        z=[point1[2], point2[2]],
                        mode='lines',
                        line=dict(color='black', width=8)
                    )
                    graphs.append(line_trace)

    # Create the figure and add the traces
    fig = go.Figure(data=graphs, layout=layout)
    fig.update_layout(updatemenus=updatemenus)
    # Show the interactive plot
    iplot(fig)


def find_nearest_points(target, binders, radius, is_ligand
                        ):
    # Convert XYZ positions of target DataFrame to a list of tuples
    target_points = target[['X', 'Y', 'Z']].values.tolist()
    if is_ligand:
        positions = list(target['atom_name'].values)
    else:
        positions = list(target['residue_index'].values)

    # Convert XYZ positions of binders DataFrame to a list of tuples
    binder_points = binders[['X', 'Y', 'Z']].values.tolist()

    # Create KDTree for binder_points
    tree = KDTree(binder_points)

    nearest_points_indices = []
    for i, target_point in enumerate(target_points):
        # Query binders DataFrame within the specified radius of target_point
        indices = tree.query_ball_point(target_point, r=radius)
        if not indices:
            continue
        nearest_points_indices.append([indices, positions[i]])

    # Create a new DataFrame with nearest points
    nearest_points_df = pd.DataFrame()
    for indices in nearest_points_indices:
        temp_df = binders.iloc[indices[0]]
        nearest_points_df = pd.concat([nearest_points_df, temp_df], ignore_index=True)

    return nearest_points_df


# Function to calculate arrow position dynamically
def calculate_arrow_position(subplot_index):
    arrow_x = (subplot_index - 1) * 0.25 + 0.1
    arrow_tail_x = arrow_x + 0.1
    return arrow_x, arrow_tail_x



def transform_to_1_letter_code(amino_acids_3_letter):
    # Mapping dictionary for 3-letter to 1-letter code
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

    amino_acids_1_letter = [aa_mapping[aa] for aa in amino_acids_3_letter]
    return amino_acids_1_letter


def calculate_frequency(character, lst):
    count = lst.count(character)

    frequency = count / len(lst)

    return frequency


def calculate_bits(list_of_AA, sequence_list):
    list_of_frequencies = []
    for AA in list_of_AA:
        list_of_frequencies.append(calculate_frequency(AA, sequence_list))
    S = math.log2(20)
    H = 0

    for f in list_of_frequencies:
        if f == 0:
            continue
        H = H + (-f) * (math.log2(f))
    H = H + 19 / (np.log(2) * 2 * len(sequence_list))
    R = S - H

    heights = []
    for f in list_of_frequencies:
        heights.append(np.abs(f * 100))
    return heights


def create_sequence_logo(df, target):
    # Create a Logo object

    logo = logomaker.Logo(df,
                          color_scheme='NajafabadiEtAl2017')
    logo.ax.set_ylabel('Frequency')
    positions = [i for i in range(len(target))]
    logo.ax.set_xticklabels(target)
    logo.ax.set_xticks(positions)
    logo.ax.set_title(f'Residues Indexes: {"-".join(map(str,target))} ')
    return logo

def create_sequence_logo_list(df_list):
    # Calculate the number of columns based on the number of logos
    num_logos = len(df_list)
    num_cols = min(num_logos, 3)  # Change this value to control the number of columns

    # Calculate the number of rows needed to display all logos in a grid
    num_rows = (num_logos + num_cols - 1) // num_cols

    # Set the figure size based on the number of rows and columns
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(5 * num_cols, 5 * num_rows),
                             gridspec_kw={'width_ratios': [1] * num_cols})

    # Flatten the axes to handle both single row and multiple row cases
    axes_flat = axes.flat if isinstance(axes, np.ndarray) else [axes]

    # Draw each sequence logo on its respective subplot
    for i, (df, ax) in enumerate(zip(df_list, axes_flat)):
        logo = logomaker.Logo(df[0], ax=ax, color_scheme='NajafabadiEtAl2017', shade_below=0.5)
        logo.ax.set_ylabel('Frequency')
        positions = [i for i in range(len(df[1]))]
        logo.ax.set_xticklabels(df[1])
        logo.ax.set_xticks(positions)

        ax.set_title(f'Residues Indexes: {"-".join(map(str, df[1]))}')

    # Hide any unused subplots if there are fewer logos than the number of axes
    for ax in axes_flat[num_logos:]:
        ax.axis('off')

    # Adjust the spacing between subplots
    plt.subplots_adjust(wspace=0.5, hspace=0.5)  # Modify wspace and hspace as needed
    
    plt.show()
def plot_sequence_logo(df, filename=None):
    # Calculate the height of each letter for each position
    stacked_df = df.apply(lambda row: pd.Series(row.sort_values(ascending=False).values), axis=1)
    bottom = stacked_df.cumsum(axis=1).shift(1, axis=1).fillna(0)
    top = bottom + df.values

    # Plotting the sequence logo
    fig, ax = plt.subplots(figsize=(10, 5))

    for base in df.columns:
        ax.fill_between(df.index, bottom[base], top[base], alpha=0.8, linewidth=0.4)

    ax.set_xticks(range(len(df.index)))
    ax.set_xticklabels(df.index)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(['0', '1'])
    ax.set_xlabel('Position')
    ax.set_ylabel('Probability')
    ax.set_title('')

    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')  # Save the figure as an image file
    else:
        plt.show()


def plot(list_of_paths, target_chain, binder, is_ligand,to_show, distance, threads_num):
    
    if is_ligand:
        target_chain_cordinates = helper_functions.extract_info_ligand(list_of_paths[0], target_chain)
        ligand_bonds = helper_functions.extract_connections(list_of_paths[0])
    else:
        target_chain_cordinates = helper_functions.extract_info_pdb(list_of_paths[0], target_chain)
    
    data_frame_target = pd.DataFrame(target_chain_cordinates)
    binder_chain_cordinates = []
    
    if threads_num == 1:
        for file in list_of_paths:        
            binder_chain_cordinates += helper_functions.extract_info_pdb(file, binder)
    else:
        chunks = helper_functions.make_chunks(list_of_paths,threads_num)
        result_queue = multiprocess.Queue()  
        threads = []
        for thread_num in range(threads_num):
            current_thread = multiprocess.Process(target=helper_functions.extract_list_pdb,args=(
            chunks[thread_num], 
            binder,
            result_queue))
            threads.append(current_thread)
            current_thread.start()
        for t in threads:
            t.join()        
        for _ in threads:
            thread_num, result = result_queue.get()
            binder_chain_cordinates += result 
            
    data_frame_binders = pd.DataFrame(binder_chain_cordinates)
    if to_show == "all":
        nearest_neighbors_df = find_nearest_points(data_frame_target, data_frame_binders, distance,is_ligand)
    else:
        if is_ligand:
            to_show_df =data_frame_target[data_frame_target['atom_name'].isin(to_show)]
        else:
            to_show_df =data_frame_target.loc[data_frame_target['residue_index'].isin(to_show)]
        nearest_neighbors_df = find_nearest_points(to_show_df, data_frame_binders, distance,is_ligand)
        
    if is_ligand:   
        create_3d_graph(nearest_neighbors_df,data_frame_target, is_ligand, ligand_bonds)
    else:
        create_3d_graph(nearest_neighbors_df,data_frame_target, is_ligand)
    
    return data_frame_target,data_frame_binders
    
def sequence_logos(data_frame_target, data_frame_binder, sequence_logo_targets, is_ligand, only_combined_logo, distance):

    warnings.filterwarnings("ignore")
    model = logomaker.get_example_matrix('ww_information_matrix',
                                         print_description=False)
    list_of_AA = model.columns.to_list()
    rows_bits= []
    residues = []
    plots = []


    for i,target in enumerate(sequence_logo_targets):
        if is_ligand:
            current_df = data_frame_target.loc[data_frame_target['atom_name'] == target]
        else:
            current_df = data_frame_target.loc[data_frame_target['residue_index'] == target ]
        near_neighbor_current = find_nearest_points(current_df,data_frame_binder,distance, is_ligand)
        if near_neighbor_current.empty:
            continue
        AA_sq = transform_to_1_letter_code(near_neighbor_current['AA'].values.tolist())
        residues.append(f' {target}, n = {len(AA_sq)} ')
        bits = calculate_bits(list_of_AA, AA_sq)
        rows_bits.append(bits)
        df = pd.DataFrame(columns=model.columns)
        df = pd.concat([df, pd.DataFrame([bits], columns=df.columns)], ignore_index=True)
        plots.append([df, [target]])
    if not len(residues) == 1:
        df = pd.DataFrame(columns=model.columns)
        df = pd.concat([df, pd.DataFrame(rows_bits, columns=df.columns)], ignore_index=True)
        plots.append([df, residues])
    if  only_combined_logo:
        df = pd.DataFrame(columns=model.columns)
        df = pd.concat([df, pd.DataFrame(rows_bits, columns=df.columns)], ignore_index=True)
        plots = [[df, residues]]
    create_sequence_logo_list(plots)
    


    
