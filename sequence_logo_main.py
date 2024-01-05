# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 17:18:58 2023

@author: camlo
"""
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
    color_hb_sc = df1['H-Bond SC'].values.tolist()
    color_hb_bb = df1['H-Bond BB'].values.tolist()


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
    buttons.append(dict(label='Shapely Colours', method='restyle',  args=[{'marker.color': [color_shapely]}, [0]]))
    buttons.append(dict(label='Amino Colours', method='restyle', args=[{'marker.color': [color_polar]}, [0]]))

    if is_ligand:
        buttons.append(dict(label='H-Bond', method='restyle', args=[{'marker.color': [color_hb_sc]}, [0]]))
    else:
        buttons.append(dict(label='H-Bond BB', method='restyle', args=[{'marker.color': [color_hb_bb]}, [0]]))
        buttons.append(dict(label='H-Bond SC', method='restyle', args=[{'marker.color': [color_hb_sc]}, [0]]))
        buttons.append(dict(label='Polar Interactions', method='restyle', args=[{'marker.color': [color_salt]}, [0]]))

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

    
    bonds_made = []

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
                        line=dict(color='black', width=8),
                        hoverinfo='skip',
                        showlegend = False
                    )
                    
                    graphs.append(line_trace)

    # Create the figure and add the traces
    fig = go.Figure(data=graphs, layout=layout)
    fig.update_layout(updatemenus=updatemenus)
    # Show the interactive plot
    iplot(fig)



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

def create_sequence_logo_list(data,only_combined, is_ligand):

    # Titles for each type of graph
    axis_label_fontsize = 28
    title_fontsize = 32
    xtick_label_fontsize = 24
    titles = ["Residues in Contact", "H-Bond Backbone", "H-Bond Side-Chain", "Polar Contact"]
    if is_ligand:
        titles[2] = "H-Bonds"
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
                    logo.ax.set_ylabel('Frequency',fontsize=axis_label_fontsize)

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
def plot(list_of_paths, target_id_chain, binder_id_chain, is_ligand, distance):


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
                exit(0)
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

    ligand_bonds ={}
    if is_ligand:
        ligand_bonds = helper_functions.extract_connections(list_of_paths[0])

    target_in_contact, binder_in_contact = find_nearest_points(target_chain_data_frame, binder_chain_data_frame,
                                                               distance)
    residues_to_plot = []
    given_ph = 7
    ph_range = 1.5

    for i in range (len(target_in_contact)): #Every file
        residues_in_contact_binder = []
        residues_in_contact_target = []
        residues_salt_bridges_binder = []
        residues_salt_bridges_target = []
        for k,residue in enumerate(chains_binder[i]):
            if k in binder_in_contact[i]:
                residues_in_contact_binder.append([hbonds_saltbridges.hydrogen_bond_donor_optimal_positions_and_directions(residue),hbonds_saltbridges.hydrogen_bond_acceptor_optimal_positions_and_directions(residue)])
                residues_salt_bridges_binder.append(residue)

        if is_ligand:
            residue = chains_target[i].get_residues().__next__()
            for k,atom in enumerate(residue):
                if k in target_in_contact[i]:
                    residues_in_contact_target.append(
                        [hbonds_saltbridges.hydrogen_bond_donor_optimal_positions_and_directions_atom(atom,residue),
                         hbonds_saltbridges.hydrogen_bond_acceptor_optimal_positions_and_directions_atom(atom,residue)])
        else:
            for k,residue in enumerate(chains_target[i]):
                if k in target_in_contact[i]:
                    residues_in_contact_target.append([hbonds_saltbridges.hydrogen_bond_donor_optimal_positions_and_directions(residue),hbonds_saltbridges.hydrogen_bond_acceptor_optimal_positions_and_directions(residue)])
                    residues_salt_bridges_target.append(residue)
                    residues_salt_bridges_target.append(residue)
        for k in range (len(residues_in_contact_target)):
            binder_h_bond_sc = False
            binder_h_bond_bb = False
            found_salt_bridge = False
            for l in range (len(residues_in_contact_binder)):
                if not is_ligand and not found_salt_bridge:
                    salt_bridges_binder = residues_salt_bridges_binder[l]
                    salt_bridges_target = residues_salt_bridges_target[k]
                    charged_binder = hbonds_saltbridges.determine_charge(salt_bridges_binder,given_ph, ph_range)
                    charged_target = hbonds_saltbridges.determine_charge(salt_bridges_target,given_ph, ph_range)
                    if charged_binder*charged_target == -1 and charged_target +charged_binder == 0:
                        found_salt_bridge = hbonds_saltbridges.compare_residue_distances(salt_bridges_target, salt_bridges_binder,4)
                for residue in residues_in_contact_binder[l][0]:
                    for residue_2 in residues_in_contact_target[k][1]:
                        if residue_2[2] and not binder_h_bond_bb:
                            binder_h_bond_bb = hbonds_saltbridges.check_hydrogen_bond(residue,residue_2)
                        if not residue_2[2]  and not binder_h_bond_sc:
                            binder_h_bond_sc = hbonds_saltbridges.check_hydrogen_bond(residue,residue_2)
                for residue in residues_in_contact_binder[l][1]:
                    for residue_2 in residues_in_contact_target[k][0]:
                        if residue_2[2]  and not binder_h_bond_bb:
                            binder_h_bond_bb = hbonds_saltbridges.check_hydrogen_bond(residue,residue_2)
                        if not residue_2[2]  and not binder_h_bond_sc:
                            binder_h_bond_sc = hbonds_saltbridges.check_hydrogen_bond(residue,residue_2)

                if k == len(residues_in_contact_target)-1:

                    residue_found = binder_chain_data_frame.iloc[i][binder_in_contact[i][l]]
                    if  binder_h_bond_sc:
                        residue_found["H-Bond SC"] = '#FF0000'
                    else:
                        residue_found["H-Bond SC"] = '#FFFFFF'

                    if binder_h_bond_bb:
                        residue_found["H-Bond BB"] = '#FF0000'
                    else:
                        residue_found["H-Bond BB"] = '#FFFFFF'
                    if not is_ligand and found_salt_bridge:
                        residue_found["Salt Bridge"] = "#FF0000"
                    if not is_ligand and not found_salt_bridge:
                        residue_found["Salt Bridge"] = "#FFFFFF"
                    residues_to_plot.append(residue_found)

                if binder_h_bond_sc and binder_h_bond_bb and found_salt_bridge:
                    residue_found = binder_chain_data_frame.iloc[i][binder_in_contact[i][l]]
                    residue_found["H-Bond SC"] = '#FF0000'
                    residue_found["H-Bond BB"] = '#FF0000'
                    if not is_ligand:
                        residue_found["Salt Bridge"] = "#FF0000"
                    residues_to_plot.append(residue_found)
                    break

    residue_found_df   = pd.DataFrame(residues_to_plot)

    target_to_to_plot = []
    for x in target_chain_data_frame.iloc[0]:
        if x is not None:
            target_to_to_plot.append(x)
    #
    create_3d_graph(residue_found_df,pd.DataFrame(target_to_to_plot), is_ligand, ligand_bonds)

    return residue_found_df,pd.DataFrame(target_chain_ca_coords[reference_id])

def sequence_logos(residues_found, target_residues, sequence_logo_targets, is_ligand, only_combined_logo, radius ):

    warnings.filterwarnings("ignore")
    model = logomaker.get_example_matrix('ww_information_matrix',
                                         print_description=False)
    list_of_AA = model.columns.to_list()
    rows_bits_all_sq= []
    rows_bits_bb= []
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
            bb_hydrogen_bonds = all_contacts.loc[
                 (residues_found['H-Bond BB'] == '#FF0000')]
            sc_hydrogen_bonds = all_contacts.loc[
                (residues_found['H-Bond SC'] == '#FF0000')]
        else:
            point = target_residues[target_residues['residue_index'] == target]
            res_nem = point["AA"].values

            target = f"{aa_mapping[res_nem[0]]}{target}"
            all_contacts = find_points_within_radius(point, residues_found, radius)
            bb_hydrogen_bonds = all_contacts.loc[
                 (residues_found['H-Bond BB'] == '#FF0000')]
            sc_hydrogen_bonds = all_contacts.loc[
                (residues_found['H-Bond SC'] == '#FF0000')]
            polar_contacts = all_contacts.loc[
                (residues_found['Salt Bridge'] == '#FF0000')]

        AA_all_contacts = transform_to_1_letter_code(all_contacts['AA'].values.tolist())
        AA_bb = transform_to_1_letter_code(bb_hydrogen_bonds['AA'].values.tolist())
        AA_sc = transform_to_1_letter_code(sc_hydrogen_bonds['AA'].values.tolist())
        if not is_ligand:
            AA_pc = transform_to_1_letter_code(polar_contacts['AA'].values.tolist())
        else:
            AA_pc = []

        # Check if any transformed lists are empty and handle accordingly
        if len( AA_all_contacts) == 0:
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
        if len( AA_bb) == 0:
            plots_rows.append([None,None])
        else:
            bb_residue_num = f' {target} \n n = {len(AA_bb)} '
            resi_combined_bb_hb.append(bb_residue_num)
            bb_bits = calculate_bits(list_of_AA, AA_bb)
            rows_bits_bb.append(bb_bits)
            df_bb = pd.DataFrame(columns=model.columns)
            df_bb = pd.concat([df_bb, pd.DataFrame([bb_bits], columns=df_bb.columns)], ignore_index=True)
            plots_rows.append([df_bb, [bb_residue_num]])

        if len(AA_sc) == 0:
            plots_rows.append([None,None])
        else:
            sc_residue_num = f' {target} \n n = {len(AA_sc)} '
            resi_combined_sc_bb.append(sc_residue_num)
            sc_bits = calculate_bits(list_of_AA, AA_sc)
            rows_bits_sc.append(sc_bits)
            df_sc = pd.DataFrame(columns=model.columns)
            df_sc = pd.concat([df_sc, pd.DataFrame([sc_bits], columns=df_sc.columns)], ignore_index=True)
            plots_rows.append([df_sc, [sc_residue_num]])

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
        if len(rows_bits_bb) == 0 :
            plots_rows.append([None,None])
        else:
            df_bb= pd.DataFrame(columns=model.columns)
            df_bb = pd.concat([df_bb, pd.DataFrame(rows_bits_bb, columns=df_bb.columns)],
                                       ignore_index=True)
            plots_rows.append([df_bb, resi_combined_bb_hb])
        if len(rows_bits_sc) == 0 :
            plots_rows.append([None,None])
        else:
            df_sc= pd.DataFrame(columns=model.columns)
            df_sc = pd.concat([df_sc, pd.DataFrame(rows_bits_sc, columns=df_sc.columns)],
                                       ignore_index=True)
            plots_rows.append([df_sc, resi_combined_sc_bb])

        if len(rows_bits_pc) == 0:
            plots_rows.append([None,None])
        else:
            df_pc= pd.DataFrame(columns=model.columns)
            df_pc = pd.concat([df_pc, pd.DataFrame(rows_bits_pc, columns=df_pc.columns)],
                                       ignore_index=True)
            plots_rows.append([df_pc, resi_combined_pc])

    plots_by_rows = []
    for i,plot in enumerate(plots):
        plots_by_rows.append([plot, plots_rows[0+i*3], plots_rows[1+i*3], plots_rows[2+i*3]])
    plots_by_rows.insert(0, plots_by_rows.pop())
    create_sequence_logo_list(plots_by_rows,only_combined_logo, is_ligand)




