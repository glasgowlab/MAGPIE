from Bio.PDB import Residue
import numpy as np
from Bio.PDB.vectors import Vector

def get_bonded_hydrogens(atom, residue, a,threshold=1.6):
    """ Return a list of atoms bonded to a given atom within a specified distance threshold. """
    bonded_atoms = []
    for other_atom in residue:
        if other_atom is not atom and other_atom.element == "H":
            distance = (atom.get_vector() - other_atom.get_vector()).norm()
            if distance < threshold:
                bonded_atoms.append(other_atom)
    return bonded_atoms

def scale_vector(vec, scalar):
    """ Scale a Bio.PDB Vector by a scalar. """
    return Vector(vec[0] * scalar, vec[1] * scalar, vec[2] * scalar)


def is_backbone_hydrogen(atom, residue):
    """ Check if the hydrogen atom is bonded to the backbone nitrogen atom. """
    if atom.element != 'H':
        return False

    nitrogen = residue['N'] if 'N' in residue else None
    try:
        if nitrogen:
            distance = atom - nitrogen
            # Assuming a threshold distance for covalent H-N bond (e.g., 1.2 Å)
            if distance < 1.2:
                return True
        return False
    except:
        return False




def calculate_distance(vec1, vec2):
    """ Calculate the Euclidean distance between two vectors. """
    return np.linalg.norm(vec1 - vec2)

def calculate_angle(vec1, vec2):
    """ Calculate the angle in degrees between vectors 'vec1' and 'vec2'"""
    cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    angle = np.arccos(cos_angle)
    return np.degrees(angle)


def find_acceptor_and_donors_from_residue_sc(residue):
    acceptors = []
    donors = []
    for atom in residue:
        if atom.element in acceptor_names:
            if atom.name in backbone_atoms:
                continue
            acceptor_coords = atom.get_coord()
            neighboors = get_bonded_hydrogens(atom, residue, 1.2)
            neighboors_H_and_direction = []
            for neighboor in neighboors:
                H_cords = neighboor.get_coord()
                direction = acceptor_coords - H_cords
                neighboors_H_and_direction.append([H_cords,direction])
            acceptors.append(atom)
            donors += neighboors_H_and_direction
    return [acceptors,donors]
def find_acceptor_and_donors_from_residue_bb(residue):
    acceptors = []
    donors = []
    for atom in residue:
        if atom.element in acceptor_names:
            if not atom.name in backbone_atoms:
                continue
            acceptor_coords = atom.get_coord()
            neighboors = get_bonded_hydrogens(atom, residue, 1.2)
            neighboors_H_and_direction = []
            for neighboor in neighboors:
                H_cords = neighboor.get_coord()
                direction = acceptor_coords - H_cords
                neighboors_H_and_direction.append([H_cords,direction])
            acceptors.append(atom)
            donors += neighboors_H_and_direction
    return [acceptors,donors]
def find_all_acceptors_and_donors(residue):
    acceptors = []
    donors = []
    for atom in residue:
        if atom.element in acceptor_names:
            neighboors = get_bonded_hydrogens(atom, residue, 1.2)
            neighboors_H_and_direction = []
            acceptor_coords = atom.get_coord()
            for neighboor in neighboors:
                H_cords = neighboor.get_coord()
                direction = acceptor_coords - H_cords
                neighboors_H_and_direction.append([H_cords,direction])
            acceptors.append(atom)
            donors += neighboors_H_and_direction
    return [acceptors,donors]
def find_hydrogen_bond(residue1, residue2,is_ligand, bb_flag):
    if bb_flag:
        acceptor_donors_1 = find_acceptor_and_donors_from_residue_bb(residue1)
    else:
        acceptor_donors_1 = find_acceptor_and_donors_from_residue_sc(residue1)
    acceptor_donors_2_all =  find_all_acceptors_and_donors(residue2)
    found_h_bond = False
    for acceptor in acceptor_donors_1[0]:
        for donor in acceptor_donors_2_all[1]:
            if found_h_bond:
                break
            found_h_bond = check_h_bond(donor,acceptor)
        if found_h_bond:
            break
    if found_h_bond:
        return True
    for donor in acceptor_donors_1[1]:
        for acceptor in acceptor_donors_2_all[0]:
            if found_h_bond:
                break
            found_h_bond = check_h_bond(donor, acceptor)
        if found_h_bond:
            break
    if found_h_bond:
        return True
    else:
        return False
def check_h_bond (donor, acceptor):
    donor_direction = donor[1]
    donor_coords = donor[0]
    acceptor_coords = acceptor.get_coord()
    donor_acceptor_direction =  acceptor_coords - donor_coords
    distance = calculate_distance(acceptor_coords,donor_coords)
    angle = calculate_angle(donor_direction, donor_acceptor_direction)
    if 1.5 <= distance <= 4 and 120 < angle < 300:
        return True
    else:
        return False

def find_salt_bridge(residue1, residue2):

    res1_name = residue1.get_resname()
    res2_name = residue2.get_resname()
    if res1_name not in dict_of_pka.keys() or res2_name not in dict_of_pka.keys():
        return False
    for atom1 in residue1:
        for atom2 in residue2:
            if atom1.name in charged_atoms[res1_name] and atom2.name in charged_atoms[res2_name]:
                distance = calculate_distance(atom1.get_coord(), atom2.get_coord())
                charge_1 = dict_of_pka[res1_name]
                charge_2 = dict_of_pka[res2_name]
                if charge_1 * charge_2 == -1 and 1 < distance < 4.5:
                    return True

    return False


acceptor_names = ["F","N","O"]
backbone_atoms = ["N","O"]
dict_of_pka = {
    'GLU': -1,  # Glutamate
    'ASP': -1,  # Aspartate
    'CYS': -1,  # Cysteine
    'TYR': -1,  # Tyrosine
    'LYS': 1,  # Lysine
    'ARG': 1,  # Arginine
    'HIS': 1  # Histidine
}
charged_atoms = {
    'GLU': ['OE1', 'OE2'],  # Glutamate
    'ASP': ['OD1', 'OD2'],  # Aspartate
    'LYS': ['NZ'],          # Lysine
    'ARG': ['NH1', 'NH2', 'NE'],  # Arginine
    'HIS': ['ND1', 'NE2'],  # Histidine
    'CYS': ['SG'],          # Cysteine
    'TYR': ['OH']           # Tyrosine
}