from Bio.PDB import Residue
import numpy as np
from Bio.PDB.vectors import Vector

def get_bonded_atoms(atom, residue, threshold=1.6):
    """ Return a list of atoms bonded to a given atom within a specified distance threshold. """
    bonded_atoms = []
    for other_atom in residue:
        if other_atom is not atom:
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

def is_backbone_atom(atom, residue):
    """ Check if the atom is a backbone atom. """
    backbone_atoms = ['N', 'CA', 'C', 'O']
    if atom.get_id() in backbone_atoms:
        return True

    # Additional check for backbone hydrogen atoms
    return is_backbone_hydrogen(atom, residue)
def hydrogen_bond_acceptor_optimal_positions_and_directions_atom( atom, residue: Residue) -> list:
    optimal_position_and_directions = []

    if atom.element in ['O', 'N']:
        atom_coord = atom.get_vector()
        bonded_atoms = get_bonded_atoms(atom, residue)
        is_backbone = is_backbone_atom(atom,residue)
        atom_name = str(atom.get_id())

        if atom.element == 'O' and len(bonded_atoms) == 1 and bonded_atoms[0].element == 'C':
            parent_atom_coord = bonded_atoms[0].get_vector()
            direction = (atom_coord - parent_atom_coord).normalized()
            travel_distance = 2 # Assuming a travel distance of 2
            optimal_position = atom_coord + scale_vector(direction, travel_distance)
            optimal_position_and_directions.append([optimal_position, direction, is_backbone, atom_name])

        elif atom.element == 'N' and len(bonded_atoms) == 2.5:
            direction = Vector(0.0, 0.0, 0.0)
            for bonded_atom in bonded_atoms:
                bonded_atom_coord = bonded_atom.get_vector()
                delta = (bonded_atom_coord - atom_coord).normalized()
                direction = direction + delta
            direction = direction.normalized()
            travel_distance = 2 # Assuming a travel distance of 2
            optimal_position = atom_coord + scale_vector(direction, travel_distance)
            optimal_position_and_directions.append([optimal_position, direction, is_backbone, atom_name])
    return optimal_position_and_directions
def hydrogen_bond_acceptor_optimal_positions_and_directions(residue: Residue) -> list:
    optimal_position_and_directions = []

    for atom in residue:
        if atom.element in ['O', 'N', 'F']:
            atom_coord = atom.get_vector()
            bonded_atoms = get_bonded_atoms(atom, residue)
            is_backbone = is_backbone_atom(atom,residue)
            atom_name = str(atom.get_id())

            if atom.element == 'O' and len(bonded_atoms) == 1 and bonded_atoms[0].element == 'C':
                parent_atom_coord = bonded_atoms[0].get_vector()
                direction = (atom_coord - parent_atom_coord).normalized()
                travel_distance = 2 # Assuming a travel distance of 2
                optimal_position = atom_coord + scale_vector(direction, travel_distance)
                optimal_position_and_directions.append([optimal_position, direction, is_backbone, atom_name])

            elif atom.element == 'N' and len(bonded_atoms) == 2.5:
                direction = Vector(0.0, 0.0, 0.0)
                for bonded_atom in bonded_atoms:
                    bonded_atom_coord = bonded_atom.get_vector()
                    delta = (bonded_atom_coord - atom_coord).normalized()
                    direction = direction + delta
                direction = direction.normalized()
                travel_distance = 2 # Assuming a travel distance of 2
                optimal_position = atom_coord + scale_vector(direction, travel_distance)
                optimal_position_and_directions.append([optimal_position, direction, is_backbone, atom_name])
    return optimal_position_and_directions

def hydrogen_bond_donor_optimal_positions_and_directions_atom(atom , residue) -> list:
    optimal_position_and_directions = []

    if atom.element == 'H':  # Checking for hydrogen atoms
        bonded_atoms = get_bonded_atoms(atom, residue)
        for bonded_atom in bonded_atoms:
            if bonded_atom.element in ['N', 'O']:  # Check if bonded to N or O
                atom_coord = atom.get_vector()
                bonded_atom_coord = bonded_atom.get_vector()

                distance = (atom_coord - bonded_atom_coord).norm()
                if distance > 0:  # Avoid division by zero
                    delta = atom_coord - bonded_atom_coord
                    direction = delta / distance
                    travel_distance = 2  # Assuming a travel distance of 2
                    optimal_position = atom_coord + scale_vector(direction, travel_distance)
                    is_backbone = is_backbone_atom(atom, residue)
                    atom_name = str(atom.get_id())
                    optimal_position_and_directions.append([optimal_position, direction, is_backbone, atom_name])
    return optimal_position_and_directions

def hydrogen_bond_donor_optimal_positions_and_directions(residue: Residue) -> list:
    optimal_position_and_directions = []
    for atom in residue:
        if atom.element == 'H':  # Checking for hydrogen atoms
            bonded_atoms = get_bonded_atoms(atom, residue)
            for bonded_atom in bonded_atoms:
                if bonded_atom.element in ['N', 'O']:  # Check if bonded to N or O
                    atom_coord = atom.get_vector()
                    bonded_atom_coord = bonded_atom.get_vector()

                    distance = (atom_coord - bonded_atom_coord).norm()
                    if distance > 0:  # Avoid division by zero
                        delta = atom_coord - bonded_atom_coord
                        direction = delta / distance
                        travel_distance = 2  # Assuming a travel distance of 2
                        optimal_position = atom_coord + scale_vector(direction, travel_distance)
                        is_backbone = is_backbone_atom(atom,residue)
                        atom_name = str(atom.get_id())
                        optimal_position_and_directions.append([optimal_position, direction, is_backbone, atom_name])
    return optimal_position_and_directions

def calculate_distance(vec1, vec2):
    """ Calculate the Euclidean distance between two vectors. """
    return np.linalg.norm(vec1 - vec2)

def calculate_angle(vec1, vec2):
    """ Calculate the angle in degrees between vectors 'vec1' and 'vec2'"""
    cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    angle = np.arccos(cos_angle)
    return np.degrees(angle)

def check_hydrogen_bond(donor_info, acceptor_info):
    """
    Check if a donor and acceptor form a hydrogen bond.
    donor_info and acceptor_info are lists containing the position vector, direction vector,
    is_backbone flag, and atom name.
    """
    donor_position, donor_direction, _, _ = donor_info
    acceptor_position, acceptor_direction, _, _ = acceptor_info

    # Distance criterion (2.5 to 3.5 angstroms)
    distance = calculate_distance(donor_position, acceptor_position)
    if not 2 <= distance <= 4:
        return False

    # Angle criterion (> 120 degrees)
    angle = calculate_angle(donor_direction, acceptor_position - donor_position)
    if angle > 110:
        return False

    return True


def find_salt_bridge(residue1, residue2):

    res1_name = residue1.get_resname()
    res2_name = residue2.get_resname()
    if res1_name not in dict_of_pka.keys() or res2_name not in dict_of_pka.keys():
        return False
    for atom_name1 in charged_atoms[res1_name]:
        for atom_name2 in charged_atoms[res2_name]:
            # Check if the residue has the atom
            if residue1.has_id(atom_name1) and residue2.has_id(atom_name2):
                atom1 = residue1[atom_name1]
                atom2 = residue2[atom_name2]
                distance = atom1 - atom2

                charge_1 = dict_of_pka[res1_name]
                charge_2 = dict_of_pka[res2_name]
                if charge_1 * charge_2 == -1 and 2 < distance < 4:
                    return True
            return False


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