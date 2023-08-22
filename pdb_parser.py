# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 17:32:03 2023

@author: kyle
"""

from __future__ import annotations
import math
import numpy as np
import numpy
import os
import gzip



class Residue:
    def __init__(self, CA : numpy.ndarray = None, N : numpy.ndarray= None, C : numpy.ndarray= None, side_chains: numpy.ndarray= None, previous_residue : Residue = None,post_residue: Residue = None):
        self.CA = CA
        self.N = N
        self.C = C
        self.side_chains = side_chains
        self.post_residue = post_residue
        self.previous_residue = previous_residue
        self.others = []
    def input_atom_info_from_pasred_line(self,parsedLine: PDBLineParser):
        if parsedLine.atom_name == "CA" and self.CA is None:
            self.CA = numpy.array([parsedLine.x_cord,parsedLine.y_cord,parsedLine.z_cord])
        elif parsedLine.atom_name == "C" and self.C is None:
            self.C = numpy.array([parsedLine.x_cord, parsedLine.y_cord, parsedLine.z_cord])
        elif parsedLine.atom_name == "N" and self.N is None:
            self.N = numpy.array([parsedLine.x_cord, parsedLine.y_cord, parsedLine.z_cord])
        else:
            self.others.append([parsedLine.atom_name,numpy.array([parsedLine.x_cord, parsedLine.y_cord, parsedLine.z_cord])])
    def get_phi_angle(self) -> float:
        if self.previous_residue is None:
            return 0.0
        plane1 = numpy.array([self.previous_residue.C, self.N, self.CA])
        plane2 = numpy.array([self.N, self.CA, self.C])
        plane1_norm = numpy.cross(numpy.subtract(plane1[0], plane1[1]), numpy.subtract(plane1[2], plane1[1]))
        plane2_norm = numpy.cross(numpy.subtract(plane2[0], plane2[1]), numpy.subtract(plane2[2], plane2[1]))
        ab = numpy.sum(numpy.multiply(plane1_norm, plane2_norm))
        anorm = numpy.linalg.norm(plane1_norm)
        bnorm = numpy.linalg.norm(plane2_norm)
        return math.degrees(math.acos(ab / (anorm * bnorm)) * -1)
    def get_psi_angle(self) -> float:
        if self.post_residue is None:
            return 0.0
        plane1 = numpy.array([self.N, self.CA, self.C])
        plane2 = numpy.array([self.CA, self.C, self.post_residue.N])
        plane1_norm = numpy.cross(numpy.subtract(plane1[0], plane1[1]), numpy.subtract(plane1[2], plane1[1]))
        plane2_norm = numpy.cross(numpy.subtract(plane2[0], plane2[1]), numpy.subtract(plane2[2], plane2[1]))
        ab = numpy.sum(numpy.multiply(plane1_norm, plane2_norm))
        anorm = numpy.linalg.norm(plane1_norm)
        bnorm = numpy.linalg.norm(plane2_norm)
        return math.degrees(math.acos(ab / (anorm * bnorm)))
    @staticmethod
    def get_phi_angle_from_residues(previous_residue: Residue, current_residue: Residue) -> float:
        plane1 = numpy.array([previous_residue.C,current_residue.N,current_residue.CA])
        plane2 = numpy.array([current_residue.N,current_residue.CA,current_residue.C])
        plane1_norm = numpy.cross(numpy.subtract(plane1[0], plane1[1]), numpy.subtract(plane1[2], plane1[1]))
        plane2_norm = numpy.cross(numpy.subtract(plane2[0], plane2[1]), numpy.subtract(plane2[2], plane2[1]))
        ab = numpy.sum(numpy.multiply(plane1_norm, plane2_norm))
        anorm = numpy.linalg.norm(plane1_norm)
        bnorm = numpy.linalg.norm(plane2_norm)
        return math.degrees(math.acos(ab / (anorm * bnorm))*-1)
    @staticmethod
    def get_psi_angle_from_residues(current_residue: Residue, post_residue: Residue) -> float:
        plane1 = numpy.array([current_residue.N, current_residue.CA,current_residue.C])
        plane2 = numpy.array([current_residue.CA, current_residue.C, post_residue.N])
        plane1_norm = numpy.cross(numpy.subtract(plane1[0], plane1[1]), numpy.subtract(plane1[2], plane1[1]))
        plane2_norm = numpy.cross(numpy.subtract(plane2[0], plane2[1]), numpy.subtract(plane2[2], plane2[1]))
        ab = numpy.sum(numpy.multiply(plane1_norm, plane2_norm))
        anorm = numpy.linalg.norm(plane1_norm)
        bnorm = numpy.linalg.norm(plane2_norm)
        return math.degrees(math.acos(ab / (anorm * bnorm)))
class PDBLineParser():

    def __init__(self, line:str = None):
        self.atom_name: None
        self.element_symbol = None
        self.segment_identifier = None
        self.temp = None
        self.occupany = None
        self.z_cord = None
        self.y_cord = None
        self.x_cord = None
        self.code_for_insertions = None
        self.atom_name = None
        self.residue_sequence_number = None
        self.chain_identifier = None
        self.residue_name = None
        self.alternate_location = None
        self.atom_serial_number = None
        self.atom_type = None
        if line is not None:
            self.line = line
            self.parse_line()
    def parse_line(self):

        if self.line[0:4].replace("\n","").replace(" ","") in "HETATM":
            self.atom_type = self.line[0:6].replace("\n", "").replace(" ", "")
            self.atom_serial_number = int(self.line[6:11].replace("\n","").replace(" ",""))
            self.atom_name = self.line[12:16].replace("\n","").replace(" ","")
            self.alternate_location = self.line[16].replace("\n","").replace(" ","")
            self.residue_name = self.line[17:20].replace("\n","").replace(" ","")
            self.chain_identifier = self.line[21].replace("\n","").replace(" ","")
            self.residue_sequence_number = int(self.line[22:26].replace("\n","").replace(" ",""))
            self.code_for_insertions = self.line[26].replace("\n","").replace(" ","")
            self.x_cord = float(self.line[30:38].replace("\n","").replace(" ",""))
            self.y_cord = float(self.line[38:46].replace("\n","").replace(" ",""))
            self.z_cord = float(self.line[46:54].replace("\n","").replace(" ",""))
            self.occupany = self.line[54:60].replace("\n","").replace(" ","")
            self.temp = float(self.line[60:66].replace("\n","").replace(" ",""))
            self.segment_identifier = self.line[72:76].replace("\n","").replace(" ","")
            self.element_symbol = self.line[76:78].replace("\n","").replace(" ","")
        else:
            self.atom_type = self.line[0:4].replace("\n","").replace(" ","")
            self.atom_serial_number = int(self.line[6:11].replace("\n","").replace(" ",""))
            self.atom_name = self.line[12:16].replace("\n","").replace(" ","")
            self.alternate_location = self.line[16].replace("\n","").replace(" ","")
            self.residue_name = self.line[17:20].replace("\n","").replace(" ","")
            self.chain_identifier = self.line[21].replace("\n","").replace(" ","")
            self.residue_sequence_number = int(self.line[22:26].replace("\n","").replace(" ",""))
            self.code_for_insertions = self.line[26].replace("\n","").replace(" ","")
            self.x_cord = float(self.line[30:38].replace("\n","").replace(" ",""))
            self.y_cord = float(self.line[38:46].replace("\n","").replace(" ",""))
            self.z_cord = float(self.line[46:54].replace("\n","").replace(" ",""))
            self.occupany = self.line[54:60].replace("\n","").replace(" ","")
            self.temp = float(self.line[60:66].replace("\n","").replace(" ",""))
            self.segment_identifier = self.line[72:76].replace("\n","").replace(" ","")
            self.element_symbol = self.line[76:78].replace("\n","").replace(" ","")
class ResidueBuilder():
    allowed_3 = ["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
    residues = []
    def __init__(self, line:str = None):
        self.current_residue = Residue()
        self.previous_residue = None
        self.pdb_line = None
        self.previous_chain = None
        self.current_chain = None
        self.previous_residue_number = None
        self.current_residue_number = None
        self.residues.append(self.current_residue)
        if line is not None:
            self.pdb_line = PDBLineParser(line)
            self.current_chain = self.pdb_line.chain_identifier
            self.current_residue_number = self.pdb_line.residue_sequence_number
    def is_con_amino_acid(self) -> bool:
        if self.pdb_line.residue_name is None:
            return False
        if self.pdb_line.residue_name.upper() in self.allowed_3:
            return True
        return False
    def is_same_chain(self) -> bool:
        if self.previous_chain is None:
            return True
        elif self.previous_chain == self.current_chain:
            return True
        return False
    def is_same_residue(self) -> bool:
        if self.previous_residue_number is None:
            return True
        elif self.previous_residue_number == self.current_residue_number:
            return True
        return False
    def process_new_line(self,line:str):
        self.pdb_line = PDBLineParser(line)
        self.previous_chain = self.current_chain
        self.current_chain = self.pdb_line.chain_identifier
        self.previous_residue_number = self.current_residue_number
        self.current_residue_number = self.pdb_line.residue_sequence_number
        if not self.is_same_residue() or not self.is_same_chain():
            if self.previous_residue is not None:
                if self.is_same_chain():
                    self.previous_residue.post_residue = self.current_residue
            self.previous_residue = self.current_residue
            self.current_residue = Residue()
            self.residues.append(self.current_residue)
            if self.is_same_chain():
                self.current_residue.previous_residue = self.previous_residue
        if not self.is_con_amino_acid():
            return
        self.current_residue.input_atom_info_from_pasred_line(self.pdb_line)
def parse_pdb_to_residues(file):
    resBuilder = ResidueBuilder()
    with open(file,"r") as inputFile:
        residue_number = 1
        chain = "A"
        for line in inputFile:
            if "ATOM " not in line:
                continue
            resBuilder.process_new_line(line.replace("\n",""))