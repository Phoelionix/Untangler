#===========================================================================
# Untangler: Free ensemble models from local minima with the wrong altlocs 
# Copyright (C)  2025 Spencer Passmore (spencerpassmore@swin.edu.au)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#===========================================================================

# Purpose: Get geometric badness for all possible "groups" of atoms that 
# we have geometric measures for, to pass to LinearOptimizer.Solver.
#
# As of July 2025, badness is the sum of the (weighted... possibly
# incorrectly) squares of the angle and bond deviations from ideal. 
# 
# TODO  
# - Change to proper chi-squared (missing denominator presently). 
# - implement all the other wonderful measures we have available. -- Noth Input and Solver
# - Clean up chunk ID system. -- Both Input and Solver
# 

from Bio.PDB import PDBParser,Structure
from Bio.PDB.Atom import Atom,DisorderedAtom
from Bio.PDB.Residue import Residue # Note we don't want the disorderedresidue here, it refers to different residue types for same res seq num.
from Bio.PDB.StructureBuilder import StructureBuilder
from pulp import *
import UntangleFunctions 
from LinearOptimizer.Swapper import Swapper
import os
import random
import numpy as np
import itertools

class OrderedAtomLookup: #TODO pandas?
    def __init__(self,atoms:list[DisorderedAtom],protein=True,waters=False): # TODO type hint for sequence?
        self.disordered_waters=[]
        self.ordered_atoms:list[Atom] = []
        self.serial_num_to_disordered_num_dict={}
        self.residue_nums=[]
        self.residue_sources=[]
        self.altlocs=[]
        self.res_names:dict[int,str]={}
        for disorderedAtom in atoms:
            is_water= UntangleFunctions.res_is_water(disorderedAtom.get_parent())
            if is_water:
                self.disordered_waters.append(disorderedAtom)
            if is_water and not waters:
                continue
            if not is_water and not protein:
                continue
            
            assert type(disorderedAtom)==DisorderedAtom, type(disorderedAtom)  # Not sure if a single altloc atom will still be stored as a disorderedatom by Bio.PDB
            for orderedAtom in disorderedAtom:
                self.ordered_atoms.append(orderedAtom)
                self.serial_num_to_disordered_num_dict[orderedAtom.get_serial_number()]=disorderedAtom.get_serial_number()
                altloc = orderedAtom.get_altloc()
                if altloc not in self.altlocs:
                    self.altlocs.append(altloc)
            res_num=OrderedAtomLookup.atom_res_seq_num(disorderedAtom)
            if res_num not in self.residue_nums:
                assert len(self.residue_nums)== 0 or self.residue_nums[-1] == res_num-1
                self.residue_nums.append(res_num)
                self.res_names[res_num]=disorderedAtom.get_parent().get_resname()
                self.residue_sources.append(disorderedAtom.get_parent())
    @staticmethod
    def atom_res_seq_num(atom:Atom)->int:
        return atom.get_parent().get_id()[1]
    @staticmethod
    def atom_res_name(atom:Atom)->str:
        return atom.get_parent().get_resname()
    def get_disordered_serial_number(self,atom:Atom)->int:
        return self.serial_num_to_disordered_num_dict[atom.get_serial_number()] 
    def get_residue_nums(self)->list[int]:
        return self.residue_nums
    def get_residue_sources(self)->list[Residue]:
        return self.residue_sources
    def get_altlocs(self)->list[str]:
        return self.altlocs
    def select_atoms_by(self,res_nums=None,serial_numbers=None,names=None,
                        disordered_serial_numbers=None,altlocs=None,
                        waters=True, protein=True, exclude_H=False)->list[Atom]:
        atom_selection = []
        for atom in self.ordered_atoms:
            is_water = UntangleFunctions.res_is_water(atom.get_parent())
            if not waters and is_water:
                continue
            if not protein and not is_water:
                continue
            if exclude_H and atom.element == "H":
                continue
            _res_num = OrderedAtomLookup.atom_res_seq_num(atom)
            _name = atom.get_name()
            _altloc = atom.get_altloc()
            _serial_number=atom.get_serial_number()
            _disordered_serial_number = self.get_disordered_serial_number(atom)
            
            for allowed_values,atom_param in zip(
            (altlocs,res_nums,serial_numbers,names,disordered_serial_numbers),
            (_altloc,_res_num,_serial_number,_name,_disordered_serial_number)):
                if allowed_values is not None and (atom_param not in allowed_values):
                    break
            else:
                atom_selection.append(atom)

        return atom_selection

# Want to calculate wE for chunks of atoms and residues and side chains and main chains 
# we then get the optimizer to minimize when choosing 
# E.g. calculate chunk over two halves of residue
# So if route includes those halves, the cost includes those.
class Chunk():
    def __init__(self,depth_tag,site_num,altloc,atoms:list[Atom],start_resnum,end_resnum,start_resname,end_resname):
        self.altloc=altloc
        self.site_num_at_depth=site_num
        self.atoms_dict:dict[int,list[Atom]]={}
        for atom in atoms:
            res_num = OrderedAtomLookup.atom_res_seq_num(atom)
            if res_num not in atoms:
                self.atoms_dict[res_num] = []
            self.atoms_dict[res_num].append(atom)
        self.start_resname=start_resname
        self.end_resname=end_resname
        self.start_resnum = start_resnum
        self.end_resnum = end_resnum
        self.depth_tag=depth_tag
    def get_site_num(self):
        return self.site_num_at_depth
    def get_altloc(self):
        return self.altloc
    def get_disordered_tag(self):
        pass
        # Making it clear this code isn't run because I made this so confusing. But might be useful idea e.g. for van der waals
        #return  f"{self.start_resnum}.{self.start_resname}x{self.end_resnum}.{self.end_resname}"
    def unique_id(self):
        pass
        #return f"{self.depth_tag}&{self.get_site_num()}.{self.get_altloc()}&{self.get_disordered_tag()}"
    
    def get_atoms_dict(self):
      return self.atoms_dict


class OrderedResidue(Chunk):
    def __init__(self,altloc,resnum,referenceResidue,atoms,site_num_override=None):
        self.resnum=resnum
        self._res=referenceResidue
        site_num = self.get_resnum()
        if site_num_override is not None:
            site_num =site_num_override
        super().__init__("RESI",site_num,altloc,atoms,
                         self.get_resnum(),self.get_resnum(),
                         self.get_resname(),self.get_resname())
    def get_resnum(self):
        return self.resnum
    def get_resname(self):
        return self._res.get_resname()
    def get_id(self):
        return self._res.get_id()
    def get_atoms(self):
        return self.atoms_dict[self.get_resname()]
    def get_disordered_tag(self):
        assert False 
        return  f"{self.get_resnum()}.{self.get_resname()}"
    


    
class ConstraintsHandler:
    class Constraint():
        def __init__(self,atom_ids,ideal,weight):
            self.atom_id = atom_ids
            self.ideal = ideal
            self.weight = weight
            self.kind=None # Why do this when can just check type??
        def num_atoms(self):
            return len(self.atom_id)
        def atom_names(self):
            return [pdb[:4].strip() for pdb in self.atom_id]
        def residues(self):
            return [int(pdb.strip().split()[-1]) for pdb in self.atom_id]
        def get_distance(self,atoms:list[Atom]):
            raise Exception("abstract method")
        def get_ordered_atoms(self,candidate_atoms:list[Atom])->list[Atom]:
            other_name_and_resnum = [(a.get_name(),OrderedAtomLookup.atom_res_seq_num(a)) for a in candidate_atoms] 
            #print(self.num_atoms(),other_name_and_resnum)
            #print(set(other_name_and_resnum))
            assert len(other_name_and_resnum)==len(set(other_name_and_resnum))==self.num_atoms()
            ordered_atoms = [] # corresponding to order of self.atom_id. Important for e.g. bond angle
            for name_num in zip(self.atom_names(),self.residues()):
                    for i, other_name_num in enumerate(other_name_and_resnum):
                        if name_num==other_name_num:
                            ordered_atoms.append(candidate_atoms[i])
            if len(ordered_atoms)!=self.num_atoms():
                return None
            assert len(np.unique(ordered_atoms))==self.num_atoms()
            return ordered_atoms

        
    class BondConstraint(Constraint):
        def __init__(self,atom_ids,ideal,weight):
            assert len (atom_ids)==2
            super().__init__(atom_ids,ideal,weight)
            self.kind="Bond"
        @staticmethod
        def separation(a:Atom,b:Atom):
            return np.sqrt(np.sum((a.get_coord()-b.get_coord())**2))
        def get_distance(self,atoms:list[Atom]):
            ordered_atoms = self.get_ordered_atoms(atoms)
            if ordered_atoms is None:
                return None
            a,b = ordered_atoms            
            sqr_deviation = (self.ideal-self.separation(a,b))**2
            return sqr_deviation *self.weight
        
    class AngleConstraint(Constraint):
        def __init__(self,atom_ids,ideal,weight):
            assert len (atom_ids)==3
            super().__init__(atom_ids,ideal,weight)
            self.kind="Angle"
        @staticmethod
        def angle(a:Atom,b:Atom,c:Atom):
            v1 = a.get_coord() - b.get_coord()
            v2 = c.get_coord()-b.get_coord()
            def unit_vector(vector):
                return vector / np.linalg.norm(vector)
            v1_u = unit_vector(v1)
            v2_u = unit_vector(v2)
            return np.arccos(np.dot(v1_u, v2_u))*180/np.pi
            #return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
        def get_distance(self,atoms:list[Atom]):
            ordered_atoms = self.get_ordered_atoms(atoms)
            if ordered_atoms is None:
                return None
            a,b,c = ordered_atoms
            sqr_deviation = (self.ideal-self.angle(a,b,c))**2

            return sqr_deviation*self.weight
        
    class ClashConstraint(Constraint):
        default_weight=1
        def __init__(self,atom_ids,altlocs,badness,weight=None):
            self.altlocs = altlocs
            self.badness = badness
            if weight is None:
                weight = ConstraintsHandler.NonbondConstraint.default_weight
            super().__init__(atom_ids,None,weight)
            self.kind="Clash" 
        def get_distance(self,atoms:list[Atom]):
            ordered_atoms = self.get_ordered_atoms(atoms)
            if ordered_atoms is None:
                return None
            a,b = ordered_atoms
            if a.get_altloc() == self.altlocs[0] and b.get_altloc() == self.altlocs[1]:
                return self.badness
            return 0 


    class NonbondConstraint(Constraint):
        default_weight=1
        def __init__(self,atom_ids,ideal_separation,weight=None):
            if weight is None:
                weight = ConstraintsHandler.NonbondConstraint.default_weight
            super().__init__(atom_ids,ideal_separation,weight)
            self.kind="Nonbond"
        
        @staticmethod
        def lennard_jones(r,r0):
            #obs=$1;ideal=$2;sigma=1;energy=lj(obs,ideal)
            # function lj0(r,r0) {if(r==0)return 1e40;return 4*((r0*2^(-1./6)/r)^12-(r0*2^(-1./6)/r)^6)}\
            # function lj(r,r0) {return lj0(r,r0)-lj0(6,r0)}'
            def lj0(r,r0):
                if r == 0:
                    return 1e40
                return 4*((r0*2**(-1/6)/r)**12-(r0*2**(-1/6)/r)**6)
            return lj0(r,r0)-lj0(6,r0)  
        @staticmethod
        def badness(r,r0,neg_badness_limit=-0.1): # TODO implement properly as in untangle_score.csh. 
            # negative when separation is greater than ideal.
            assert neg_badness_limit < 0 
            return abs(max(neg_badness_limit,ConstraintsHandler.NonbondConstraint.lennard_jones(r,r0)))
        def get_distance(self,atoms:list[Atom]):
            ordered_atoms = self.get_ordered_atoms(atoms)
            if ordered_atoms is None:
                return None
            a,b = ordered_atoms            
            r = ConstraintsHandler.BondConstraint.separation(a,b)
            r0 = self.ideal
            energy = ConstraintsHandler.NonbondConstraint.badness(r,r0)
            return energy * self.weight
        

    def __init__(self,constraints:list[Constraint]=[]):
        self.constraints:list[ConstraintsHandler.Constraint]=constraints

    def atom_in_constraints(self,atom_name,res_num):
        for constraint in self.constraints:
            for other_name,other_res_num in zip(constraint.atom_names(), constraint.residues()):
                if atom_name== other_name and res_num==other_res_num:
                    return True
        return False

    def get_constraints(self,atom_name,res_num):
        atom_constraints:list[ConstraintsHandler.Constraint] = [] 
        for constraint in self.constraints:
            for other_atom_name, other_res_num in zip(constraint.atom_names(),constraint.residues()):
                # print("---")
                # print(other_atom_name)
                # print(other_res_num)
                # print(atom_name)
                # print(res_num)
                if other_atom_name == atom_name and other_res_num==res_num:
                    atom_constraints.append(constraint)
        return ConstraintsHandler(atom_constraints)

            
    def load_all_constraints(self,constraints_file,nonbond_scores_path,nonbond_water_flipped_scores_path,water_clashes,ordered_atom_lookup:OrderedAtomLookup,water_water_nonbond:bool):
        print(f"Loading constraints from {constraints_file}")

        self.constraints: list[ConstraintsHandler.Constraint]=[]
        with open(constraints_file,"r") as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if line.startswith("bond"):
                    constraint = lines[i:i+4]
                    pdb1=constraint[0].strip().split("\"")[1]
                    pdb2=constraint[1].strip().split("\"")[1]
                    ideal,  _,  _, _,  weight, _ = constraint[3].strip().split()
                    altloc = pdb1.strip()[3]
                    if altloc!="A": # just look at A altloc to get constraint. TODO check constraint isn't in self.constraints instead.
                        continue
                    ideal,weight = float(ideal), float(weight)
                    #print(pdb1,"|","|",pdb2,"|",ideal,"|",weight)
                    self.constraints.append(ConstraintsHandler.BondConstraint((pdb1,pdb2),ideal,weight))  
                if line.startswith("angle"):
                    constraint = lines[i:i+5]
                    pdb1=constraint[0].strip().split("\"")[1]
                    pdb2=constraint[1].strip().split("\"")[1]
                    pdb3=constraint[2].strip().split("\"")[1]
                    ideal,  _,  _, _,  weight, _ = constraint[4].strip().split()
                    altloc = pdb1.strip()[3]
                    if altloc!="A":
                        continue
                    ideal,weight = float(ideal), float(weight)
                    #print(pdb1,"|","|",pdb2,"|",ideal,"|",weight)
                    self.constraints.append(ConstraintsHandler.AngleConstraint((pdb1,pdb2,pdb3),ideal,weight))  
       
       
        # Add nonbonds that are flagged as issues for current structure AND when waters are swapped
        print("WARNING: assuming residue numbers are all unique")
        pdb_ids_added = []
        for file,flipped_water in zip((nonbond_scores_path, nonbond_water_flipped_scores_path),(False,True)):
        

            if file is None:
                continue
            with open(file) as f:

                lines = f.readlines()
                for i, line in enumerate(lines):
                    if line.startswith("nonbonded"):
                        constraint = lines[i:i+4]
                        #     pdb1 = f"{name}     ARES     A      {res_num}"
                        pdb1=constraint[0].strip().split("\"")[1]
                        pdb2=constraint[1].strip().split("\"")[1]
                        name1 = pdb1.strip()[0]
                        name2 = pdb2.strip()[0]
                        resnum1 = int(pdb1.strip().split()[-1])
                        resnum2 = int(pdb2.strip().split()[-1])
                        resname1 = pdb1.strip().split()[1][1:]
                        resname2 = pdb2.strip().split()[1][1:]
                        both_waters = resname1 == resname2=="HOH"
                        if both_waters and not water_water_nonbond:
                            continue
                            
                        
                        # I mean, they should be bonded. But just in case. 
                        if resnum1 == resnum2:
                            if name1[0] == "H":
                                nonH_in_res = ordered_atom_lookup.select_atoms_by(res_nums=[resnum1],exclude_H=True)
                                if nonH_in_res == name2:
                                    continue


                        separation = constraint[3].strip().split()[0]
                        ideal = constraint[3].strip().split()[1]
                        if separation > ideal:
                            continue 
                        altloc = pdb1.strip()[3]
                        if altloc!="A": 
                            pass
                            #continue
                        ideal = float(ideal)
                        pdb_ids = (pdb1,pdb2)
                        pdb_ids_flipped = (pdb2,pdb1)
                        if pdb_ids not in pdb_ids_added and pdb_ids_flipped not in pdb_ids_added:
                            self.constraints.append(ConstraintsHandler.NonbondConstraint(pdb_ids,ideal))
                            pdb_ids_added.append(pdb_ids)          


                #for line in f:
                    #ideal_separation = float(line.strip().split()[4])
                    #atom_a,atom_b = line.strip().split()[-2:]
                    # for atom in atom_a,atom_b:
                    #     name, res_num = atom.split("_")[0],int(atom.split("_")[1])
                    #     if name[0] == "H":
                    #         nonH_in_res = ordered_atom_lookup.select_atoms_by(res_nums=[res_num],exclude_H=True)
                    #         name = UntangleFunctions.H_get_parent_fullname(name,set([nonH.get_fullname() for nonH in nonH_in_res]))
                    #     # if name == "O":
                    #     #     # assume water... # FIXME # distinguish between water 
                    #     #     continue 
                    #     # TODO we are assuming no duplicates
                    #     pdb1 = f"{name}     ARES     A      {res_num}"
                    #     self.constraints.append(ConstraintsHandler.NonbondConstraint([pdb1,pdb1],file,flipped,badness))
                    #pdb_ids = []
                    # for atom in atom_a,atom_b:
                    #     name, res_num = atom.split("_")[0],int(atom.split("_")[1])
                    #     if name[0] == "H":
                    #         nonH_in_res = ordered_atom_lookup.select_atoms_by(res_nums=[res_num],exclude_H=True)
                    #         name = UntangleFunctions.H_get_parent_fullname(name,set([nonH.get_fullname() for nonH in nonH_in_res]))
                    #     # if name == "O":
                    #     #     # assume water... # FIXME # distinguish between water 
                    #     #     continue 
                    #     # TODO we are assuming no duplicates
                    #     pdb_ids.append(f"{name}     ARES     A      {res_num}")
                    # if pdb_ids not in pdb_ids_added:
                    #     self.constraints.append(ConstraintsHandler.NonbondConstraint(pdb_ids,ideal_separation))
                    #     pdb_ids_added.append(pdb_ids)
        num_nonbonded_from_geo = len(pdb_ids_added)
        waters = ordered_atom_lookup.select_atoms_by(protein=False)
        ideal_water_separation=2.200
        for atom in waters:
            for other_atom in waters:
                if other_atom == atom:
                    continue
                if ConstraintsHandler.BondConstraint.separation(atom,other_atom) < 8:
                    pdb1 = f"{atom.name}     ARES     A      {OrderedAtomLookup.atom_res_seq_num(atom)}"
                    pdb2 = f"{other_atom.name}     ARES     A      {OrderedAtomLookup.atom_res_seq_num(other_atom)}"
                    pdb_ids = (pdb1,pdb2)
                    pdb_ids_flipped = (pdb2,pdb1)
                    if pdb_ids not in pdb_ids_added and pdb_ids_flipped not in pdb_ids_added:
                        self.constraints.append(ConstraintsHandler.NonbondConstraint(pdb_ids,ideal_water_separation))
                        pdb_ids_added.append(pdb_ids)     
        num_nonbonded_extra = len(pdb_ids_added)-num_nonbonded_from_geo

        print(f"Nonbonded from geo: {num_nonbonded_from_geo}, extra: {num_nonbonded_extra}")

        for name, res_num, badness, altloc in water_clashes:
            pdb_ids = [f"{n}     ARES     A      {r}" for (n,r) in zip(name,res_num)]
            self.constraints.append(ConstraintsHandler.ClashConstraint(pdb_ids,altloc,badness))

# Ugh never do inheritance kids. TODO refactor to composition.
class AtomChunk(OrderedResidue):
    # Just an atom c:
    def __init__(self,site_num,altloc,resnum,referenceResidue,atom:Atom,is_water:bool,constraints_handler):
        self.name=atom.get_name()
        self.is_water = is_water
        self.element = atom.element
        self.coord = atom.get_coord()
        super().__init__(altloc,resnum,referenceResidue,[atom],site_num)
        self.generate_bonds(constraints_handler)
    def generate_bonds(self,constraints_handler:ConstraintsHandler):
        self.constraints_holder:ConstraintsHandler = \
            constraints_handler.get_constraints(self.name,self.resnum)
    def get_coord(self):
        return self.coord
    def get_disordered_tag(self):
        return  f"{self.get_resnum()}.{self.name}"
    def unique_id(self):
        return f"{self.depth_tag}&{self.get_site_num()}.{self.get_altloc()}&{self.get_disordered_tag()}"


class MTSP_Solver:
    ### This commented out code computes wE for combinations of larger groups of atoms.  
    # Could be useful in future as a quick coarse step... but would probably be better 
    # to build up from the atom-site approach, coding the wE measure directly.   
    ### 

    # class ChunkConnection():
    #     def __init__(self,A:Chunk,B:Chunk):
    #         self.start = A.unique_id()
    #         self.end = B.unique_id()
    #         self.altlocs=[A.get_altloc(),B.get_altloc()]
    #         self.res_nums=[A.start_resnum,B.end_resnum]
    #         assert A.start_resnum <= B.end_resnum
    #         #self.res_init_args = [[res.get_resname(),*res.get_id()] for res in (A,B)]
    #         self.ts_distance=None  # NOTE As in the travelling salesman problem sense
    #         self.residue_atoms:list[dict[int,list[Atom]]]=[A.get_atoms_dict(),B.get_atoms_dict()]
    #         self.dry_calc=False
    #     def calculated_dry(self):
    #         return self.dry_calc
    #     def calculate_distance(self,parent_structure_path,tmp_out_folder_path,disordered_waters:list[DisorderedAtom],dry=False,quick_wE=False):
    #         # NOTE disatnces in the travelling salesman problem sense
    #         print(f"Computing wE of chunks {self.start} & {self.end}")
            
            
    #         assert parent_structure_path[-4:]==".pdb"
    #         parent_model_handle = os.path.basename(parent_structure_path)[:-4]
    #         builder = StructureBuilder()
    #         builder.init_structure("duo")
    #         builder.init_model("M")
    #         builder.init_chain("D")
            
    #         for residue_atoms_dict in self.residue_atoms:
    #             for res_num, residue_atoms in residue_atoms_dict.items(): # NOT necessarily all residue atoms
    #                 init_res_args = (residue_atoms[0].get_parent().get_resname(),
    #                                  *residue_atoms[0].get_parent().get_id())
    #                 builder.init_seg(res_num)
    #                 builder.init_residue(*init_res_args)
    #                 for atom in residue_atoms:
    #                     builder.init_atom(
    #                         atom.get_name(),atom.get_coord(),atom.get_bfactor(),
    #                         occupancy=atom.occupancy,altloc=atom.get_altloc(),fullname=atom.get_fullname(),
    #                         element=atom.element
    #                     )
    #         #builder.init_chain("S")
    #         for disordered_water in disordered_waters:
    #             builder.init_seg(disordered_water.get_parent().get_id()[1])
    #             #builder.init_residue(disordered_water.get_parent().get_resname(),*disordered_water.get_parent().get_id())
    #             builder.init_residue(disordered_water.get_parent().get_resname(),"W",*disordered_water.get_parent().get_id()[1:])
    #             for atom in disordered_water:
    #                 builder.init_atom(
    #                 atom.get_name(),atom.get_coord(),atom.get_bfactor(),
    #                 occupancy=atom.occupancy,altloc=atom.get_altloc(),fullname=atom.get_fullname(),
    #                 element=atom.element
    #                 )                   

    #         connection_structure_save_path =tmp_out_folder_path+"tmpXtion_"+parent_model_handle+".pdb"
    #         UntangleFunctions.save_structure(builder.get_structure(),parent_structure_path,connection_structure_save_path)
    #         if dry:
    #             self.ts_distance = random.random()
    #             self.dry_calc=True
    #             return
    #         _,self.ts_distance,_ = UntangleFunctions.assess_geometry_wE(tmp_out_folder_path,connection_structure_save_path,phenixgeometry_only=quick_wE) 

    class AtomChunkConnection():
        def __init__(self, atom_chunks:list[AtomChunk],ts_distance,connection_type,hydrogen_names):
            self.atom_chunks = atom_chunks
            self.from_altlocs=[a.get_altloc() for a in atom_chunks]
            self.res_nums=[a.get_resnum() for a in atom_chunks]
            self.connection_type=connection_type
            self.ts_distance=ts_distance  # NOTE As in the travelling salesman problem sense
            self.hydrogen_tag=""
            if len(hydrogen_names)>0:
                self.hydrogen_tag = "_"+''.join(hydrogen_names)
            self.hydrogen_name_set = set(hydrogen_names)
        def get_disordered_connection_id(self):
            return f"{self.connection_type}{self.hydrogen_tag}_{'_'.join([a_chunk.get_disordered_tag() for a_chunk in self.atom_chunks])}"
    class AtomChunkNonBondConnection(AtomChunkConnection):
        def __init__(self, atom_chunks:list[AtomChunk],ts_distance,connection_type,hydrogen_names):
            assert len(atom_chunks) == 2
            super().__init__(atom_chunks,ts_distance,connection_type,hydrogen_names)
            self.to_altloc = self.from_altlocs[1]
            self.from_altlocs = [self.from_altlocs[0]]
            self.from_altloc = self.from_altlocs[0]
            self.site_name = atom_chunks[0].get_disordered_tag()
            assert self.site_name == atom_chunks[1].get_disordered_tag()


    def __init__(self,pdb_file_path:str, align_uncertainty=False,ignore_waters=False):
        # Note if we ignore waters then we aren't considering nonbond clashes between macromolecule and water.
        self.model_path = pdb_file_path
        original_structure = PDBParser().get_structure("struct",pdb_file_path)
        if align_uncertainty:
            self.align_uncertainty(original_structure)
        self.ordered_atom_lookup = OrderedAtomLookup(original_structure.get_atoms(),protein=True,waters=not ignore_waters)   
    def align_uncertainty(self,structure:Structure.Structure):
        # in x-ray data and geom.
        # Experimental and doesn't help.
        for atom in structure.get_atoms():
            a,b = atom.__iter__()
            coord_diff = np.sqrt(np.sum((a.get_coord()-b.get_coord())**2))
            mean_coord = (a.get_coord()+b.get_coord())/2
            if coord_diff < 0.08: #TODO resolution based #TODO properly do the rotation and find minimum method
                a.coord = b.coord =  mean_coord
                atom.coord = mean_coord


    def calculate_paths(self,quick_wE=False,dry_run=False,atoms_only=True,
                        clash_punish_thing=False,nonbonds=True,water_water_nonbond=True)->tuple[list[Chunk],dict[str,list[AtomChunkConnection]]]: #disorderedResidues:list[Residue]
        print("Calculating geometric costs for all possible connections between chunks of atoms (pairs for bonds, triplets for angles, etc.)")
        
        tmp_out_folder_path=UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+"LinearOptimizer/tmp_out/"
        os.makedirs(tmp_out_folder_path,exist_ok=True)
        
        assert atoms_only, "solving with wEs of chunks of residues not supported right now"

        chunk_sets = []
        connection_types=[]
        if not atoms_only:
            assert False
            orderedResidues: list[OrderedResidue]= [] # residues for each conformation
            self.dry_run=dry_run
            for res_num,referenceResidue in zip(self.ordered_atom_lookup.get_residue_nums(),self.ordered_atom_lookup.get_residue_sources()):
                assert not UntangleFunctions.res_is_water(referenceResidue)
                for altloc in self.ordered_atom_lookup.get_altlocs():
                    atoms = self.ordered_atom_lookup.select_atoms_by(res_nums=[res_num],altlocs=[altloc])
                    orderedResidues.append(OrderedResidue(altloc,res_num,referenceResidue,atoms))
            chunk_sets.append(orderedResidues)
            connection_types.append(MTSP_Solver.ChunkConnection)


            #Triplets of residues
            tripletResidues: list[Chunk]= []
            n_size=3
            for i,start_resnum in enumerate(self.ordered_atom_lookup.get_residue_nums()[::n_size]):
                end_resnum=min(start_resnum+n_size-1,max(self.ordered_atom_lookup.get_residue_nums()))
                start_resname = self.ordered_atom_lookup.res_names[start_resnum]
                end_resname = self.ordered_atom_lookup.res_names[end_resnum]
                for altloc in self.ordered_atom_lookup.get_altlocs():
                    atoms = self.ordered_atom_lookup.select_atoms_by(res_nums=range(start_resnum,end_resnum+1),altlocs=[altloc])
                    triplet = Chunk("TRPL",int((i+2)/2),altloc,atoms,start_resnum,end_resnum,start_resname,end_resname)
                    tripletResidues.append(triplet)
            # half-halfs of residues
            chunk_sets.append(tripletResidues)
            connection_types.append(MTSP_Solver.ChunkConnection)

            tripletAtoms = "TODO"

        def atom_id_from_params(atom_name,atom_altloc,atom_res_seq_num):
            return f"{atom_name,}.{atom_altloc}{atom_res_seq_num}"
        atom_chunks:dict[str,AtomChunk] = {}

        def atom_id(atom:Atom):
            #return f"{atom.get_name()}.{atom.get_altloc()}{OrderedAtomLookup.atom_res_seq_num(atom)}"
            return atom_id_from_params(atom.get_name(),atom.get_altloc(),OrderedAtomLookup.atom_res_seq_num(atom))
        
        nonbond_scores_path=nonbond_water_flipped_scores_path=None
        assert self.model_path[-4:]==".pdb"
        model_handle = os.path.basename(self.model_path)[:-4]

        if nonbonds:
            # generate nonbond files
            debug_no_wE = False
            geo_log_out_folder = UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+"StructureGeneration/HoltonOutputs/"
            if not debug_no_wE:
                UntangleFunctions.assess_geometry_wE(geo_log_out_folder,self.model_path) 
            #nonbond_scores_path = UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+f"StructureGeneration/HoltonOutputs/{model_handle}_scorednonbond.txt"
            nonbond_scores_path = f"{UntangleFunctions.UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/HoltonOutputs/{model_handle}.geo"
            


            # CLASH 26.6682 0.79 | A  62 AARG  HD2| S 128 AHOH  O 
            # Terrible code. XXX
            def get_water_clashes(handle,waters_flipped:bool)->list:
                clashes = []
                nonbond_path = UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+f"StructureGeneration/HoltonOutputs/{handle}_scorednonbond.txt"
                clashes_path = UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+f"StructureGeneration/HoltonOutputs/{handle}_clashes.txt"
                with open(clashes_path) as clashes_file:
                    for line in clashes_file:
                        #pdb1 = f"{name}     ARES     A      {res_num}"
                        res_name = [entry.split()[-2][1:] for entry in line.split("|")[1:]] 
                        if not res_name.count("HOH")==1:
                            continue
                        res_num = [int(entry.split()[1]) for entry in line.split("|")[1:]] 
                        altloc = [entry.split()[-2][0] for entry in line.split("|")[1:]] 
                        name = [entry.strip().split()[-1] for entry in line.split("|")[1:]] 
                        if waters_flipped:
                            water_idx = res_name.index("HOH")
                            altloc[water_idx] = {"A":"B","B":"A"}[altloc[water_idx]]
                        badness = float(line.split()[1])
                        #new_line = f"CLASH   {badness} XXXXX XXXXX XXXXX X |  {name[0]}_{res_num[0]} {name[1]}_{res_num[1]}"
                        clashes.append((name, res_num, badness, altloc))
                return clashes
            # if clash_punish_thing:
            #     assert False
            #     nonbond_scores_path = add_clashes_to_nonbond(model_handle)
            
            model_water_swapped_handle=model_handle+"WS"
            model_water_swapped_path=UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+f"StructureGeneration/HoltonOutputs/{model_water_swapped_handle}.pdb"
            if not debug_no_wE:
                Swapper.MakeSwapWaterFile(self.model_path,model_water_swapped_path)
                UntangleFunctions.assess_geometry_wE(geo_log_out_folder,model_water_swapped_path) 

            water_clashes =  get_water_clashes(model_handle,False) + get_water_clashes(model_water_swapped_handle,True)
            
            #nonbond_water_flipped_scores_path = UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+f"StructureGeneration/HoltonOutputs/{model_water_swapped_handle}_scorednonbond.txt"
            nonbond_water_flipped_scores_path = f"{UntangleFunctions.UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/HoltonOutputs/{model_water_swapped_handle}.geo"


        constraints_file = f"{UntangleFunctions.UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/HoltonOutputs/{model_handle}.geo" # NOTE we only read ideal and weights.
        constraints_handler=ConstraintsHandler()
        constraints_handler.load_all_constraints(constraints_file,nonbond_scores_path,nonbond_water_flipped_scores_path,water_clashes,self.ordered_atom_lookup,water_water_nonbond=water_water_nonbond)
        #for n,atom in enumerate(self.ordered_atom_lookup.select_atoms_by(names=["CA","C","N"])):
        for n,atom in enumerate(self.ordered_atom_lookup.select_atoms_by()):
            if atom.element=="H":
                continue
            is_water = UntangleFunctions.res_is_water(atom.get_parent())
            atom_chunks[atom_id(atom)]=AtomChunk(
                                         int((n+2)/2),
                                         atom.get_altloc(),
                                         OrderedAtomLookup.atom_res_seq_num(atom),
                                         atom.get_parent(),
                                         atom,
                                         is_water,
                                         constraints_handler
                                         )
        chunk_sets.append(atom_chunks)
        connection_types.append(MTSP_Solver.AtomChunkConnection)


        ############ Calculate wE for different combinations (don't need known constraints).
        '''
        possible_connections:dict[str,dict[str,MTSP_Solver.ChunkConnection]]={}
        
        # Add chunk connections
        for chunk_set,connection_type in zip(chunk_sets,connection_types):
            for n,A in enumerate(chunk_set):
                possible_connections[A.unique_id()]={}
                for m,B in enumerate(chunk_set):
                    if B==A:
                        continue
                    if connection_type is MTSP_Solver.AtomChunkConnection:
                        if n>=m:
                            continue
                        connection = self.AtomChunkConnection(A,B)
                        connection.calculate_distance()
                        if connection.ts_distance is not None:
                            #print(connection.ts_distance)
                            possible_connections[A.unique_id()][B.unique_id()]=connection
                    elif connection_type is MTSP_Solver.ChunkConnection:
                        if B.start_resnum-A.end_resnum!=1:
                            continue
                        connection = self.ChunkConnection(A,B)
                        connection.calculate_distance(self.model_path,tmp_out_folder_path,self.ordered_atom_lookup.disordered_waters,quick_wE=quick_wE,dry=dry_run)
                        possible_connections[A.unique_id()][B.unique_id()]=connection
                        #print(connection.ts_distance)
                    else: assert False, connection_type
                for connection in possible_connections[A.unique_id()].values():
                    print(connection.A.get_disordered_tag(),connection.A.altloc,connection.B.get_disordered_tag(),connection.B.altloc,connection.ts_distance)
        '''
        ###########

        possible_connections:list[MTSP_Solver.AtomChunkConnection]=[]
        for constraint in constraints_handler.constraints:
            assert constraint.kind is not None
            constraints_that_include_H = ["Angle","non"]
            atoms_for_constraint = self.ordered_atom_lookup.select_atoms_by(
                names=constraint.atom_names(),
                res_nums=constraint.residues(),
                exclude_H=constraint.kind not in constraints_that_include_H
            )
            # Generate each combination of n atoms
            combinations_iterator:list[list[Atom]] = itertools.combinations(atoms_for_constraint, constraint.num_atoms())
            for atoms in combinations_iterator:
                #if constraint.kind!="Nonbond":
                # Don't have two alt locs of same atom in group
                if len(set([(a.get_name(),OrderedAtomLookup.atom_res_seq_num(a)) for a in atoms]))!= len(atoms): 
                    continue
                distance = constraint.get_distance(atoms)
                if distance is not None:
                    # print(constraint.atom_id)
                    # print(constraint.atom_names(),constraint.residues())
                    
                    atom_chunks_selection:list[AtomChunk] = []
                    for a in atoms:
                        if a.get_name()[0] != "H":
                            atom_chunks_selection.append(atom_chunks[atom_id(a)])
                        else:
                            # Convert restraints on riding hydrogens to restraints on parent atoms
                            assert constraint.kind in constraints_that_include_H, (constraint.kind,[a.get_name() for a in atoms])
                            res_num, altloc = OrderedAtomLookup.atom_res_seq_num(a), a.get_altloc()
                            possible_parents: list[Atom] = self.ordered_atom_lookup.select_atoms_by(
                                res_nums=[res_num],
                                altlocs=[altloc],
                                exclude_H=True,
                            )
                            parent_name = UntangleFunctions.H_get_parent_fullname(
                                a.get_name(),
                                [p.get_fullname() for p in possible_parents]
                            ).strip()
                            parent_atom = self.ordered_atom_lookup.select_atoms_by(
                                res_nums=[res_num],
                                altlocs=[altloc],
                                names=[parent_name],
                            )
                            assert len(parent_atom) == 1, (parent_atom,parent_name,res_num,altloc)
                            parent_atom = parent_atom[0]
                            atom_chunks_selection.append(atom_chunks[atom_id(parent_atom)])

                    hydrogens:list[str] = [a.get_name() for a in atoms if a.get_name()[0]=="H"]
                    #atom_chunks_selection:list[AtomChunk]= [atom_chunks[atom_id(a)] for a in atoms]
                    for ach in atom_chunks_selection:
                        assert constraint in ach.constraints_holder.constraints
                    #print([ch.name for ch in atom_chunks_selection])

                    if constraint.kind == "Nonbond":
                        connection = self.AtomChunkNonBondConnection(atom_chunks_selection,distance,constraint.kind,hydrogens)
                    else:
                        connection = self.AtomChunkConnection(atom_chunks_selection,distance,constraint.kind,hydrogens)
                    possible_connections.append(connection)

                        #print([(ch.name,ch.resnum,distance) for ch in connection.atom_chunks])
        # for connection in possible_connections:
        #     print(connection.ts_distance)
        #     print("_".join(ch.unique_id() for ch in connection.atom_chunks))

            

        disordered_connections:dict[str,list[MTSP_Solver.AtomChunkConnection]] ={} # options for each alt connection
        for connection in possible_connections:
            connection_id = connection.get_disordered_connection_id()
            if connection_id not in disordered_connections:
                disordered_connections[connection_id]=[]
            assert connection not in disordered_connections[connection_id] 
            disordered_connections[connection_id].append(connection)

        #finest_depth_chunks=orderedResidues
        finest_depth_chunks=atom_chunks
        return finest_depth_chunks,disordered_connections
