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
# - Consider turning sigmas back off? Don't actually make much sense if following unrestrained refinement. Should divide by ideal instead.
# - If only some atoms are allowed, need to include all restraints that involve the atoms, even those that do not involve allowed atoms
# - # Swaps can create nonbond issues that are not recorded due to not being present in geo file?
 

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
from multiprocessing import Pool
import scipy.stats as st
from statistics import NormalDist
from LinearOptimizer.VariableID import *


NO_INDIV_WEIGHTS=False  # (Edit: Also, could be good for avoiding undue weight placed on genuine outliers) Idea being that when we do unrestrained refinement, all geometry is ignored, and all atoms are treated "equally" in fitting to xray data. So it makes little sense to weight connections of the same type differently.
APPLY_TENSION_MOD=False
TURN_OFF_MIN_SIGMAS=False # Since there are genuine outliers this could be harming things.
ALLOW_OUTLIERS_FROM_POTENTIAL_OVERRIDES=False # Will not do min sigmas with potential genuine outliers.
#weight_mod_for_allowed_outliers=1e-1
weight_mod_for_allowed_outliers=1
#high_tension_penalty=1.75  # Penalty for current connections (i.e. geometries? groups?) that correspond to a high tension (above min_tension_where_anything_goes). TODO should just be a penalty for not having any geometries involving the site change.
#high_tension_penalty=10  # Penalty for current connections (i.e. geometries? groups?) that correspond to a high tension (above min_tension_where_anything_goes). TODO should just be a penalty for not having any geometries involving the site change.
high_tension_penalty=1  # Penalty for current connections (i.e. geometries? groups?) that correspond to a high tension (above min_tension_where_anything_goes). TODO should just be a penalty for not having any geometries involving the site change.

#TODO could try reverting model back to preswap model.

class OrderedAtomLookup: #TODO pandas?
    def __init__(self,atoms:list[DisorderedAtom],protein=True,waters=False,altloc_subset=None,allowed_resnums=None,allowed_resnames=None): # TODO type hint for sequence?
        assert allowed_resnums is None, "Not working"
        
        self.disordered_waters=[]
        self.ordered_atoms:list[Atom] = []
        self.serial_num_to_disordered_num_dict={}
        self.residue_nums=[]
        self.residue_sources=[]
        self.altlocs=[]
        self.res_names:dict[int,str]={}
        self.better_dict:dict[str,dict[str,dict[str,Atom]]] = {}
        self.water_residue_nums:list[int]=[]

        self.waters_allowed=waters

        atoms = list(atoms)


        #TODO check whether zero occ affects wE/geom stuff. If not can skip zero occ safely. Otherwise need to deal with it.
        num_zero_occ_skip=0
        last_skipped_res_num=None
        skip_zero_occ = False
        for disorderedAtom in atoms:
            res_num=OrderedAtomLookup.atom_res_seq_num(disorderedAtom)
            res_name=OrderedAtomLookup.atom_res_name(disorderedAtom)
            if allowed_resnums is not None and res_num not in allowed_resnums:
                continue
            if allowed_resnames is not None and res_name not in allowed_resnames:
                continue
            
            if disorderedAtom.get_occupancy()==0:
                if skip_zero_occ:
                    num_zero_occ_skip+=1
                    last_skipped_res_num=OrderedAtomLookup.atom_res_seq_num(disorderedAtom)
                    continue
            is_water= UntangleFunctions.res_is_water(disorderedAtom.get_parent())
            if is_water:
                self.disordered_waters.append(disorderedAtom)
            if is_water and not waters:
                continue
            if not is_water and not protein:
                continue
            
            if type(disorderedAtom)!=DisorderedAtom:
                assert type(disorderedAtom)==Atom, type(disorderedAtom)
                tmp = DisorderedAtom(disorderedAtom.name)
                tmp.set_parent(disorderedAtom.get_parent())
                tmp.disordered_add(disorderedAtom)
                disorderedAtom=tmp
                #"Checking can get res num"
                try:
                    disorderedAtom.get_parent().get_id()[1]
                except:
                    raise Exception("Something went wrong when making ordered atom disordered")
            if res_num not in self.better_dict:
                self.better_dict[res_num] = {}
                if allowed_resnums is not None and allowed_resnames is not None:
                    assert len(self.residue_nums)== 0 or res_num-1 in [self.residue_nums[-1],last_skipped_res_num], (res_num,self.residue_nums)
                self.residue_nums.append(res_num)
                self.res_names[res_num]=disorderedAtom.get_parent().get_resname()
                self.residue_sources.append(disorderedAtom.get_parent())
                if is_water:
                    self.water_residue_nums.append(res_num)
            else: 
                if is_water:
                    assert res_num in self.water_residue_nums, "Waters cannot reuse residue sequence numbers used in protein"

            for orderedAtom in disorderedAtom:
                altloc = orderedAtom.get_altloc()
                if (altloc_subset is not None) and (altloc not in altloc_subset):
                    continue
                
                assert type(orderedAtom) == Atom
                self.ordered_atoms.append(orderedAtom)
                self.serial_num_to_disordered_num_dict[orderedAtom.get_serial_number()]=disorderedAtom.get_serial_number()
                
                if altloc not in self.altlocs:
                    self.altlocs.append(altloc)

                

                ## BETTER APPROACH                    
                if disorderedAtom.name not in self.better_dict[res_num]:
                    self.better_dict[res_num][disorderedAtom.name] = {}
                #self.better_dict[res_num][disorderedAtom.name].append(orderedAtom)
                self.better_dict[res_num][disorderedAtom.name][altloc] = orderedAtom
        if num_zero_occ_skip>0:
            print(f"Skipped {num_zero_occ_skip} zero occupancy ordered atoms")
    def output_as_pdb_file(self, reference_pdb_file,out_path):
        # Only outputs based on self.better_dict at present!!
        with open(reference_pdb_file) as R, open(out_path,'w') as O:
            start_strs_considered = ["ATOM","HETATM"]
            for line in R:
                if line == "TER\n":
                    continue
                for s in start_strs_considered:
                    if line.startswith(s):
                        break
                else: # Not an atom entry
                    O.write(line)
                    continue
                name = line[12:16].strip()
                altloc = line[16]
                resnum = int(line[22:26])
                if resnum in self.better_dict and name in self.better_dict[resnum] and altloc in self.better_dict[resnum][name]:
                    O.write(line)
            
            


    @staticmethod
    def atom_res_seq_num(atom:Atom)->int:
        try:
            return atom.get_parent().get_id()[1]
        except:
            print(atom)
            print(atom.get_parent())
            print(atom.get_parent().get_id())
            raise Exception("AAAA")
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
    

    def select_atoms_for_sites(self,sites:list['DisorderedTag'],altlocs=None,**kwargs):
        atom_semi_selection = []
        if altlocs is None:
            for site in sites:
                if site.resnum() not in self.better_dict:
                    print(f"Warning: res num {site.resnum()} skipped")
                    continue
                atom_semi_selection.extend(self.better_dict[site.resnum()][site.atom_name()].values())
        else:
            for altloc in altlocs:
                for site in sites:
                    atom_semi_selection.extend(self.better_dict[site.resnum()][site.atom_name()][altloc])
        return self.select_atoms_from(atom_semi_selection,**kwargs)

    def select_atoms_by(self,**kwargs)->list[Atom]:
        

        return self.select_atoms_from(self.ordered_atoms, **kwargs)

    def select_atoms_from(self,ordered_atoms:list[Atom], res_nums=None,serial_numbers=None,names=None,
                        disordered_serial_numbers=None,altlocs=None,
                        waters=True, protein=True, exclude_H=False)->list[Atom]:
        atom_selection = []
        

        allowed_values = {} 
        for key, v in zip(("altloc","num","sn","name","dsn"),(altlocs,res_nums,serial_numbers,names,disordered_serial_numbers)):
            if v is not None:
                allowed_values[key] = set(v)
        for atom in ordered_atoms:
            is_water = UntangleFunctions.res_is_water(atom.get_parent())
            if not waters and is_water:
                continue
            if not protein and not is_water:
                continue
            if exclude_H and atom.element == "H":
                continue
            
            for key, v in allowed_values.items():
                if key == "altloc" and atom.get_altloc() not in v:
                    break
                if key == "num" and OrderedAtomLookup.atom_res_seq_num(atom) not in v:
                    break
                if key == "sn" and atom.get_serial_number() not in v:
                    break
                if key == "name" and atom.get_name() not in v:
                    break
                if key == "dsn" and self.get_disordered_serial_number(atom) not in v:
                    break
            else:
                assert atom.get_altloc() in self.altlocs
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
    @staticmethod
    def scoring_function(dev,sigma,ideal,num_bound_e):  # TODO how to put on same scale as nonbond and clashes?
        assert False, "not set" 
    @staticmethod
    def prob_weighted_stat(dev,sigma,ideal,num_bound_e):
        z_score=abs(dev/sigma)
        stat_energy=z_score**2
        prob = NormalDist().cdf(z_score)*2-1 # probability not noise
        assert 0<=prob <= 1,(z_score,NormalDist().cdf(z_score))
        return prob*stat_energy
    @staticmethod
    def chi(dev,sigma,ideal,num_bound_e):
        return  dev**2/ideal
    @staticmethod
    def scaled_dev(dev,sigma,ideal,num_bound_e):
        return abs(dev)/ideal*10
    @staticmethod
    def e_density_scaled_dev(dev,sigma,ideal,num_bound_e):
        return abs(dev)/ideal*num_bound_e
    @staticmethod
    def dev_sqr(dev,sigma,ideal,num_bound_e):
        return  dev**2


    class Constraint():
        def __init__(self,atom_ids:list[str],outlier_ok:bool,ideal:float,weight:float,sigma:float):
            self.site_tags = [DisorderedTag(pdb.strip().split()[-1], pdb[:4].strip()) for pdb in atom_ids]
            self.ideal = ideal
            self.weight = weight
            if NO_INDIV_WEIGHTS:
                self.weight = 1
            self.sigma=sigma
            self.outlier_ok=outlier_ok
            #self.sigma=1
            self.num_bound_e= sum([site.num_bound_e() for site in self.site_tags])
            self.tension = 0
            self.max_site_tension = 0
        def get_tension(self):
            return self.tension
        def get_max_site_tension(self):
            return self.max_site_tension
        def set_tension(self,tension):
            self.tension = tension
        def set_max_site_tension(self,max_site_tension):
            self.max_site_tension = max_site_tension
        def num_atoms(self):
            return len(self.site_tags)
        def atom_names(self):
            return [site.atom_name() for site in self.site_tags]
        def residues(self):
            return [site.resnum() for site in self.site_tags]
        def get_distance(self,atoms:list[Atom],scoring_function)->tuple[float,float]:
            raise Exception("abstract method")
        def __repr__(self):
            return f"({ConstraintsHandler.Constraint.kind(type(self))} : {self.site_tags})"
        def __eq__(self, other:'ConstraintsHandler.Constraint'):
            return (type(self),self.site_tags) == (type(other),other.site_tags)
        @staticmethod
        def kind(constraint_type:Type):
            return {
                ConstraintsHandler.BondConstraint:"Bond",
                ConstraintsHandler.AngleConstraint:"Angle",
                ConstraintsHandler.NonbondConstraint:"Nonbond",
                ConstraintsHandler.ClashConstraint:"Clash",
                }[constraint_type]

        DEBUG=False
        def get_ordered_atoms(self,candidate_atoms:list[Atom])->list[Atom]:
            other_name_and_resnum = [(a.get_name(),OrderedAtomLookup.atom_res_seq_num(a)) for a in candidate_atoms] 
            #print(self.num_atoms(),other_name_and_resnum)
            #print(set(other_name_and_resnum))

            ordered_atoms = [] # corresponding to order of self.atom_id. Important for e.g. bond angle
            for name_num in zip(self.atom_names(),self.residues()):
                    for i, other_name_num in enumerate(other_name_and_resnum):
                        if name_num==other_name_num:
                            ordered_atoms.append(candidate_atoms[i])
            if len(ordered_atoms)!=self.num_atoms():
                return None
            if self.DEBUG:
                assert len(other_name_and_resnum)==len(set(other_name_and_resnum))==self.num_atoms()
                assert len(np.unique(ordered_atoms))==self.num_atoms()
            return ordered_atoms

        
    class BondConstraint(Constraint):
        def __init__(self,atom_ids,outlier_ok,ideal,weight,sigma):
            assert len (atom_ids)==2
            super().__init__(atom_ids,outlier_ok,ideal,weight,sigma)
        @staticmethod
        def separation(a:Atom,b:Atom):
            return np.sqrt(np.sum((a.get_coord()-b.get_coord())**2))
        def get_distance(self,atoms:list[Atom],scoring_function)->tuple[float,float]:
            ordered_atoms = self.get_ordered_atoms(atoms)
            if ordered_atoms is None:
                return None
            a,b = ordered_atoms      
            dev =  self.ideal-self.separation(a,b)     
            z_score=abs(dev/self.sigma) # XXX
            return z_score, scoring_function(dev,self.sigma,self.ideal,self.num_bound_e) * self.weight
            # badness = (self.ideal-self.separation(a,b))**2/self.ideal
            # return 0, badness


        
    class AngleConstraint(Constraint):
        def __init__(self,atom_ids,outlier_ok,ideal,weight,sigma):
            assert len (atom_ids)==3
            super().__init__(atom_ids,outlier_ok,ideal,weight,sigma)
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
        def get_distance(self,atoms:list[Atom],scoring_function)->tuple[float,float]:
            ordered_atoms = self.get_ordered_atoms(atoms)
            if ordered_atoms is None:
                return None
            a,b,c = ordered_atoms
            dev = (self.ideal-self.angle(a,b,c))
            z_score=abs(dev/self.sigma) # XXX
            return z_score, scoring_function(dev,self.sigma,self.ideal,self.num_bound_e) * self.weight
        
            # badness = (self.ideal-self.angle(a,b,c))**2/self.ideal
            # return 0, badness
        
    class ClashConstraint(Constraint):
        default_weight=1
        def __init__(self,atom_ids,outlier_ok,altlocs,badness,weight=None,sigma=None):
            self.altlocs = altlocs
            self.badness = badness
            if weight is None:
                weight = ConstraintsHandler.ClashConstraint.default_weight
            super().__init__(atom_ids,outlier_ok,None,weight,None)
        def get_distance(self,atoms:list[Atom],scoring_function)->tuple[float,float]:
            ordered_atoms = self.get_ordered_atoms(atoms)
            if ordered_atoms is None:
                return None
            a,b = ordered_atoms
            if a.get_altloc() == self.altlocs[0] and b.get_altloc() == self.altlocs[1]:
                return 0, self.badness*self.weight
            return None 


    class NonbondConstraint(Constraint):
        # TODO currently just looks at the smallest separation of all symmetries.
        default_weight=1
        def __init__(self,atom_ids,outlier_ok,ideal_separation,symmetries,weight=None):
            if weight is None:
                weight = ConstraintsHandler.NonbondConstraint.default_weight
            super().__init__(atom_ids,outlier_ok,ideal_separation,weight,None)
            # if (DisorderedTag(17,"H") in self.site_tags) and (DisorderedTag(81,"O") in self.site_tags):
            self.symmetries=symmetries
        @staticmethod
        def symm_min_separation(a:Atom,b:Atom,symmetries):
            coord_dict={"a":[],"b":[]}
            for atom, key in zip([a,b],coord_dict):
                for symm in symmetries:
                    coord_dict[key].append(
                        UntangleFunctions.get_sym_xfmed_point(
                            atom.get_coord(),
                            symm)
                    )
            min_separation = np.inf
            for coord_a in coord_dict["a"]:
                for coord_b in coord_dict["b"]:
                    min_separation = min(min_separation,np.sqrt(np.sum((coord_a-coord_b)**2)))
            assert min_separation != np.inf, (coord_dict,symmetries)
            # if min_separation < ConstraintsHandler.BondConstraint.separation(a,b):
            #     print(f"Symmetry nonbond found for {a,b}")

            return min_separation
        
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
        def get_distance(self,atoms:list[Atom],scoring_function)->tuple[float,float]:
            ordered_atoms = self.get_ordered_atoms(atoms)
            if ordered_atoms is None:
                return None
            a,b = ordered_atoms            
            r = ConstraintsHandler.NonbondConstraint.symm_min_separation(a,b,self.symmetries)
            r0 = self.ideal
            energy = ConstraintsHandler.NonbondConstraint.badness(r,r0)
            return 0, energy * self.weight
        

    def __init__(self,constraints:list[Constraint]=[]):
        self.constraints:list[ConstraintsHandler.Constraint]=constraints
        self.atom_constraints:dict[DisorderedTag,list[ConstraintsHandler.Constraint]]={}
        self.atom_residuals:dict[DisorderedTag,float] = {}

    # def atom_in_constraints(self,atom_name,res_num):
    #     for constraint in self.constraints:
    #         for other_name,other_res_num in zip(constraint.atom_names(), constraint.residues()):
    #             if atom_name== other_name and res_num==other_res_num:
    #                 return True
    #     return False

    def get_constraints(self,site:'DisorderedTag',no_constraints_ok=False):
        if site not in self.atom_constraints:
            if no_constraints_ok:
                return ConstraintsHandler([])

            raise Exception(f"No constraints found for {site}")
        return ConstraintsHandler(self.atom_constraints[site])

    def add(self,constraint:Constraint,residual):
        if constraint not in self.constraints:
            self.constraints.append(constraint)
        # For each disordered atom site, track constraints it must use
            for site in constraint.site_tags:
                #site_id = VariableID(site,VariableKind.Atom)
                if site not in self.atom_constraints:
                    self.atom_constraints[site]=[]
                    self.atom_residuals[site]=0
                self.atom_constraints[site].append(constraint)
        if residual is not None:
            for site in constraint.site_tags:
                #site_id = VariableID(site,VariableKind.Atom)
                self.atom_residuals[site]+=residual
            # if site == DisorderedTag(17,"N"):
            #     print(self.atom_constraints[site])

    def load_all_constraints(self,constraints_file,nonbond_scores_files,water_clashes,ordered_atom_lookup:OrderedAtomLookup,symmetries, water_water_nonbond:bool,constraints_to_skip=[],potential_overrides_file=None):
        print(f"Parsing constraints in {constraints_file}")


        outliers_to_ignore:list[tuple[str,list[str]]]=[]
        def outlier_ok(kind:str,pdb_strings:list[str]):
            return False
        if ALLOW_OUTLIERS_FROM_POTENTIAL_OVERRIDES and potential_overrides_file is not None:
            with open(potential_overrides_file) as f:
                for line in f:
                    kind = line.split()[0]
                    atom_ids = line.split("|")[-1].strip().split()
                    outliers_to_ignore.append((kind,atom_ids))
            def outlier_ok(kind:str,pdb_strings:list[str]):
                atom_strings = []
                for p in pdb_strings:
                    name = p[0:4].strip()
                    resnum = p.strip().split()[-1]
                    atom_strings.append(f"{name}_{resnum}")
                return (kind,atom_strings) in outliers_to_ignore

        self.constraints: list[ConstraintsHandler.Constraint]=[]
        with open(constraints_file,"r") as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                ideal= _= sigma=  weight = None
                if line.startswith("bond"):
                    if ConstraintsHandler.BondConstraint in constraints_to_skip:
                        continue
                    constraint = lines[i:i+4]
                    pdb1=constraint[0].strip().split("\"")[1]
                    pdb2=constraint[1].strip().split("\"")[1]
                    ideal,  _,  _, sigma,  weight, residual = [float(v) for v in constraint[3].strip().split()]
                    altloc = pdb1.strip()[3]
                    # if altloc!=ordered_atom_lookup.altlocs[0]: # just look at one altloc to get constraint. 
                    #     continue
                    #print(pdb1,"|","|",pdb2,"|",ideal,"|",weight)
                    pdb_ids = (pdb1,pdb2)
                    self.add(ConstraintsHandler.BondConstraint(pdb_ids,outlier_ok("BOND",pdb_ids),ideal,weight,sigma),residual)
                if line.startswith("angle"):
                    if ConstraintsHandler.AngleConstraint in constraints_to_skip:
                        continue
                    constraint = lines[i:i+5]
                    pdb1=constraint[0].strip().split("\"")[1]
                    pdb2=constraint[1].strip().split("\"")[1]
                    pdb3=constraint[2].strip().split("\"")[1]
                    ideal,  _,  _, sigma,  weight, residual = [float(v) for v in constraint[4].strip().split()]
                    altloc = pdb1.strip()[3]
                    # if altloc!=ordered_atom_lookup.altlocs[0]:
                    #     continue
                    #print(pdb1,"|","|",pdb2,"|",ideal,"|",weight)
                    pdb_ids = (pdb1,pdb2,pdb3)
                    self.add(ConstraintsHandler.AngleConstraint(pdb_ids,outlier_ok("ANGLE",pdb_ids),ideal,weight,sigma),residual)
                
       
        # Add nonbonds that are flagged as issues for current structure AND when waters are swapped

        print("WARNING: assuming residue numbers are all unique")
        print("WARNING: assuming elements all single character")
        NB_pdb_ids_added = []
        for file in (nonbond_scores_files):
            if ConstraintsHandler.NonbondConstraint in constraints_to_skip:
                break

            with open(file) as f:

                lines = f.readlines()
                for i, line in enumerate(lines):
                    if line.startswith("nonbonded"):
                        constraint = lines[i:i+4]
                        #     pdb1 = f"{name}     ARES     A      {res_num}"
                        pdb1=constraint[0].strip().split("\"")[1]
                        pdb2=constraint[1].strip().split("\"")[1]
                        name1 = pdb1[0:4].strip()
                        name2 = pdb2[0:4].strip()
                        resnum1 = int(pdb1.strip().split()[-1])
                        resnum2 = int(pdb2.strip().split()[-1])
                        resname1 = pdb1.strip().split()[1][1:]
                        resname2 = pdb2.strip().split()[1][1:]
                        both_waters = resname1 == resname2=="HOH"
                        if both_waters and not water_water_nonbond:
                            continue
                            
                        
                        # Don't put constraint between H and its bonded atom. I mean, it shouldn't be a possible constraint. But just in case. 
                        if resnum1 == resnum2:
                            if name1[0] == "H":
                                nonH_in_res = ordered_atom_lookup.select_atoms_by(res_nums=[resnum1],exclude_H=True)
                                if nonH_in_res == name2:
                                    continue


                        separation = constraint[3].strip().split()[0]
                        ideal = constraint[3].strip().split()[1]
                        if separation > ideal:
                            continue 
                        altloc1 = pdb1.strip()[3]
                        altloc2 = pdb2.strip()[3]
                        broke=False
                        for a,r,n in zip((altloc1,altloc2),(resnum1,resnum2),(name1,name2)):
                            if r not in ordered_atom_lookup.better_dict or n not in ordered_atom_lookup.better_dict[r]: 
                                break
                        if broke:
                            continue
                        
                        #  Use pdb id string that is blind to altloc. TODO this is silly and confusing.
                        pdb1 = f"{name1}     ARES     A      {resnum1}"
                        pdb2 = f"{name2}     ARES     A      {resnum2}"

                            #continue
                        ideal = float(ideal)
                        pdb_ids = (pdb1,pdb2)
                        pdb_ids_flipped = (pdb2,pdb1)
                        if pdb_ids not in NB_pdb_ids_added and pdb_ids_flipped not in NB_pdb_ids_added:
                            self.add(ConstraintsHandler.NonbondConstraint(pdb_ids,outlier_ok("NONBOND",pdb_ids),ideal,symmetries),None)  
                            NB_pdb_ids_added.append(pdb_ids)  
        num_nonbonded_from_geo = len(NB_pdb_ids_added)
        if water_water_nonbond:  # TODO Symmetry clashes
            waters = ordered_atom_lookup.select_atoms_by(protein=False)
            ideal_water_separation=2.200
            for atom in waters:
                for other_atom in waters:
                    if other_atom == atom:
                        continue
                    if ConstraintsHandler.NonbondConstraint.symm_min_separation(atom,other_atom,symmetries) < 8:
                        pdb1 = f"{atom.name}     ARES     A      {OrderedAtomLookup.atom_res_seq_num(atom)}"
                        pdb2 = f"{other_atom.name}     ARES     A      {OrderedAtomLookup.atom_res_seq_num(other_atom)}"
                        pdb_ids = (pdb1,pdb2)
                        pdb_ids_flipped = (pdb2,pdb1)
                        if pdb_ids not in NB_pdb_ids_added and pdb_ids_flipped not in NB_pdb_ids_added:
                            self.add(ConstraintsHandler.NonbondConstraint(pdb_ids,outlier_ok("NONBOND",pdb_ids),ideal_water_separation,symmetries),None)
                            NB_pdb_ids_added.append(pdb_ids)     
        num_nonbonded_extra = len(NB_pdb_ids_added)-num_nonbonded_from_geo

        print(f"Nonbonded from geo: {num_nonbonded_from_geo}, extra: {num_nonbonded_extra}")

        if ordered_atom_lookup.waters_allowed and ConstraintsHandler.ClashConstraint not in constraints_to_skip:
            for name, res_num, badness, altloc in water_clashes:
                for n, r, a in zip(name,res_num,altloc):
                    if r not in ordered_atom_lookup.better_dict or n not in ordered_atom_lookup.better_dict[r]: 
                        break
                else:
                    pdb_ids = [f"{n}     ARES     A      {r}" for (n,r) in zip(name,res_num)]
                    self.add(ConstraintsHandler.ClashConstraint(pdb_ids,outlier_ok("CLASH",pdb_ids),altloc,badness))




# Ugh, never do inheritance. TODO refactor to composition.
class AtomChunk(OrderedResidue):
    # Just an atom c:
    def __init__(self,site_num,altloc,resnum,referenceResidue,atom:Atom,is_water:bool,constraints_handler):
        self.name=atom.get_name()
        self.is_water = is_water
        self.element = atom.element
        self.coord = atom.get_coord()
        super().__init__(altloc,resnum,referenceResidue,[atom],site_num)
        self.generate_bonds(constraints_handler)
        self.depth_tag="ATOM"
    def generate_bonds(self,constraints_handler:ConstraintsHandler):
        self.constraints_holder:ConstraintsHandler = \
            constraints_handler.get_constraints(self.get_disordered_tag(),no_constraints_ok=self.is_water)
            
    def get_coord(self):
        return self.coord
    def get_disordered_tag(self):
        return DisorderedTag(self.get_resnum(),self.name)
    def unique_id(self):
        return f"{self.depth_tag}&{self.get_site_num()}.{self.get_altloc()}&{self.get_disordered_tag()}"

class DisorderedTag():
    def __init__(self,res_num,atom_name):
        self._resnum, self._name = int(res_num),str(atom_name)
    def resnum(self):
        return self._resnum
    def atom_name(self):
        return self._name 
    def element(self):
        return self.atom_name()[0]   # FIXME !!!!!!!!!!
    def num_bound_e(self):
        return UntangleFunctions.NUM_E[self.element()]
    def __repr__(self):
        return f"{self.resnum()}.{self.atom_name()}"
    # def __format__(self,format_spec):
    #     return str(self)
    # def __str__(self):
    #     return f"{self.resnum()}.{self.atom_name()}"
    
    def __hash__(self):
        return hash((self._resnum, self.atom_name()))

    def __eq__(self, other:'DisorderedTag'):
        return (self._resnum, self.atom_name()) == (other._resnum, other.atom_name())
    def __ne__(self, other):
        return not(self == other)


class MTSP_Solver:
    MODE= "NO_RESTRICTIONS" #"HIGH_TOL" #"NO_RESTRICTIONS" # PHENIX REFMAC

    if MODE=="NO_RESTRICTIONS":
        max_sigmas=min_sigmas_where_anything_goes=min_tension_where_anything_goes={}
    if MODE=="HIGH_TOL":
        max_sigmas={
            # ConstraintsHandler.BondConstraint:4,
            # ConstraintsHandler.AngleConstraint:4,
            #ConstraintsHandler.BondConstraint:6,
            #ConstraintsHandler.AngleConstraint:3,
            ConstraintsHandler.BondConstraint:3,
            ConstraintsHandler.AngleConstraint:3,
            # ConstraintsHandler.BondConstraint:2,
            # ConstraintsHandler.AngleConstraint:2,
            # ConstraintsHandler.BondConstraint:99,
            # ConstraintsHandler.AngleConstraint:99,
        }    
        min_sigmas_where_anything_goes={
            ConstraintsHandler.BondConstraint:2,
            ConstraintsHandler.AngleConstraint:2,
            # ConstraintsHandler.BondConstraint:4,
            # ConstraintsHandler.AngleConstraint:4,
            # ConstraintsHandler.BondConstraint:99,
            # ConstraintsHandler.AngleConstraint:99,
        } 
        min_tension_where_anything_goes={
            ConstraintsHandler.BondConstraint:5,
            ConstraintsHandler.AngleConstraint:5,
            ConstraintsHandler.NonbondConstraint:5,
            ConstraintsHandler.ClashConstraint:5,
            # ConstraintsHandler.BondConstraint:7,
            # ConstraintsHandler.AngleConstraint:7,
            # ConstraintsHandler.NonbondConstraint:7,
            # ConstraintsHandler.ClashConstraint:7,
        } 
    if MODE=="PHENIX":
        max_sigmas={
            # ConstraintsHandler.BondConstraint:4,
            # ConstraintsHandler.AngleConstraint:4,
            #ConstraintsHandler.BondConstraint:6,
            #ConstraintsHandler.AngleConstraint:3,
            ConstraintsHandler.BondConstraint:4,
            ConstraintsHandler.AngleConstraint:3,
            # ConstraintsHandler.BondConstraint:2,
            # ConstraintsHandler.AngleConstraint:2,
            # ConstraintsHandler.BondConstraint:99,
            # ConstraintsHandler.AngleConstraint:99,
        } 
        min_sigmas_where_anything_goes={
            ConstraintsHandler.BondConstraint:2.5,
            ConstraintsHandler.AngleConstraint:2.5,
            # ConstraintsHandler.BondConstraint:4,
            # ConstraintsHandler.AngleConstraint:4,
            # ConstraintsHandler.BondConstraint:99,
            # ConstraintsHandler.AngleConstraint:99,
        } 
        # min_tension_where_anything_goes={
        #     ConstraintsHandler.BondConstraint:4,
        #     ConstraintsHandler.AngleConstraint:5,
        #     ConstraintsHandler.NonbondConstraint:4,
        #     ConstraintsHandler.ClashConstraint:4,
        # } 
        min_tension_where_anything_goes={
            ConstraintsHandler.BondConstraint:6,
            ConstraintsHandler.AngleConstraint:6,
            ConstraintsHandler.NonbondConstraint:6,
            ConstraintsHandler.ClashConstraint:6,
            # ConstraintsHandler.BondConstraint:7,
            # ConstraintsHandler.AngleConstraint:7,
            # ConstraintsHandler.NonbondConstraint:7,
            # ConstraintsHandler.ClashConstraint:7,
        } 
    if MODE == "REFMAC":
        max_sigmas={
            # ConstraintsHandler.BondConstraint:4,
            # ConstraintsHandler.AngleConstraint:4,
            # ConstraintsHandler.BondConstraint:2.5,
            # ConstraintsHandler.AngleConstraint:2.5,
            ConstraintsHandler.BondConstraint:4,
            ConstraintsHandler.AngleConstraint:4,
            # ConstraintsHandler.BondConstraint:99,
            # ConstraintsHandler.AngleConstraint:99,
        } 
        min_sigmas_where_anything_goes={
            # ConstraintsHandler.BondConstraint:2.5,
            # ConstraintsHandler.AngleConstraint:2.5,
            ConstraintsHandler.BondConstraint:2,
             ConstraintsHandler.AngleConstraint:3,
            # ConstraintsHandler.BondConstraint:99,
            # ConstraintsHandler.AngleConstraint:99,
        } 
        # min_tension_where_anything_goes={
        #     ConstraintsHandler.BondConstraint:4,
        #     ConstraintsHandler.AngleConstraint:5,
        #     ConstraintsHandler.NonbondConstraint:4,
        #     ConstraintsHandler.ClashConstraint:4,
        # } 
        min_tension_where_anything_goes={
            ConstraintsHandler.BondConstraint:6,
            ConstraintsHandler.AngleConstraint:6,
            ConstraintsHandler.NonbondConstraint:6,
            ConstraintsHandler.ClashConstraint:6,
            # ConstraintsHandler.BondConstraint:5,
            # ConstraintsHandler.AngleConstraint:5,
            # ConstraintsHandler.NonbondConstraint:5,
            # ConstraintsHandler.ClashConstraint:5,
        } 
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
    #         _,self.ts_distance,_ = UntangleFunctions.assess_geometry_wE(connection_structure_save_path,tmp_out_folder_path,phenixgeometry_only=quick_wE) 

    class AtomChunkConnection(): # TODO really need to have a disordered connection class to have some of these properties (e.g. max site tension, outlier_ok)

        def __init__(self, atom_chunks:list[AtomChunk],ts_distance,connection_type,hydrogen_names,z_score,max_site_tension,outlier_ok):
            assert hydrogen_names==None or len(hydrogen_names)==len(atom_chunks)
            self.atom_chunks = atom_chunks
            self.from_altlocs=[a.get_altloc() for a in atom_chunks]
            self.res_nums=[a.get_resnum() for a in atom_chunks]
            self.connection_type=connection_type
            self.ts_distance=ts_distance  # NOTE As in the travelling salesman problem sense
            self.hydrogen_tag=""
            self.hydrogen_name_set=set([])
            self.z_score=z_score
            ###
            self.max_site_tension=max_site_tension
            self.outlier_ok=outlier_ok
            ###
            self.forbidden= (connection_type in MTSP_Solver.max_sigmas) and (self.z_score > MTSP_Solver.max_sigmas[self.connection_type]) 
            if hydrogen_names is not None:
                self.hydrogen_tag = "_"+''.join(hydrogen_names)
                self.hydrogen_name_set = set(hydrogen_names)
        def single_altloc(self):
            return len({ch.altloc for ch in self.atom_chunks})==1
        def get_disordered_connection_id(self):
            kind = ConstraintsHandler.Constraint.kind(self.connection_type)
            return f"{kind}{self.hydrogen_tag}_{'_'.join([str(a_chunk.get_disordered_tag()) for a_chunk in self.atom_chunks])}"


    def __init__(self,pdb_file_path:str,APPLY_TENSION_MOD:dict[DisorderedTag,float], symmetries, align_uncertainty=False,ignore_waters=False,altloc_subset=None,resnums=None,resnames=None): 
        # TODO when subset size > 2, employ fragmentation/partitioning.
         
        # Note if we ignore waters then we aren't considering nonbond clashes between macromolecule and water.
        original_structure = PDBParser().get_structure("struct",pdb_file_path)
        # if align_uncertainty:
        #     self.align_uncertainty(original_structure)
        
        # debug_quick=False
        # if debug_quick:
        #     resnums = range(64)


        # Disabling because phenix seems to be deleting remark 290 :C
        '''
        self.symmetries = UntangleFunctions.parse_symmetries_from_pdb(pdb_file_path)
        '''
        self.symmetries = symmetries

        self.ordered_atom_lookup = OrderedAtomLookup(original_structure.get_atoms(),
                                                     protein=True,waters=not ignore_waters,
                                                     altloc_subset=altloc_subset,
                                                     allowed_resnums=resnums,allowed_resnames=resnames)   
        self.model_path=MTSP_Solver.subset_model_path(pdb_file_path,altloc_subset)
        self.APPLY_TENSION_MOD=APPLY_TENSION_MOD
        
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

    @staticmethod 
    def subset_model_path(pdb_file_path,altloc_subset:list[str]):
        tag=""
        if altloc_subset is not None:
            if len(altloc_subset)==1:
                tag = f"_conformer-{''.join(altloc_subset)}"
            else:
                tag = f"_subset-{''.join(altloc_subset)}"
        assert pdb_file_path[-4:]==".pdb",pdb_file_path
        return os.path.join(UntangleFunctions.separated_conformer_pdb_dir(), os.path.basename(pdb_file_path)[:-4]+tag+".pdb")
    @staticmethod
    def prepare_geom_files(base_model_path,all_altloc_subsets,num_threads=10,water_swaps=True,allowed_resnums=None,allowed_resnames=None,waters=True):
        if all_altloc_subsets==[None]:
            print("Considering full set of altlocs")

        original_structure = PDBParser().get_structure("struct",base_model_path)
        subset_model_paths=[]
        for altloc_subset in all_altloc_subsets:
            ordered_atom_lookup = OrderedAtomLookup(original_structure.get_atoms(),
                                                     protein=True,waters=waters,
                                                     altloc_subset=altloc_subset,
                                                     allowed_resnums=allowed_resnums,allowed_resnames=allowed_resnames)
                    
            subset_model = MTSP_Solver.subset_model_path(base_model_path,altloc_subset)
            assert subset_model not in subset_model_paths
            subset_model_paths.append(subset_model)
            ordered_atom_lookup.output_as_pdb_file(reference_pdb_file=base_model_path,out_path=subset_model)


        global pooled_method # not sure if this is a good idea. Did this because it tries to pickle but fails if local. Try replacing with line: multiprocessing.set_start_method(‘fork’)
        def pooled_method(i):
            MTSP_Solver.prepare_geom_files_for_one_subset(subset_model_paths[i],water_swaps=water_swaps)

        with Pool(num_threads) as p:
            p.map(pooled_method,range(len(subset_model_paths)))
            

    '''
    def geo_model_paths(pdb_file_path,altloc_subset:list[str]):
        assert pdb_file_path[-4:]==".pdb",pdb_file_path
        assert altloc_subset is not None
        files=[]
        for altloc in altloc_subset:
            tag = f"_conformer-{altloc}"
            files.append(os.path.join(
                UntangleFunctions.separated_conformer_pdb_dir(), os.path.basename(pdb_file_path)[:-4]+tag+".pdb"
            ))
        return files
            

    @staticmethod
    def prepare_geom_files(base_model_path,all_altloc_subsets,num_threads=10,water_swaps=True):

        original_structure = PDBParser().get_structure("struct",base_model_path)

        # Create subset pdb files
        for altloc_subset in all_altloc_subsets:
            ordered_atom_lookup = OrderedAtomLookup(original_structure.get_atoms(),
                                                     protein=True,waters=True,
                                                     altloc_subset=altloc_subset)
                    
            subset_model = MTSP_Solver.subset_model_path(base_model_path,altloc_subset)
            ordered_atom_lookup.output_as_pdb_file(reference_pdb_file=base_model_path,out_path=subset_model)

        altlocs = []
        for altloc_subset in all_altloc_subsets:
            for altloc in altloc_subset:
                if altloc not in altlocs:
                    altlocs.append(altloc)

        # Create geo files for each conformation
        geo_model_paths=[]
        for altloc in altlocs:
            ordered_atom_lookup = OrderedAtomLookup(original_structure.get_atoms(),
                                                     protein=True,waters=True,
                                                     altloc_subset=[altloc])
                    
            conformation = MTSP_Solver.geo_model_paths(base_model_path,[altloc])
            assert len(conformation)==1 and type(conformation)==list
            conformation=conformation[0]
            assert conformation not in geo_model_paths
            geo_model_paths.append(conformation)
            ordered_atom_lookup.output_as_pdb_file(reference_pdb_file=base_model_path,out_path=conformation)
        

        global pooled_method # not sure if this is a good idea. Did this because it tries to pickle but fails if local. Try replacing with line: multiprocessing.set_start_method(‘fork’)
        def pooled_method(i):
            MTSP_Solver.prepare_geom_files_for_one_subset(subset_model_paths[i],water_swaps=water_swaps)

        with Pool(num_threads) as p:
            p.map(pooled_method,range(len(subset_model_paths)))
    '''


    @staticmethod
    def geo_log_out_folder():
        return UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+"StructureGeneration/HoltonOutputs/"
    
    @staticmethod
    def prepare_geom_files_for_one_subset(model_path,water_swaps=True):
        UntangleFunctions.assess_geometry_wE(model_path,MTSP_Solver.geo_log_out_folder()) 
        model_water_swapped_path=UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+f"StructureGeneration/HoltonOutputs/{MTSP_Solver.water_swapped_handle(model_path)}.pdb"
        if water_swaps:
            Swapper.MakeSwapWaterFile(model_path,model_water_swapped_path)
            UntangleFunctions.assess_geometry_wE(model_water_swapped_path,MTSP_Solver.geo_log_out_folder()) 

    @staticmethod 
    def water_swapped_handle(model_path):
        return UntangleFunctions.model_handle(model_path)+"_WaSw"
        
    def calculate_paths(self,scoring_function,quick_wE=False, dry_run=False,atoms_only=True,
                        clash_punish_thing=False,nonbonds=True,water_water_nonbond=None,
                        constraint_weights:dict[Type,float]=None,
                        )->tuple[list[Chunk],dict[str,list[AtomChunkConnection]]]: #disorderedResidues:list[Residue]
        print("Calculating geometric costs for all possible connections between chunks of atoms (pairs for bonds, triplets for angles, etc.)")
        
        if constraint_weights is None:
            constraint_weights = {
                ConstraintsHandler.BondConstraint: 1,
                ConstraintsHandler.AngleConstraint: 1,
                ConstraintsHandler.NonbondConstraint: 1,
                ConstraintsHandler.ClashConstraint: 1,
            }

        if water_water_nonbond is None:
            water_water_nonbond = self.ordered_atom_lookup.has_water()
        if not self.ordered_atom_lookup.waters_allowed and (water_water_nonbond):
            assert False

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
        model_handle = UntangleFunctions.model_handle(self.model_path)

        needToFixWaterAltlocsDebugging=False
        if needToFixWaterAltlocsDebugging:
            nonbonds = False; water_water_nonbond=False ###################### XXX TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!
            print("WARNING: nonbonds force-disabled for debugging")

        # Generate geo file
        nonbond_scores_files = []
        water_clashes=[]
        if nonbonds:

            # CLASH 26.6682 0.79 | A  62 AARG  HD2| S 128 AHOH  O 
            # Terrible code. XXX
            def get_water_clashes(handle,from_altloc_dict:bool)->list: # from_altloc_dict, with keys being the to_altlocs.
                clashes = []
                #nonbond_path = UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+f"StructureGeneration/HoltonOutputs/{handle}_scorednonbond.txt"
                clashes_path = UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+f"StructureGeneration/HoltonOutputs/{handle}_clashes.txt"
                with open(clashes_path) as clashes_file:
                    for line in clashes_file:
                        #pdb1 = f"{name}     ARES     A      {res_num}"
                        res_name = [entry.split()[-2][1:] for entry in line.split("|")[1:]] 
                        if not res_name.count("HOH")==1:
                            continue
                        res_num = [int(entry.split()[1]) for entry in line.split("|")[1:]] 
                        #  NOTE At the moment altlocs will all be the same. But hope to 
                        # one day read from file that does clashes between all altlocs!!!!
                        to_altloc = [entry.split()[-2][0] for entry in line.split("|")[1:]] 
                        name = [entry.strip().split()[-1] for entry in line.split("|")[1:]] 
                        from_altloc = [a for a in to_altloc]
                        water_idx = res_name.index("HOH") # Because we continue above if don't have exactly one HOH involved in clash.
                        from_altloc[water_idx] = from_altloc_dict[to_altloc[water_idx]]
                        badness = float(line.split()[1])
                        #new_line = f"CLASH   {badness} XXXXX XXXXX XXXXX X |  {name[0]}_{res_num[0]} {name[1]}_{res_num[1]}"
                        clashes.append((name, res_num, badness, from_altloc))
                return clashes

            

            nonbond_scores_path = f"{UntangleFunctions.UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/HoltonOutputs/{model_handle}.geo"

            nonbond_scores_files = [nonbond_scores_path]

            model_water_swapped_handle=self.water_swapped_handle(self.model_path)

            unflipped_water_dict = {}
            for altloc in self.ordered_atom_lookup.get_altlocs():
                unflipped_water_dict[altloc]=altloc
            water_clashes =  get_water_clashes(model_handle,unflipped_water_dict) 
            
            more_water_swaps=True
            if water_water_nonbond and more_water_swaps and len(self.ordered_atom_lookup.get_altlocs())==2:
                flipped_water_dict = {}
                altlocs = self.ordered_atom_lookup.get_altlocs()
                for i in range(2):
                    flipped_water_dict[altlocs[i]]=altlocs[-i-1]

                water_clashes+= get_water_clashes(model_water_swapped_handle,flipped_water_dict)

            
                #nonbond_water_flipped_scores_path = UntangleFunctions.UNTANGLER_WORKING_DIRECTORY+f"StructureGeneration/HoltonOutputs/{model_water_swapped_handle}_scorednonbond.txt"
                nonbond_water_flipped_scores_path = f"{UntangleFunctions.UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/HoltonOutputs/{model_water_swapped_handle}.geo"
                nonbond_scores_files.append(nonbond_water_flipped_scores_path)

        constraints_file = f"{UntangleFunctions.UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/HoltonOutputs/{model_handle}.geo" # NOTE we only read ideal and weights.
        #post_unrestrained_constraints_file = f"{UntangleFunctions.UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/HoltonOutputs/{model_handle}.geo" # NOTE we only read ideal and weights.
        #pre_unrestrained_constraints_file = f"{UntangleFunctions.UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/HoltonOutputs/{model_handle}_init.geo" # NOTE we only read ideal and weights.
        #constraints_file = create_delta_constraints_file.create_delta_constraints_file(pre_unrestrained_constraints_file,post_unrestrained_constraints_file) 
        
        constraints_handler=ConstraintsHandler()
        constraints_to_skip=[]
        constraints_to_skip = [kind for kind,value in constraint_weights.items() if value <= 0]

        potential_overrides = f"{UntangleFunctions.UNTANGLER_WORKING_DIRECTORY}/StructureGeneration/HoltonOutputs/{model_handle}_potential_overrides.txt"
        constraints_handler.load_all_constraints(constraints_file,nonbond_scores_files,water_clashes,self.ordered_atom_lookup,symmetries=self.symmetries,water_water_nonbond=water_water_nonbond,constraints_to_skip=constraints_to_skip,
                                                 potential_overrides_file=potential_overrides)
        print("Nonordered constraint properties loaded.")
        #for n,atom in enumerate(self.ordered_atom_lookup.select_atoms_by(names=["CA","C","N"])):

        print("Creating ordered atom chunks")
        ordered_atoms = self.ordered_atom_lookup.select_atoms_by()
        altloc_counts={}
        last_site=None
        for n,atom in enumerate(ordered_atoms):
            if atom.get_altloc() not in altloc_counts:
                altloc_counts[atom.get_altloc()]=0
            altloc_counts[atom.get_altloc()]+=1
            
            if n%1000==0:
                print(f"Creating chunk {n} / {len(ordered_atoms)} ")
            if atom.element=="H":
                continue
            is_water = UntangleFunctions.res_is_water(atom.get_parent())
            atom_chunks[atom_id(atom)]=AtomChunk(
                                         altloc_counts[atom.get_altloc()]+1,
                                         atom.get_altloc(),
                                         OrderedAtomLookup.atom_res_seq_num(atom),
                                         atom.get_parent(),
                                         atom,
                                         is_water,
                                         constraints_handler
                                         )
            
            
        del ordered_atoms
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

        #################################################################
        print("Computing connection costs")
        possible_connections:list[MTSP_Solver.AtomChunkConnection]=[]
        for c, constraint in enumerate(constraints_handler.constraints):
            if c%1000 == 0: 
                print(f"Calculating constraint {c} / {len(constraints_handler.constraints)} ({constraint}) ")
            constraints_that_include_H = [ConstraintsHandler.AngleConstraint,ConstraintsHandler.NonbondConstraint]
            # atoms_for_constraint = self.ordered_atom_lookup.select_atoms_by(
            #     names=constraint.atom_names(),
            #     res_nums=constraint.residues(),
            #     exclude_H=type(constraint) not in constraints_that_include_H
            # )
            
            # contains_H = False
            # for site_tag in constraint.site_tags:
            #     if site_tag.atom_name()[0]=="H":
            #         contains_H=True
            # if contains_H:
            #     continue
            atoms_for_constraint = self.ordered_atom_lookup.select_atoms_for_sites(
                constraint.site_tags,
                exclude_H=type(constraint) not in constraints_that_include_H
            )
            res_name_num_dict={}
            for a in atoms_for_constraint:
                res_name_num_dict[a] = (a.name,OrderedAtomLookup.atom_res_seq_num(a) ) 
            # Generate each combination of n atoms
            combinations_iterator:itertools.combinations[tuple[Atom]] = itertools.combinations(atoms_for_constraint, constraint.num_atoms())
            
            debug_print=False
            if debug_print:
                combinations_iterator = list(combinations_iterator)
                print(len(list(combinations_iterator)))
                old_num_connections=len(list(possible_connections))
            # print(set([res_name_num_dict[a] for a in atoms]))
            # print(atoms)
            # assert False
            # NOTE that the tension is contributed by ALL geometric restraints affecting the sites, not just the one represented by 'constraint'
            tension_mod=None
            if self.APPLY_TENSION_MOD is not None:
                constraint.set_tension(np.sum([self.APPLY_TENSION_MOD[tag] for tag in constraint.site_tags])/len(constraint.site_tags))
                constraint.set_max_site_tension(max([self.APPLY_TENSION_MOD[tag] for tag in constraint.site_tags]))
                if APPLY_TENSION_MOD:
                    #tension_mod = (constraint.get_tension()+1)**2
                    #tension_mod = (constraint.get_max_site_tension()+1)**2
                    tension_mod = (constraint.get_max_site_tension()+1)**2
                    assert tension_mod>=1


            for atoms in combinations_iterator:
                # Don't have two alt locs of same atom in group
                if len(set([res_name_num_dict[a] for a in atoms]))!= len(atoms):                     
                    continue
                output = constraint.get_distance(atoms,scoring_function)
                if output is not None:
                    z_score,distance = output
                    distance*=constraint_weights[type(constraint)]

                    # print(constraint.atom_id)
                    # print(constraint.atom_names(),constraint.residues())
                    
                    atom_chunks_selection:list[AtomChunk] = []
                    for a in atoms:
                        if a.get_name()[0] != "H":
                            atom_chunks_selection.append(atom_chunks[atom_id(a)])
                        else:
                            # Convert restraints on riding hydrogens to restraints on parent atoms
                            assert type(constraint) in constraints_that_include_H, (type(constraint),[a.get_name() for a in atoms])
                            res_num, altloc = OrderedAtomLookup.atom_res_seq_num(a), a.get_altloc()
                            
                            res_atoms = []
                            for name in self.ordered_atom_lookup.better_dict[res_num]:
                                res_atoms.extend(self.ordered_atom_lookup.better_dict[res_num][name].values())
                            possible_parents: list[Atom] = self.ordered_atom_lookup.select_atoms_from(
                                res_atoms,
                                altlocs=[altloc],
                                exclude_H=True,
                            )
                            parent_name = UntangleFunctions.H_get_parent_fullname(
                                a.get_name(),
                                [p.get_fullname() for p in possible_parents]
                            ).strip()
                            parent_atom = self.ordered_atom_lookup.better_dict[res_num][parent_name][altloc]

                            #assert len(parent_atom) == 1, (parent_atom,parent_name,res_num,altloc)
                            #parent_atom = parent_atom[0]
                            atom_chunks[atom_id(parent_atom)].constraints_holder.add(constraint,None) # NOTE XXX
                            atom_chunks_selection.append(atom_chunks[atom_id(parent_atom)])

                    num_hydrogens = len([a.get_name() for a in atoms if a.element=="H"])
                    hydrogens=None
                    if num_hydrogens > 0:
                        hydrogens:list[str] = [a.get_name() if a.element=="H" else "x" for a in atoms]
                    
                    for ach in atom_chunks_selection:
                        assert constraint in ach.constraints_holder.constraints, (constraint, ach.get_disordered_tag(), [c for c in ach.constraints_holder.constraints])

                    if tension_mod is not None:
                        distance*=tension_mod
                    if constraint.outlier_ok:
                        distance*=weight_mod_for_allowed_outliers

                    connection = self.AtomChunkConnection(atom_chunks_selection,distance,type(constraint),hydrogens,z_score,constraint.get_max_site_tension(),constraint.outlier_ok)  # XXX putting max site tension in here is bad
                    possible_connections.append(connection)
            if debug_print:
                print("added:",len(list(possible_connections))-old_num_connections)
                print("total:",len(list(possible_connections)))
                        #print([(ch.name,ch.resnum,distance) for ch in connection.atom_chunks])
        #################################################################

        

        disordered_connections:dict[str,list[MTSP_Solver.AtomChunkConnection]] ={} # options for each alt connection
        
        
        for connection in possible_connections:
            connection_id = connection.get_disordered_connection_id()
            if connection_id not in disordered_connections:
                disordered_connections[connection_id]=[]
            assert connection not in disordered_connections[connection_id] 
            for other in disordered_connections[connection_id]:
                assert (connection.atom_chunks!=other.atom_chunks) or (connection.hydrogen_tag!=other.hydrogen_tag) or (connection.connection_type!=other.connection_type), f"duplicate connections of kind {connection.connection_type}, involving ordered atoms {[ch.unique_id() for ch in connection.atom_chunks]}, with badness {connection.ts_distance} and {other.ts_distance}!"
            disordered_connections[connection_id].append(connection)
        

        # If the current connections have a tension above a certain value, let solver consider all alternatives.
        print(f"Re-enabling connections that are alternatives to current connections with tension > {MTSP_Solver.min_tension_where_anything_goes}")
        num_bad_current_disordered_connections={k:0 for k in MTSP_Solver.min_tension_where_anything_goes}
        num_connections_re_enabled={k:0 for k in MTSP_Solver.min_tension_where_anything_goes}
        for disordered_connection_id, ordered_connections in disordered_connections.items():

            conn_type = ordered_connections[0].connection_type
            if conn_type not in MTSP_Solver.min_tension_where_anything_goes:
                continue
            if ordered_connections[0].max_site_tension >= MTSP_Solver.min_tension_where_anything_goes[conn_type]:
                num_bad_current_disordered_connections[conn_type]+=1
            else: continue

            for conn in ordered_connections:
                if conn.forbidden:
                    num_connections_re_enabled[conn.connection_type]+=1
                conn.forbidden=False
                if conn.single_altloc():
                    if conn_type == ConstraintsHandler.AngleConstraint:
                        conn.ts_distance*=high_tension_penalty
        print(f"Number of high tension disordered connections detected: {num_bad_current_disordered_connections}") 
        print(f"Ordered connections re-enabled: {num_connections_re_enabled}")
        # If the current connections have a sigma above a certain value, let solver consider all alternatives.

        if not TURN_OFF_MIN_SIGMAS:
            print(f"Re-enabling connections that are alternatives to current connections with sigma > {MTSP_Solver.min_sigmas_where_anything_goes}")
            num_bad_current_disordered_connections={k:0 for k in MTSP_Solver.min_sigmas_where_anything_goes}
            num_connections_re_enabled={k:0 for k in MTSP_Solver.min_sigmas_where_anything_goes}
            for disordered_connection_id, ordered_connections in disordered_connections.items():
                if ordered_connections[0].outlier_ok:  # XXX represents disordered connection
                    continue
                
                if ordered_connections[0].connection_type not in MTSP_Solver.min_sigmas_where_anything_goes:
                    continue
                for conn in ordered_connections:
                    if conn.single_altloc() and conn.z_score >= MTSP_Solver.min_sigmas_where_anything_goes[conn.connection_type]:
                        num_bad_current_disordered_connections[conn.connection_type]+=1
                        break
                else: continue
                for conn in ordered_connections:
                    if conn.forbidden:
                        num_connections_re_enabled[conn.connection_type]+=1
                    conn.forbidden=False
            print(f"Number of high sigma disordered connections detected: {num_bad_current_disordered_connections}") 
            print(f"Ordered connections re-enabled: {num_connections_re_enabled}")






                


        #finest_depth_chunks=orderedResidues
        finest_depth_chunks=atom_chunks
        return finest_depth_chunks,disordered_connections
