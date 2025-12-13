import math
from typing import Type
from Bio.PDB.Atom import Atom,DisorderedAtom
import UntangleFunctions 
import numpy as np
from multiprocessing import Pool
import scipy.stats as st
from statistics import NormalDist
from LinearOptimizer.VariableID import *
from LinearOptimizer.OrderedAtomLookup import OrderedAtomLookup
from LinearOptimizer import mon_lib_read



#from PhenixEnvScripts.cross_conformation_nonbonds import get_cross_conf_nonbonds
from LinearOptimizer.get_cross_conf_nonbonds_wrapper import get_cross_conf_nonbonds


MON_LIB_NONBOND=False # must be false at present.
SUPPRESS_CB_ANGLE=False
NO_INDIV_WEIGHTS=False  # (Edit: Also, could be good for avoiding undue weight placed on genuine outliers) Idea being that when we do unrestrained refinement, all geometry is ignored, and all atoms are treated "equally" in fitting to xray data. So it makes little sense to weight connections of the same type differently.
CLIP_NEG_LJ=True # try False
IGNORE_HYDROGEN_CLASHES=True # Does not apply to TwoAtomPenalty
VDW_BUFFER=0 
CLASH_OVERLAP_THRESHOLD=0.6 #0.8+0.4 # 0.4 0.6


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
    def log_chi(dev,sigma,ideal,num_bound_e):
        return math.log(1+ConstraintsHandler.chi(dev,sigma,ideal,num_bound_e))
    @staticmethod
    def scaled_dev(dev,sigma,ideal,num_bound_e):
        return abs(dev)/ideal*10
    @staticmethod
    def e_density_scaled_dev(dev,sigma,ideal,num_bound_e):
        return abs(dev)/ideal*num_bound_e
    @staticmethod
    def dev_sqr(dev,sigma,ideal,num_bound_e):
        return  dev**2
    #@staticmethod
    # def chi_piecewise(dev,sigma,ideal,num_bound_e):
    #     return  


    class Constraint():
        def specific_weight_mod(self,atom_names):
            # An additional weight that is dependent on the specific atoms
            raise NotImplementedError("abstract method")
        def __init__(self,pdb_ids:list[str],outlier_ok:bool,ideal:float,weight:float,sigma:float):
            if type(pdb_ids[0])==OrderedTag:
                self.site_tags = [ot.disordered_tag() for ot in pdb_ids]
            elif type(pdb_ids[0])!=DisorderedTag:
                self.site_tags = ConstraintsHandler.Constraint.site_tags_from_pdb_ids(pdb_ids)
            self.ideal = ideal
            self.weight = weight
            if NO_INDIV_WEIGHTS:
                self.weight = 1
            self.weight*=self.specific_weight_mod(self.atom_names())
            self.sigma=sigma
            self.outlier_ok=outlier_ok
            #self.sigma=1
            self.num_bound_e= sum([site.num_bound_e() for site in self.site_tags])
            self.tension = 0
            self.max_site_tension = 0

        @staticmethod
        def site_tags_from_pdb_ids(pdb_ids):
            return [DisorderedTag(pdb.strip()[-4:], pdb[:4].strip()) for pdb in pdb_ids]
        @staticmethod
        def ORDERED_site_tags_from_pdb_ids(pdb_ids):
            return [OrderedTag(pdb.strip()[-4:],pdb[:4].strip(),pdb[4]) for pdb in pdb_ids]
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
        def get_cost(self,atoms:list[Atom],scoring_function)->tuple[float,float,float]:  # TODO this returns an ideal value, z_score, and a cost, but should make a class that holds this info and return that
            raise NotImplementedError("Abstract method")
        def __repr__(self):
            return f"({ConstraintsHandler.Constraint.kind(type(self))} : {self.site_tags})"
        def __eq__(self, other:'ConstraintsHandler.Constraint'):
            return (type(self),self.site_tags) == (type(other),other.site_tags)
        @staticmethod
        def kind(constraint_type:Type):
            return {
                ConstraintsHandler.BondConstraint:VariableKind.Bond.value,
                ConstraintsHandler.AngleConstraint:VariableKind.Angle.value,
                ConstraintsHandler.NonbondConstraint:VariableKind.Nonbond.value,
                ConstraintsHandler.ClashConstraint:VariableKind.Clash.value,
                ConstraintsHandler.TwoAtomPenalty:VariableKind.Penalty.value,
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
        def specific_weight_mod(self,atom_names):
            return 1
        def __init__(self,atom_ids,outlier_ok,ideal,weight,sigma):
            assert len (atom_ids)==2
            super().__init__(atom_ids,outlier_ok,ideal,weight,sigma)
        @staticmethod
        def separation(a:Atom,b:Atom):
            return np.sqrt(np.sum((a.get_coord()-b.get_coord())**2))
        def get_cost(self,atoms:list[Atom],scoring_function)->tuple[float,float,float]:
            ordered_atoms = self.get_ordered_atoms(atoms)
            if ordered_atoms is None:
                return None
            a,b = ordered_atoms      
            dev =  self.ideal-self.separation(a,b)     
            z_score=abs(dev/self.sigma) # XXX
            return self.ideal, z_score, scoring_function(dev,self.sigma,self.ideal,self.num_bound_e) * self.weight
            # badness = (self.ideal-self.separation(a,b))**2/self.ideal
            # return 0, badness
    

        
    class AngleConstraint(Constraint):
        def specific_weight_mod(self,atom_names):
            if SUPPRESS_CB_ANGLE and ("CB" in atom_names and ("C" in atom_names or "N" in atom_names)):
                return 0.2
            return 1
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
        def get_cost(self,atoms:list[Atom],scoring_function)->tuple[float,float]:
            ordered_atoms = self.get_ordered_atoms(atoms)
            if ordered_atoms is None:
                return None
            a,b,c = ordered_atoms
            dev = (self.ideal-self.angle(a,b,c))
            z_score=abs(dev/self.sigma) # XXX

            #if atoms_in_LO_variable_string("Angle_51.C_B|51.CA_A|51.CB_A",ordered_atoms):
            # if atoms_in_LO_variable_string("Angle_51.C_B|51.CA_B|51.CB_B",ordered_atoms):
            #     print(dev, self.ideal, z_score, scoring_function(dev,self.sigma,self.ideal,self.num_bound_e), self.weight)
            #     print(a.get_coord(),b.get_coord(),c.get_coord())
            #     assert False
            return self.ideal, z_score, scoring_function(dev,self.sigma,self.ideal,self.num_bound_e) * self.weight
        
            # badness = (self.ideal-self.angle(a,b,c))**2/self.ideal
            # return 0, badness
        
    class ClashConstraint(Constraint):
        def specific_weight_mod(self,atom_names):
            return 1
        def __init__(self,atom_ids,outlier_ok,symmetries,weight=1):
            self.altlocs_vdw_dict={}
            super().__init__(atom_ids,outlier_ok,None,weight,None)
            self.symmetries=symmetries
        def add_ordered(self,altlocs,vdw,is_symm):
            assert len(altlocs)==2
            assert tuple(altlocs) not in self.altlocs_vdw_dict
            
            key = tuple(altlocs)
            if key not in self.altlocs_vdw_dict:
                self.altlocs_vdw_dict[key]=[None,None] # (same-asu contacts, crystal-packing contacts) 
            sep_idx = 1 if is_symm else 0
            self.altlocs_vdw_dict[key][sep_idx]=vdw
        def get_cost(self,atoms:list[Atom],scoring_function)->tuple[float,float]:
            ordered_atoms = self.get_ordered_atoms(atoms)
            if ordered_atoms is None:
                return None
            a,b = ordered_atoms
            altlocs=(a.get_altloc(),b.get_altloc())  # NOTE order matters!!
            if (altlocs not in self.altlocs_vdw_dict):
                return -1, 0, 0
            
            r0,r0_symm = self.altlocs_vdw_dict[altlocs]
            r, r_symm = ConstraintsHandler.NonbondConstraint.symm_min_separation(a,b,self.symmetries)
            assert any([r is not None for r in (r0,r0_symm)])
            dev = 0
            if r0 is not None:
                 dev += max(0,r0-r)
            if r0_symm is not None:
                dev += max(0,r0_symm-r_symm)
            #TODO consider making it a sum of scoring_function(dev,...) + scoring_function(dev_sym,...)

            sigma=1
            multiplicity_correction = 0.5 if len(set(altlocs))==1 else 1
            return r0, 0, scoring_function(dev,sigma,r0,self.num_bound_e) * self.weight * multiplicity_correction # note a z-score doesnt make sense here.

        def badness(r, vdw_gap):
            #if r > (vdw_gap-CLASH_OVERLAP_THRESHOLD):
            if r > (vdw_gap-0.4):
                return 0
            return 3*((vdw_gap - r)/0.4)**2  # TODO how Holton calculates it?
        
    class TwoAtomPenalty(Constraint):
        def specific_weight_mod(self,atom_names):
            return 1
        def __init__(self,atom_ids,outlier_ok,weight=1):
            self.altlocs_clash_dict={}
            super().__init__(atom_ids,outlier_ok,None,weight,None)
        def add_ordered(self,altlocs,badness):
            assert len(altlocs)==2
            assert tuple(altlocs) not in self.altlocs_clash_dict
            self.altlocs_clash_dict[tuple(altlocs)]=badness
        def get_cost(self,atoms:list[Atom],scoring_function):
            ordered_atoms = self.get_ordered_atoms(atoms)
            if ordered_atoms is None:
                return None
            a,b = ordered_atoms
            altlocs=(a.get_altloc(),b.get_altloc())  # NOTE order matters!!
            if (altlocs not in self.altlocs_clash_dict):
                return -1, 0, 0 #  TODO return None once LinearOptimizer.Solver behaves fine with it
            return -1, 0, self.altlocs_clash_dict[altlocs]*self.weight


    class NonbondConstraint(Constraint):
        def specific_weight_mod(self,atom_names):
            return 1
        # TODO currently just looks at the smallest separation of all symmetries.
        def __init__(self,atom_ids,outlier_ok,symmetries,weight=1):
            self.altlocs_vdw_dict={}
            super().__init__(atom_ids,outlier_ok,None,weight,None)
            # if (DisorderedTag(17,"H") in self.site_tags) and (DisorderedTag(81,"O") in self.site_tags):
            self.symmetries=symmetries
        def add_ordered(self,altlocs,vdw,is_symm):
            assert len(altlocs)==2
            assert tuple(altlocs) not in self.altlocs_vdw_dict
            
            key = tuple(altlocs)
            if key not in self.altlocs_vdw_dict:
                self.altlocs_vdw_dict[key]=[None,None] # (same-asu contacts, crystal-packing contacts) 
            idx = 1 if is_symm else 0
            self.altlocs_vdw_dict[key][idx]=vdw
        @staticmethod
        def symm_min_separation(a:Atom,b:Atom,symmetries):
            coord_dict={"a":[],"b":[]}
            assert (np.array([2,2,2])==UntangleFunctions.get_sym_xfmed_point(np.array([2,2,2]),symmetries[0])).all()
            symmetries = symmetries[1:] 
            #TODO FIXME neighbouring unit cells
            for atom, key in zip([a,b],coord_dict):
                for symm in symmetries:
                    coord_dict[key].append(
                        UntangleFunctions.get_sym_xfmed_point(
                            atom.get_coord(),
                            symm)
                    )
            min_separation_symm = np.inf
            #for coord_a in coord_dict["a"]:
                # for coord_b in coord_dict["b"]:
                #     min_separation_symm = min(min_separation_symm,np.sqrt(np.sum((coord_a-coord_b)**2)))
            for coord_a in coord_dict["a"]:
                min_separation_symm = min(min_separation_symm,np.sqrt(np.sum((coord_a-b.get_coord())**2)))
            for coord_b in coord_dict["b"]:
                min_separation_symm = min(min_separation_symm,np.sqrt(np.sum((a.get_coord()-coord_b)**2)))
            assert min_separation_symm != np.inf, (coord_dict,symmetries)
            separation_nonsymm=np.sqrt(np.sum((a.get_coord()-b.get_coord())**2))
            # if min_separation < ConstraintsHandler.BondConstraint.separation(a,b):
            #     print(f"Symmetry nonbond found for {a,b}")

            return separation_nonsymm,min_separation_symm
        
        @staticmethod
        def lennard_jones(r,r0):
            # NOTE: E_min not used.

            #obs=$1;ideal=$2;sigma=1;energy=lj(obs,ideal)
            # function lj0(r,r0) {if(r==0)return 1e40;return 4*((r0*2^(-1./6)/r)^12-(r0*2^(-1./6)/r)^6)}\
            # function lj(r,r0) {return lj0(r,r0)-lj0(6,r0)}'
            #def lj0(r,r0):
            if r == 0:
                return 1e40
            return 4*((r0*2**(-1/6)/r)**12-(r0*2**(-1/6)/r)**6)
            #return lj0(r,r0)-lj0(6,r0)
        @staticmethod  
        def badness(r,r0): # TODO implement properly as in untangle_score.csh. 
            lj_energy = ConstraintsHandler.NonbondConstraint.lennard_jones(r,r0)
            if CLIP_NEG_LJ and r>r0:
                return 0
            if r > r0:
                assert lj_energy <= 0, (r, r0, lj_energy)
                # Logic: Expect atoms to be near bottom of potential if within a certain range. 
                if lj_energy>=-0.0163169:  # R=2.5 sigma TODO look for a good value based on structures generated by phenix geometry refine etc. 
                    #return None
                    return 0
            return lj_energy+1
            #return abs(max(neg_badness_limit,ConstraintsHandler.NonbondConstraint.lennard_jones(r,r0,E_min)))
        #@staticmethod
        # def lennard_jones(Rij,Rmin,EPSij):
        #     # ENERGY =  EPSij * ( (Rmin/Rij)**12 - 2 * (Rmin/Rij)**6 )  NB: Sign is flipped here since EPSij is negative... (ENERGY = -EPSij when Rij=Rmin) 
        #     if Rij == 0:
        #         return 1e40
        #     return -1 * EPSij * ( (Rmin/Rij)**12 - 2 * (Rmin/Rij)**6 )
        #@staticmethod
        # def badness(r,r0,E_min): # TODO implement properly as in untangle_score.csh. 
        #     # negative when separation is greater than ideal.
        #     assert E_min < 0
        #     lj_energy = ConstraintsHandler.NonbondConstraint.lennard_jones(r,r0,E_min)
        #     assert lj_energy <= 0
        #     # Logic: Expect atoms to be near bottom of potential if within a certain range. 
        #     if lj_energy>=0.0163169*E_min:  # R=2.5 sigma TODO look for a good value based on structures generated by phenix geometry refine etc. 
        #         return 0
        #     return lj_energy-E_min
        #     #return abs(max(neg_badness_limit,ConstraintsHandler.NonbondConstraint.lennard_jones(r,r0,E_min)))
        def get_cost(self,atoms:list[Atom],scoring_function)->tuple[float,float]:
            ordered_atoms = self.get_ordered_atoms(atoms)
            if ordered_atoms is None:
                return None
            a,b = ordered_atoms            
            altlocs=(a.get_altloc(),b.get_altloc()) # NOTE order matters!!
            if (altlocs not in self.altlocs_vdw_dict):
                return -1, 0, 0
            
            r0,r0_symm = self.altlocs_vdw_dict[altlocs]
            r, r_min_symm = ConstraintsHandler.NonbondConstraint.symm_min_separation(a,b,self.symmetries)
            assert any([r is not None for r in (r0,r0_symm)])
            energy = 0
            if r0 is not None:
                 energy += ConstraintsHandler.NonbondConstraint.badness(r,r0)
            if r0_symm is not None:
                energy+= ConstraintsHandler.NonbondConstraint.badness(r_min_symm,r0_symm)
            # if energy is None:
            #     return None

            # if is_atom(a,62,"CD","A") and is_atom(b,128,"O","A"):
            #     print(r,r0,energy)
            #     asdads
            # if is_atom(a,62,"NH1","B") and is_atom(b,102,"O","B"):
            #     print(r,r0,energy)
            #     asdads

            #return 0, energy * self.weight
            #return np.sqrt(energy), energy * self.weight
            sigma=1
            z_score = dev = np.sqrt(energy) # Not a distance!
            #return z_score, scoring_function(dev,sigma,r0,self.num_bound_e) * self.weight
            multiplicity_correction = 0.5 if len(set(altlocs))==1 else 1
            return r0, np.sqrt(energy), energy * self.weight *multiplicity_correction
        

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

    def add_two_atom_penalty(self,constraint:TwoAtomPenalty,residual,altlocs,badness): # TODO same function as add_nonbond_constraint
        self.add(constraint,residual)
        found=0
        first_site = constraint.site_tags[0]
        for held_constraint in self.atom_constraints[first_site]:
            if held_constraint == constraint:
                constraint.add_ordered(altlocs,badness)
                found+=1
        assert found == 1
    def add_nonbond_constraint(self,constraint:NonbondConstraint,residual,altlocs,vdw_sum:float,is_symm:bool): # or clash constraint.
        self.add(constraint,residual)
        found=0
        first_site = constraint.site_tags[0]
        for held_constraint in self.atom_constraints[first_site]:
            if held_constraint == constraint:
                constraint.add_ordered(altlocs,vdw_sum,is_symm)
                found+=1
        assert found == 1
    def scale_constraint_weight(self,pdb_ids:list[str],constraint_type:Type,weight_factor:float):
        site_tags = ConstraintsHandler.Constraint.site_tags_from_pdb_ids(pdb_ids)
        site = site_tags[0]
        found=0
        for constraint in self.atom_constraints[site]:
            if (type(constraint) == constraint_type) and (set(constraint.site_tags) == set(site_tags)):
                found+=1
                constraint.weight*=weight_factor
        assert found<=1
        if found < 1:
            print(f"Warning: could not scale non-existent constraint: {constraint_type} {site_tags}")
        
    # For constraints that have same ideal and weight regardless of order of from conformer labels. (maybe this should never be used due to conformation-dependent library) 
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

    def load_all_constraints(self,pdb_file,ordered_atom_lookup:OrderedAtomLookup,symmetries:list, calc_nonbonds:bool=True,water_water_nonbond:bool=None,constraints_to_skip=[],two_atom_penalty_tuples:list[tuple[tuple[str,str],tuple[int,int],float,tuple[str,str]]]=[],outliers_to_ignore_file=None,turn_off_cdl=False):
        # two_atom_penalty_tuples: list of tuples like ( (CA,O), (12,110), 3, (A,B) ) --> their names, their res nums, the badness, their altlocs 

        constraints_file = UntangleFunctions.geo_file_name(pdb_file,turn_off_cdl=turn_off_cdl) # NOTE we only read ideal and weights.
        
        print(f"Parsing constraints in {constraints_file}")

        if calc_nonbonds:
            assert symmetries is not None, "Cannot set symmetries to `None` if doing nonbonds"
            assert water_water_nonbond is not None, "Must specify whether water-water nonbonds should be considered"

        outliers_to_ignore:list[tuple[str,list[str]]]=[]
        def outlier_ok(kind:str,pdb_strings:list[str]):
            return False
        if outliers_to_ignore_file is not None:
            with open(outliers_to_ignore_file) as f:
                for line in f:
                    kind = line.split()[0]
                    atom_ids = line.split("|")[-1].strip().split()
                    outliers_to_ignore.append((kind,atom_ids))
            def outlier_ok(kind:str,pdb_strings:list[str]):
                assert False # need to make work with OrderedSiteTag for pdb_strings variable (for reference, see how clashes and nonbond constraints are added )
                atom_strings = []
                for p in pdb_strings:
                    name = p[0:4].strip()
                    resnum = p.strip().split()[-1]
                    atom_strings.append(f"{name}_{resnum}")
                return (kind,atom_strings) in outliers_to_ignore

        self.constraints: list[ConstraintsHandler.Constraint]=[]
        bonds_added:list[DisorderedTag]=[]
        AngleEnds_added=[]
        end_point_angle_scale_factor=1
        local_nb_scale_factor=1
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
                    name1 = pdb1[0:4].strip()
                    name2 = pdb2[0:4].strip()
                    resnum1 = int(pdb1.strip().split()[-1])
                    resnum2 = int(pdb2.strip().split()[-1])

                    sites =  ConstraintsHandler.Constraint.site_tags_from_pdb_ids((pdb1,pdb2))
                    bonds_added.append(frozenset(tuple(sites)))
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

                    name1 = pdb1[0:4].strip()
                    name2 = pdb2[0:4].strip()
                    name3 = pdb3[0:4].strip()
                    resnum1 = int(pdb1.strip().split()[-1])
                    resnum3 = int(pdb3.strip().split()[-1])
                    AngleEnds_added.append( ((name1,resnum1),(name3,resnum3)) )
                    # TODO try reducing any angles that involve a "tip/end-point/dead-end" atom like O
                    if name1 == "O" or name2 == "O" or name3=="O":
                        weight*=end_point_angle_scale_factor
                    self.add(ConstraintsHandler.AngleConstraint(pdb_ids,outlier_ok("ANGLE",pdb_ids),ideal,weight,sigma),residual)
                
       
        # Add nonbonds that are flagged as issues for current structure AND when waters are swapped

        print("WARNING: assuming residue numbers are all unique")
        print("WARNING: assuming elements all single character")
        NB_pairs_added = []





        phenix_vdw_distances_table:dict[tuple[OrderedTag,OrderedTag,bool],float]={}        
        vdw_mon_lib_energy_name_dict:dict[dict[str]] = {}
        lj_mon_lib_energy_name_dict:dict[dict[str]] = {}

        skip_nonbonds = not calc_nonbonds or all([constr in constraints_to_skip for constr in (ConstraintsHandler.NonbondConstraint,ConstraintsHandler.ClashConstraint)])
        if skip_nonbonds:
            pass
        elif MON_LIB_NONBOND:
            assert False
            for res in "HOH, ALA, ARG, ASN, ASP, CYS, GLU, GLN, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL".split(", "):
                vdw_mon_lib_energy_name_dict[res],lj_mon_lib_energy_name_dict[res]=mon_lib_read.get_mon_lib_names(res)
            
            vdw_radii = mon_lib_read.read_vdw_radii(with_H=not IGNORE_HYDROGEN_CLASHES)
            lj_params = mon_lib_read.read_lj_parameters()
        else:
            # TODO don't need to be storing info for every conformation. All that matters is atom "energy types" and whether it's a crystal-packing contact or same-ASU nonbond 
            for pdb1, pdb2, vdw_sum,is_symm in get_cross_conf_nonbonds(pdb_file):
                conformer_tags = ConstraintsHandler.Constraint.ORDERED_site_tags_from_pdb_ids((pdb1,pdb2)) 
                key = tuple(conformer_tags)+(is_symm,) # NOTE Order matters
                # if IGNORE_HYDROGEN_CLASHES and any([t.element()=="H" for t in conformer_tags]):
                #     continue
                if key not in phenix_vdw_distances_table:
                    phenix_vdw_distances_table[key]=vdw_sum
                else:
                    assert phenix_vdw_distances_table[key]==vdw_sum, f"vdw sums differ for {key}, {phenix_vdw_distances_table[key]} != {vdw_sum}"
        


        

        general_water_nonbond=True
        water_water_nonbond=water_water_nonbond
        protein_protein_nonbonds=True
        cross_conformation_clashes = True
        num_clashes_found=0
        if water_water_nonbond: 
            assert general_water_nonbond
        if not skip_nonbonds and (general_water_nonbond or protein_protein_nonbonds): 
            waters_outer_loop=water_water_nonbond
            # atoms = ordered_atom_lookup.select_atoms_by(protein=True, waters=waters_outer_loop,exclude_H=IGNORE_HYDROGEN_CLASHES)
            # other_atoms = ordered_atom_lookup.select_atoms_by(protein=protein_protein_nonbonds,waters=general_water_nonbond, exclude_H=IGNORE_HYDROGEN_CLASHES
            #                                                   )
                                                              #exclude_atom_names=["C","N","CA","CB"])
            
            max_nonbond_sep=4
            #TODO can speed up a lot by generating symmetric grid and assigning atoms to voxels and checking if in same or neighbouring voxels,
            # rather than calculating distances between every single atom and checking if within max_nonbond_sep...



            def phenix_vdw_from_pdb_ids(pdb1:str,pdb2:str):
                conformer_tags = ConstraintsHandler.Constraint.ORDERED_site_tags_from_pdb_ids((pdb1,pdb2)) 
                key = tuple(conformer_tags)
                if key not in phenix_vdw_distances_table:
                    return None
                return phenix_vdw_distances_table[key]
            def phenix_vdw(confA:OrderedTag,confB:OrderedTag,is_symm:bool):
                key = (confA,confB,is_symm)
                if key not in phenix_vdw_distances_table:
                    return None
                return phenix_vdw_distances_table[key]
            

            last_num_found_LJ=last_num_found_clashes=0
            # TODO  XXX XXX CRITICAL !!!!!!!!!! Loop through the phenix_vdw_distances_table instead.
            for i, (confA,confB,is_symm) in enumerate(phenix_vdw_distances_table): # iterate over conformers
                if i%10000==0:
                    if i !=0:
                        if ConstraintsHandler.NonbondConstraint not in constraints_to_skip:
                            print(f"LJ potentials found: {len(NB_pairs_added)-last_num_found_LJ}")
                        if ConstraintsHandler.ClashConstraint not in constraints_to_skip:
                            print(f"Clashes found: {num_clashes_found-last_num_found_clashes}")
                    last_num_found_LJ=len(NB_pairs_added)
                    last_num_found_clashes=num_clashes_found
                    #print(f"Flagging cross-conformation nonbonds for{' protein' if not waters_outer_loop else ''} conformer {i}/{len(phenix_vdw_distances_table)}")
                    print(f"Flagging cross-conformation nonbonds {i}/{len(phenix_vdw_distances_table)}")

                if IGNORE_HYDROGEN_CLASHES and (confA.element()=="H" or confB.element=="H"):
                    continue

                if (((confA.element()=="H") and (confB.element() == "H")) # Do not consider H-H clashes (expect to be more harmful than helpful due to H positions being poor.)
                    or frozenset((confA.disordered_tag(),confB.disordered_tag())) in bonds_added): # Nonbonded only!
                    continue
                if (confA.atom_name() == confB.atom_name()) and (confA.resnum()==confB.resnum()):
                    continue
                #protein=True, waters=waters_outer_loop,exclude_H=IGNORE_HYDROGEN_CLASHES
                atomA = confA.lookup_atom(ordered_atom_lookup)
                if not waters_outer_loop\
                and UntangleFunctions.res_is_water(atomA.get_parent()):
                    continue
                #protein=protein_protein_nonbonds,waters=general_water_nonbond, exclude_H=IGNORE_HYDROGEN_CLASHES
                atomB = confB.lookup_atom(ordered_atom_lookup)
                if not general_water_nonbond\
                and UntangleFunctions.res_is_water(atomB.get_parent()):
                    continue
                if not protein_protein_nonbonds\
                and not UntangleFunctions.res_is_water(atomB.get_parent()):
                    continue

                    
                sep_idx = 1 if is_symm else 0
                min_separation = ConstraintsHandler.NonbondConstraint.symm_min_separation(
                    atomA,
                    atomB,
                    symmetries)[sep_idx]
                
                
                #if atoms_in_LO_variable_string("Nonbond_30.C_B|31.CA_B",(atomA,atomB)):
                # if atoms_in_LO_variable_string("Nonbond_1.C_A|2.CA_A",(atomA,atomB)):
                #     print(min_separation)
                #     print(phenix_vdw(confA,confB))
                
                phenix_vdw_sum = phenix_vdw(confA,confB,is_symm)
                if phenix_vdw_sum is None:
                    assert False

                nb_weight=1
                if abs((confA.resnum()-confB.resnum()))<=2: # XXX quite arbitrary
                    nb_weight*=local_nb_scale_factor

                
        
                if min_separation > phenix_vdw_sum-VDW_BUFFER:
                    continue
                    
                conf_pair=(confA,confB)
                # VDW overlap (clashes)
                altlocs = (confA.altloc(),confB.altloc())
                if (cross_conformation_clashes and (ConstraintsHandler.ClashConstraint not in constraints_to_skip)):
                    vdw_gap = phenix_vdw_sum
                    
                    if vdw_gap - min_separation >= CLASH_OVERLAP_THRESHOLD:  
                        num_clashes_found+=1                            
                        self.add_nonbond_constraint(ConstraintsHandler.ClashConstraint(conf_pair,outlier_ok("CLASH",conf_pair),symmetries,weight=nb_weight),
                                                    residual=None,altlocs=altlocs,vdw_sum=vdw_gap,is_symm=is_symm)

                # Lennard-Jones (nonbonded)
                if ConstraintsHandler.NonbondConstraint not in constraints_to_skip:
                    #if ((((B_A_check in AngleEnds_added) or (B_A_check_flipped in AngleEnds_added)) and other_atom.name in ["C","N","CA","CB","O"])
                    if ((confA.element() == "H" or confB.element() == "H") # Do not consider any LJ involving H.
                    #or (B_A_check in AngleEnds_added) or (B_A_check_flipped in AngleEnds_added)): 
                    ):
                        continue    
                    R_min = phenix_vdw_sum
                    self.add_nonbond_constraint(ConstraintsHandler.NonbondConstraint(conf_pair,outlier_ok("NONBOND",conf_pair),symmetries,weight=nb_weight),
                                                residual=None,altlocs=altlocs,vdw_sum=R_min,is_symm=is_symm)
                    NB_pairs_added.append(conf_pair)
        print(f"Added {num_clashes_found} clashes")
     
        num_nonbonded_general = len(NB_pairs_added)
        print(f"Added {num_nonbonded_general} nonbond constraints")
        

        # If clash was present at end of last loop, prioritise finding a solution without a clash and/or punish solution that does not assign differing altlocs to the conformers.
        if ConstraintsHandler.TwoAtomPenalty not in constraints_to_skip:
            for name, res_num, badness, altloc in two_atom_penalty_tuples:
                assert len(name)==len(res_num)==len(altloc)==2
                for n, r, a in zip(name,res_num,altloc):
                    if r not in ordered_atom_lookup.better_dict or n not in ordered_atom_lookup.better_dict[r]: 
                        break
                else:
                    pdb_ids = [f"{n}     ARES     A      {r}" for (n,r) in zip(name,res_num)]
                    #self.scale_constraint_weight(pdb_ids,ConstraintsHandler.ClashConstraint,10*badness)
                    self.add_two_atom_penalty(ConstraintsHandler.TwoAtomPenalty(pdb_ids,outlier_ok("PENALTY",pdb_ids)),None,altloc,100*badness)


class DisorderedTag():
    def __init__(self,res_num,atom_name):
        self._resnum, self._name = int(res_num),str(atom_name)
    def is_entry(self,pdb_entry:UntangleFunctions.PDB_Atom_Entry)->bool:
        if not pdb_entry.valid:
            return False
        return DisorderedTag(pdb_entry.res_num,pdb_entry.atom_name)==self
    def is_riding_H_entry(self,H_pdb_entry:UntangleFunctions.PDB_Atom_Entry)->bool:
        if not H_pdb_entry.valid:
            return False
        if not H_pdb_entry.atom_name[0]=="H":
            return False
        if not int(H_pdb_entry.res_num)==self.resnum():
            return False
        name_length = len(self.atom_name())
        assert 1 <= name_length <= 4 
        if name_length<4:
            full_name = ' '+self.atom_name() + ' '*(3-name_length)
        else:
            full_name=self.atom_name()
        if UntangleFunctions.H_get_parent_fullname(H_pdb_entry.atom_name,[full_name],debug_none_return=False)==full_name:
            return True
        return False
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

    def ordered_tag(self,altloc:str):
        return OrderedTag(self.resnum(),self.atom_name(),altloc)

    def __eq__(self, other:'DisorderedTag'):
        return (self._resnum, self.atom_name()) == (other._resnum, other.atom_name())
    def __ne__(self, other):
        return not(self == other)
    @staticmethod 
    def from_atom(a:Atom):
        return DisorderedTag(OrderedAtomLookup.atom_res_seq_num(a),a.get_name())

class OrderedTag(DisorderedTag):
    def __init__(self,res_num,atom_name,altloc):
        super().__init__(res_num,atom_name)
        self._altloc=altloc
    def altloc(self):
        return self._altloc
    def lookup_atom(self,ordered_atom_lookup:OrderedAtomLookup)->Atom:
        return ordered_atom_lookup.better_dict[self.resnum()][self.atom_name()][self.altloc()]
    def disordered_tag(self):
        return DisorderedTag(self.resnum(),self.atom_name())
    def __repr__(self):
        return f"{self.resnum()}.{self.atom_name()}.{self.altloc()}"
    def __hash__(self):
        return hash((self._resnum, self.atom_name(),self.altloc()))

    def __eq__(self, other:'OrderedTag'):
        return (self._resnum, self.atom_name(),self.altloc()) == (other._resnum, other.atom_name(),other.altloc())
    def __ne__(self, other):
        return not(self == other)