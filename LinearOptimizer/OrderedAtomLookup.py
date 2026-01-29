from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom,DisorderedAtom
from Bio.PDB.Residue import Residue # Note we don't want the disorderedresidue here, it refers to different residue types for same res seq num.
import UntangleFunctions 
from LinearOptimizer.Tag import *
from LinearOptimizer.VariableID import *
import os
from typing import Union


class OrderedAtomLookup: #TODO pandas?
    def __init__(self,atoms:Union[list[DisorderedAtom],str],protein=True,waters=False,altloc_subset=None,allowed_resnums=None,allowed_resnames=None,excluded_resnames=None,alternate_atoms:list[DisorderedAtom]=[]): # TODO type hint for sequence?
        # atoms can be a path to a pdb file to get atoms from 

        if type(atoms)==str:
            assert os.path.isfile(atoms), f"Missing file: {atoms}" 
            assert atoms[-4:]==".pdb", f"Invalid file extension: {atoms}" 
            atoms=PDBParser().get_structure("struct",atoms).get_atoms()
            
        
        assert allowed_resnums is None, "Not working"
        
        self.disordered_waters=[]
        self.ordered_atoms:list[Atom] = []
        self.serial_num_to_disordered_num_dict={}
        self.residue_nums=[]
        self.res_names:dict[int,str]={} # key: res num
        # self.residue_sources:Residue=[]
        self.altlocs=[]
        self.protein_altlocs=[]
        self.better_dict:dict[str,dict[str,dict[str,Atom]]] = {}  # res_num, atom name, altloc
        self.alt_pos_options:dict[DisorderedTag,list[dict[str,Atom]]]={}
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
            if excluded_resnames is not None and res_name in excluded_resnames:
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
                #if allowed_resnums is not None and allowed_resnames is not None:
                if allowed_resnums is None and allowed_resnames is None and excluded_resnames is None:
                    assert len(self.residue_nums)== 0 or res_num-1 in [self.residue_nums[-1],last_skipped_res_num], (res_num,self.residue_nums)
                self.residue_nums.append(res_num)
                self.res_names[res_num]=disorderedAtom.get_parent().get_resname()
                #self.residue_sources.append(disorderedAtom.get_parent())
                if is_water:
                    self.water_residue_nums.append(res_num)
            else: 
                if is_water:
                    assert res_num in self.water_residue_nums, f"Waters cannot reuse residue sequence numbers used in protein (res {res_num})"

            for orderedAtom in disorderedAtom:
                altloc = orderedAtom.get_altloc()
                if (altloc_subset is not None) and (altloc not in altloc_subset):
                    continue
                
                assert type(orderedAtom) == Atom
                self.ordered_atoms.append(orderedAtom)
                self.serial_num_to_disordered_num_dict[orderedAtom.get_serial_number()]=disorderedAtom.get_serial_number()
                
                if altloc not in self.altlocs:
                    self.altlocs.append(altloc)
                    if not is_water:
                        assert altloc not in self.protein_altlocs
                        self.protein_altlocs.append(altloc)

                

                ## BETTER APPROACH                    
                if disorderedAtom.name not in self.better_dict[res_num]:
                    self.better_dict[res_num][disorderedAtom.name] = {}
                #self.better_dict[res_num][disorderedAtom.name].append(orderedAtom)
                self.better_dict[res_num][disorderedAtom.name][altloc] = orderedAtom
        if num_zero_occ_skip>0:
            print(f"Skipped {num_zero_occ_skip} zero occupancy ordered atoms")

        for disordered_other in alternate_atoms:
            disordered_tag = DisorderedTag.from_atom(disordered_other)
            if disordered_tag not in self.alt_pos_options:
                self.alt_pos_options[disordered_tag]=[]
            ordered_dict={ordered_other.get_altloc(): ordered_other for ordered_other in disordered_other}
            self.alt_pos_options[disordered_tag].append(ordered_dict)
                    
    
    def from_tag(self,ordered_tag:OrderedTag,debug=True)->Atom:
        if debug:
            assert  ordered_tag.resnum() in self.better_dict, ordered_tag 
            assert  ordered_tag.atom_name() in self.better_dict[ordered_tag.resnum()], ordered_tag 
            assert  ordered_tag.altloc() in self.better_dict[ordered_tag.resnum()][ordered_tag.atom_name()], ordered_tag 
        return self.better_dict[ordered_tag.resnum()][ordered_tag.atom_name()][ordered_tag.altloc()] # NOTE If error here during constraint loading, might be that waters weren't loaded in OrderedAtomLookup.

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
    # def get_residue_sources(self)->list[Residue]:
    #     return self.residue_sources
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
                        waters=True, protein=True, exclude_H=False,exclude_atom_names=[],
                        only_protein_altlocs=True)->list[Atom]:
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
            if atom.get_name() in exclude_atom_names:
                continue
            if only_protein_altlocs and (atom.altloc not in self.protein_altlocs): 
                 # Skip solvent atoms that aren't assigned a conformer label of the protein. 
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