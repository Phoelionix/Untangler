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

import numpy as np
import os
from UntangleFunctions import H_get_parent_fullname


class Swapper():
    class SwapGroup:
        def __init__(self,badness):
            self.badness =badness
            self.swaps: list[tuple[str]] = [] 
        def add(self,res_seq_num,atom_name):
            self.swaps.append((res_seq_num,atom_name))
        def get_atom_names(self):
            return [swap[1] for swap in self.swaps]
        def __repr__(self):
            return str(self.swaps)  
        
    def __init__(self):
        self.num_times_swapped:dict[str,int]={}
        self.sol_idx=None
        self.swap_groups_sorted:list[Swapper.SwapGroup] = []

    def get_priority(self,swap_group:SwapGroup):
        #TODO reduce priority substantially if a group is obtained by combining multiple other groups (sharing no atoms between each other) that have less badness
        priority = -swap_group.badness
        atoms_to_swap=swap_group.get_atom_names()
        badness = swap_group.badness
        if len(atoms_to_swap)==0:
            return -9999
        
        pct = badness/100

        # Prioritise exploring big changes that don't score much worse
        interesting_swaps = [s for s in swap_group.swaps if s[1][0] not in ["O","H"]]
        num_separated_swaps = 0
        last_swap_start_res_num=-1
        for swap in interesting_swaps:
            res_num = int(swap[0])
            if res_num > last_swap_start_res_num+1:
                last_swap_start_res_num = res_num
                num_separated_swaps+=1

        priority += 10*pct*max(0,min( num_separated_swaps-1, 3)) 

        num_very_interesting = swap_group.get_atom_names().count("CB")  #TODO test
        priority += 5*pct*min(num_very_interesting,2) 
    
        num_interesting = len(interesting_swaps)
        priority+=1*pct*min(num_interesting,5)

        string_of_swaps =str(swap_group.swaps) 
        if string_of_swaps  in self.num_times_swapped:           
            priority-=15*pct*self.num_times_swapped[string_of_swaps]
            #priority-=5*pct*self.num_times_swapped[swap]

        return priority
    
    def flag_swapping(self,swap_group:SwapGroup):
        # for swap in swap_group.swaps:
        string_of_swaps =str(swap_group.swaps) 
        print(f"Swapping {string_of_swaps}")
        if string_of_swaps not in self.num_times_swapped:
            self.num_times_swapped[string_of_swaps]=0
        self.num_times_swapped[string_of_swaps]+=1

    class PDB_Params:
        # Functions return atom entry with changes applied. Does not change in place 
        def __init__(self,atom_entry:str):
            self.valid=False
            self.atom_entry=atom_entry
            valid_record_types=["ATOM","HETATM"]
            if np.any([atom_entry.startswith(k) for k in valid_record_types]):
                self.valid=True
                self.res_name = atom_entry[17:20]
                self.atom_name_unstripped = atom_entry[12:16] 
                self.atom_name=self.atom_name_unstripped.strip()
                self.altloc=atom_entry[16]
                self.res_num =atom_entry[22:26].strip()
        def swap_altloc(self):
            new_altloc = dict(A="B",B="A")[self.altloc]   
            return self.atom_entry[:16]+new_altloc+self.atom_entry[17:]
        
    def clear_candidates(self):
        self.swap_groups_sorted = []
    def add_candidates(self,swaps_file_path):
        self.sol_idx=None
        swap_group_candidates :list[Swapper.SwapGroup] = []
        with open(swaps_file_path,'r') as f:
            for line in f.readlines():
                if line.startswith(" "):
                    continue
                if line.startswith("solution"):
                    badness = float(line.strip().split()[-1])
                    swap_group_candidates.append(Swapper.SwapGroup(badness))
                if line.startswith("Flagging"):
                    swap = line.split('---')[-1]
                    _res_num, _atom_name = swap.split('.')
                    _atom_name = _atom_name.strip("\n") 
                    swap_group_candidates[-1].add(_res_num, _atom_name)

        assert len(swap_group_candidates)>0, f"nothing ({swap_group_candidates}) parsed from {swaps_file_path}"


        # Remove groups that are compositions of other non-intersecting groups
        # NOTE these should be forbidden in the linear optimizer. So this should be removed.
        for composition_candidate in swap_group_candidates:
            S = set(composition_candidate.swaps)
            if len(S)==0: # don't exclude the no-swaps solution
                continue
            for swap_group in swap_group_candidates:
                A = set(swap_group.swaps)
                broken=False
                for other_group in swap_group_candidates:
                    B= set(other_group.swaps)
                    if len(A & B)!=0: 
                        # sets intersect
                        continue
                    if S==A or A == B:
                        continue
                    if (S == A | B) and (composition_candidate.badness > swap_group.badness) and (composition_candidate.badness > other_group.badness):
                        #swap_group_candidates.remove(composition_candidate)
                        composition_candidate.badness*=100
                        broken=True
                        break
                if broken:
                    break
        # Decrease priority of groups that are contain less bad groups
        for composition_candidate in swap_group_candidates:
            S = set(composition_candidate.swaps)
            for swap_group in swap_group_candidates:
                # if swap_group.badness > composition_candidate.badness: # TODO check if this is better
                #     continue
                A = set(swap_group.swaps)
                if A==S or len(A)==0:
                    continue
                if len(S & A)== len(A):
                    composition_candidate.badness*=10
                    print(composition_candidate.swaps)
                
        try:
            sorted_groups = sorted(swap_group_candidates,key=self.get_priority,reverse=True)
            print("Swap group | priority")
            for sg, priority in zip(sorted_groups,[self.get_priority(sg) for sg in sorted_groups]):
                print(f"{sg} | {priority}")
            self.swap_groups_sorted.extend(sorted_groups)
        except:
            print(swap_group_candidates)
            raise Exception()
        else:
            self.sol_idx=0
            print(f"Loaded {len(sorted_groups)} candidate solutions successfully ({len(self.swap_groups_sorted)} total)")
        
    def solutions_remaining(self):
        return len(self.swap_groups_sorted) - self.sol_idx
    def run(self,model_path) -> tuple[str,SwapGroup]:
        assert model_path[-4:]==".pdb", model_path
        model_handle = os.path.basename(model_path)[:-4]
        out_path = f"{os.path.abspath(os.getcwd())}/output/{model_handle}_lpSwapped.pdb"
        out_path_water_swapped = f"{os.path.abspath(os.getcwd())}/output/{model_handle}_lpSwapped_water_swapped.pdb"       

        swap_group: Swapper.SwapGroup = self.swap_groups_sorted[self.sol_idx]
        self.flag_swapping(swap_group)
        self.sol_idx+=1 

        with open(model_path,'r') as f:
            new_lines = ""
            polarity_backbone=dict(A=1,B=1) # Track polarity for each (current) conformation
            polarity_O=dict(A=1,B=1)
            polarity_side_chain=dict(A=1,B=1)
            H_polarity_dict=dict()
            last_resnum=None
            for line in f.readlines():
                # Don't change irrelevant lines
            
                P = Swapper.PDB_Params(line)
                if not P.valid:
                    new_lines+=line
                    continue
                if P.res_name == "HOH":
                    if (P.res_num,P.atom_name) in swap_group.swaps:
                        new_lines += P.swap_altloc()
                    else:
                        new_lines += line 
                    continue
                if P.res_num!=last_resnum:
                    H_polarity_dict=dict()
                    last_resnum=P.res_num
                polarity_dict = None

                # Swap according to swaps list
                if P.atom_name[0]=="H":
                    polarity_dict = H_polarity_dict[H_get_parent_fullname(P.atom_name_unstripped,H_polarity_dict.keys())]
                else:

                    if P.atom_name in ["CA","C","N"]:
                        polarity_dict = polarity_backbone
                    elif P.atom_name == "O":
                        polarity_dict = polarity_O
                    else:
                        polarity_dict = polarity_side_chain
                    assert P.res_num.isnumeric()
                    
                    if (P.res_num,P.atom_name) in swap_group.swaps:  
                        polarity_dict[P.altloc]*=-1

                # Alter line if polarity is negative 1
                new_line = line
                # if P.res_num=="37" and P.altloc == "A":
                #     print(polarity_dict[P.altloc],P.atom_name)
                if polarity_dict[P.altloc]==-1:
                    new_line = P.swap_altloc()
                elif polarity_dict[P.altloc]!=1:
                    assert False,polarity_dict[P.altloc]
                new_lines+=new_line

                if P.atom_name[0]!="H":
                    # Set branches to polarity
                    branches = dict(
                        CA=polarity_side_chain,
                        C=polarity_O,
                    )
                    if P.atom_name in branches:
                        assert polarity_dict is polarity_backbone
                        branches[P.atom_name][P.altloc]=polarity_dict[P.altloc]
                    if P.atom_name_unstripped not in H_polarity_dict: # Sorry for this
                        H_polarity_dict[P.atom_name_unstripped] = {"A":None,"B":None}
                    H_polarity_dict[P.atom_name_unstripped][P.altloc] =polarity_dict[P.altloc]
            
            

        #print(new_lines)
        with open(out_path,'w') as f:
            f.writelines(new_lines)
        
        #self.MakeSwapWaterFileByLines(new_lines.split("\n"),out_path_water_swapped)

        return out_path,swap_group
                    
    @staticmethod
    def MakeSwapWaterFileByLines(lines,out_path):
        # swap waters
        water_swapped_lines = ""
        for line in lines:
            P = Swapper.PDB_Params(line)
            if (not P.valid) or (not P.res_name=="HOH"):
                water_swapped_lines+=line+"\n"
            else:
                water_swapped_lines+=P.swap_altloc()+"\n"
        #print(water_swapped_lines)
        with open(out_path,'w') as f:
            f.writelines(water_swapped_lines)

    @staticmethod
    def MakeSwapWaterFile(in_pdb_file, out_path):
        # swap waters
        lines = []
        with open(in_pdb_file) as f:
            for line in f.readlines():
                lines.append(line)
        Swapper.MakeSwapWaterFileByLines(lines,out_path)
