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
from UntangleFunctions import H_get_parent_fullname, PDB_Atom_Entry
import json
import copy
import gc

class Swapper():
    class Swap():
        def __init__(self,res_num:int,atom_name:str,from_altloc:str,to_altloc:str,new_pos=None):
            self.res_num = int(res_num)
            self.atom_name = atom_name 
            self.from_altloc = from_altloc
            self.to_altloc = to_altloc
            self.new_pos=new_pos
        def atom_id(self):
            return (self.res_num,self.atom_name,self.from_altloc)
        def matches_id(self,res_num,atom_name:str,from_altloc:str):
            return self.atom_id() == (int(res_num),atom_name,from_altloc)
            #return (self.res_num,self.atom_name,self.from_altloc) == (int(res_num),atom_name,from_altloc)
    
    class SwapGroup:
        def __init__(self,badness):
            self.badness =badness
            self.swaps: list[Swapper.Swap] = [] 
        def add(self,res_seq_num,atom_name:str,from_altloc,to_altloc:str,new_pos=None):
            self.swaps.append(Swapper.Swap(int(res_seq_num),atom_name,from_altloc=from_altloc,to_altloc=to_altloc,new_pos=new_pos))
        def get_atom_names(self):
            return [swap.atom_name for swap in self.swaps]
        def get_atom_ids(self):
            return [swap.atom_id() for swap in self.swaps]
        # Why not just use a dict?
        def get_to_altloc(self,res_num,atom_name:str,from_altloc:str):
            for swap in self.swaps:
                if swap.matches_id(res_num,atom_name,from_altloc):
                    return swap.to_altloc
            return from_altloc # No swap!
        def get_new_coord(self,res_num,atom_name:str,from_altloc:str):
            for swap in self.swaps:
                if swap.matches_id(res_num,atom_name,from_altloc):
                    return swap.new_pos
            return None # No swap!
        def get_line(self,P:PDB_Atom_Entry,name_for_altloc_override=None):
            atom_name_for_altloc=P.atom_name if name_for_altloc_override is None else name_for_altloc_override  #XXX
            to_altloc = self.get_to_altloc(P.res_num,atom_name_for_altloc, P.altloc)
            new_coord = self.get_new_coord(P.res_num,P.atom_name, P.altloc)
            new_line = P.new_altloc(to_altloc)
            if new_coord is not None:
                new_line = P.new_coord(new_coord)
            return new_line

        def __repr__(self):
            return str(self.swaps)  
        
    def __init__(self):
        self.num_times_swapped:dict[str,int]={}
        self.sol_idx=None
        self.swap_groups_sorted:list[Swapper.SwapGroup] = []

    
    def flag_swapping(self,swap_group:SwapGroup):
        raise Exception("Unimplemented")
        # for swap in swap_group.swaps:
        string_of_swaps =str(swap_group.swaps) 
        print(f"Swapping {string_of_swaps}")
        if string_of_swaps not in self.num_times_swapped:
            self.num_times_swapped[string_of_swaps]=0
        self.num_times_swapped[string_of_swaps]+=1


        
    def clear_candidates(self):
        self.swap_groups_sorted = []

    #TODO this is an unnecessary conversion of the solutions dict to a new data type (Swapper.SwapGroup).
    def add_candidates(self,swaps_file_path):
        self.sol_idx=None # XXX
        swap_group_candidates :list[Swapper.SwapGroup] = []
        assert swaps_file_path[-5:]==".json"
        with open(swaps_file_path,'r') as f:
            #solutions =  json.load(f)["solutions"]
            solutions =  copy.deepcopy(json.load(f))["solutions"]  # May avoid memory leak? TODO test. https://stackoverflow.com/questions/49369778/python-memory-increases-despite-garbage-collection
            for solution_name, solution_dict in solutions.items():
                badness:float = solution_dict["badness"]
                moves: dict[str,dict[str,str]] = solution_dict["moves"]
                swap_group = Swapper.SwapGroup(badness)
                # LinearOptimizer.Solver: solution_dict[site_key][from_altloc] = to_altloc
                swap_group_candidates.append(swap_group)
                for site_key in moves:
                    res_num, atom_name = site_key.split()[-1].split('.')
                    #_atom_name = _atom_name.strip("\n") 
                    for from_altloc, to_altloc in moves[site_key].items():
                        new_pos=None
                        if len(to_altloc.split())!=1:
                            to_altloc, new_pos_str = to_altloc.split(maxsplit=1)
                            i1,i2=new_pos_str.index("[")+1,new_pos_str.index("]")
                            new_pos= np.fromstring(new_pos_str[i1:i2],dtype=float,sep=' ')
                            #new_pos= np.fromstring(new_pos_str[i1:i2],dtype=float,sep=' ')

                        swap_group.add(res_num, atom_name,from_altloc,to_altloc,new_pos)

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
            #sorted_groups = sorted(swap_group_candidates,key=self.get_priority,reverse=True)
            sorted_groups = swap_group_candidates
            #TODO reimplement
            '''print("Swap group | priority")
            for sg, priority in zip(sorted_groups,[self.get_priority(sg) for sg in sorted_groups]):
                print(f"{sg} | {priority}")
            '''
            self.swap_groups_sorted.extend(sorted_groups)
        except:
            print(swap_group_candidates)
            raise Exception()
        else:
            self.sol_idx=0
            print(f"Loaded {len(sorted_groups)} candidate solutions successfully ({len(self.swap_groups_sorted)} total)")
        gc.collect()
        
    def solutions_remaining(self):
        return len(self.swap_groups_sorted) - self.sol_idx
    def run(self,model_path) -> tuple[str,SwapGroup]:
        assert model_path[-4:]==".pdb", model_path
        model_handle = os.path.basename(model_path)[:-4]
        out_path = f"{os.path.abspath(os.getcwd())}/output/{model_handle}_lpSwapped.pdb"
        assert os.path.abspath(out_path)!=os.path.abspath(model_path)

        swap_group: Swapper.SwapGroup = self.swap_groups_sorted[self.sol_idx]
        assert isinstance(swap_group,Swapper.SwapGroup)
         #NOTE Do not delete, plan to reimplement 
        '''
        self.flag_swapping(swap_group)
        '''
        self.sol_idx+=1 

        with open(model_path,'r') as f:
            new_lines = ""
            last_resnum=None
            for line in f.readlines():
                # Don't change irrelevant lines
            
                P = PDB_Atom_Entry(line)
                if not P.valid:
                    new_lines+=line
                    continue
                if P.res_name == "HOH":
                    new_lines += swap_group.get_line(P)
                    continue
                if P.res_num!=last_resnum:
                    nonHatom_full_names=[]
                    last_resnum=P.res_num

                # Swap according to swaps list
                anchored_atom_name = P.atom_name 
                if P.atom_name[0]=="H":
                    anchored_atom_name = H_get_parent_fullname(P.atom_name_unstripped,nonHatom_full_names).strip()
                assert P.res_num.isnumeric()


                

                # Alter line if polarity is negative 1
                new_lines+= swap_group.get_line(P,name_for_altloc_override=anchored_atom_name)

                # Sorry for this
                if P.atom_name[0]!="H" and P.atom_name_unstripped not in nonHatom_full_names:
                    nonHatom_full_names.append(P.atom_name_unstripped)
            
            

        #print(new_lines)
        with open(out_path,'w') as f:
            f.writelines(new_lines)
        
        #self.MakeSwapWaterFileByLines(new_lines.split("\n"),out_path_water_swapped)

        return out_path,swap_group

    # def get_priority(self,swap_group:SwapGroup):
    #     #TODO reduce priority substantially if a group is obtained by combining multiple other groups (sharing no atoms between each other) that have less badness
    #     priority = -swap_group.badness
    #     atoms_to_swap=swap_group.get_atom_names()
    #     badness = swap_group.badness
    #     if len(atoms_to_swap)==0:
    #         return -9999
        
    #     pct = badness/100

    #     # Prioritise exploring big changes that don't score much worse
    #     #TODO Single or odd number of swaps should be flagged as interesting

    #     # TODO Reimplement
    #     ''' 
    #     interesting_swaps = [s for s in swap_group.swaps if s.atom_name() not in ["O","H"]]
    #     num_separated_swaps = 0
    #     last_swap_start_res_num=-1
    #     for swap in interesting_swaps:
    #         res_num = int(swap[0])
    #         if res_num > last_swap_start_res_num+1:
    #             last_swap_start_res_num = res_num
    #             num_separated_swaps+=1

    #     priority += 10*pct*max(0,min( num_separated_swaps-1, 3)) 

    #     num_very_interesting = swap_group.get_atom_names().count("CB")  #TODO test
    #     priority += 5*pct*min(num_very_interesting,2) 
    
    #     num_interesting = len(interesting_swaps)
    #     priority+=1*pct*min(num_interesting,5)

    #     string_of_swaps =str(swap_group.swaps) 
    #     if string_of_swaps  in self.num_times_swapped:           
    #         priority-=15*pct*self.num_times_swapped[string_of_swaps]
    #         #priority-=5*pct*self.num_times_swapped[swap]
    #     '''
        

    #     return priority


    # @staticmethod
    # def MakeSwapWaterFileByLines(lines,out_path):
    #     # swap waters
    #     water_swapped_lines = ""
    #     altlocs = []
    #     for line in lines:
    #         P = Swapper.PDB_Atom_Entry(line)
    #         if (not P.valid):
    #             continue
    #         if P.altloc not in altlocs:
    #             altlocs.append(P.altloc)

    #     assert len(altlocs)==2, (altlocs, out_path)
    #     for line in lines:
    #         P = Swapper.PDB_Atom_Entry(line)
    #         if (not P.valid) or (not P.res_name=="HOH"):
    #             water_swapped_lines+=line
    #         else:
    #             altloc_idx = altlocs.index(P.altloc)
    #             water_swapped_lines+=P.new_altloc(altlocs[(altloc_idx+1)%2])
    #     #print(water_swapped_lines)
    #     with open(out_path,'w') as f:
    #         f.writelines(water_swapped_lines)

    # @staticmethod
    # def MakeSwapWaterFile(in_pdb_file, out_path):
    #     # swap waters
    #     lines = []
    #     with open(in_pdb_file) as f:
    #         for line in f.readlines():
    #             lines.append(line)
    #     Swapper.MakeSwapWaterFileByLines(lines,out_path)
