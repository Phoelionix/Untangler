#===========================================================================
# Untangler: Free ensemble models from local minima with the wrong altlocs 
# Copyright (C)  2026 Spencer Passmore (spencerpassmore@swin.edu.au)

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

# Purpose: Construct "highway" ILP variables for geometries that involve atoms separated by more than 2 covalent bonds. In particular, nonbonds/clashes, cross-link bonds, and rings.

from LinearOptimizer.Tag import DisorderedTag,OrderedTag
import pulp as pl
import LinearOptimizer.Input
from LinearOptimizer.Input import LP_Input, get_altlocs_key
from LinearOptimizer.VariableID import VariableID
from LinearOptimizer.RestraintsHandler import RestraintsHandler


def angle_to_highways(disordered_angles:dict[str,tuple[LP_Input.Geomection,pl.LpVariable]]):
    assert len(disordered_angles)>0

    highway_dict: dict[OrderedTag,dict[OrderedTag, list[pl.LpVariable]]] = {}
    LH_dtag = RH_dtag = None
    for angle_geomection, var in disordered_angles.values():
        LH_altloc,RH_altloc = (angle_geomection.from_altlocs[0],
                       angle_geomection.from_altlocs[-1]) 
        # end_pos_option_indices = None
        # if angle_geomection.involves_position_changes():
        #     end_pos_option_indices=(angle_geomection.position_option_indices[0],
        #                         angle_geomection.position_option_indices[-1])
        # altlocs_key = get_altlocs_key(end_altlocs,end_pos_option_indices)
        # if altlocs_key not in highway_dict:
        #     highway_dict[altlocs_key] = []
        # highway_dict[altlocs_key].append(var)
        if LH_dtag is None:
            LH_dtag=DisorderedTag(angle_geomection.res_nums[0],angle_geomection.atom_names[0])
            RH_dtag=DisorderedTag(angle_geomection.res_nums[-1],angle_geomection.atom_names[-1])
        else:
            assert LH_dtag == DisorderedTag(angle_geomection.res_nums[0],angle_geomection.atom_names[0])
            assert RH_dtag == DisorderedTag(angle_geomection.res_nums[-1],angle_geomection.atom_names[-1])

        LH_tag = LH_dtag.ordered_tag(LH_altloc)
        RH_tag = RH_dtag.ordered_tag(RH_altloc)
        if LH_tag not in highway_dict:
            highway_dict[LH_tag]={}
        if RH_tag not in highway_dict[LH_tag]:
            highway_dict[LH_tag][RH_tag]=[]
        highway_dict[LH_tag][RH_tag].append(var)

    #return {altlocs_key: pl.lpSum(ILP_vars) for altlocs_key, ILP_vars in highway_dict.items()}
    #highway_dict: dict[str, list[pl.LpAffineExpression]]= {altlocs_key: pl.lpSum(ILP_vars) for altlocs_key, ILP_vars in highway_dict.items()}
    
    highway_var_lpConstr_dict:dict[OrderedTag,dict[OrderedTag,tuple[pl.LpVariable,pl.LpConstraint]]] = {LH_tag:{} for LH_tag in highway_dict}
    for LH_tag in highway_dict:
        for RH_tag,ILP_vars in highway_dict[LH_tag].items():
            tag=f"{LH_tag};{RH_tag}"
            var_active = pl.LpVariable(f"Highway_{tag}",  #TODO cat=pl.LpBinary
                                lowBound=0,upBound=1,cat=pl.LpBinary)
            
            var_active.setInitialValue(int(LH_tag.altloc()==RH_tag.altloc()))
            
            constr = (var_active==pl.lpSum(ILP_vars),
                    f"HighwayConstraint_{tag}")
            highway_var_lpConstr_dict[LH_tag][RH_tag]=(var_active,constr)

    return highway_var_lpConstr_dict



def highway_from_highways(LH_highways:dict[OrderedTag,dict[OrderedTag,tuple[pl.LpVariable,pl.LpConstraint]]],
                          RH_highways:dict[OrderedTag,dict[OrderedTag,tuple[pl.LpVariable,pl.LpConstraint]]],
                          already_constructed_highway_dict:dict[OrderedTag,dict[OrderedTag,DisorderedTag]]={}):
    composite_highway_joint_dict:dict[OrderedTag,dict[OrderedTag,DisorderedTag]] = {LH_tag:{} for LH_tag in LH_highways} # gives the joining (disordered) atom. The composite highways are then [LH_tag][ordered_center_tag] and [ordered_center_tag][RH_tag] for every ordered_center_tag of the disordered atom
    highway_var_lpConstr_dict:dict[OrderedTag,dict[OrderedTag,tuple[pl.LpVariable,pl.LpConstraint]]] = {LH_tag:{} for LH_tag in LH_highways}
    
    for LH_tag in LH_highways:
        for center_tag, RH_tags in RH_highways.items():
            for RH_tag in RH_tags:
                if RH_tag in highway_var_lpConstr_dict[LH_tag]:
                    continue 
                already_constructed_highway = get_highway(LH_tag,RH_tag,already_constructed_highway_dict)
                if already_constructed_highway is not None:
                    highway_var_lpConstr_dict[LH_tag][RH_tag]=already_constructed_highway
                    continue

                # Add the new composite highway variable
                tag=f"{LH_tag};{RH_tag}"
                var_active = pl.LpVariable(f"Highway_{tag}",  #TODO cat=pl.LpBinary
                                    lowBound=0,upBound=1,cat=pl.LpBinary)

                var_active.setInitialValue(int(LH_tag.altloc()==RH_tag.altloc()))
                highway_var_lpConstr_dict[LH_tag][RH_tag]=(var_active,[])

        for center_tag, (LH_highway_var,_) in LH_highways[LH_tag].items():
            if center_tag not in RH_highways:
                continue
            for RH_tag, (RH_highway_var,_) in RH_highways[center_tag].items():
                if RH_tag not in composite_highway_joint_dict[LH_tag]: 
                    composite_highway_joint_dict[LH_tag][RH_tag]= center_tag.disordered_tag()

                if get_highway(LH_tag,RH_tag,already_constructed_highway_dict) is not None:
                    continue

                ILP_vars = (LH_highway_var, RH_highway_var)
                #composite_highway_joint_dict[LH_tag][RH_tag]= ( (LH_tag,center_tag, RH_tag) )

        
                # Now add the constraint corresponding to the shared center.
                tag=f"{LH_tag};{center_tag};{RH_tag}" 
                var_active=highway_var_lpConstr_dict[LH_tag][RH_tag][0]
                constr = (var_active>=pl.lpSum(ILP_vars) - len(ILP_vars)+1,
                        f"HighwayConstraint_{tag}")
                highway_var_lpConstr_dict[LH_tag][RH_tag][1].append(constr)
    

    return highway_var_lpConstr_dict,composite_highway_joint_dict


# TODO Put in a class for holding highways
def get_highway(left_tag,right_tag,composite_highways)->tuple:
    if left_tag in composite_highways and right_tag in composite_highways[left_tag]:
        return composite_highways[left_tag][right_tag]
    elif  right_tag in composite_highways and  left_tag  in composite_highways[right_tag]:
        return composite_highways[right_tag][left_tag]
    return None



def add_CA_highway(left_num,right_num,composite_highways, composition_dictionary,CA_altlocs,left_altlocs=None,right_altlocs=None):
    # NB: Whole thing with specifying left and right altlocs is redundant and makes things much more complex. We are adding all. (Though the approach could be useful when connecting to sidechains. )
    same_args=(composite_highways, composition_dictionary,CA_altlocs)
    
    assert right_num > left_num+1, (right_num,left_num)
    if left_altlocs is None:
        left_altlocs = CA_altlocs[left_num]
    if right_altlocs is None:
        right_altlocs = CA_altlocs[right_num]
    for altloc_left in left_altlocs:
      for altloc_right in right_altlocs:
        left_tag=OrderedTag(left_num,"CA",altloc_left)
        right_tag=OrderedTag(right_num,"CA",altloc_right)
        if  get_highway(left_tag,right_tag,composite_highways) is not None:
            continue
        else:
            start_with_left_resnum= [tag.resnum() for tag in composite_highways[left_tag] if tag.atom_name()=="CA"]
            # Find the longest highway constructed between left num and right num.   left_num ---- longest_left -- longest_right ---- right_num
            try:
                longest_left= max(num for num in start_with_left_resnum if num < right_num)
            except Exception as e:
                print(start_with_left_resnum)
                print(left_tag)
                print([ tag for tag in composite_highways[left_tag] if tag.atom_name()=="CA"])
                raise(e)

            # and from right to left (but past longest left):
            end_with_right_resnum= []
            for longest_right in range(longest_left,right_num):
                if get_highway(OrderedTag(longest_right, "CA",CA_altlocs[longest_right][0]), right_tag,composite_highways) is not None:
                    break
            else:
                print(left_tag,right_tag,longest_left)
                print(composite_highways[OrderedTag(right_num-1,"CA",CA_altlocs[right_num-1][0])])
                assert False 

            #longest_right=  min(num for num in end_with_right_resnum if num > longest_left) 

            # XXX Disgusting!
            left_highways= {
                left_tag:{OrderedTag(longest_left,"CA",alt):get_highway(left_tag,  OrderedTag(longest_left,"CA",alt),composite_highways) for alt in CA_altlocs[longest_left]}
            }
            right_highways={
                OrderedTag(longest_right,"CA",alt): {right_tag: get_highway(OrderedTag(longest_right,"CA",alt),right_tag,composite_highways) } 
                for alt in CA_altlocs[longest_right] }
            ###
            if longest_right==longest_left:
                # Join! 
                CA_to_CA,CA_to_CA_composition=highway_from_highways(left_highways,right_highways)
                #print(f"Joined {left_num} - {right_num} at {longest_right}")
                
                for key in CA_to_CA:
                    composite_highways[key] |= CA_to_CA[key]
                    composition_dictionary[key] |= CA_to_CA_composition[key]
            else:
                # Try to join up
                if abs(left_num-longest_right) <= abs(right_num-longest_left):  # left_num ----- longest_right | longest_right ---- right_num
                    add_CA_highway(left_num,longest_right,*same_args,left_altlocs,CA_altlocs[longest_right])
                else:
                    add_CA_highway(longest_left,right_num,*same_args,CA_altlocs[longest_left],right_altlocs) # left_num ----- longest_left  | longest_left ---- right_num
                print(f"Joining {left_num} - {right_num}")
                add_CA_highway(left_num,right_num,*same_args,(altloc_left,),(altloc_right,))
                
        #return CA_to_CA
        #return CA_to_CA
        

        
        #add_CA_highway(left_num,right_num)
        #highway_from_highways()


def construct_backbone_highways(LD_geomections:list[LP_Input.Geomection],constr_var_dict:dict[str,dict[str,tuple[LP_Input.Geomection,pl.LpVariable]]]): # LD = long-distance.
    # bond_var_dict={}
    # angle_var_dict={}
    # for constr_var_subdict in constr_var_dict.values():
    #     bond_var_dict |= {k:v for k, v in constr_var_subdict.items() if v[0].connection_type is RestraintsHandler.BondRestraint}
    #     angle_var_dict |= {k:v for k, v in constr_var_subdict.items() if v[0].connection_type is RestraintsHandler.AngleRestraint}
    bond_var_dict = {k:v for k, v in constr_var_dict.items() if list(v.values())[0][0].connection_type is RestraintsHandler.BondRestraint}
    angle_var_dict = {k:v for k, v in constr_var_dict.items() if list(v.values())[0][0].connection_type is RestraintsHandler.AngleRestraint}

    max_resnum= max(max(list(v.values())[0][0].res_nums) for v in bond_var_dict.values())
    
    
    ### Construct highways so as to link all C alphas ###

    ###-Ca--C--N--Ca-
    ### |         |
    ### Cb        Cb
    N_to_CA_bonds:dict[int,dict[OrderedTag, dict[OrderedTag, tuple[pl.LpVariable,None]]]]={} # [resnum][LH_from_altloc][RH_from_altloc]
    composite_highways:dict[OrderedTag,dict[OrderedTag,tuple[pl.LpVariable,pl.LpConstraint]]]={} 
    composition_dictionary:dict[OrderedTag,dict[OrderedTag,DisorderedTag]]= {}
    
    protein_residues=[1,]  # TODO exclude ligands
    print("TODO exclude ligands")
    CA_altlocs:dict[int,list[str]]={1:[]}
    for i in range(1,max_resnum): # NOTE 1 iteration for each N-C bond. 
        # 1. Create the 3-atom highways. 

        CA_to_N_id = LP_Input.Geomection.construct_disordered_connection_id(
            RestraintsHandler.AngleRestraint,
            [DisorderedTag(i,"CA"),DisorderedTag(i,"C"),DisorderedTag(i+1,"N")]
            )
        N_to_CA_id = LP_Input.Geomection.construct_disordered_connection_id(
            RestraintsHandler.BondRestraint,
            [DisorderedTag(i+1,"N"),DisorderedTag(i+1,"CA")]
            )
        
        if not CA_to_N_id in angle_var_dict:
            assert i>3, (CA_to_N_id,list(angle_var_dict.keys()))
            # #TODO assert next residue is separate ligand OR not N_to_CA_id in bond_var_dict, (CA_to_N_id,N_to_CA_id,list(angle_var_dict.keys()))
            #assert not N_to_CA_id in bond_var_dict, 
            continue 
        if not N_to_CA_id in bond_var_dict:
            assert i>2, (N_to_CA_id,list(bond_var_dict.keys()))
            continue
        if i not in protein_residues:
            protein_residues.append(i)
        if i+1 not in protein_residues:
            protein_residues.append(i+1)
        
        N_to_CA_bonds[i]={}
        CA_altlocs[i+1]=[]
        for bond_geomection, var in bond_var_dict[N_to_CA_id].values():
            LH_tag=bond_geomection.atom_chunks[0].get_ordered_tag()
            RH_tag=bond_geomection.atom_chunks[-1].get_ordered_tag()
            CA_altlocs[i+1].append(RH_tag.altloc())
            if i==1: #XXX
                CA_altlocs[i].append(LH_tag.altloc())
                assert LH_tag.altloc() in [bond_var[0].from_altlocs[0] for bond_var in bond_var_dict[LP_Input.Geomection.construct_disordered_connection_id(
                RestraintsHandler.BondRestraint,
                [DisorderedTag(i,"CA"),DisorderedTag(i,"C")]
                )].values()], f"{LH_tag.altloc()} does not appear to be present for first CA"

            if LH_tag not in N_to_CA_bonds[i]:
                N_to_CA_bonds[i][LH_tag]={}
            
            N_to_CA_bonds[i][LH_tag][RH_tag]=(var,None)
            
        # CB_to_C=CB_to_C_highways[i]  = angle_to_highways(angle_var_dict[CB_to_C_id])
        # C_to_N=C_to_N_bonds[i]
        # N_to_CB=N_to_CB_highways[i+1]= angle_to_highways(angle_var_dict[N_to_CB_id])      
        
        # # 2. Create 4-atom highways from Cb-Ca-C highway, C-N bond (equiv. to 2-atom highway).
        # CB_to_N,CB_to_N_composition=highway_from_highways(CB_to_C,C_to_N)
        # # 3. Create 6-atom highways from Cb-Ca-C-N highway, N-Ca-Cb highways
        # CB_to_CB,CB_to_CB_composition=highway_from_highways(CB_to_N,N_to_CB)

        CA_to_N=angle_to_highways(angle_var_dict[CA_to_N_id])
        N_to_CA=N_to_CA_bonds[i]
        
        # 2. Create 4-atom highways from Ca-C-N highway, N-Ca bond (equiv. to 2-atom highway).
        CA_to_CA,CA_to_CA_composition=highway_from_highways(CA_to_N,N_to_CA)
        
        composition_dictionary |= CA_to_CA_composition

        composite_highways = composite_highways | CA_to_N| CA_to_CA
        

    # 4 For each residue, determine which residues they are connected to by a LD_geomection. 
    linked_residues:dict[int,list[int]]={}
    for geomection in LD_geomections:
        resnums = [ch.resnum for ch in geomection.atom_chunks]
        if not all(resnum in protein_residues for resnum in resnums):
            continue
        
        if min(resnums) not in linked_residues:
            linked_residues[min(resnums)]=[]
        if max(resnums) not in linked_residues[min(resnums)]: 
            linked_residues[min(resnums)].append(max(resnums))
    print(linked_residues)

        

    # 5. Construct highways composed from CB_to_CB highways. 
    



    # "It'll do for now" binary search method
    # Create binary search tree
    # Construct highways for each connection as needed
    def make_tree(parent_node, segment,tree=None):
        if tree is None:
            tree = {}
        if parent_node not in tree:
            tree[parent_node]=[]
        if len(segment)==1:
            tree[parent_node].append(segment[0])
            return tree
        split=int(len(segment)/2)
        split_point = segment[split]
        child_segments = segment[:split], segment[split+1:]
        tree[parent_node].append(split_point)
        assert len(tree[parent_node])<=2, tree
        for child_segment in child_segments:
            if len(child_segment)>0:
                tree = make_tree(split_point,child_segment,tree)
    
        return tree
    
    tree={}
    segment=list(range(1, max_resnum+1))
    try:
        tree = make_tree(None,segment,tree) 
    except Exception as e:
        print(tree)
        raise e

    
    def find_smallest_subtree(val1,val2,x=None):
        if x is None:
            assert len(tree[None])==1
            x=tree[None][0]
        if x <val1 and x < val2:
            return find_smallest_subtree(val1,val2,max(tree[x]))
        elif x > val1 and x> val2:
            return find_smallest_subtree(val1,val2,min(tree[x]))
        else:
            return make_subtree(x,val1),make_subtree(x,val2)
    def make_subtree(x,val,tree_path=None):
        if tree_path is None:
            tree_path = []
        tree_path.append(x)
        if x == val:
            return tree_path
        if x < val:
            return make_subtree(max(tree[x]),val,tree_path) 
        if x > val:
            return make_subtree(min(tree[x]),val,tree_path) 





    for resnum, linked_resnums in linked_residues.items():
        print(f"Adding CA-CA highways to {resnum}")
        for linked_resnum in linked_resnums: 
            if abs(linked_resnum-resnum)>1:
                tree_paths = find_smallest_subtree(linked_resnum,resnum)
                for tree_path in tree_paths:
                    for i in reversed(range(len(tree_path)-1)):
                        num_prior=tree_path[i]; num_next = tree_path[i+1]
                        if not all(r in  CA_altlocs for r in [num_prior,num_next]):
                            continue
                        if num_next>(num_prior+1):
                            add_CA_highway(num_prior,num_next,composite_highways,composition_dictionary,CA_altlocs)
                        elif num_prior>(num_next+1):
                            add_CA_highway(num_next,num_prior,composite_highways,composition_dictionary,CA_altlocs)
    #print("Non-adjacent CA highways:")
    min_gap_print=2
    print(f"CA highways longer than {min_gap_print+1}:")
    print("-------------------------")
    total_adjacent_disordered=0
    total_non_adjacent_disordered=0
    for resnum in CA_altlocs:
        key=OrderedTag(resnum,"CA",CA_altlocs[resnum][0])
        if key not in composition_dictionary:
            continue
        for other_key in composition_dictionary[key].keys():
            if other_key.altloc()==CA_altlocs[other_key.resnum()][0]:
                if abs(other_key.resnum()-key.resnum())>1:
                    total_non_adjacent_disordered+=1
                else:
                    total_adjacent_disordered+=1
                    
                if abs(other_key.resnum()-key.resnum())>min_gap_print:
                    print(key.disordered_tag(),"<-->",other_key.disordered_tag())
    print(f"Total non-adjacent (single altloc): {total_non_adjacent_disordered}, (adjacent: {total_adjacent_disordered})")
    print("-------------------------")



    # Recursive-Tree-Search(x, key)
    #     if x = NIL or key = x.key then
    #         return x
    #     if key < x.key then
    #         return Recursive-Tree-Search(x.left, key)
    #     else
    #         return Recursive-Tree-Search(x.right, key)
    #     end if



    ### Remove highways that do not connect residues linked by a LD_geomection AND are not used to compose a highway that IS. ###
    # (This is effectively done by tracking the composite variables. We will choose to add all variables/constraints that a necessary highway is composed of, 
    # where they are not already in the solution. Therefore, any variables that aren't necessary won't be added.)
    return composite_highways,composition_dictionary, CA_altlocs # NOTE returning CA_altlocs is temporary


def construct_easy_highways(LD_geomections:list[LP_Input.Geomection],constr_var_dict:dict[str,dict[str,tuple[LP_Input.Geomection,pl.LpVariable]]],
                            composite_highways:dict[OrderedTag,dict[OrderedTag,tuple[pl.LpVariable,pl.LpConstraint]]],
                            composition_dictionary:dict[OrderedTag,dict[OrderedTag,DisorderedTag]],
                            CA_altlocs): # NOTE will remove CA_altlocs once bug fixed
    # Connect each ordered atom involved in LD geomection where each atom separated by 1-2 bonds from Ca.
    # Temporary/prototype method for sake of testing concept.
    
    connected_to_CA=[]
    for dID in constr_var_dict:
        is_bond=False
        if dID.startswith(RestraintsHandler.Constraint.kind(RestraintsHandler.BondRestraint)): # XXX
            is_bond=True
        elif dID.startswith(RestraintsHandler.Constraint.kind(RestraintsHandler.AngleRestraint)):
            pass
        else:
            continue
        bonded_end_tags = dID.split('_')[1],dID.split('_')[-1]  # XXX
        not_CA = tuple(v for v in bonded_end_tags if v.split('.')[1]!="CA") # XXX
        is_CA = tuple(v for v in bonded_end_tags if v.split('.')[1]=="CA") # XXX
        if len(not_CA)!=1:
            continue
        assert len(is_CA)==1
        connected_to_CA.append((DisorderedTag(not_CA[0].split('.')[0],not_CA[0].split('.')[1]),is_bond,dID)) #XXX
    
    easy_geoms_found=[]
    for geom in LD_geomections:
        d_end_tags=(geom.atom_chunks[0].get_disordered_tag(),geom.atom_chunks[-1].get_disordered_tag())
        is_bond_list=[]
        dID_list=[]
        for t in d_end_tags:
            for tag,is_bond,dID in connected_to_CA:
                if tag == t:
                    if not is_bond and (tag.atom_name() in "C","N"): # Always have a bond to Ca - don't unnecessarily make a highway variable for these.
                        continue
                    is_bond_list.append(is_bond)
                    dID_list.append(dID)
                    break
        if len(dID_list)!=2:
            continue   

        easy_geoms_found.append(geom.get_disordered_connection_id())
        
        end_tags=(geom.atom_chunks[0].get_ordered_tag(),geom.atom_chunks[-1].get_ordered_tag())
        if get_highway(*end_tags,composite_highways) is not None:
            continue
        
        altlocs_of_CAs:tuple[list[str],list[str]]=([],[]) 
        end_highways:tuple[dict[OrderedTag, dict[OrderedTag, tuple[pl.LpVariable,]]]]=({},{})
        for i, (tag,is_bond,dID) in enumerate(zip(end_tags,is_bond_list,dID_list)):
            bonded_end_tags = dID.split('_')[1],dID.split('_')[-1]
            for CA_idx_in_dID, bonded_end_tag in zip((0,-1), bonded_end_tags):
                if bonded_end_tag.split('.')[1]=="CA":
                    break
            else:
                assert False,(dID,tag)
            if is_bond:
                for altlocs, (bond,lp_var) in constr_var_dict[dID].items():
                    if tag == bond.atom_chunks[1+CA_idx_in_dID].get_ordered_tag():
                        altlocs_of_CAs[i].append(altlocs[CA_idx_in_dID])
                        LH_tag,RH_tag = bond.atom_chunks[0].get_ordered_tag(), bond.atom_chunks[1].get_ordered_tag()
                        if (LH_tag.atom_name()=="CA") ==  (i==0):
                            assert (LH_tag.atom_name()=="CA") + (RH_tag.atom_name()=="CA") == 1
                            LH_tag,RH_tag=RH_tag,LH_tag

                        if LH_tag not in end_highways[i]:
                            end_highways[i][LH_tag]={}
                        end_highways[i][LH_tag][RH_tag]=(lp_var,None)
            else:
                angles_to_make_highway:dict[str,tuple[LP_Input.Geomection,pl.LpVariable]] = {}
                already_made_highways={}
                for altlocs, (angle,lp_var) in constr_var_dict[dID].items():
                    if tag in [ch.get_ordered_tag() for ch in angle.atom_chunks]:
                        altlocs_of_CAs[i].append(altlocs[CA_idx_in_dID])
                        LH_tag=angle.atom_chunks[0].get_ordered_tag()
                        RH_tag=angle.atom_chunks[1].get_ordered_tag()
                        lookedup_highway=get_highway(LH_tag,RH_tag,composite_highways)
                        if lookedup_highway is None:
                            angles_to_make_highway[altlocs]= (angle,lp_var)
                        else:
                            if LH_tag not in already_made_highways:
                                already_made_highways[LH_tag]={}
                            already_made_highways[LH_tag][RH_tag]=lookedup_highway

                assert len(angles_to_make_highway)>1,(angles_to_make_highway,constr_var_dict[dID])
                altloc_to_CA=angle_to_highways(angles_to_make_highway)

                for key in altloc_to_CA:
                    if key not in composite_highways:
                        composite_highways[key]={}
                    composite_highways[key]|=altloc_to_CA[key]

                end_highways[i]=altloc_to_CA

                # TODO function for adding highway dicts like this
                for key in already_made_highways:
                    if key not in end_highways[i]:
                        end_highways[i][key]={}
                    end_highways[key]|=already_made_highways[key]

        # Combine into one highway
        assert (end_tags[1].resnum() > end_tags[0].resnum()) or (end_tags[0].atom_name()=="N") or (end_tags[1].atom_name()=="C") , end_tags
        central_highways={}
        assert len(altlocs_of_CAs)==2
        if end_tags[0].resnum()==end_tags[1].resnum() :
            # No need to use CA-CA highway, these connect at the same CA!
            left_to_centre=end_highways[0]
            for key in left_to_centre:
                if key not in composite_highways:
                    composite_highways[key]={}
                    composition_dictionary[key]={}
        else:
            for LH_altloc in altlocs_of_CAs[0]:
                LH_tag=OrderedTag(end_tags[0].resnum(),"CA",LH_altloc)
                central_highways[LH_tag]={}
                for RH_altloc in altlocs_of_CAs[1]:
                    RH_tag=OrderedTag(end_tags[1].resnum(),"CA",RH_altloc)
                    central_highways[LH_tag][RH_tag] = get_highway(LH_tag,RH_tag,composite_highways)
                    if central_highways[LH_tag][RH_tag] is None:
                        print(f"ERROR: missing highway between {LH_tag} and {RH_tag}. (Needed for {geom.connection_type}, {[ch.get_ordered_tag() for ch in geom.atom_chunks]}). Constructing ad hoc.")
                        add_CA_highway(LH_tag.resnum(),RH_tag.resnum(),composite_highways,composition_dictionary,CA_altlocs)
                        central_highways[LH_tag][RH_tag] = get_highway(LH_tag,RH_tag,composite_highways)
                        assert central_highways[LH_tag][RH_tag] is not None

                        
                        
                
            # for lhaltlocrhaltloc in xxx:
            #     if gethighwayisnone:
            #         left_to_centre,left_to_centre_dict = highway_from_highways(end_highways[0],central_highways)
            #     else:
            #         left_to_centre=get_highwaything
            left_to_centre,left_to_centre_dict = highway_from_highways(end_highways[0],central_highways,
                                                                       already_constructed_highway_dict=composite_highways)
            for key in left_to_centre:
                if key not in composite_highways:
                    composite_highways[key]={}
                    composition_dictionary[key]={}
                composite_highways[key]|=left_to_centre[key]
                composition_dictionary[key]|=left_to_centre_dict[key]
        left_highway=left_to_centre
        right_highway=end_highways[-1]
        left_to_right,left_to_right_dict = highway_from_highways(left_highway,right_highway)
        for key in left_to_right:
            composite_highways[key]|=left_to_right[key]
            composition_dictionary[key]|=left_to_right_dict[key]


        
    
    print(easy_geoms_found)
    print(f"Found {len(easy_geoms_found)} easy LD geoms")
    return composite_highways,composition_dictionary



def construct_CA_to_other_highways():
        raise NotImplementedError()


# def temporary_add_highway(disordered_connections:dict[str,list[LP_Input.Geomection]],constr_var_dict:dict[str,dict[str,tuple[LP_Input.Geomection,pl.LpVariable]]]):
#  ##  TEMPORARY
#     CA_to_CG_id = LP_Input.Geomection.construct_disordered_connection_id(
#         RestraintsHandler.AngleRestraint,
#         [DisorderedTag(147,"CA"),DisorderedTag(147,"CB"),DisorderedTag(147,"CG")]
#         )
def construct_highways(disordered_connections:dict[str,list[LP_Input.Geomection]],constr_var_dict:dict[str,dict[str,tuple[LP_Input.Geomection,pl.LpVariable]]]):
    LD_geomections = []
    for key, conns in disordered_connections.items():
        # TODO also include proline, aromatic rings, etc.
        if any( key.startswith(RestraintsHandler.Constraint.kind(x)) for x in (RestraintsHandler.NonbondRestraint,RestraintsHandler.ClashRestraint) ):
            LD_geomections.extend(conns)
        
    composite_highways, composition_dictionary,CA_altlocs = construct_backbone_highways(LD_geomections,constr_var_dict)
    composite_highways, composition_dictionary = construct_easy_highways(LD_geomections,constr_var_dict,composite_highways, composition_dictionary,CA_altlocs)    


    testing = True
    if testing:
        print("WARNING: TESTING WITHOUT FULL HIGHWAYS!")
        impossible_LD=[]
    else:
        construct_CA_to_other_highways()

    return composite_highways,impossible_LD





