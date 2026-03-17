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

from Tag import DisorderedTag,OrderedTag
import pulp as pl
import Input
from Input import LP_Input, get_altlocs_key
from VariableID import VariableID
from ConstraintsHandler import ConstraintsHandler


def angle_to_highways(disordered_angles:dict[str,tuple[LP_Input.Geomection,pl.LpVariable]]):
    assert len(disordered_angles)>0

    highway_dict: dict[str,dict[str, list[pl.LpVariable]]] = {}
    LH_tag = RH_tag = None
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
        if LH_altloc not in highway_dict:
            highway_dict[LH_altloc]={}
        if RH_altloc not in highway_dict[LH_altloc]:
            highway_dict[LH_altloc][RH_altloc]=[]
        highway_dict[LH_altloc][RH_altloc].append(var)

        if LH_tag is None:
            LH_tag=DisorderedTag(angle_geomection.res_nums[0],angle_geomection.atom_names[0])
            RH_tag=DisorderedTag(angle_geomection.res_nums[0],angle_geomection.atom_names[0])
        else:
            assert LH_tag == DisorderedTag(angle_geomection.res_nums[0],angle_geomection.atom_names[0])
            assert RH_tag == DisorderedTag(angle_geomection.res_nums[0],angle_geomection.atom_names[0])
    #return {altlocs_key: pl.lpSum(ILP_vars) for altlocs_key, ILP_vars in highway_dict.items()}
    #highway_dict: dict[str, list[pl.LpAffineExpression]]= {altlocs_key: pl.lpSum(ILP_vars) for altlocs_key, ILP_vars in highway_dict.items()}
    
    highway_var_lpConstr_dict:dict[str,dict[str,tuple[pl.LpVariable,pl.LpConstraint]]] = {LH_altloc:{} for LH_altloc in highway_dict}
    for LH_altloc in highway_dict:
        for RH_altloc,ILP_vars in highway_dict[LH_altloc].items():
            tag={LH_tag.ordered_tag(LH_altloc)};{RH_tag.ordered_tag(RH_altloc)}
            var_active = pl.LpVariable(f"Highway_{tag}",  #TODO cat=pl.LpBinary
                                lowBound=0,upBound=1,cat=pl.LpBinary)
            
            constr = (var_active==pl.lpSum(ILP_vars),
                    f"HighwayConstraint_{tag}")
            highway_var_lpConstr_dict[LH_altloc][RH_altloc]=(var_active,constr)

    return highway_var_lpConstr_dict

def highway_from_highways(LH_tag:DisorderedTag,RH_tag:DisorderedTag,
                          LH_highways:dict[str,dict[str,tuple[pl.LpVariable,pl.LpConstraint]]],
                          RH_highways:dict[str,dict[str,tuple[pl.LpVariable,pl.LpConstraint]]]):
    composite_highway_var_dict:dict[str,dict[str,list[pl.LpVariable]]] = {LH_altloc:{} for LH_altloc in LH_highways}
    highway_var_lpConstr_dict:dict[str,dict[str,tuple[pl.LpVariable,pl.LpConstraint]]] = {LH_altloc:{} for LH_altloc in LH_highways}
    
    for LH_altloc in LH_highways:
        for centre_altloc, (LH_highway_var,_) in LH_highways[LH_altloc].items():
            for RH_altloc, (RH_highway_var,_) in RH_highways[centre_altloc].items():
                if RH_altloc not in  composite_highway_var_dict[LH_altloc]:
                    composite_highway_var_dict[LH_altloc][RH_altloc]=[LH_highway_var]
                composite_highway_var_dict[LH_altloc][RH_altloc].append(RH_highway_var)
        
        # Now add the new composite highway variable
        for RH_altloc,ILP_vars in composite_highway_var_dict[LH_altloc].items():
            tag={LH_tag.ordered_tag(LH_altloc)};{RH_tag.ordered_tag(RH_altloc)}
            var_active = pl.LpVariable(f"Highway_{tag}",  #TODO cat=pl.LpBinary
                                lowBound=0,upBound=1,cat=pl.LpBinary)
            
            constr = (var_active==pl.lpSum(ILP_vars),
                    f"HighwayConstraint_{tag}")
            highway_var_lpConstr_dict[LH_altloc][RH_altloc]=(var_active,constr)
    

    return highway_Var_lpConstr_dict,composite_highway_var_lpConstr_dict

def construct_backbone_highways(LD_geomection_tags:DisorderedTag,constr_var_dict:dict[str,dict[str,tuple[LP_Input.Geomection,pl.LpVariable]]]): # LD = long-distance.
    # For each residue, determine which residues they are connected to by a LD_geomection. 
    C_beta_link_dict:dict[int,int] 
    bond_var_dict = {k:v for k, v in constr_var_dict.items() if v[0].connection_type is ConstraintsHandler.BondConstraint}
    angle_var_dict = {k:v for k, v in constr_var_dict.items() if v[0].connection_type is ConstraintsHandler.AngleConstraint}

    max_resnum= max(max(v[0].res_nums) for v in constr_var_dict.values())
    
    
    ### Construct highways so as to link all C betas ###
    CB_to_C_highways:dict[int,dict[str, dict[str, tuple[pl.LpVariable,pl.LpConstraint]]]]={} # [resnum][LH_from_altloc][RH_from_altloc]
    N_to_CB_highways:dict[int,dict[str, dict[str, tuple[pl.LpVariable,pl.LpConstraint]]]]={} # [resnum][LH_from_altloc][RH_from_altloc]
    C_to_N_bonds:dict[int,dict[str, dict[str, tuple[pl.LpVariable,None]]]]={} # [resnum][LH_from_altloc][RH_from_altloc]
    for i in range(1,max_resnum): # NOTE 1 iteration for each N-C bond. 
        # 1. Create the 3-atom highways. 

        CB_to_C_id = LP_Input.Geomection.construct_disordered_connection_id(
            ConstraintsHandler.AngleConstraint,
            [DisorderedTag(i,"CB"),DisorderedTag(i,"CA"),DisorderedTag(i,"C")]
            )
        C_to_N_id = LP_Input.Geomection.construct_disordered_connection_id(
            ConstraintsHandler.BondConstraint,
            [DisorderedTag(i,"C"),DisorderedTag(i+1,"N")]
            )
        N_to_CB_id = LP_Input.Geomection.construct_disordered_connection_id(
            ConstraintsHandler.AngleConstraint,
            [DisorderedTag(i+1,"N"),DisorderedTag(i+1,"CA"),DisorderedTag(i+1,"CB")]
            )
        
        for bond_geomection, var in bond_var_dict[C_to_N_id].values():
            LH_altloc,RH_altloc = (bond_geomection.from_altlocs[0],
                        bond_geomection.from_altlocs[-1])
            if LH_altloc not in C_to_N_bonds[i]:
                C_to_N_bonds[i][LH_altloc]={}
            
            C_to_N_bonds[i][LH_altloc][RH_altloc]=(var,None)
            
        CB_to_C=CB_to_C_highways[i]  = angle_to_highways(angle_var_dict[CB_to_C_id])
        C_to_N=C_to_N_bonds[i]
        N_to_CB=N_to_CB_highways[i+1]= angle_to_highways(angle_var_dict[N_to_CB_id])      
        
        # 2. Create 4-atom highways from Cb-Ca-C highway, C-N bond (equiv. to 2-atom highway).
        CB_to_N,CB_to_N_composition=highway_from_highways(CB_to_C,C_to_N)
        # 3. Create 6-atom highways from Cb-Ca-C-N highway, N-Ca-Cb highways
        CB_to_CB,CB_to_CB_composition=highway_from_highways(CB_to_N,N_to_CB)

    # 4. Construct highways composed from CB_to_CB highways. 
    composite_CB_highways:dict[int,dict[str, dict[str, list[tuple[pl.LpVariable,pl.LpConstraint]]]]]={} # List of all variable/constraint tuples that the highway is composed of.  



    ### Remove highways that do not connect residues linked by a LD_geomection AND are not used to compose a highway that IS. ###
    # (This is effectively done by tracking the composite variables. We will choose to add all variables/constraints that a necessary highway is composed of, 
    # where they are not already in the solution. Therefore, any variables that aren't necessary won't be added.)
    return composite_CB_highways
    

    


