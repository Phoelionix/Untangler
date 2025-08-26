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

# Finds assignment of altlocs that minimizes the geometric "badness" defined in LinearOptimizer.Input. 
# It may be easier to think about it as connecting up the atoms in a way that minimizes the total energy.


import pulp as pl
import numpy as np
from LinearOptimizer.Input import *
import itertools
import UntangleFunctions
#import pulp as pl
# just for residues for now



#def solve(chunk_sites: list[Chunk],connections:dict[str,dict[str,MTSP_Solver.ChunkConnection]],out_handle:str): # 
def solve(chunk_sites: dict[str,AtomChunk],connections:list[MTSP_Solver.AtomChunkConnection],out_dir,out_handle:str,force_no_flips=False,num_solutions=20,force_sulfur_bridge_swap_solutions=True,
          protein_sites:bool=True,water_sites:bool=True): # 
    # protein_sites, water_sites: Whether these can be swapped.

    # This sucks but is fixed in more_than_two branch.
    #first_atom_id = "1.N"
    first_atom_id = None
    lowest_site_num = np.inf
    for chunk_site in chunk_sites.values():
        if chunk_site.get_site_num() < lowest_site_num:
            lowest_site_num = chunk_site.get_site_num()
            first_atom_id = chunk_site.get_disordered_tag()

    nodes = [chunk.unique_id() for chunk in chunk_sites.values()]

    def get_variables(unique_id)->dict:
        d={}
        d["depth"] = unique_id.split("&")[0]
        d["site_num"],d["altloc"] = unique_id.split("&")[1].split(".")
        d["name"] = unique_id.split("&")[2]
        return d
    
    def get(unique_id,key):
        return get_variables(unique_id)[key]


    forced_swap_solutions=[]




    lp_problem = pl.LpProblem("Untangling_Problem", LpMinimize)

    

   

    ##### Constraints ####

    def pair_name(unique_id):
        return f"{get(unique_id,'name')}"
    pair_names = []
    for n in nodes:
        if pair_name(n) not in pair_names:
            pair_names.append(pair_name(n))
    
    #TODO generalise to VariableID, with any number of unique  d-Atoms (disordered atoms)
    # def BondID(a,b): # id of d-Bond (disordered bond) 
    #     ordered = sorted((a,b),key=lambda x: get(x,"site_num"))
    #     return (pair_name(ordered[0])+"_"+pair_name(ordered[1]))  #TODO make unique regardless of order of a,b
    def ConnectionID(connection:MTSP_Solver.AtomChunkConnection): # id of d-Bond (disordered bond) 
        #ordered = sorted((a,b),key=lambda x: get(x,"site_num"))
        #return (pair_name(ordered[0])+"_"+pair_name(ordered[1]))  #TODO make unique regardless of order of a,b
        return connection.get_disordered_connection_id()
    
    alt_connections:dict[str,list[MTSP_Solver.AtomChunkConnection]] ={} # options for each alt connection



    for connection in connections:
        if ConnectionID(connection) not in alt_connections:
            #assert BondID(a,b)==BondID(b,a)  
            #alt_bonds[BondID(a,b)]=[]
            alt_connections[ConnectionID(connection)]=[]
        assert connection not in alt_connections[ConnectionID(connection)] 
        alt_connections[ConnectionID(connection)].append(connection)
        



    class VariableID:
        Atom="Atom"
        Bond = "Bond"
        Nonbond="Nonbond"
        Clash="Clash"
        Angle = "Angle"

        def __init__():
            pass
    #bond_choices = {} # Each must sum to 1
    distance_vars = [] 
    var_dict:dict[VariableID,dict[str,pl.LpVariable]] = {}

    # Setup whether atoms flipped
    for i, (atom_id, chunk) in enumerate(chunk_sites.items()):
        # Only iterate through each disordered atom once
        if get(chunk.unique_id(),"altloc")=="A":
            disordered_tag = chunk.get_disordered_tag()

            # var_unflipped = pl.LpVariable(f"unflipped_{disordered_tag}",
            #                     lowBound=0,upBound=1,cat=pl.LpInteger)
            var_flip = pl.LpVariable(f"flippedAtom_{disordered_tag}",
                                lowBound=0,upBound=1,cat=pl.LpInteger)     
            var_unflipped = 1-var_flip
             
            var_dict[disordered_tag]=dict(
                unflipped=var_unflipped,
                flipped = var_flip
            )
            lp_problem += (
                var_flip+var_unflipped==1,
                f"polarity_{disordered_tag}"
            )

            
            if (force_no_flips 
                or  (not protein_sites and not chunk.is_water) 
                or (not water_sites and chunk.is_water )) :
                lp_problem += (
                    var_flip==0,
                    f"forceNoFlips_{disordered_tag}"
                )

            # Remove equivalent solutins - set solution to one where altlocs for first (disordered) atom in protein sequence is left unchanged.
            if protein_sites:
                if (chunk.get_disordered_tag()==first_atom_id):
                    lp_problem += (
                        var_flip==0,
                        f"anchorSymmetry_{disordered_tag}"
                    )
             

            # #TODO bad idea if we aren't allowing all waters to flip... 
            # #elif i == 0:
            # # XXX NOTE  "fixed" by making random...
            # elif i == random.randint(0,len(chunk_sites)-1):
            #     lp_problem += (
            #         var_flip==0,
            #         f"break_symmetry_by_fixing_{disordered_tag}_altloc"
            #     )   

    for connection_id,connection_choices in alt_connections.items():

        # NOTE: This is way overcomplicated but it's fixed in the more_than_two branch
        constraint_type = connection_id.split('_')[0]
        if constraint_type in (VariableID.Bond, VariableID.Nonbond,VariableID.Clash): #TODO make LinearOptimizer/Input.py use this too.
            var_flip = pl.LpVariable(f"flipped_{connection_id}",
                                lowBound=0,upBound=1,cat=pl.LpInteger)
            var_unflipped = pl.LpVariable(f"unflipped_{connection_id}",
                                lowBound=0,upBound=1,cat=pl.LpInteger)
            #var_unflipped = 1 - var_flip # TODO
            

            var_dict[connection_id]=dict(
                unflipped=var_unflipped,
                flipped = var_flip
            )

            distance_flipped = 0
            distance_unflipped = 0
            disordered_atom_ids=None
            for connection in connection_choices:
                a,b = [ch.unique_id() for ch in connection.atom_chunks]
                if disordered_atom_ids is None:
                    disordered_atom_ids = [ch.get_disordered_tag() for ch in connection.atom_chunks]
                else:
                    assert disordered_atom_ids == [ch.get_disordered_tag() for ch in connection.atom_chunks]
                flipped = get(a,"altloc")!=get(b,"altloc")
                if flipped:
                    distance_flipped += connection.ts_distance
                else:
                    distance_unflipped += connection.ts_distance

            distance_vars.append(var_flip * distance_flipped)
            distance_vars.append(var_unflipped * distance_unflipped)
            #assert False, (distance_vars[-1],distance_vars[-2]) 
            
            # choose one polarity
            lp_problem += (
                var_flip+var_unflipped==1,
                f"polarity_{connection_id}"
            )
            ### Constrain atom flips
            # NB:var_flip and var_unflipped are whether the bond **between** atoms is flipped
            # flipped constraint

            # This is getting it backwards. The bonds, angles, etc. should be set to be active when all the atoms they contain are active. See the more_than_two branch. 
            atom_a,atom_b = disordered_atom_ids
            lp_problem += (
                 (var_dict[atom_a]["flipped"]+  var_dict[atom_b]["flipped"])>=var_flip,
                f"atomFlippedA_{connection_id}")
            lp_problem += (
                 (var_dict[atom_a]["unflipped"]+  var_dict[atom_b]["unflipped"])>=var_flip,
                f"atomFlippedB_{connection_id}")
            lp_problem += (
                 (var_dict[atom_a]["flipped"]+  var_dict[atom_b]["unflipped"])>=var_unflipped,
                f"atomUnflippedA_{connection_id}")
            lp_problem += (
                 (var_dict[atom_a]["unflipped"]+  var_dict[atom_b]["flipped"])>=var_unflipped,
                f"atomUnflippedB_{connection_id}")
            
            if  force_sulfur_bridge_swap_solutions \
                and constraint_type == VariableID.Bond \
                and [ch.element for ch in connection.atom_chunks]==["S","S"]:
                forced_swap_solutions.append(
                    var_dict[atom_a]["flipped"]==var_dict[atom_b]["flipped"]==1 - var_dict[first_atom_id]["flipped"]
                )


    for connection_id,connection_choices in alt_connections.items():
        constraint_type = connection_id.split('_')[0]
        if constraint_type == VariableID.Angle:

            vars = []
            for connection in connection_choices:
                #tag = "_".join(ch.unique_id() for ch in connection.atom_chunks)
                tag = "_".join(ch.unique_id() for ch in connection.atom_chunks)
                if connection.hydrogen_tag!="":
                    tag = connection.hydrogen_tag+"_"+tag
                var_active = pl.LpVariable(f"{constraint_type}_{tag}",
                                    lowBound=0,upBound=1,cat=pl.LpInteger)
                vars.append(var_active)
                distance_vars.append(var_active*connection.ts_distance)
            lp_problem += (
                lpSum([var for var in vars])==2,
                f"mirroredPolarity_{connection_id}"
            )

            # 1 of each polarity
            polarity_pairs = []
            var_pairs:list[LpVariable] = []
            for connectionA,var_a in zip(connection_choices,vars):
                found_partner=False
                for connectionB,var_b in zip(connection_choices,vars):
                    if np.all([a.altloc!=b.altloc for a,b in zip(connectionA.atom_chunks,connectionB.atom_chunks)]):
                        assert not found_partner
                        if connectionB.hydrogen_name_set != connectionA.hydrogen_name_set:
                            continue
                        if (connectionB,connectionA) not in polarity_pairs:
                            polarity_pairs.append((connectionA,connectionB))
                            var_pairs.append((var_a,var_b))
                        found_partner=True
                assert found_partner


            i=0
            for var_a,var_b in var_pairs:
                #print(var_a.name,var_b.name)
                lp_problem += (
                    var_a==var_b,
                    f"mustEqual_{var_a.name}_{var_b.name}_{i}"
                )
                i+=1        
            
            


            def get_sub_vars(connection:MTSP_Solver.AtomChunkConnection)->dict[str,LpVariable]:
                subset_bond_ids = []
                is_flipped={}
                combinations_iterator:list[list[AtomChunk]] = itertools.combinations(connection.atom_chunks, 2)
                for atom_chunks in combinations_iterator:
                    bond_id = MTSP_Solver.AtomChunkConnection(atom_chunks,None,"Bond",[]).get_disordered_connection_id()
                    if bond_id not in subset_bond_ids:
                        subset_bond_ids.append(bond_id)  # HACK THIS SUCKS!!!!!
                        is_flipped[bond_id] = atom_chunks[0].altloc!=atom_chunks[1].altloc
                sub_vars = {} 
                for bond_id in var_dict.keys():   # (Note to self: probably more pythonic way to get a subset of dict according to condition like this)
                    if bond_id in subset_bond_ids:
                        subkey = "flipped" if is_flipped[bond_id] else "unflipped"
                        sub_vars[bond_id] = var_dict[bond_id][subkey]
                return sub_vars
            

            for connection,var_active in zip(connection_choices,vars): 
        
                #if subset of a constraint, must be 0 if other is 0.
                sub_vars:dict[str,LpVariable] = get_sub_vars(connection)
                print(sub_vars,var_active.name)
                print()
                for subset_var in sub_vars.values():
                    lp_problem += (
                        subset_var >= var_active,
                        f"alignParent_{subset_var.name}_{var_active.name}"
                    )
    lp_problem += (
    lpSum([dist for dist in distance_vars]),
    "badness",)

    # for constraint, vars in var_dict.items():

    #TODO if subset of a constraint, must be 0 if other is 0.
    # def constrain_subsets():
    #     for constraint, vars in var_dict.items():
    #         for polarity_var,polarity in zip(vars,("flipped","unflipped")):
    #             less_constraints_version = True
    #             if less_constraints_version:
    #                 p = LOConstraint.get_smallest_parent_var(constraint,polarity)

    #                 lp_problem += (
    #                     polarity_var <= p,
    #                     f"align_parent_{constraint}_{polarity}"
    #                 )
    #             else:
    #                 for p,parent_name in LOConstraint.get_parent_vars(constraint,polarity):
    #                     lp_problem += (
    #                         polarity_var <= p,
    #                         f"align_parent_{constraint}_{polarity}_{parent_name}",
    #                     )              
    
    # TODO: Suspect we want the two badness to be close to equal. Because if they aren't similar badness, we wouldn't expec to see both?




    # # Must have one connection 
    # for pair_name in pair_names:

    # for d in departures:
    #     arrival_sites = [a for (o,a) in routes if o ==d]

    #     #d = LpAffineExpression(dict(x_0=1, x_1=-3, x_2=4))

    #     lp_problem += (
    #     lpSum([vars[o][d] * distances[o][d] for (o, d) in routes]),
    #     "Soma_das_distancias_percorridas", )

    
    ######################
    

    assert num_solutions >= len(forced_swap_solutions) 

    bonds_to_flip_arrays=[]
    distances=[]
    for l in range(num_solutions):
        if l > 0 and l <= len(forced_swap_solutions):
            lp_problem.constraints.pop("forcedSwap")
        if l < len(forced_swap_solutions):
            lp_problem += (
                forced_swap_solutions[l],
                f"forcedSwap"
            )
        
        print()
        print(f"------- Start loop {l+1} of {num_solutions} ---------")
        if l == 0:
            print("Solving best solution")
        elif l == 1:
            print("Solving next-best solution")
        else:
            print(f"Solving {l}th best solution")
        print(f"-----------------------------------------------------")
        print()

        lp_problem.writeLP(f"{out_dir}/TSP_{out_handle}.lp")

        lp_problem.solve()

        def get_status(verbose=False):
            print("Status:", LpStatus[lp_problem.status])
            

            if verbose:
                for v in lp_problem.variables():
                    if v.varValue > 0:
                        print(v.name, "=", v.varValue)

            print(f"Target: {out_handle}")
            print("Total distance = ", value(lp_problem.objective))
            #plt.scatter()
        get_status(verbose=False)

        # dry = None
        # for val in connections.values():
        #     for v in val.values():
        #         dry=v.calculated_dry()
        #         break
        # tag = "_dry" if dry else ""
        with open(f"{out_dir}/TSP-Out{out_handle}.txt",'w') as f:
            f.write(f"Status: {LpStatus[lp_problem.status]}\n")

            for v in lp_problem.variables():
                try:
                    if v.varValue > 0:
                        f.write(f"{v.name} = {v.varValue}\n")
                except:
                    raise Exception(v,v.varValue)
            f.write(f"Total distance = {value(lp_problem.objective)}")
            

        nonzero_variables = {}
        for v in lp_problem.variables():
            n=3 # best_route_6_A_5_A or Flow_6_A_5_A, n corresponds to number of underscores from back
            #TODO wrong now that we have hydrogen tag
            var_prefix= "_".join(v.name.split("_")[:1])
            if var_prefix not in nonzero_variables:
                nonzero_variables[var_prefix]=[]
            if v.varValue > 0: #and ends_in_connection_id(v.name):
                variable_name="var"
                if len(v.name.split("_"))>1:
                    variable_name = v.name.split("_")[1:]
                nonzero_variables[var_prefix].append("_".join(variable_name))

        
        bonds_to_flip=[]
        bonds_to_flip_arrays.append(bonds_to_flip)


        assert f"flippedAtom_{first_atom_id}" in [v.name for v in lp_problem.variables()] 
        if "flippedAtom" in nonzero_variables:
            # Start with opposite polarity
            if first_atom_id in nonzero_variables["flippedAtom"]:
                bonds_to_flip.append(f"START---{first_atom_id}")
                print(f"Flagging swap of {first_atom_id}")
            for key in nonzero_variables["flippedAtom"]:
                water_start_res_num = 65 # XXX
                if int(key.split('.')[0]) >= water_start_res_num:
                    bonds_to_flip.append(f"WATER---{key}")

        if "flipped" not in nonzero_variables:
            #TODO assert that we are only looking at waters
            for chunk_site in chunk_sites.values():
                assert chunk_site.name == "O"
            nonzero_variables["flipped"] = []
        for variable in nonzero_variables["flipped"]:
            variable_type,N,M = variable.split("_")
            if variable_type != "Bond":
                continue
            N_res_num, N_name  = N.split('.') 
            M_res_num, M_name = M.split('.') 

            
            bond_to_flip = f"{N}---{M}"
            assert bond_to_flip not in bonds_to_flip
            bonds_to_flip.append(bond_to_flip)
            print(f"Flagging swap of {variable_type} {bond_to_flip}")
            
       
        solution = [var.value() for var in lp_problem.variables()] 
        lp_variables = lp_problem.variables()
        lp_problem += pulp.lpSum(solution[i]*lp_variables[i] for i in range(len(solution))) <= len(solution)-1

        distances.append(value(lp_problem.objective))

        flipped_flip_variables = []
        flip_variables=[]
        for chunk in chunk_sites.values():
            if protein_sites and chunk.is_water: # not interested in different solutions for water. 
                continue # This means water atoms may or may not swap for the single solution where no protein atoms swap.
            var = var_dict[chunk.get_disordered_tag()]["flipped"]
            val = var.varValue
            flip_variables.append(var)
            if val == 1:
                flipped_flip_variables.append(1-var)
        if len(flipped_flip_variables) == 0:
            # Require one variable to be flipped
            lp_problem += pulp.lpSum(flip_variables) >= 1, f"force_next_best_solution_{l}"
        else:
            # require at least one flip to be different
            lp_problem += pulp.lpSum(flipped_flip_variables) >= 1, f"force_next_best_solution_{l}"

        ##################

    lines:list[str] = []
    lines.append(f"Target: {out_handle}")
    for i, (distance, bonds_to_flip) in enumerate(zip(distances,bonds_to_flip_arrays)):
        lines.append("")
        lines.append("--------")
        lines.append(f"solution {i+1} with distance {distance}")
        for bond_to_flip in bonds_to_flip:
            lines.append(f"Flagging swap of {bond_to_flip}")
    
    swaps_file =  f"{out_dir}/TSP-toFlip_{out_handle}.txt"
    with open(swaps_file,'w') as f: 
        f.writelines([l +"\n" for l in lines])
    
    return swaps_file


if __name__=="__main__":
    handle = sys.argv[1]
    out_dir = f"{os.path.abspath(os.getcwd())}/output/{handle}/"
    os.makedirs(out_dir,exist_ok=True)
    pdbPath = f"{UntangleFunctions.pdb_data_dir()}/{handle}.pdb"
    finest_chunks,connections= MTSP_Solver(pdbPath).calculate_paths(
        atoms_only=True,
    )
    solve(finest_chunks,connections,handle)


# %%
