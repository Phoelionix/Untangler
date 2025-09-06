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


# TODO 
# We will quickly reach an intractable number of variables as we consider more conformations.
# (for each angle, there will be 3^n connections.)
# So it will be important to break down the problem in some way. Possibly a stochastic way.
# Also need to look at exploring different branches of solutions in an intelligent way, and with parallel processing.  


import pulp as pl
import numpy as np
from enum import Enum
from LinearOptimizer.Input import *
import itertools
import UntangleFunctions
import json
#import pulp as pl
# just for residues for now



#def solve(chunk_sites: list[Chunk],connections:dict[str,dict[str,MTSP_Solver.ChunkConnection]],out_handle:str): # 
def solve(chunk_sites: dict[str,AtomChunk],disordered_connections:dict[str,list[MTSP_Solver.AtomChunkConnection]],out_dir,out_handle:str,force_no_flips=False,num_solutions=20,force_sulfur_bridge_swap_solutions=True,
          protein_sites:bool=True,water_sites:bool=True): # 
    # protein_sites, water_sites: Whether these can be swapped.

    """
    # This sucks but is fixed in more_than_two branch. HOW
    #first_atom_id = "1.N"
    first_atom_id = None
    lowest_site_num = np.inf
    for chunk_site in chunk_sites.values():
        if chunk_site.get_site_num() < lowest_site_num:
            lowest_site_num = chunk_site.get_site_num()
            first_atom_id = chunk_site.get_disordered_tag()
    """

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


    class VariableKind(Enum):
        Atom="Atom"
        Bond = "Bond"
        Nonbond="Nonbond"
        Clash="Clash"
        Angle = "Angle"

    #Disordered variable id
    class VariableID:
        @staticmethod
        def Atom(chunk:AtomChunk):
            return VariableID(chunk.get_disordered_tag(),VariableKind.Atom)
        def __init__(self,name:str,kind:VariableKind):
            self.name=str(name)
            self.kind = kind
        def __repr__(self):
            return self.name
        def __hash__(self):
            return hash((self.name, self.kind))
        def __eq__(self,other):
            return (self.name,self.kind) == (other.name,other.kind) 
        def __ne__(self,other):
            return not(self == other) 
    #bond_choices = {} # Each must sum to 1
    distance_vars = [] 
    # NOTE When assigning LP variable names, the "class" of variable should follow format of variableName_other_stuff   (class of variable is variable_type.split("_")[0]) 
    constraint_var_dict:dict[VariableID,dict[str,pl.LpVariable]] = {} 
    site_var_dict:dict[VariableID,dict[str,dict[str,pl.LpVariable]]] ={}
    var_dictionaries = dict(
        constraints = constraint_var_dict,
        sites = site_var_dict,
    )

    # Setup where atoms are assigned.
    # Track which are which based on original altloc.

    altlocs:list[str]=[]
    for ch in chunk_sites.values():
        if ch.altloc not in altlocs:
            altlocs.append(ch.altloc)
    # def site_variable_name(site: VariableID, from_altloc:str, to_altloc:str):
    #     return f"atomState_{site}_{from_altloc}.{to_altloc}"
    disordered_atom_sites:list[str] = []
    
    # Setup atom swap variables
    for i, (atom_id, chunk) in enumerate(chunk_sites.items()):
        site = VariableID.Atom(chunk)
        if site not in disordered_atom_sites:
            disordered_atom_sites.append(site)

        from_altloc = get(chunk.unique_id(),"altloc")
        
        if site not in site_var_dict:
            site_var_dict[site] = {}
        if from_altloc not in site_var_dict[site]:
            site_var_dict[site][from_altloc]={}

        # Create variable for each possible swap to other swap
        for to_altloc in altlocs:
            var_atom_assignment =  pl.LpVariable(
                f"atomState_{site}_{from_altloc}.{to_altloc}",
                lowBound=0,upBound=1,cat=pl.LpInteger
            )
            site_var_dict[site][from_altloc][to_altloc]=var_atom_assignment

        # Each ordered atom is assigned to one conformer from:to == n:1
        lp_problem += (  
            lpSum(site_var_dict[site][from_altloc])==1,
            f"fromAltLoc_{site}_{from_altloc}"
        )

    for site in disordered_atom_sites:
        # Each conformer is assigned one ordered atom. from:to == 1:n
        for to_altloc in altlocs:
            to_altloc_vars:list[LpVariable] = []
            for from_alt_loc_dict in site_var_dict[site].values():
                to_altloc_vars.append(from_alt_loc_dict[to_altloc])
        lp_problem += (  
            lpSum(to_altloc_vars)==1,
            f"toAltLoc_{site}"
        )

        if force_no_flips:
            for altloc in altlocs:
                lp_problem += (
                    site_var_dict[site][altloc][altloc]==1,
                    f"forceNoFlips_{site}_{altloc}"
                )  


    # TODO improve terminology
    # A "connection" just refers to a group of atoms with a constraint assigned by LinearOptimizer.Input.
    # A disordered connection refers to all the *possible* groupings of these atoms.

    def add_constraints_from_disordered_connection(constraint_type:VariableKind,disordered_connection: list[MTSP_Solver.AtomChunkConnection]):
        # Rule: If all atom assignments corresponding to a connection are active,
        # then all those atoms must be swapped to the same assignment.
        nonlocal lp_problem  


        # if constraint_type == VariableKind.Nonbond:
        #     for connection in disordered_connection:
        #         connection:MTSP_Solver.AtomChunkNonBondConnection = connection
        #         assert connection.altlocs[0]!=connection.altlocs[1], connection.altlocs
        #         to_altloc = connection.to_altloc
        #         from_altloc = connection.from_altloc
        #         site = VariableID(connection.site_name,VariableKind.Atom)
        #         var = site_var_dict[site][from_altloc][to_altloc]
        #         distance_vars.append(var*connection.ts_distance)
        #         #print(var*connection.ts_distance)
        #     return
        
        for ordered_connection_option in disordered_connection:
            # Variable
            tag = "_".join([ch.unique_id() for ch in ordered_connection_option.atom_chunks])
            if ordered_connection_option.hydrogen_tag!="":
                tag = "Htag["+ordered_connection_option.hydrogen_tag+"]_"+tag
            var_active = pl.LpVariable(f"{constraint_type}_{tag}",  #TODO cat=pl.LpBinary
                                lowBound=0,upBound=1,cat=pl.LpInteger)
                
            # Connections contains ordered atoms from any and likely multiplke altlocs that 
            # *are to be assigned to the same altlocs*
            assignment_options:dict[str,list[LpVariable]]={}
            for to_altloc in altlocs:
                assignment_options[to_altloc]=[]
                for ch in ordered_connection_option.atom_chunks:
                    from_altloc = ch.get_altloc()
                    site = VariableID.Atom(ch)
                    assignment_options[to_altloc].append(site_var_dict[site][from_altloc][to_altloc])
            # if variable is inactive, cannot have all atoms assigned to the same altloc.
            # Note that since every connection option is looped through, this also means 
            # that if variable is active, all atoms will be assigned to the same altloc.
            try:
                for to_altloc, assignment_vars in assignment_options.items():
                    swaps = '|'.join([f"{assignment.name}" for assignment in assignment_vars])
                    num_assignments = len(assignment_vars)
                    lp_problem += (
                    lpSum(assignment_vars) <= num_assignments*var_active,
                    f"{var_active.name}_{swaps}"
                    )
            except:
                for to_altloc, assignment_vars in assignment_options.items():
                    print(to_altloc)
                    print(assignment_vars)
                    print('|'.join([f"{assignment.name}" for assignment in assignment_vars]))
                    print()
                for other in disordered_connection:
                    if ordered_connection_option == other:
                        print("!!!!")
                    if ordered_connection_option.atom_chunks == other.atom_chunks and ordered_connection_option.hydrogen_tag==other.hydrogen_tag:
                        print("!",ordered_connection_option)
                raise(Exception(""))
            distance_vars.append(ordered_connection_option.ts_distance*var_active)

            
            
                    
    for connection_id, ordered_connection_choices in disordered_connections.items():
        constraint_type = VariableKind[connection_id.split('_')[0]] 
        add_constraints_from_disordered_connection(constraint_type,ordered_connection_choices)

    # if  force_sulfur_bridge_swap_solutions \
    #     and [ch.element for ch in connection.atom_chunks]==["S","S"]:
    #     forced_swap_solutions.append(
    #         var_dict[atom_a]["flipped"]==var_dict[atom_b]["flipped"]==1
    #     )

    if  force_sulfur_bridge_swap_solutions:
        raise Exception("force swap solutions unimplemented")








    lp_problem += (
    lpSum([dist for dist in distance_vars]),
    "badness",)
        
    
    # TODO: Suspect we want the two badness to be close to equal. Because if they aren't similar badness, we wouldn't expec to see both?

    
    ######################
    

    assert num_solutions >= len(forced_swap_solutions) 

    site_assignment_arrays:list[dict[VariableID,dict[str,str]]]=[]
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
            



        
        site_assignments:dict[VariableID,dict[str,dict[str,int]]] = {}
        site_assignment_arrays.append(site_assignments)

        # Determine which atom has been assigned where.
        for site in site_var_dict:
            site_assignments[site]={}
            for from_altloc in site_var_dict[site]:
                to_altloc_found = False
                for to_altloc in site_var_dict[site][from_altloc]:
                    if site_var_dict[site][from_altloc][to_altloc].value()==1:
                        assert not to_altloc_found
                        to_altloc_found=True
                        site_assignments[site][from_altloc]=to_altloc
                assert to_altloc_found


        # print("KEYS")
        # for key in nonzero_variables:
        #     print(key)
        
            
       
        solution = [var.value() for var in lp_problem.variables()] 
        lp_variables = lp_problem.variables()
        lp_problem += pulp.lpSum(solution[i]*lp_variables[i] for i in range(len(solution))) <= len(solution)-1

        distances.append(value(lp_problem.objective))

        flipped_flip_variables = []
        flip_variables=[]
        for chunk in chunk_sites.values():
            if protein_sites and chunk.is_water: # not interested in different solutions for water. 
                continue # This means water atoms may or may not swap for the single solution where no protein atoms swap.
            site = VariableID.Atom(chunk)
            from_altloc = chunk.get_altloc()
            for to_altloc, var in site_var_dict[site][from_altloc].items():
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

    # Create json file that lists all the site *changes* that are required to meet the solution. 
    out_dict = {"target": out_handle,"solutions":{}}
    verbose=True
    for i, (distance, atom_assignments) in enumerate(zip(distances,site_assignment_arrays)):
        solution_dict = {"badness": distance}
        out_dict["solutions"][f"solution {i+1}"] = solution_dict
        moves:dict[str,dict[str]] = {}
        solution_dict["moves"]=moves
        for site in atom_assignments:
            site_key = f"site {site}"
            for from_altloc,to_altloc in atom_assignments[site].items():
                if from_altloc == to_altloc:
                    # No change needed
                    continue
                if verbose:
                    print(f"Flagging assignment of {site} {from_altloc} to {to_altloc}")
                if site_key not in moves:
                    moves[site_key] = {}
                moves[site_key][from_altloc] = to_altloc
    
    swaps_file =  f"{out_dir}/LO-toFlip_{out_handle}.json"
    with open(swaps_file,'w') as f: 
        json.dump(out_dict,swaps_file,indent=4)
    
    return swaps_file


if __name__=="__main__":
    handle = sys.argv[1]
    out_dir = f"{os.path.abspath(os.getcwd())}/output/{handle}/"
    os.makedirs(out_dir,exist_ok=True)
    pdbPath = f"{UntangleFunctions.pdb_data_dir()}/{handle}.pdb"
    finest_chunks,disordered_connections= MTSP_Solver(pdbPath).calculate_paths(
        atoms_only=True,
    )
    solve(finest_chunks,disordered_connections,handle)


# %%
