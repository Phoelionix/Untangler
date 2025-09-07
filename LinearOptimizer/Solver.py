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
          protein_sites:bool=True,water_sites:bool=True,max_mins_start=60,mins_extra_per_loop=0): # 
    # protein_sites, water_sites: Whether these can be swapped.

    #first_atom_id = "1.N"
    lowest_site_num = np.inf
    for chunk_site in chunk_sites.values():
        if chunk_site.get_site_num() < lowest_site_num:
            lowest_site_num = chunk_site.get_site_num()


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
            return VariableID(chunk.get_disordered_tag(),VariableKind.Atom,chunk.get_site_num())
        def __init__(self,name:str,kind:VariableKind,site_num=None):
            self.name=str(name)
            self.kind = kind
            self.site_num=site_num
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

    all_altlocs:list[str]=[]
    for ch in chunk_sites.values():
        if ch.altloc not in all_altlocs:
            all_altlocs.append(ch.altloc)
    # def site_variable_name(site: VariableID, from_altloc:str, to_altloc:str):
    #     return f"atomState_{site}_{from_altloc}.{to_altloc}"
    disordered_atom_sites:list[VariableID] = []
    
    # Setup atom swap variables
    #Crucial TODO: Too many swap variables. E.g. with two altlocs we will have four variables per atom site when we only need one.
    for i, (atom_id, chunk) in enumerate(chunk_sites.items()):
        site = VariableID.Atom(chunk)
        if site not in disordered_atom_sites:
            disordered_atom_sites.append(site)

        from_altloc = get(chunk.unique_id(),"altloc")
        
        if site not in site_var_dict:
            site_var_dict[site] = {}
        if from_altloc not in site_var_dict[site]:
            site_var_dict[site][from_altloc]={}


    for site in disordered_atom_sites:
        site_altlocs = []
        for possible_altloc in site_var_dict[site]:
            site_altlocs.append(possible_altloc)

        for from_altloc in site_altlocs:
            # Create variable for each possible swap to other swap
            for to_altloc in site_altlocs:
                var_atom_assignment =  pl.LpVariable(
                    f"atomState_{site}_{from_altloc}.{to_altloc}",
                    lowBound=0,upBound=1,cat=pl.LpBinary #TODO pl.LpBinary
                )
                # For warm start.
                var_atom_assignment.setInitialValue(0)
                if to_altloc == from_altloc: 
                    var_atom_assignment.setInitialValue(1)
                site_var_dict[site][from_altloc][to_altloc]=var_atom_assignment

            # Each ordered atom is assigned to one conformer from:to == n:1
            lp_problem += (  
                #lpSum(site_var_dict[site][from_altloc])==1,
                lpSum(site_var_dict[site][from_altloc].values())==1,
                f"fromAltLoc_{site}_{from_altloc}"
            )

        # Each conformer is assigned one ordered atom. from:to == 1:n
        for to_altloc in  site_altlocs:
            to_altloc_vars:list[LpVariable] = []
            for from_alt_loc_dict in site_var_dict[site].values():
                to_altloc_vars.append(from_alt_loc_dict[to_altloc])
            lp_problem += (  
                lpSum(to_altloc_vars)==1,
                f"toAltLoc_{site}.{to_altloc}"
            )

        if (site.site_num == lowest_site_num # Anchor solution to one where first disordered atom is unchanged.  
        or force_no_flips) :
            for altloc in site_altlocs:
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
        altlocs = None
        for ch in disordered_connection[0].atom_chunks:
            site = VariableID.Atom(chunk)
            if altlocs is None:
                altlocs = site_var_dict[site].keys()
            else:
                if set(altlocs) != set(site_var_dict[site].keys()):
                    assert False
                    #print(f"Warning: altlocs don't match, skipping {disordered_connection[0].get_disordered_connection_id()}")
                    #return
        altlocs = set(altlocs)
        for ordered_connection_option in disordered_connection:
            #continue
            # Variable
            tag = "_".join([ch.unique_id() for ch in ordered_connection_option.atom_chunks])
            if ordered_connection_option.hydrogen_tag!="":
                tag = "Htag["+ordered_connection_option.hydrogen_tag+"]_"+tag
            var_active = pl.LpVariable(f"{constraint_type}_{tag}",  #TODO cat=pl.LpBinary
                                lowBound=0,upBound=1,cat=pl.LpBinary)
            var_active.setInitialValue(0)
            if len(set([ch.get_altloc() for ch in ordered_connection_option.atom_chunks])) == 1:
                var_active.setInitialValue(1)
                
            # Connections contains ordered atoms from any and likely multiplke altlocs that 
            # *are to be assigned to the same altlocs*
            assignment_options:dict[str,list[LpVariable]]={}
            for to_altloc in altlocs:
                assignment_options[to_altloc]=[]
                for ch in ordered_connection_option.atom_chunks:
                    from_altloc = ch.get_altloc()
                    site = VariableID.Atom(ch)

                    if to_altloc not in site_var_dict[site][from_altloc]:
                        print(f"Warning: to_altloc {to_altloc} not found in {site_var_dict[site][from_altloc]}, skipping {disordered_connection[0].get_disordered_connection_id()}")
                        return

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
    swaps_file =  f"{out_dir}/LO-toFlip_{out_handle}.json"
    def update_swaps_file(distances, site_assignment_arrays):
    # Create json file that lists all the site *changes* that are required to meet the solution. 
        out_dict = {"target": out_handle,"solutions":{}}
        verbose=False
        for i, (distance, atom_assignments) in enumerate(zip(distances,site_assignment_arrays)):
            solution_dict = {"badness": distance}
            print(solution_dict)
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
                        # TODO Flag only when swaps have **changed**, as did in 2 altloc version.
                        print(f"Flagging assignment of {site} {from_altloc} to {to_altloc}")
                    if site_key not in moves:
                        moves[site_key] = {}
                    moves[site_key][from_altloc] = to_altloc
        with open(swaps_file,'w') as f: 
            json.dump(out_dict,f,indent=4)


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


        class Solver(Enum):
            COIN="COIN",
            CPLX_PY="CPLX_PY",
            CPLX_CMD="CPLX_CMD"

        extra_time_per_loop=60*mins_extra_per_loop
        timeLimitStart=60*max_mins_start
        timeLimit=timeLimitStart+l*extra_time_per_loop
        threads=26
        logPath=out_dir+out_handle+"-LP.log"
        pulp_solver = Solver.CPLX_PY
        warmStart=True
        


        solver_class=PULP_CBC_CMD
        solver_options=[]
        if pulp_solver == Solver.COIN: 
            solver_class = PULP_CBC_CMD
        elif pulp_solver == Solver.CPLX_CMD:   # cant get working
            assert False
            pulp_solver = CPLEX_CMD
        # CPLEX parameters: https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.6.0/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/tutorials/InteractiveOptimizer/settingParams.html
        # CPLEX status: https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.10.0/ilog.odms.cplex.help/refcallablelibrary/macros/Solution_status_codes.html
            solver_options.append("set parallel -1")
            #path='~/ibm/ILOG/CPLEX_STUDIO2211/cplex'
        elif pulp_solver == Solver.CPLX_PY:
            solver_class = CPLEX_PY
        else:
            raise Exception("not implemented")
        #solver = solver_class(timeLimit=timeLimit,threads=threads,logPath=logPath,warmStart=warmStart,path=path)
        solver = solver_class(timeLimit=timeLimit,threads=threads,warmStart=warmStart,logPath=logPath)
        lp_problem.solve(solver)
        print(solver)
        # if pulp_solver == Solver.COIN:
        #     # solver = PULP_CBC_CMD(threads=NTHREADS)
        #     # lp_problem.solve(solver=solver)
        #     maxNodes=None
        #     #pulpTestAll()
        #     #asdds
        #     solver = PULP_CBC_CMD(logPath=,threads=threads,timeLimit=timeLimit,maxNodes=maxNodes)
        # elif pulp_solver == Solver.COIN:
        #     solver = CPLEX_PY(timeLimit=timeLimit,threads=threads)
        # lp_problem.solve(solver)



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
        distances.append(value(lp_problem.objective))
        
        update_swaps_file(distances,site_assignment_arrays)
      
       
        solution = [var.value() for var in lp_problem.variables()] 
        lp_variables = lp_problem.variables()
        lp_problem += pulp.lpSum(solution[i]*lp_variables[i] for i in range(len(solution))) <= len(solution)-1


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
