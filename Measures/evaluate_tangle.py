import os, sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent))
from LinearOptimizer.Input import *
from LinearOptimizer import Solver 
from UntangleFunctions import parse_symmetries_from_pdb, prepare_pdb, UNTANGLER_WORKING_DIRECTORY
from LinearOptimizer.Swapper import Swapper

# How many bonds need to change to get it right?

def evaluate_tangle(model, ground_truth):
    def get_out_path(model_handle,out_tag):
        output_dir = os.path.join(UNTANGLER_WORKING_DIRECTORY,"output","")
        return f"{output_dir}{model_handle}_{out_tag}.pdb"
    
    # Find ground-truth conformers nearest to model conformers
    model_lookup=OrderedAtomLookup(
         PDBParser().get_structure("struct",model).get_atoms()
        )
    truth_lookup=OrderedAtomLookup(
         PDBParser().get_structure("struct",ground_truth).get_atoms()
        )
    force_solution_reference:dict[OrderedTag,str]={}
    def dist(a:Atom,b:Atom,translation_a,translation_b):
            a_shifted=a.get_coord()+translation_a
            b_shifted=b.get_coord()+translation_b
            #print(a_shifted,b_shifted)
            return np.sum((a_shifted-b_shifted)**2)
    
    translations_model={}; translations_truth={}

    def get_avg_coord_dict(lookup:OrderedAtomLookup): # average position of conformers for an atom
        avg_coord_dict={}
        for res_num in lookup.better_dict:
            avg_coord_dict[res_num]={}
            for name in lookup.better_dict[res_num]:
                avg=0
                for altloc, atom in lookup.better_dict[res_num][name].items():
                    avg+=atom.get_coord()
                avg/=len(lookup.better_dict[res_num][name])
                avg_coord_dict[res_num][name]=avg
        return avg_coord_dict
    avg_dict_model=get_avg_coord_dict(model_lookup)
    avg_dict_truth=get_avg_coord_dict(truth_lookup)

    for res_num in truth_lookup.better_dict:
        for name in truth_lookup.better_dict[res_num]:
            truth_conformers_assigned=[]
            #print(res_num,name)
            for model_altloc, model_atom in model_lookup.better_dict[res_num][name].items():
                min_distance=np.inf
                closest_truth_conformer_altloc=None
                for truth_altloc, truth_atom in truth_lookup.better_dict[res_num][name].items():
                    this_distance=dist(model_atom,truth_atom,-avg_dict_model[res_num][name],-avg_dict_truth[res_num][name])
                    if this_distance<min_distance:
                         closest_truth_conformer_altloc=truth_altloc
                         min_distance=this_distance
                #print(min_distance,closest_truth_conformer_altloc)
                # We are assuming that there is no case where a ground-truth conformer is the closest conformer for two model atoms. TODO handle this by finding whatever minimizes RMSD
                assert closest_truth_conformer_altloc not in truth_conformers_assigned, "Error, unhandled case" # TODO
                truth_conformers_assigned.append(closest_truth_conformer_altloc)
                force_solution_reference[OrderedTag(res_num,name,model_altloc)]=closest_truth_conformer_altloc
                if is_atom(model_atom,7,"CA",model_altloc):
                    print(f"{model_altloc}>>{closest_truth_conformer_altloc}")
            
 
    # Pass model conformers into LinearOptimizer.Input, with ground-truth conformers as reference (labelled by serial number) so it can forbid any other changes.
    weight_factors = {
        ConstraintsHandler.BondConstraint: 0.1,
        ConstraintsHandler.AngleConstraint: 80,#1,
        ConstraintsHandler.NonbondConstraint: 0, # 0.1
        ConstraintsHandler.ClashConstraint: 0, 
        ConstraintsHandler.TwoAtomPenalty: 0,
    }
    
    symmetries=parse_symmetries_from_pdb(model)
    model_handle=os.path.basename(model)[:-4]
    untangle_fmted_model = get_out_path(model_handle,"fmted")
    prepare_pdb(model,untangle_fmted_model)
    model = untangle_fmted_model
    
    LP_Input.prepare_geom_files(model,None)
    atoms, connections = LP_Input(model, model, None, symmetries).calculate_paths(
        scoring_function=ConstraintsHandler.log_chi,
        constraint_weights=weight_factors,
        force_solution_reference=force_solution_reference
    )

    # Turn off symmetry breaking thing in LinearOptimizer.Solver, and instead add a cost for the number of changes.
    out_handle="cheat_"+model_handle
    change_punish_factor=0
    if change_punish_factor<=0:
        print("WARNING: not doing minimal changes. Returned tangle level likely to be much higher than it is.")
    swaps_file_path,bonds_replaced_each_loop=Solver.solve(atoms,connections,out_dir=os.path.join(UNTANGLER_WORKING_DIRECTORY,"output",""),
                out_handle=out_handle,
                num_solutions=1,
                modify_forbid_conditions=False,
                change_punish_factor=change_punish_factor#0.01 # Need to make non-zero 
                )
    swapper=Swapper()
    swapper.add_candidates(swaps_file_path) #
    working_model, swapGroup = swapper.run(model)
    # Read ChangedConnections to see where bond changes occur.

    bonds_replaced=bonds_replaced_each_loop[0]
    out_str=""
    all_terms=[]
    for bond in bonds_replaced:
        assert len(bond.atom_chunks)==2
        tmp = zip(bond.atom_names,bond.res_nums,bond.from_altlocs)
        terms = [f"{resnum}.{name}.{altloc}" for name,resnum,altloc in tmp]
        out_str+= ' '.join(terms)+"\n"
        all_terms.extend(terms)

    bond_changes_file="tangled_bonds.txt"
    with open(bond_changes_file,"w") as f:
        f.write(out_str)
    print(f"Necessary bond changes written to {bond_changes_file}")
    tangled_bonds_path=get_out_path(model_handle,"tangled_bonds")
    write_to_b_factors(model,all_terms,tangled_bonds_path)
    print(f"Written bad bonds to B factors in {tangled_bonds_path}")

    return len(bonds_replaced)



def write_to_b_factors(pdb_path,terms,out_path): # 0/1 if in/not in terms, 1 i.
            
        def replace_B_factor(line,b):
            # b = str(b)
            # b=f"{b:.3f}"
            assert type(b)==int
            b = str(b)
            b = ' '*(6-len(b))+b
            return line[:60] + b + line[66:]
        

        with open(pdb_path) as I:

            out_lines=[]
            for line in I:
                if line.startswith("TER") or line.startswith("ANISOU"):
                    continue
                start_strs_considered = ["ATOM","HETATM"]
                for s in start_strs_considered:
                    if line.startswith(s):
                        break
                else: 
                    out_lines.append(line)
                    continue 

                # Atom entries
                altloc = line[16]
                name = line[12:16].strip()
                resnum = int(line[22:26])
               
                if f"{resnum}.{name}.{altloc}" in terms:
                    out_lines.append(replace_B_factor(line,1))
                else:
                    out_lines.append(replace_B_factor(line,0))
                    
        with open(out_path,'w') as O:
            O.writelines(out_lines)

if __name__ == "__main__":
     evaluate_tangle(*sys.argv[1:])

