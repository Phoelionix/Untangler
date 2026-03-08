# Prototype of conformation tree system, creates a second "layer" of altlocs (splits are only made off the "main" altlocs) 
#%%
import sys,pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent))
from LinearOptimizer.Tag import *
from LinearOptimizer.ConstraintsHandler import ConstraintsHandler
from LinearOptimizer.Input import LP_Input
from LinearOptimizer.OrderedAtomLookup import OrderedAtomLookup
from Bio.PDB import PDBParser,Structure,PDBIO
from UntangleFunctions import parse_symmetries_from_pdb
import os

# refinement {
#   geometry_restraints.edits {
#     bond {
#       action = *add
#       atom_selection_1 = "name O and resseq 1 and chain S and altid A"
#       atom_selection_2 = "name O and resseq 1 and chain S and altid B"
#       distance_ideal = 0.4436
#       sigma = 0.2
#     }
#   }
# }


# phenix-2.0-5793/lib/python3.9/site-packages/cctbx/geometry_restraints/bond.h
# //! weight * delta_slack**2.
# /*! See also: Hendrickson, W.A. (1985). Meth. Enzym. 115, 252-270.
# */
# double
# residual() const {
# // unlike the dihedral angle restraint, the harmonic potential is
# // always used if distance_model < distance_ideal, to compensate for
# // the lack of a nonbonded restraint for the bonded atoms.
# if ((top_out) && (delta_slack < 0)) {
#     double top = weight * limit * limit;
#     //top*(1-exp(-weight*x**2/top))
#     return top * (1.0-std::exp(-weight*delta_slack*delta_slack/top));
# } else {
#     return weight * scitbx::fn::pow2(delta_slack);
# }
# }


# TODO Dont create restraints if they already exist.

def create_all_child_restraints(model_path,child_parent_altlocs_dict:dict,include_nonbonds=True):

    print("creating sub-conformation restraints")


    # child_parent_altlocs_dict. Dictionary of altlocs to have interactions with in addition to itself. e.g. {"A":"","B":"","C":"A"}. <- Note here the A and B keys could be left out.
    # If there are two altlocs, it should mean that one of the altlocs is a child of another. e.g. {"C":"A","D":"CA"}
    def get_out_path(model_handle,out_tag):
        output_dir = os.path.join(UntangleFunctions.UNTANGLER_WORKING_DIRECTORY,"output","")
        return f"{output_dir}{model_handle}_{out_tag}.pdb"
    fmted_model = get_out_path(UntangleFunctions.model_handle(model_path),"fmtd")
    UntangleFunctions.prepare_pdb(model_path,fmted_model,
                ring_name_grouping=False)


    parent_altlocs=[]
    for val in child_parent_altlocs_dict.values():
        parent_altlocs+=val
    parent_altlocs = list(set(parent_altlocs)) 

    subset_model_path=LP_Input.create_altloc_subset_model(fmted_model,parent_altlocs)
    LP_Input.prepare_geom_files(subset_model_path,None)
    struct=PDBParser().get_structure("struct",subset_model_path)
    ordered_atom_lookup=OrderedAtomLookup(struct.get_atoms(),waters=True,altloc_subset=parent_altlocs)


    
    child_atom_tags=list(set(
        [DisorderedTag.from_atom(a) for a in OrderedAtomLookup(fmted_model, waters=False).ordered_atoms if a.get_altloc() in child_parent_altlocs_dict]
        ))
    print(child_atom_tags[:100])
    all_ordered_tags = [OrderedTag.from_atom(a) for a in OrderedAtomLookup(fmted_model, waters=False).ordered_atoms]

    
    
    assert len(child_atom_tags)>0
    assert len(all_ordered_tags)>0

    constraints_handler=ConstraintsHandler()
    constraints_to_skip=[ConstraintsHandler.ClashConstraint,ConstraintsHandler.TwoAtomPenalty]
    if not include_nonbonds:
        constraints_to_skip.append(ConstraintsHandler.NonbondConstraint)
    symmetries=parse_symmetries_from_pdb(model_path)
    constraints_handler.load_all_constraints(subset_model_path,ordered_atom_lookup,
                                             symmetries=symmetries,
                                             water_water_nonbond=False,
                                             constraints_to_skip=constraints_to_skip,
                                             all_restraints_mode=True)
    
    all_atom_lookup=OrderedAtomLookup(fmted_model,waters=True)
    print("Writing constraints")
    text="""refinement {
  geometry_restraints.edits {\n"""
    for child_altloc,parent_altlocs in child_parent_altlocs_dict.items():
        if len(parent_altlocs)==0:
            continue 
        # Create all restraints for single child altloc
        text+=create_child_restraints(child_altloc,parent_altlocs,child_atom_tags,all_ordered_tags,constraints_handler,
                                      ordered_atom_lookup.chain_dict,all_atom_lookup,symmetries)
    text+="""  }
}\n"""
    return text

def create_child_restraints(child_altloc,parent_altlocs,child_atom_tags:list[DisorderedTag],all_ordered_tags:list[OrderedTag],
                            constraints_handler:ConstraintsHandler, chain_dict,
                            all_atom_lookup:OrderedAtomLookup,symmetries):

    allowed_constraints = [
        ConstraintsHandler.BondConstraint,
        ConstraintsHandler.AngleConstraint,
        ConstraintsHandler.Dihedral,
        ConstraintsHandler.Planarity,
        ConstraintsHandler.NonbondConstraint,
    ]
    # Create all geometry restraints for child atoms to mimic their parents.
    text=""
    processed_constraints:list[ConstraintsHandler.Constraint]=[]
    for disordered_tag, constraints in constraints_handler.atom_constraints.items():
        if not disordered_tag in child_atom_tags:
            continue
        for constraint in constraints:
            if constraint in processed_constraints:
                continue
            if type(constraint) not in allowed_constraints:
                continue
            parent_site_tags= [site_tag for site_tag in constraint.site_tags 
                               if (site_tag not in child_atom_tags 
                                   and site_tag.ordered_tag(child_altloc) not in all_ordered_tags)]
            if len(parent_site_tags)==0:
                continue # Constraint does not involve a parent altloc
            for parent_altloc in parent_altlocs:
                atom_selection_lines=[]
                for i, site_tag in enumerate(constraint.site_tags):
                    if site_tag in parent_site_tags:
                        line_altloc=parent_altloc
                    else:
                        line_altloc=child_altloc
                    atom_selection_lines.append(
                        f"      atom_selection_{i+1} = name {site_tag.atom_name()} and resseq {site_tag.resnum()} and chain {chain_dict[site_tag.resnum()]} and altid {line_altloc}"
                    )
                if type(constraint) == ConstraintsHandler.Planarity:
                    atom_selection_lines = [line[line.index('=')+1:] for line in atom_selection_lines ]
                    atom_selection_lines='      atom_selection = (' + ') \\\n      or ('.join(atom_selection_lines)+ ')'
                else:
                    atom_selection_lines='\n'.join(atom_selection_lines)
                parameter_scope_name =constraint.get_str_rep_kind().lower()
                if type(constraint) != ConstraintsHandler.NonbondConstraint:
                    ideal = constraint.ideal
                else:
                    def include_nonbond(crystal_packing=False):
                        assert not crystal_packing, "Unimplemented"
                        atoms = [all_atom_lookup.better_dict[site_tag.resnum()][site_tag.atom_name()][
                                parent_altloc if site_tag in parent_site_tags else child_altloc] 
                                for site_tag in constraint.site_tags] 
                        r0,r0_sym = constraint.altlocs_vdw_dict[None] # XXX TODO
                        r, r_sym_min = ConstraintsHandler.NonbondConstraint.symm_min_separation(atoms[0],atoms[1],symmetries)
                        threshold_factor=1.5
                        return r is not None and r/r0<=threshold_factor
                    if not include_nonbond():
                        continue
                            
                    #  ###
                    #  involves_water=any([st.atom_name()=="O" for st in constraint.site_tags])
                    #  is_H = [st.element()=="H" for st in constraint.site_tags]
                    #  if any(is_H) and not all(is_H) and not involves_water:
                    #      continue  # Creating bonds between H and non-H messes with riding H logic in phenix.  # Does it?
                    #  ####
                    if (parent_altloc,parent_altloc) not in constraint.altlocs_vdw_dict:
                        if None in constraint.altlocs_vdw_dict:
                            ideal_same_asu,ideal_crystal_packing = constraint.altlocs_vdw_dict[None]
                        else:
                            continue
                    else:
                     # TODO crystal-packing?
                        ideal_same_asu,ideal_crystal_packing = constraint.altlocs_vdw_dict[(parent_altloc,parent_altloc)]  
                    ideal = ideal_same_asu  
                    parameter_scope_name="bond" # FIXME

                if type(constraint) in [ConstraintsHandler.AngleConstraint,ConstraintsHandler.Dihedral]:
                    ideal_variable_name="angle_ideal"
                else:
                    ideal_variable_name="distance_ideal"

                text += (f"    {parameter_scope_name}"+" {\n"
                +"      action = *add\n"
                +f"{atom_selection_lines}\n")
                if type(constraint)!=ConstraintsHandler.Planarity:
                    text+=f"      {ideal_variable_name} = {ideal:.4f}\n"
                    
                text+=((f"      sigma = {constraint.sigma:.4f}\n" if constraint.sigma is not None else "      sigma = 1\n      "+f"limit = 0\n      "+f"top_out = True"+"\n") 
                + "    }\n")
                processed_constraints.append(constraint)
    return text

if __name__ == "__main__":
    model_path="/home/speno/Untangler/data/4PSS_6conf18conf_reduced_unrefined.pdb"
    
    include_nonbonds=True
    text=create_all_child_restraints(model_path,
                                     child_parent_altlocs_dict=UntangleFunctions.TEMP_TEST_CHILD_PARENT,
                                     include_nonbonds=include_nonbonds)
    out_path=os.path.join(UntangleFunctions.UNTANGLER_WORKING_DIRECTORY,"ConformationTree","output",f"split_conformations_restraints{'' if include_nonbonds else '_noNB'}.eff")
    with open(out_path,"w") as f:
        f.write(text)
# %%
