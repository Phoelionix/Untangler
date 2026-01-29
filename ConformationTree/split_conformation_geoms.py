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


# TODO Dont create restraints if they already exist.

def create_all_child_restraints(model_path,altloc_parents_dict:dict,child_atom_tags:list[DisorderedTag],all_ordered_tags:list[OrderedTag]):

    print("creating sub-conformation restraints")


    # altloc_parents_dict. Dictionary of altlocs to have interactions with in addition to itself. e.g. {"A":"","B":"","C":"A"}. <- Note here the A and B keys could be left out.
    # If there are two altlocs, it should mean that one of the altlocs is a child of another. e.g. {"C":"A","D":"CA"}

    LP_Input.prepare_geom_files(model_path,None)
    struct=PDBParser().get_structure("struct",model_path)
    ordered_atom_lookup=OrderedAtomLookup(struct.get_atoms(),waters=True)
    constraints_handler=ConstraintsHandler()
    constraints_handler.load_all_constraints(model_path,ordered_atom_lookup,symmetries=parse_symmetries_from_pdb(model_path),water_water_nonbond=False,
                                             constraints_to_skip=[ConstraintsHandler.ClashConstraint,ConstraintsHandler.TwoAtomPenalty])
    text="""refinement {
  geometry_restraints.edits {\n"""
    for child_altloc,parent_altlocs in altloc_parents_dict.items():
        if len(parent_altlocs)==0:
            continue 
        # Create all restraints for single child altloc
        text+=create_child_restraints(child_altloc,parent_altlocs,child_atom_tags,all_ordered_tags,constraints_handler)
    text+="""  }
}\n"""
    return text

def create_child_restraints(child_altloc,parent_altlocs,child_atom_tags:list[DisorderedTag],all_ordered_tags:list[OrderedTag],
                            constraints_handler:ConstraintsHandler, chain="A"):

    allowed_constraints = [
        ConstraintsHandler.BondConstraint,
        ConstraintsHandler.AngleConstraint,
        ConstraintsHandler.NonbondConstraint,
    ]
    # Create all geometry restraints for child atoms to mimic their parents.
    text=""
    for disordered_tag, constraints in constraints_handler.atom_constraints.items():
        if not disordered_tag in child_atom_tags:
            continue
        for constraint in constraints:
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
                    atom_selection_lines.append(f"      atom_selection_{i+1} = name {site_tag.atom_name()} and resseq {site_tag.resnum()} and chain {chain} and altid {line_altloc}")
                atom_selection_lines='\n'.join(atom_selection_lines)
                parameter_scope_name =constraint.get_str_rep_kind().lower()
                if type(constraint) != ConstraintsHandler.NonbondConstraint:
                    ideal = constraint.ideal
                else:
                     if (parent_altloc,parent_altloc) not in constraint.altlocs_vdw_dict:
                         continue
                     # TODO crystal-packing?
                     ideal_same_asu,ideal_crystal_packing = constraint.altlocs_vdw_dict[(parent_altloc,parent_altloc)]  
                     ideal = ideal_same_asu  
                     parameter_scope_name="bond" # FIXME
                     nonbond_limit_param=ideal # TODO check this does what expect... 
                if type(constraint) == ConstraintsHandler.AngleConstraint:
                    ideal_variable_name="angle_ideal"
                else:
                    ideal_variable_name="distance_ideal"

                text += (f"    {parameter_scope_name}"+" {\n"
                + f"""      action = *add
{atom_selection_lines}
      {ideal_variable_name} = {ideal:.4f}\n"""
+(f"      sigma = {constraint.sigma:.4f}\n" if constraint.sigma is not None else "      sigma = None\n      "+f"limit = {nonbond_limit_param}"+"\n") 
                + "    }\n")
    return text

if __name__ == "__main__":
    model_path="/home/speno/Untangler/output/2conf_start.pdb"
    child_atom_tags=[DisorderedTag(51,"CB"),]
    parent_altlocs={"C":"A"}
    text=create_all_child_restraints(model_path,parent_altlocs,child_atom_tags)
    with open("split_conformations_restraints.eff","w") as f:
        f.write(text)
# %%
