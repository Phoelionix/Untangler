# "It's a twister! It's a twister!"
import numpy as np
import re
import os, sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent))
from UntangleFunctions import H_get_parent_fullname, UNTANGLER_WORKING_DIRECTORY, PDB_Atom_Entry
import shutil
from LinearOptimizer.ConstraintsHandler  import DisorderedTag 
import numpy as np
def apply_untwists(model_path, untwist_file):

    def create_untwisted(model_path,untwist_move:tuple[str,DisorderedTag,str,np.typing.NDArray],out_path,allow_overwrite=False):
        _,site,altlocs,coords = untwist_move 

        def line_coords(line:str):
            return np.array([ float(v) for v in (line[30:38],line[38:46],line[46:54]) ])

        if not allow_overwrite:
            assert os.path.abspath(out_path)!=os.path.abspath(model_path)

        with open(model_path,'r') as f:
            new_lines = ""
            parent_coords_delta={a:None for a in altlocs}
            found_riding_H=False
            found_resname=None
            for line in f.readlines():
                # Don't change irrelevant lines
            
                P = PDB_Atom_Entry(line)
                #if not P.valid:
                if not P.valid or P.altloc not in altlocs:
                    pass
                elif site.is_entry(P):
                    if found_resname is not None:
                        assert found_resname==P.res_name
                    else:
                        found_resname=P.res_name
                    idx = altlocs.index(P.altloc)
                    new_lines+=P.new_line(new_coord=coords[idx])
                    parent_coords_delta[P.altloc] = coords[idx]-line_coords(line)
                    continue
                elif site.is_riding_H_entry(P):
                    idx = altlocs.index(P.altloc)
                    assert parent_coords_delta[P.altloc] is not None, "H's ridden atom hasn't been parsed"
                    #assert parent_coords_delta[P.altloc] != "USED", f"{site} is considered riding H, but {riding_H_found} was previously identified as the riding H!"
                    new_lines+=P.new_line(new_coord=line_coords(line)+parent_coords_delta[P.altloc])
                    #parent_coords_delta[P.altloc] = "USED"; riding_H_found=site
                    found_riding_H=True
                    continue

                # Do not modify
                new_lines+=line
                continue 
        assert found_riding_H or site.atom_name()=="C" or site.atom_name()=="O" \
        or (site.atom_name()=="N" and found_resname=="PRO"), site

        with open(out_path,'w') as f:
            f.writelines(new_lines)


    def get_candidate_model_dir(tag=""): # NOTE Set tag to the loop number
        assert model_path[-4:]==".pdb", model_path
        tag = tag if tag in ["",None] else f"_{tag}"
        model_handle = os.path.basename(model_path)[:-4]
        output_dir = os.path.join(UNTANGLER_WORKING_DIRECTORY,"output","")
        return os.path.join(output_dir,"twistOptions",f"{model_handle}{tag}","")
    
    def changed_entries_only(original,changed,out_path): # Very brittle code
        with open(original,'r') as I, open(changed,'r') as F:
            new_lines = ""
            for lineI,lineF in zip(I.readlines(),F.readlines()):
            
                P = PDB_Atom_Entry(lineI)
                #if not P.valid:
                if not P.valid or P.altloc not in altlocs:
                    assert lineI==lineF
                elif lineI!=lineF:
                    pass 
                else:
                    continue # Same lines
                new_lines+=lineF
        with open(out_path,'w') as f:
            f.writelines(new_lines)
    

    print("Applying untwists")
    untwist_candidates = {}
    with open(untwist_file) as f:
        for line in f:
            if line[0]=="#":
                continue
            site_string,altlocs = line.split('[')[0].split()
            coord_strings =  re.findall(r"\[([^\]]*)\]", line) #re.findall(r"\[[^\]]*\]", line)
            requires_altloc_optimization = line.split(']')[-1].strip() == "Y"

            assert len(altlocs)==len(coord_strings), (altlocs,coord_strings)

            site = DisorderedTag(*site_string.split("."))
            coords = [np.fromstring(cs,dtype=float,sep=' ') for cs in coord_strings]
            untwist_name = f"{site}{altlocs}"
            if untwist_name not in untwist_candidates:
                untwist_candidates[untwist_name]=[]
            untwist_candidates[untwist_name].append((untwist_name,site,altlocs,coords))
    # the maximum number of untwists for a single atom site has.
    num_models_needed = max(len(cndts) for cndts in untwist_candidates.items())

    if os.path.isdir(get_candidate_model_dir()):
        shutil.rmtree(get_candidate_model_dir())
    os.makedirs(get_candidate_model_dir())
    # for untwist_move in untwist_candidates:
    #     untwist_name = untwist_move[0]
    #     out_path =  f"{get_candidate_model_dir()}{untwist_move[0]}.pdb"
    #     create_untwisted(model_path, untwist_move,out_path)

    all_applied_models=[]
    changes_only_models=[]
    for model_idx in range(num_models_needed):

        all_applied_out_path =  f"{get_candidate_model_dir()}all_untwists-{model_idx}.pdb"
        shutil.copy(model_path,all_applied_out_path)
        for untwist_moves in untwist_candidates.values():
            if model_idx < len(untwist_moves):
                move = untwist_moves[model_idx]
                create_untwisted(all_applied_out_path, move,all_applied_out_path,allow_overwrite=True)
        all_applied_models.append(all_applied_out_path)
        changes_only_path =  f"{get_candidate_model_dir()}untwists_isolated-{model_idx}.pdb"
        changed_entries_only(model_path,all_applied_out_path,changes_only_path)
        changes_only_models.append(changes_only_path)



    return get_candidate_model_dir(), all_applied_models, changes_only_models

if __name__=="__main__":
    apply_untwists(*sys.argv[1:])


