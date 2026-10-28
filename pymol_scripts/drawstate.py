#%%
import sys,os,pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent))
import pymol
import shutil
from pymol import cmd
from pymol_scripts.confmode import *
from pymol_scripts.highlightconnections import *
from pymol_scripts.loadBfacts import *
#from pymol_scripts.pairconformations import pairconformations
#from pymol_scripts.colorbadbonds import *
from UntangleFunctions import UNTANGLER_WORKING_DIRECTORY,model_handle
from Measures.evaluate_tangle import evaluate_tangle,get_tangled_bonds_path
from LinearOptimizer.RestraintsHandler import RestraintsHandler




DEFAULT_FIGURE_FOLDER=os.path.join(UNTANGLER_WORKING_DIRECTORY,"_figures","")

def drawstate(
    model_path:str,
    image_out:str,
    focus:str = None, # PyMOL selection
    altlocs:list[str]=None, # altlocs to show
    reference_ground_truth_pdb_for_bad_bonds:str=None,
    #color_by_diff_from:str=None,
    bad_bonds:str=None, # pdb file that indicates bad bonds with B factors. # Output of Measures/evaluate_tangle.py
    changed_bonds_file:str=None, # ChangedBonds file in linear optimizer log folder
    changed_bonds:str=None, # pymol selection string for bonds that changed.
    old_file_for_old_bonds:str=None, # pdb file for refrence to draw old bonds if have changed_bonds.
    atom_worst_geometry_devs:str=None, # pdb file where B factors correspond to the z value of the worst geometry. 
    atom_worst_xray_devs:str=None,  # pdb file where B factors correspond to things like B factor deviations, and whether it is supported by the maps.
    geom_devs_highlight_threshold:float=0.2,
    xray_devs_highlight_threshold:float=0.2,
    save_to_most_recent_image_path:bool=True,
    linewidth:float=1,
    auto_focus:bool=False,
    transparent_bg:bool=False,
        
):


    try:
        import pymol
    except:
        raise Exception("Could not import pymol. Please try running: pip install pymol-open-source")
    pymol.finish_launching(['pymol', '-cq'])
    cmd.delete("all")

    if not os.path.exists(DEFAULT_FIGURE_FOLDER):
        os.mkdir(DEFAULT_FIGURE_FOLDER)

    #assert (reference_ground_truth_pdb_for_bad_bonds,bad_bonds).count(None) >= 1

    mol=model_handle(model_path)
    cmd.load(model_path,mol)

    confmode()

    cmd.set("line_width",linewidth)



    # if color_by_diff_from is not None:
    #    difference_color_reference_model=model_handle(color_by_diff_from)
    #    cmd.load(color_by_diff_from,difference_color_reference_model)
    #    pairconformations(mol,difference_color_reference_model) 
    #    cmd.hide(f"not {mol}")

    if reference_ground_truth_pdb_for_bad_bonds is not None:
        assert bad_bonds is None
        #evaluate_tangle(model_path,reference_ground_truth_pdb_for_bad_bonds)
        bad_bonds=get_tangled_bonds_path(model_handle(model_path))
    if bad_bonds is not None:
        bad_bonds_sele=highlightconnections(bad_bonds,"bad",stick_radius=0.2,color="red")
        if focus is None and auto_focus:
            focus=bad_bonds_sele

    if changed_bonds_file is not None:
        assert changed_bonds is None
        changed_bonds=changed_bonds_file_to_pymol_sele(changed_bonds_file)

    if changed_bonds is not None:
        highlightconnections(model_path,"new",stick_radius=0.06,color="white",sele=changed_bonds)
        
    if old_file_for_old_bonds is not None and changed_bonds is not None:
        highlightconnections(old_file_for_old_bonds,"old",stick_radius=0.05,color="gray",sele=changed_bonds)
        
    if altlocs is not None:
        cmd.hide(f"(not altloc {'+'.join(altlocs)})")

    if atom_worst_geometry_devs is not None:
        print("Warning: Implementation of worst atom geometry devs not tested") # TODO
        loadBfacts(mol, atom_worst_geometry_devs)
        cmd.select(f"{mol} and b >= {geom_devs_highlight_threshold}")
        cmd.spectrum(b, "blue yellow red", "sele")
        cmd.hide("spheres", "mol and not sele")


    if focus is not None:
        zoom_buffer=5
    else:
        focus=mol
        zoom_buffer=-4
    cmd.center(focus)
    cmd.zoom(focus,buffer=zoom_buffer)
    
    #cmd.capture()
    #cmd.save(image_out)
    if transparent_bg:
        cmd.set("ray_opaque_background",0)
    else:
        cmd.set("ray_opaque_background",1)
    cmd.png(image_out,width=1440,height=1080) 
    if save_to_most_recent_image_path:
        most_recent_image_path=os.path.join(DEFAULT_FIGURE_FOLDER,"latest_capture.png")
        shutil.copy(image_out,most_recent_image_path)

    cmd.delete("all")


def changed_bonds_file_to_pymol_sele(changed_bonds:str):
    bond_selections=[]
    with open(changed_bonds) as f:
        for line in f:
            if line.startswith("Bond"):
                bond_str,altlocs=line.split(' ')
                atom_strings = bond_str.split('_')[1:]
                atom_seles=[]
                for atom_str in atom_strings:
                    resnum,name = atom_str.split('.')
                    atom_seles.append(f"{resnum}/{name}")
                bond_selections.append(
                    f"({' or '.join(atom_seles)}) and altloc {"+".join(altlocs)}" 
                )
    return " or ".join(bond_selections)

# if __name__=="__main__":
#     model_path,image_out=sys.argv[1:3]
#     if len(sys.argv)>3:
#         most_recent_image_path=sys.argv[3]
#     assert len(sys.argv) < 4
#     drawstate(model_path)
# %% Example
if __name__=="__main__":
    
    if len(sys.argv)>3:
        most_recent_image_path=sys.argv[3]
    assert len(sys.argv) < 4
    drawstate("drawstate_example/cheat_lrt_article_ogW_loopEnd1/lrt_article_ogW_loopEnd1_tangled_bonds.pdb",DEFAULT_FIGURE_FOLDER + "lrt_test.png",
              #reference_ground_truth_pdb_for_bad_bonds="../data/synthetic/best.pdb",
              #color_by_diff_from="../data/synthetic/best.pdb",
              bad_bonds="drawstate_example/cheat_lrt_article_ogW_loopEnd1/lrt_article_ogW_loopEnd1_tangled_bonds.pdb",
              #old_file_for_old_bonds="drawstate_example/cheat_longrangetraps/longrangetraps_tangled_bonds.pdb",
              changed_bonds_file="drawstate_example/cheat_lrt_article_ogW_loopEnd1/ChangedBondsFromInitial.txt",
              )
# %%
