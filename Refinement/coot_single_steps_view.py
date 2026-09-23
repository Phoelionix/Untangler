import sys
print(sys.version)
import coot
import os
from time import sleep
import shutil


batch="cheat"
#batch="mixedup"
refinement_folders="Refinement/MovieMaking/"+batch+"/"
out_dir="Refinement/MovieMaking/frames_"+batch+"/"


def folder_order_by_end_number(string):
    # By number at end
    assert '.' not in string, string
    last_non_numeric=None
    for i, c in enumerate(string):
        if not c.isdigit():
            last_non_numeric=i

    if last_non_numeric==i: 
        return 0 # does not end with number
    return int(string[last_non_numeric+1:])

pdb_files=[]
for folder in sorted(os.listdir(refinement_folders),key=folder_order_by_end_number):
    folder_path=os.path.join(refinement_folders,folder,"",)
    for file in os.listdir(folder_path):
        if file.endswith("_999.pdb"):
            pdb_files.append(os.path.join(folder_path,file))
            break
    else:
        assert False, "Could not find pdb file in "+folder_path




# Animation
if os.path.exists(out_dir):
    shutil.rmtree(out_dir)
os.mkdir(out_dir)

sigma_2fofc=1.5
#sigma_fofc=5
sigma_fofc=5
coot.set_default_initial_contour_level_for_map(sigma_2fofc)
coot.set_default_initial_contour_level_for_difference_map(sigma_fofc)

frame_duration=0.5
rots=10
assert len(pdb_files)<=50
imol_pdb=imol_map=imol_diff_map=None
frame=0
for pdb_file in pdb_files:
    handle=pdb_file[:-4]
    print(handle)
    if imol_pdb is None:
        imol_pdb = coot.read_pdb(pdb_file)
    else:
        coot.clear_and_update_model_molecule_from_file(imol_pdb,pdb_file)
    old_imol_map=imol_map
    old_imol_diff_map=imol_diff_map
    imol_map = coot.read_mtz(handle+".mtz", "2FOFCWT", "PH2FOFCWT", "2FOFCWT", False, False)
    imol_diff_map = coot.read_mtz(handle+".mtz", "FOFCWT", "PHFOFCWT", "FOFCWT", False, True)

    if old_imol_map is not None:
        coot.close_molecule(old_imol_map)
    if old_imol_diff_map is not None:
        coot.close_molecule(old_imol_diff_map)


    coot.set_contour_level_in_sigma(imol_map, sigma_2fofc)
    coot.set_contour_level_in_sigma(imol_diff_map,sigma_fofc)

    rotate=True
    #deg_rot=0.5
    deg_rot=0.15
    for _ in range(rots):
        coot.screendump_image(os.path.join(out_dir,str(frame)+".ppm"))
        if rotate:
            frame+=1
            coot.rotate_y_scene(1,deg_rot)
            sleep(frame_duration/rots)

# Reset to original view
if rotate:
    for _ in range(frame):
        coot.rotate_y_scene(1,-deg_rot)


# ffmpeg -framerate 10 -i Refinement/MovieMaking/frames_cheat/%d.ppm -c:v libx264 -crf 25 -vf "scale=2490:1276,format=yuv420p" -movflags +faststart Refinement/MovieMaking/cheat.mp4
# ffmpeg -framerate 10 -i Refinement/MovieMaking/frames_cheat/%d.ppm -c:v libx264 -crf 25 -vf "scale=4982:2552,format=yuv420p" -movflags +faststart Refinement/MovieMaking/cheat.mp4
# ffmpeg -framerate 10 -i Refinement/MovieMaking/frames_cheat/%d.ppm -c:v libx264 -crf 25 -filter:v "crop=iw-1:ih" -pix_fmt yuv420p -movflags +faststart Refinement/MovieMaking/cheat.mp4
