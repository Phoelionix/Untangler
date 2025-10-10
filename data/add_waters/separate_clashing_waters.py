# Changes altloc

import shutil
import sys
import water_processing
fname=sys.argv[1]
fname_out=sys.argv[2]
tol_distance=1.6

#water_altloc="z"

with open(fname) as f:
    lines_out = f.readlines()
    water_processing.separate_clashing_waters(lines_out,tol_distance)

if fname_out == fname:
    shutil.move(fname,fname+"#")
with open(fname_out,"w+") as f:
    f.writelines(lines_out)