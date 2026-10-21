# Spencer Passmore 2026
#-----------------------

#from mmtbx.monomer_library import pdb_interpretation
from mmtbx.command_line import pdb_interpretation
import sys, os
import shutil

# Takes in model, writes .eff file excluding inter-unit-cell nonbond (Van der Waals) restraints.
# TODO if a water/residue has clash geomections with a residue/ligand in all conformations,
# should NOT exclude, since no arrangement of neighbors can accommodate it.

assert len(sys.argv[1:])==2, "Usage: phenix.python disable_intercell_nonbonds.py model.pdb path/to/make/exclusions.eff"

model, out_exclusions_file = sys.argv[1:]

assert out_exclusions_file.endswith(".eff")

model_handle = '.'.join(model.split('.')[:-1])


reuse_geo_file=True

geo_file=model_handle+"_geometry.geo"
args=[model,"--write_geo_files"]
print("Running phenix.pdb_interpretation "+' '.join(args))
if not reuse_geo_file or not os.path.exists(geo_file):
    pdb_interpretation.run(args=args)
    shutil.move(os.path.basename(model)+".geo",geo_file)
else:
    "Loaded previously generated file: "+geo_file

#print(os.path.realpath(geo_file))


min_overlap_to_exclude=0.6 # in Angstroms


def inter_unit_cell(sym_op:str):
    x,y,z=sym_op.split(',')
    
    for substr in [s for s in (x,y,z) if s[-1].isnumeric()]:
        # Example possibilities at this point: "x+1/2", "y-3/2", "x-2", "x-y",
        for translation in [t for t in substr.replace('+',' ').replace("-",' ').split() if t[-1].isnumeric()]:
            assert len(translation)>0, (x,y,z)
            # "1/2", "3/2", "2"
            if '/' not in translation:
                return True
            # "1/2", "3/2"
            numerator,denominator = [int(v) for v in translation.split('/')]
            if numerator>=denominator:
                return True

    return False



exclusions=[]
water_exclusions=[]
with open(geo_file,"r") as f:
    lines = f.readlines()
for i, line in enumerate(lines):
    model= vdw= None
    if line.startswith("nonbonded"):
        constraint = lines[i:i+4]
        pdb1=constraint[0].strip().split("\"")[1]
        pdb2=constraint[1].strip().split("\"")[1]
        quant_info=[v for v in constraint[3].strip().split()]
        if len(quant_info)==2:
            continue # Not a symmetry relation (intra-molecule)
        model, vdw, sym_op = quant_info
        try:
            if not inter_unit_cell(sym_op):
                continue
        except:
            assert False, "Unexpected symop string form: "+sym_op
        overlap=float(vdw)-float(model)
        if overlap<min_overlap_to_exclude:
            continue


        assert len(pdb1)==15,pdb1
        assert len(pdb2)==15,pdb2
        #" O  FGLY A 130 "
        exclusions.append(((pdb1,pdb2),overlap))
 

exclusions.sort(key=lambda x: -x[1])


# convert pdb strings to selections # TODO phenix.python probably has a method for this.
exclude_individual_atoms=False # if False, the entire residue is excluded from symmetry nonbond interactions #TODO when True group by residue and altloc label
exclusion_selections=[]
for pdb_strings,overlap in exclusions+water_exclusions:
  individual_parts=[]
  for pdb_string in pdb_strings:
    name = pdb_string[0:4].strip()
    altloc=pdb_string[4]
    resname=pdb_string[5:8].strip()
    resnum = pdb_string[-5:].strip()
    chain=pdb_string[-6]
    selection=("resname "+resname+" and resseq "+resnum+" and chain "+chain
                                + ((" and altloc "+altloc) if altloc!=" " else "")
                                + ((" and name"+name) if exclude_individual_atoms else ""))
                                #+'"'" # overlap: "+str(round(overlap,2))+" Ã…")
    individual_parts.append(selection)

  selection_union=('"('+") or (".join(individual_parts)+')'
# Hack to track the overlaps in the phenix.refine log.
                +  ' or name Ovrlp'+str(round(overlap,2)))+'"'
  if selection_union not in exclusion_selections:
    exclusion_selections.append(selection_union)

#mega_selection = ' or '.join(exclusion_selections)


prefix="  refinement.pdb_interpretation.custom_nonbonded_symmetry_exclusions = "
if len(exclusion_selections)>0:
    exclusion_selection_lines=prefix+ ('\n'+prefix).join(exclusion_selections)
else:
    exclusion_selection_lines="custom_nonbonded_symmetry_exclusions = None"

with open(out_exclusions_file,'w') as f:
    f.write(exclusion_selection_lines)
print(str(len(exclusion_selections))+" exclusions written to "+out_exclusions_file)



# phenix.python disable_intercell_nonbonds.py ../44PS.pdb exclusions_44PS.eff
# phenix.refine exclusions_44PS.eff ../44PS.pdb ../7KR0_S0053_refine_data.mtz overwrite=true