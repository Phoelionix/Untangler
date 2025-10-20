#%%
import os.path as path
import UntangleFunctions


# TODO need to group nearby ordered waters into same residue (i.e. make them the same disordered atom site!!).

pdb_file="refmacout_minRfree_current.pdb"
data_dir= UntangleFunctions.pdb_data_dir()

waters = []
with open(data_dir+pdb_file) as f:
    last_protein_res_num=0
    for line in f:
        resname = line[17:20] 
        if resname=="HOH":
            waters.append(line.strip("\n"))
        elif line.startswith("ATOM"):
            resnum=int(line[22:26].strip())
            last_protein_res_num=max(last_protein_res_num,resnum)  

for line in waters:
    old_res_num = int(line[22:26].strip())
    new_res_num = str(last_protein_res_num + old_res_num)
    new_res_num = " "*(4-len(new_res_num))+new_res_num
    new_line = line[:22] + new_res_num + line[26:]
    assert len(new_line) == len(line), len(new_res_num)
    print(new_line)
# %%

eff_constraints = """bond {
      action = *add
      atom_selection_1 = "name O and resseq 1 and chain S and altid A"
      atom_selection_2 = "name O and resseq 1 and chain S and altid B"
      distance_ideal = 0.4436
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 2 and chain S and altid A"
      atom_selection_2 = "name O and resseq 2 and chain S and altid B"
      distance_ideal = 0.3167
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 3 and chain S and altid A"
      atom_selection_2 = "name O and resseq 3 and chain S and altid B"
      distance_ideal = 0.4769
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 4 and chain S and altid A"
      atom_selection_2 = "name O and resseq 4 and chain S and altid B"
      distance_ideal = 0.3964
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 5 and chain S and altid A"
      atom_selection_2 = "name O and resseq 5 and chain S and altid B"
      distance_ideal = 0.5163
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 6 and chain S and altid A"
      atom_selection_2 = "name O and resseq 6 and chain S and altid B"
      distance_ideal = 0.7324
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 7 and chain S and altid A"
      atom_selection_2 = "name O and resseq 7 and chain S and altid B"
      distance_ideal = 0.7803
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 8 and chain S and altid A"
      atom_selection_2 = "name O and resseq 8 and chain S and altid B"
      distance_ideal = 1.3207
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 9 and chain S and altid A"
      atom_selection_2 = "name O and resseq 9 and chain S and altid B"
      distance_ideal = 0.5708
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 10 and chain S and altid A"
      atom_selection_2 = "name O and resseq 10 and chain S and altid B"
      distance_ideal = 0.6301
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 11 and chain S and altid A"
      atom_selection_2 = "name O and resseq 11 and chain S and altid B"
      distance_ideal = 0.7798
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 12 and chain S and altid A"
      atom_selection_2 = "name O and resseq 12 and chain S and altid B"
      distance_ideal = 0.7461
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 13 and chain S and altid A"
      atom_selection_2 = "name O and resseq 13 and chain S and altid B"
      distance_ideal = 0.7094
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 14 and chain S and altid A"
      atom_selection_2 = "name O and resseq 14 and chain S and altid B"
      distance_ideal = 1.1348
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 15 and chain S and altid A"
      atom_selection_2 = "name O and resseq 15 and chain S and altid B"
      distance_ideal = 0.8058
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 16 and chain S and altid A"
      atom_selection_2 = "name O and resseq 16 and chain S and altid B"
      distance_ideal = 0.7146
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 17 and chain S and altid A"
      atom_selection_2 = "name O and resseq 17 and chain S and altid B"
      distance_ideal = 0.6214
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 18 and chain S and altid A"
      atom_selection_2 = "name O and resseq 18 and chain S and altid B"
      distance_ideal = 0.7285
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 19 and chain S and altid A"
      atom_selection_2 = "name O and resseq 19 and chain S and altid B"
      distance_ideal = 1.5855
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 20 and chain S and altid A"
      atom_selection_2 = "name O and resseq 20 and chain S and altid B"
      distance_ideal = 0.8156
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 21 and chain S and altid A"
      atom_selection_2 = "name O and resseq 21 and chain S and altid B"
      distance_ideal = 0.7867
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 22 and chain S and altid A"
      atom_selection_2 = "name O and resseq 22 and chain S and altid B"
      distance_ideal = 0.8393
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 23 and chain S and altid A"
      atom_selection_2 = "name O and resseq 23 and chain S and altid B"
      distance_ideal = 0.8325
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 24 and chain S and altid A"
      atom_selection_2 = "name O and resseq 24 and chain S and altid B"
      distance_ideal = 0.8886
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 25 and chain S and altid A"
      atom_selection_2 = "name O and resseq 25 and chain S and altid B"
      distance_ideal = 0.7608
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 26 and chain S and altid A"
      atom_selection_2 = "name O and resseq 26 and chain S and altid B"
      distance_ideal = 0.6139
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 27 and chain S and altid A"
      atom_selection_2 = "name O and resseq 27 and chain S and altid B"
      distance_ideal = 0.7522
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 28 and chain S and altid A"
      atom_selection_2 = "name O and resseq 28 and chain S and altid B"
      distance_ideal = 0.9834
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 29 and chain S and altid A"
      atom_selection_2 = "name O and resseq 29 and chain S and altid B"
      distance_ideal = 0.7141
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 30 and chain S and altid A"
      atom_selection_2 = "name O and resseq 30 and chain S and altid B"
      distance_ideal = 0.9085
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 31 and chain S and altid A"
      atom_selection_2 = "name O and resseq 31 and chain S and altid B"
      distance_ideal = 0.8265
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 32 and chain S and altid A"
      atom_selection_2 = "name O and resseq 32 and chain S and altid B"
      distance_ideal = 2.2270
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 33 and chain S and altid A"
      atom_selection_2 = "name O and resseq 33 and chain S and altid B"
      distance_ideal = 1.2347
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 34 and chain S and altid A"
      atom_selection_2 = "name O and resseq 34 and chain S and altid B"
      distance_ideal = 0.8067
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 35 and chain S and altid A"
      atom_selection_2 = "name O and resseq 35 and chain S and altid B"
      distance_ideal = 0.9120
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 36 and chain S and altid A"
      atom_selection_2 = "name O and resseq 36 and chain S and altid B"
      distance_ideal = 6.9923
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 37 and chain S and altid A"
      atom_selection_2 = "name O and resseq 37 and chain S and altid B"
      distance_ideal = 1.1617
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 38 and chain S and altid A"
      atom_selection_2 = "name O and resseq 38 and chain S and altid B"
      distance_ideal = 1.0718
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 39 and chain S and altid A"
      atom_selection_2 = "name O and resseq 39 and chain S and altid B"
      distance_ideal = 1.1046
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 40 and chain S and altid A"
      atom_selection_2 = "name O and resseq 40 and chain S and altid B"
      distance_ideal = 1.1078
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 41 and chain S and altid A"
      atom_selection_2 = "name O and resseq 41 and chain S and altid B"
      distance_ideal = 0.9540
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 42 and chain S and altid A"
      atom_selection_2 = "name O and resseq 42 and chain S and altid B"
      distance_ideal = 0.8325
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 43 and chain S and altid A"
      atom_selection_2 = "name O and resseq 43 and chain S and altid B"
      distance_ideal = 1.0343
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 44 and chain S and altid A"
      atom_selection_2 = "name O and resseq 44 and chain S and altid B"
      distance_ideal = 0.9965
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 45 and chain S and altid A"
      atom_selection_2 = "name O and resseq 45 and chain S and altid B"
      distance_ideal = 1.0373
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 46 and chain S and altid A"
      atom_selection_2 = "name O and resseq 46 and chain S and altid B"
      distance_ideal = 1.0454
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 47 and chain S and altid A"
      atom_selection_2 = "name O and resseq 47 and chain S and altid B"
      distance_ideal = 1.2136
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 48 and chain S and altid A"
      atom_selection_2 = "name O and resseq 48 and chain S and altid B"
      distance_ideal = 1.3593
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 49 and chain S and altid A"
      atom_selection_2 = "name O and resseq 49 and chain S and altid B"
      distance_ideal = 1.3604
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 50 and chain S and altid A"
      atom_selection_2 = "name O and resseq 50 and chain S and altid B"
      distance_ideal = 1.1944
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 51 and chain S and altid A"
      atom_selection_2 = "name O and resseq 51 and chain S and altid B"
      distance_ideal = 1.0083
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 52 and chain S and altid A"
      atom_selection_2 = "name O and resseq 52 and chain S and altid B"
      distance_ideal = 1.1515
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 53 and chain S and altid A"
      atom_selection_2 = "name O and resseq 53 and chain S and altid B"
      distance_ideal = 1.0003
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 54 and chain S and altid A"
      atom_selection_2 = "name O and resseq 54 and chain S and altid B"
      distance_ideal = 0.8963
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 55 and chain S and altid A"
      atom_selection_2 = "name O and resseq 55 and chain S and altid B"
      distance_ideal = 1.0567
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 56 and chain S and altid A"
      atom_selection_2 = "name O and resseq 56 and chain S and altid B"
      distance_ideal = 1.0341
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 57 and chain S and altid A"
      atom_selection_2 = "name O and resseq 57 and chain S and altid B"
      distance_ideal = 1.4880
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 58 and chain S and altid A"
      atom_selection_2 = "name O and resseq 58 and chain S and altid B"
      distance_ideal = 2.0742
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 59 and chain S and altid A"
      atom_selection_2 = "name O and resseq 59 and chain S and altid B"
      distance_ideal = 1.0331
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 60 and chain S and altid A"
      atom_selection_2 = "name O and resseq 60 and chain S and altid B"
      distance_ideal = 1.2260
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 61 and chain S and altid A"
      atom_selection_2 = "name O and resseq 61 and chain S and altid B"
      distance_ideal = 1.1611
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 62 and chain S and altid A"
      atom_selection_2 = "name O and resseq 62 and chain S and altid B"
      distance_ideal = 1.1641
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 63 and chain S and altid A"
      atom_selection_2 = "name O and resseq 63 and chain S and altid B"
      distance_ideal = 0.9210
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 64 and chain S and altid A"
      atom_selection_2 = "name O and resseq 64 and chain S and altid B"
      distance_ideal = 1.0943
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 65 and chain S and altid A"
      atom_selection_2 = "name O and resseq 65 and chain S and altid B"
      distance_ideal = 0.8738
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 66 and chain S and altid A"
      atom_selection_2 = "name O and resseq 66 and chain S and altid B"
      distance_ideal = 1.7408
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 67 and chain S and altid A"
      atom_selection_2 = "name O and resseq 67 and chain S and altid B"
      distance_ideal = 1.1932
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 68 and chain S and altid A"
      atom_selection_2 = "name O and resseq 68 and chain S and altid B"
      distance_ideal = 1.3023
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 69 and chain S and altid A"
      atom_selection_2 = "name O and resseq 69 and chain S and altid B"
      distance_ideal = 1.0004
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 70 and chain S and altid A"
      atom_selection_2 = "name O and resseq 70 and chain S and altid B"
      distance_ideal = 1.2823
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 71 and chain S and altid A"
      atom_selection_2 = "name O and resseq 71 and chain S and altid B"
      distance_ideal = 1.1208
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 72 and chain S and altid A"
      atom_selection_2 = "name O and resseq 72 and chain S and altid B"
      distance_ideal = 1.1993
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 73 and chain S and altid A"
      atom_selection_2 = "name O and resseq 73 and chain S and altid B"
      distance_ideal = 1.0859
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 74 and chain S and altid A"
      atom_selection_2 = "name O and resseq 74 and chain S and altid B"
      distance_ideal = 1.0914
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 75 and chain S and altid A"
      atom_selection_2 = "name O and resseq 75 and chain S and altid B"
      distance_ideal = 1.3890
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 76 and chain S and altid A"
      atom_selection_2 = "name O and resseq 76 and chain S and altid B"
      distance_ideal = 1.1576
      sigma = 0.2
    }
    bond {
      action = *add
      atom_selection_1 = "name O and resseq 77 and chain S and altid A"
      atom_selection_2 = "name O and resseq 77 and chain S and altid B"
      distance_ideal = 1.9518
      sigma = 0.2
    }"""

for line in eff_constraints.split("\n"):
    if not "atom_selection" in line:
        print(line)
        continue 
    old_res_num = int(line.split()[6])
    new_res_num = str(64 + old_res_num)
    new_line = " "*6+ ' '.join(line.split()[:6]) + " " + new_res_num + " " + ' '.join(line.split()[7:])
    print(new_line)
# %%
