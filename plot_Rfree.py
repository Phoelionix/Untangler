#%%
import matplotlib.pyplot as plt 
import os


#TODO plot lines to points showing each potential Rfree at each step, 
# but which don't continue further (dead-ends) except for the one chosen
# and are fainter colour 

handle="4conf_sim" #"longrangetraps_TW"
last_num=30 #11


files = []
for n in range(0,last_num+1):
    #files.append(f"output/{handle}_loopEnd{n}.pdb")
    #files.append(f"output/{handle}_Accepted{n}.pdb")
    files.append(f"output/{handle}_current{n}.pdb")

Rfree=[]
#last_val = None
for file in files:
    # if not os.path.exists(file):
    #     if last_val is not None:
    #         Rfree.append(last_val)
    #     continue
    with open(file) as f:
        for line in f:
            if line.startswith("REMARK   3   FREE R VALUE                     :"):
                last_val=float(line.split()[-1])
                Rfree.append(last_val)
                break
        else:
            raise Exception("Rfree line missing")

plt.plot(range(len(Rfree)),Rfree,marker='o')
#plt.ylim(0.09,None)
plt.grid()
plt.show()
                

# %%
