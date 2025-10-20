#%%
import matplotlib.pyplot as plt 


#TODO plot lines to points showing each potential Rfree at each step, 
# but which don't continue further (dead-ends) except for the one chosen
# and are fainter colour 

handle="5conf"
last_num=9


files = []
for n in range(0,last_num+1):
    files.append(f"output/{handle}_loopEnd{n}.pdb")

Rfree=[]
for file in files:
    with open(file) as f:
        for line in f:
            if line.startswith("REMARK   3   FREE R VALUE                     :"):
                Rfree.append(float(line.split()[-1]))
                break
        else:
            raise Exception("Rfree line missing")

plt.plot(range(len(Rfree)),Rfree,marker='o')
#plt.ylim(0.09,None)
plt.grid()
plt.show()
                

# %%
