#%%
import matplotlib.pyplot as plt 
import matplotlib.rcsetup as rcsetup
import os
from cycler import cycler



#TODO plot lines to points showing each potential Rfree at each step, 
# but which don't continue further (dead-ends) except for the one chosen
# and are fainter colour 

include_start=True
skip_missing=True
# handle="4PSS_4conf_control" #"4conf_sim" #"longrangetraps_TW"
# handle="4PSS_4conf" 
#handles = ["4PSS_4conf_control","4PSS_4conf","4PSS_4confH","4PSS_4confH_control"]
handles = ["4PSS_4confH","4PSS_4confH_control","4PSS_4confHV2","4PSS_4confH_controlV2"]
#handles = ["4PSS_4confHV2","4PSS_4confH_controlV2"]

RVALS=0
RWORK=1
RFREE=2
BONDS=3
ANGLES=4
CLASHES=5

wu_decrease_points=[]

last_num=24
ylim=[None,None]
#####

# TODO non sidechain bond/angle rmsz. 

for mode in (RWORK,RFREE,BONDS,ANGLES,CLASHES):
    fig, ax = plt.subplots(figsize=(6,6))
    if mode is RVALS:
        ax.set_prop_cycle(cycler(color=plt.cm.tab20.colors))
    #ax.set_prop_cycle(rcsetup.cycler('color', plt.cm.tab20.colors))
    def add_trace(handle):
        files=[]
        X=[]
        if include_start:
            files.append(f"output/{handle}_start.pdb")
            X.append(0)
        for n in range(0,last_num+1):
            #files.append(f"output/{handle}_loopEnd{n}.pdb")
            #files.append(f"output/{handle}_Accepted{n}.pdb")
            files.append(f"output/{handle}_Accepted{n}.pdb")
            X.append(n+1)
            if not os.path.exists(files[-1]):
                files.pop()
                X.pop()
                if not skip_missing:
                    break

        Rwork=[]; Rfree=[]; bond_rmsz=[]; angle_rmsz=[]; clash=[]; 
        #last_val = None
        for file in files:
            # if not os.path.exists(file):
            #     if last_val is not None:
            #         Rfree.append(last_val)
            #     continue
            with open(file) as f:
                for line in f:
                    for measure, identifying_str in zip(
                        (Rwork,Rfree,bond_rmsz,angle_rmsz,clash),
                        ("R VALUE            (WORKING SET) :",
                        "FREE R VALUE                     :",
                        "BOND      : ",
                        "ANGLE     :",
                        "ALL-ATOM CLASHSCORE :",
                        #"COVALENT GEOMETRY    : BOND",
                        #"COVALENT GEOMETRY    : ANGLE",
                        )
                    ):
                        if identifying_str in line:
                            last_val=float(line.split()[-1])
                            measure.append(last_val)              

        if mode is RVALS:
            ax.plot(X,Rfree,marker='o',label=handle)
            ax.plot(X,Rwork,linestyle="--",marker='o')
            ax.set_ylabel("Rwork/Rfree")
        elif mode is RWORK:
            ax.plot(X,Rwork,linestyle="--",marker='o',label=handle)
            ax.set_ylabel("$R_{\\rm{work}}$")
        elif mode is RFREE:
            ax.plot(X,Rfree,marker='o',label=handle)
            ax.set_ylabel("$R_{\\rm{free}}$ ")
        elif mode is BONDS:
            ax.plot(X,bond_rmsz,marker='o',label=handle)
            ax.set_ylabel("$RMSZ(\\rm{bonds})$")
        elif mode is ANGLES:
            ax.plot(X,angle_rmsz,marker='o',label=handle)
            ax.set_ylabel("$RMSZ(\\rm{angles})$")
        elif mode is CLASHES:
            ax.plot(X,clash,linestyle="--",marker='o',label=handle)
            ax.set_ylabel("Clash score")
        if mode in [BONDS,ANGLES]: 
            # Deal with much lower Z value for starting model 
            measure = bond_rmsz if mode is BONDS else angle_rmsz
            #ax.set_ylim(max(measure[1]*0.99,ax.get_ylim()[0]),max(max(measure)*1.001,ax.get_ylim()[1]))
        ax.set_xlabel("Loop")

    for handle in handles:
        add_trace(handle)

    ax.legend()
    plt.ylim(*ylim)
    ax.set_xlim(0,last_num+0.6)
    plt.grid()
    plt.show()
                

# %%
