#%%
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker
import matplotlib.rcsetup as rcsetup
import numpy as np
import os
from cycler import cycler
from UntangleFunctions import assess_geometry_wE,score_file_name,get_score
import UntangleFunctions

#TODO plot lines to points showing each potential Rfree at each step, 
# but which don't continue further (dead-ends) except for the one chosen
# and are fainter colour 

include_start=True
skip_missing=True
# handle="4PSS_4conf_control" #"4conf_sim" #"longrangetraps_TW"
# handle="4PSS_4conf" 



## Angle 2x weight, refine positions until worse:
# V2 : main chain geometries weight x 10  
# V3 : Forbid side-chain-only geometry changes  
## Refine positions once:
# V4 : challenge sol weights. 
# V5 : Challenge sol weights, challenge score (chi z sqr)

RVALS=0
RWORK=1
RFREE=2
BONDS=3
ANGLES=4
CLASHES=5
WE=6

wu_decrease_points=[]

last_num=31# 20  #25
ylim=[None,None]
#####

# TODO non sidechain bond/angle rmsz. 


UntangleFunctions.TEMP_SCORE_WITH_FIRST_PROTEIN_ALTLOC_ONLY=False # Ensure we get wE score for entire structure and not a single conformation.

#include_unrestrained=mode in [RWORK,RFREE]
include_unrestrained=False
def get_data(handle,skip_wE=False):
    files=[]
    X=[]
    disjoint_x_vals = []
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
        elif include_unrestrained:
            files.insert(-1,f"output/{handle}_unrestrained{n}.pdb")
            x=n+1/2
            X.insert(-1,x)
            disjointed=False
            if disjointed:
                disjoint_x_vals.append(x)

    Rwork=[]; Rfree=[]; wE=[]; bond_rmsz=[]; angle_rmsz=[]; clash=[]; 
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
        if not skip_wE:
            if os.path.exists(score_file_name(file)):
                wE.append(get_score(score_file_name(file))[1])
            else:
                wE.append(assess_geometry_wE(file)[1])            
    
    y_data= {
        RWORK:[100* v for v in Rwork],
        RFREE:[100* v for v in Rfree],
        WE:wE if not skip_wE else None,
        BONDS:bond_rmsz,
        ANGLES:angle_rmsz,
        CLASHES:clash,
    }
    return X, y_data, disjoint_x_vals,  

def add_trace(handle,data_dict,measure_mode,label=None):
    label = handle if label is None else label

    X,measure,disjoint_x_vals = data_dict[handle]
    measure=measure[measure_mode]

    if measure_mode is RVALS:
        assert False
        ax.plot(X,Rwork,linestyle="--",marker='o')
        ax.plot(X,Rfree,marker='o',label=label)
        ax.set_ylabel("Rwork/Rfree")
    elif measure_mode is RWORK:
        ax.set_ylabel("$R_{\\rm{work}}$")
    elif measure_mode is RFREE:
        ax.set_ylabel("$R_{\\rm{free}}$ ")
    elif measure_mode is WE:
        ax.set_ylabel("$\\text{wE}$")
    elif measure_mode is BONDS:
        ax.set_ylabel("Bond length RMSZ")
    elif measure_mode is ANGLES:
        ax.set_ylabel("Angle RMSZ")
    elif measure_mode is CLASHES:
        ax.set_ylabel("Clash score")
        ax.set_ylim(0,max(max(measure)*1.01,ax.get_ylim()[1]))

    
    X_joined =[x for x in X if x not in disjoint_x_vals] 
    Y_joined = [v for x, v in zip(X,measure) if x not in disjoint_x_vals]
    ax.plot(X_joined,Y_joined, marker='o',label=label)
    ax.scatter(disjoint_x_vals,[v for x, v in zip(X,measure) if x in disjoint_x_vals])
    #ax.plot(X,measure,marker='o',label=handle)


    # if mode in [BONDS,ANGLES]: 
    #     # Deal with much lower Z value for starting model 
    #     if len(measure) > 1:
    #         ax.set_ylim(min(Y_joined[1:])*0.99,max(max(Y_joined)*1.001,ax.get_ylim()[1]))
            
    ax.set_xlabel("Loop")
    #ax.set_xlabel("Loop (10 restrained macrocycles ea.)")



def rollavg_cumsum_edges(a, n):
    #https://stackoverflow.com/questions/14313510/how-can-i-calculate-a-rolling-moving-average-using-python-numpy-scipy
    'numpy.cumsum, edge handling'
    assert n%2 == 1
    N = len(a)
    cumsum_vec = np.cumsum(np.insert(np.pad(a, (n-1, n-1), 'constant'), 0, 0))
    d = np.hstack((np.arange(n//2 + 1, n), np.ones(N-n)*n, np.arange(n, n//2, -1)))
    return (cumsum_vec[n + n//2:-n//2 + 1] - cumsum_vec[n//2:-n - n//2]) / d

def add_trace_with_error(handles,data_dict,measure_mode,label,col=None,moving_average_window_width=None):
    if len(handles)==0:
        print(f"Empty handles for {label}")
        return
    series=[]
    for handle in handles:
        X,measure,disjoint_x_vals = data_dict[handle]
        series.append(measure[measure_mode][:last_num+1])
    if moving_average_window_width is not None:
            # w=moving_average_window_width
            # cumsum_vec = np.cumsum(np.insert(series, 0, 0)) 
            # ma_vec = (cumsum_vec[w:] - cumsum_vec[:-w]) / w
            for i in range(len(series)):
                series[i] = rollavg_cumsum_edges(series[i],moving_average_window_width)                     
    X=X[:last_num+1]
    mean = np.mean(series,axis=0)
    stdv = np.std(series,axis=0)
    



    #with sns.axes_style("darkgrid"):  
    if measure_mode is RVALS:
        assert False
        ax.plot(X,Rwork,linestyle="--",marker='o')
        ax.plot(X,Rfree,marker='o',label=label)
        ax.set_ylabel("Rwork/Rfree")
    elif measure_mode is RWORK:
        ax.set_ylabel("$R_{\\rm{work}}$ (%)")
        ax.set_ylim(11.25,13.1)
        #ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
    elif measure_mode is RFREE:
        ax.set_ylabel("$R_{\\rm{free}}$ (%)")
        ax.set_ylim(12.5,14.1)
        #ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
    elif measure_mode is WE:
        ax.set_ylabel("$\\text{wE}$")
    elif measure_mode is BONDS:
        ax.set_ylabel("Bond length RMSZ")
    elif measure_mode is ANGLES:
        ax.set_ylabel("Angle RMSZ")
    elif measure_mode is CLASHES:
        ax.set_ylabel("Clash score")
        ax.set_ylim(0,max(max(mean)*1.01,ax.get_ylim()[1]))

    
    X_joined =[x for x in X if x not in disjoint_x_vals] 
    Y_joined = [v for x, v in zip(X,mean) if x not in disjoint_x_vals]
    #ax.plot(X_joined,Y_joined, marker='o',label=label,c=col)
    ax.plot(X_joined,Y_joined, marker='o',label=label,c=col,markersize=0)
    ax.scatter(disjoint_x_vals,[v for (x, v) in zip(X,mean) if x in disjoint_x_vals],c=col)
    ax.fill_between(X, mean-stdv, mean+stdv,alpha=0.3,color=col)# facecolor=clrs[i])


    #ax.plot(X,measure,marker='o',label=handle)


    # if mode in [BONDS,ANGLES]: 
    #     # Deal with much lower Z value for starting model 
    #     if len(measure) > 1:
    #         ax.set_ylim(min(Y_joined[1:])*0.99,max(max(Y_joined)*1.001,ax.get_ylim()[1]))
            
    ax.set_xlabel("Loop")
    #ax.set_xlabel("Loop (10 restrained macrocycles ea.)")



# %%



moving_average_window_width=3
# handle_sets = [
#     #["4PSS_4confHV4","altlocoptimizer3","altlocoptimizer4","altlocoptimizer5","altlocoptimizer6","altlocoptimizer7","altlocoptimizer8"],
#     ["4PSS_4confHV4","altlocoptimizer3","altlocoptimizer4","altlocoptimizer5","altlocoptimizer6"],
#     ["4PSS_4confH_controlV4","control1","control2","control3","control4",],
#     ["4PSS_4confH_controlRestrained","controlRestrained1","controlRestrained2","controlRestrained3","controlRestrained4",]
# ]
handle_sets = [
    #["4PSS_4confHV4","altlocoptimizer3","altlocoptimizer4","altlocoptimizer5","altlocoptimizer6","altlocoptimizer7","altlocoptimizer8"],
    #["FrzS_3conf_coot1",],
    ["FrzS_3conf_57To60Fix"],
]

skip_wE=True
data_dict={}
for handle_set in handle_sets:
    for handle in handle_set:
        data_dict[handle] = get_data(handle,skip_wE=skip_wE)


colors = [plt.cm.tab10.colors[j] for j in (1,2,0)]
for mode,measure_name in zip((RWORK,RFREE,WE,BONDS,ANGLES,CLASHES), 
                             ("Rwork","Rfree","wE","bonds","angles","clash_score")):
    if skip_wE and mode == WE:
        continue    
    fig, ax = plt.subplots(figsize=(2.3,2.5),)
    ax.set_prop_cycle(cycler(color=plt.cm.tab10.colors[0:]))
    if mode is RVALS:
        ax.set_prop_cycle(cycler(color=plt.cm.tab20.colors))
    #ax.set_prop_cycle(rcsetup.cycler('color', plt.cm.tab20.colors))
    #labels=["With conn. optim.", "No conn. optim.", "Restrained only"]
    labels=["FrzS"]
    # labels=["1 Unrestr. $\\rightarrow$ altloc optim. $\\rightarrow$ 10 restr.", 
    #         "1 Unrestr.  $\\rightarrow$ 10 restr.", 
    #         "10 restr."]
    # labels=["Unrestr., altloc optim., restr.",
    #         #"With altloc optim." 
    #         "Unrestr., restr.", 
    #         "Restr."]
    assert len(labels)==len(handle_sets)
    for i, label in zip(range(len(labels)), labels):
        add_trace_with_error(handle_sets[i],data_dict,mode,labels[i],col=colors[i],
                             moving_average_window_width=moving_average_window_width)


    
    ax.legend(handlelength=1,framealpha=0.7,labelspacing=0.2)
    #ax.legend(handlelength=1,framealpha=0.7,labelspacing=0.2,title="Steps each loop")
    plt.ylim(*ylim)
    #ax.set_xlim(0,last_num+0.6)
    ax.set_xlim(0,last_num)
    os.makedirs("measures_by_cycle",exist_ok=True)
    #plt.grid()
    fig.savefig(os.path.join("measures_by_cycle",f"{measure_name}.png"),bbox_inches='tight',dpi=800)#,pad_inches=0)
    plt.show()

#ax.set_ylim(0.12,0.14)
#plt.ylim()






# #%% OLD







# #handles = ["4PSS_4conf","4PSS_4conf_control"]; include_start=False; 
# #handles = ["4PSS_4confH","4PSS_4confH_control"]
# #handles = ["4PSS_4confH","4PSS_4confH_control","4PSS_4confHV2","4PSS_4confH_controlV2"]
# #handles = ["4PSS_4confH_controlRestrained","4PSS_4confHV2","4PSS_4confH_controlV2","4PSS_4confHV3","4PSS_4confHV4","4PSS_4confHV5"]
# #handles = ["4PSS_4confH_controlRestrained","4PSS_4confH_controlV2","4PSS_4confHV3"]



# #handles = ["4PSS_4confHV4","4PSS_4confH_controlRestrained","4PSS_4confH_controlV4"]
# #handles = ["altlocoptimizer2","control1","controlRestrained1"]
# #handles = ["control2","control3","control4"]
# #handles = ["controlRestrained2","controlRestrained3","controlRestrained4"]
# #handles = ["altlocoptimizer4","altlocoptimizer5","altlocoptimizer6"]
# handles = ["FrzS_3conf_coot1"]

# data_dict={}
# for handle in handles:
#     data_dict[handle]=get_data(handle,skip_wE=skip_wE)

# for mode,measure_name in zip((RWORK,RFREE,WE,BONDS,ANGLES,CLASHES), 
#                              ("Rwork","Rfree","wE","bonds","angles","clash_score")):   
#     fig, ax = plt.subplots(figsize=(3.3,3.3*9/7))
#     ax.set_prop_cycle(cycler(color=plt.cm.tab10.colors[0:]))
#     if mode is RVALS:
#         ax.set_prop_cycle(cycler(color=plt.cm.tab20.colors))
#     #ax.set_prop_cycle(rcsetup.cycler('color', plt.cm.tab20.colors))
#     #labels=["With altloc optimization", "Without altloc optimization", "Ordinary refinement"]
#     labels=["1 Unrestr. $\\rightarrow$ altloc optim. $\\rightarrow$ 10 restr.", 
#             "1 Unrestr.  $\\rightarrow$ 10 restr.", 
#             "10 restr."]
#     for handle, label in zip(handles, labels):
#         add_trace(handle,data_dict,measure_mode=mode,label=label)

    
#     ax.legend(title="Loop stages")
#     plt.ylim(*ylim)
#     ax.set_xlim(0,last_num+0.6)
#     os.makedirs("measures_by_cycle",exist_ok=True)
#     plt.grid()
#     fig.savefig(os.path.join("measures_by_cycle",f"{measure_name}.png"))
#     plt.show()
# # %%

# %%
