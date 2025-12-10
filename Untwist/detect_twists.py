#%%
import os, sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent))
from Bio.PDB.Structure import Structure#, parse_pdb_header
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom,DisorderedAtom
# from typing import List # if your version of python is older than 3.9 uncomment this and replace instances of "list[" with "List[" 
import os
import subprocess
import matplotlib
from time import sleep
import numpy as np
from numpy.typing import NDArray
#matplotlib.use('AGG')
import shutil
import itertools
from scipy.linalg import norm
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
from LinearOptimizer.ConstraintsHandler import ConstraintsHandler,DisorderedTag
from LinearOptimizer.OrderedAtomLookup import OrderedAtomLookup
import UntangleFunctions


#(I) they share the same midpoint
# (II) they are farther apart
# (III) each conformer position is related via rotation of a model conformer around the vector formed by a and b.
# 
# 
# 
# 

def get_angle(a:Atom,b:Atom,c:Atom):
    v1 = a.get_coord() - b.get_coord()
    v2 = c.get_coord()-b.get_coord()
    def unit_vector(vector):
        return vector / np.linalg.norm(vector)
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.dot(v1_u, v2_u))*180/np.pi

def separation(v1,v2):
    if type(v1)==Atom and type(v2)==Atom:
        v1, v2 = v1.get_coord(),v2.get_coord()
    assert np.array(v1).shape==np.array(v2).shape==(3,), (np.array(v1).shape, np.array(v2).shape)
    return np.sqrt(np.sum((v1-v2)**2))

# def atom_sep(a:Atom,b:Atom):
#     return separation(a.get_coord(),b.get_coord())




def get_angle_atoms()->list[tuple[Atom,Atom,Atom]]:
    assert False
def get_ideal_angle(constraints_handler:ConstraintsHandler,atoms:tuple[Atom,Atom,Atom])->float:

    site_tags = [DisorderedTag.from_atom(a) for a in atoms]

    valid_constraints = [constraint for constraint in constraints_handler.atom_constraints[site_tags[0]] \
                         if (type(constraint)==ConstraintsHandler.AngleConstraint
                         and all([
                             (constraint in constraints_handler.atom_constraints[site_tag]) for site_tag in site_tags[1:]
                             ]))]
    assert len(valid_constraints)==1
    return valid_constraints[0].ideal
def get_ideal_bond_length(constraints_handler:ConstraintsHandler,atoms:tuple[Atom,Atom])->float:
    site_tags = [DisorderedTag.from_atom(a) for a in atoms]

    valid_constraints = [constraint for constraint in constraints_handler.atom_constraints[site_tags[0]] \
                         if (type(constraint)==ConstraintsHandler.BondConstraint
                         and all([
                             (constraint in constraints_handler.atom_constraints[site_tags[i+1]]) for i in range(1)
                             ]))]
    assert len(valid_constraints)==1
    return valid_constraints[0].ideal

def get_constraints_handler(pdb_file_path,geo_file_needs_generation=True):
    struct = PDBParser().get_structure("struct",pdb_file_path)
    if geo_file_needs_generation:
        UntangleFunctions.create_score_file(pdb_file_path,turn_off_cdl=True)
    ordered_atom_lookup = OrderedAtomLookup(struct.get_atoms())
    constraints_handler=ConstraintsHandler()
    constraints_handler.load_all_constraints(pdb_file_path,ordered_atom_lookup,symmetries=None,calc_nonbonds=False,turn_off_cdl=True)
    return constraints_handler

# Get the circle of possible values for position x that satisfy an angle between start_anchor--mid_anchor--x and distance for mid_anchor--x
def get_arc(start_anchor:tuple[float],mid_anchor:tuple[float],bond_length,atom_angle,anchor_for_direction_to_arc_midpoint:NDArray,num_points=50,arc_angle=25):
    # direction_to_arc_midpoint, should just correspond to the current conformer coords.

    #https://stackoverflow.com/questions/48703275/3d-truncated-cone-in-python

    # phi = cone angle
    # c = slant_height

    phi = (180-atom_angle)*np.pi/180
    arc_angle = arc_angle*np.pi/180
    c=bond_length

    v = mid_anchor-start_anchor 
    v =  v/norm(v) # unit vector along cone axis.
    
    h = c*np.cos(phi) # cone height
    r = c*np.sin(phi) # cone radius

    
    
    # now get perpendicular vectors
    # make some vector not in the same direction as v
    not_v = np.array([1,1,0])
    if (v == not_v).all():
        not_v = np.array([0, 1, 0])
    # make vector perpendicular to v
    x = np.cross(v, not_v)
    # print n1,'\t',norm(n1)
    # normalize n1
    x /= norm(x)
    # make unit vector perpendicular to v and n1
    y = np.cross(v, x) 

    # Get angle for centre of arc
    D=anchor_for_direction_to_arc_midpoint-mid_anchor-h*v # point relative to circle centre
    D = D-np.dot(D,v)*v
    D/=norm(D)
    theta_0 = np.arctan2(np.dot(D,y),np.dot(D,x))

    theta = np.linspace(theta_0-arc_angle/2, theta_0+arc_angle/2, num_points)

    points = mid_anchor+h*v + r*(np.cos(theta[...,None])*x + np.sin(theta[...,None])*y)
    return points


    
def separation_condition(atom_pair,max_separation):
    A,B = [a for a in atom_pair]
    assert type(A) == type(B)==Atom, (A,B)
    return separation(A,B)<=max_separation

def relative_orientation(v1,v2):
    v1/=norm(v1)
    v2/=norm(v2)
    return np.arctan2(norm(np.cross(v1, v2)), np.dot(v1, v2))*180/np.pi
    

def detect_twist(
    atoms:tuple[Atom,Atom,Atom],
    constraints_handler=None,
    snap_to_ideal=True,  # Snap to ideal angle and bond length. Otherwise, uses current (potentially warped angle and bond length)
    snap_interp=0.5, # for angles and bonds: (1-snap_interp)*current + snap_interp*ideal
    plot_zoomed_in=False,
    plot_zoomed_out=False,
    plot_rot_speed=1,
    min_twist_angle=25,
    max_fake_separation=0.25, # TODO make it resolution dependent
    max_true_separation=0.25,
    midpoint_match_tol_frac=0.1,
    min_ratio_real_sep_on_fake_sep=1.1, # note that at 1, requires separation of real conformers are at least as large
    max_difference_real_fake_sep=0.1, # TODO should be replaced with a check for gaussian overlap with tolerance determined by the resolution
    take_average = False, # Takes average of all twist point pairs that satisfy the conditions
    take_closest = False, # Takes the twist point pair that satisfies the conditions that has a separation most similar to the current atoms.
    verbose=False,
):
    
    if snap_interp == 0:
        snap_to_ideal=False # NOTE
    if snap_to_ideal:
        assert 0 < snap_interp <=1
        assert constraints_handler is not None 

    assert not (take_average and take_closest)

    twist_points=[]
    
    A,B,C=atoms
    # Check C is close enough. 
    if not separation_condition(C,max_fake_separation):
        return twist_points # empty list
    
    ordered_atoms:dict[str,list[Atom]]={}
    for i, disordered_atom in enumerate(atoms):
        for atom in disordered_atom:
            if atom.get_altloc() not in ordered_atoms:
                ordered_atoms[atom.get_altloc()]=[None,None,None]
            ordered_atoms[atom.get_altloc()][i]=atom

    c_coords=[c.get_coord() for c in C]
    midpoint = np.mean(c_coords,axis=0)
    fake_separation=separation(*c_coords)

    # a✓ -- b✓ -- c? 
    #arcs = {altloc:[] for altloc in ordered_atoms.keys()}
    arcs={}
    arcs_to_plot=[]
    bond_lengths,angles = [],[]
    for altloc, (a,b,c) in ordered_atoms.items():
        bond_length,angle=separation(b,c),get_angle(a,b,c)
        if snap_to_ideal:
            bond_length = (snap_interp*get_ideal_bond_length(constraints_handler,(b,c)) + (1-snap_interp)*bond_length)
            angle = (snap_interp*get_ideal_angle(constraints_handler,(a,b,c)) + (1-snap_interp)*angle)
            if verbose:
                print(f"Snapping to bond length {bond_length} and angle {angle}")

        arc = get_arc(a.get_coord(),b.get_coord(),bond_length,angle,c.get_coord()) # TODO current phi... to speed things up
        arcs[altloc] = arc
        arcs_to_plot.append(arc)
        bond_lengths.append(bond_length)
        angles.append(angle)
    
    assert 0 < min_twist_angle <=90

    arcA_coords,arcB_coords = arcs.values()
    possible_true_coords=[]
    for arcA_coord in arcA_coords:
        for arcB_coord in arcB_coords:
            if match_midpoint(arcA_coord,arcB_coord,midpoint,tol_frac=midpoint_match_tol_frac): # May be empty list
                possible_true_coords.append((arcA_coord,arcB_coord))
    c_vector=c_coords[1]-c_coords[0]
    for p in possible_true_coords:
        p_vector=p[1]-p[0]

        # multivariate_normal.pdf()
        #TODO Condition A should change to the electron density overlap of the two atoms
        condition_A = max(fake_separation*min_ratio_real_sep_on_fake_sep,separation(*p)-max_difference_real_fake_sep) <= separation(*p) <= min(max_true_separation,separation(*p)+max_difference_real_fake_sep)
        condition_B = min_twist_angle < relative_orientation(p_vector,c_vector) < 180-min_twist_angle

        if all((condition_A,condition_B)):
            twist_points.append(p)
    if take_closest:
        smallest_diff=np.inf # difference in separation
        for pair in twist_points:
            diff = abs(separation(*C)-separation(*pair))
            if diff < smallest_diff:
                twist_points=np.array([pair])
                smallest_diff=diff


    if take_average and len(twist_points)>0:
        twist_points = np.array([[np.mean(np.array(twist_points)[:,i],axis=0) for i in range(2)]])
    if plot_zoomed_in or plot_zoomed_out:
        a_coords,b_coords,c_coords = [],[],[]
        for altloc in ordered_atoms:
            a_coords.append(ordered_atoms[altloc][0].get_coord())
            b_coords.append(ordered_atoms[altloc][1].get_coord())
            c_coords.append(ordered_atoms[altloc][2].get_coord())
        #print(bond_length,angle);print(a_coords,b_coords)
        names = [atom[0].get_name() for atom in atoms]
        if plot_zoomed_in:
            plot_it(a_coords,b_coords,c_coords,twist_points,arcs_to_plot,
                    names=names,
                    focus_twist=True,rot_speed=plot_rot_speed)
        if plot_zoomed_out:
            plot_it(a_coords,b_coords,c_coords,twist_points,arcs_to_plot,
                    names=names,
                    focus_twist=False,rot_speed=plot_rot_speed)
        
    return twist_points
        

def match_midpoint(coordA,coordB,needed_midpoint,tol_frac):
    tol = tol_frac*separation(coordA,coordB)
    midpoint = (coordA+coordB)/2
    return separation(midpoint,needed_midpoint)<=tol



def plot_it(A,B,C,twist_points,arcs,focus_twist=True,rot_speed=1,names=[]):
    ((a,b,c),(a2,b2,c2)) = zip(A,B,C)


    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    X,Y,Z = zip(*arcs[0]) 
    ax.scatter(*zip(a,b,c),color=["red","purple","blue"],s=100)
    for i, x,y,z, name in zip(range(3),*zip(a,b,c),names):
        ax.text(x,y,z,name)

    ax.scatter(X,Y,Z,s=10,color="blue")
    
    
    X,Y,Z = zip(*arcs[1])
    ax.scatter(X,Y,Z,s=10,color="cyan")
    ax.scatter(*zip(a2,b2,c2),color=["brown","orange","cyan"],s=100)
    ax.view_init(0, 0, 0)

    for pair in twist_points:
        ax.scatter(*[np.array(pair)[:,i] for i in range(3)],marker="X",s=250,color=["blue","cyan"]) # color="maroon",


    buff=0.2
    box_length = max([max([x[i] for x in C]) - min([x[i] for x in C]) for i in range(3)])+buff
    centre = np.mean(C,axis=0)
    if focus_twist:
        for i, func in zip(range(2),(ax.axes.set_ylim3d,ax.axes.set_zlim3d)):
            func(bottom=centre[i+1]-box_length/2, top=centre[i+1]+box_length/2)
        ax.axes.set_xlim3d(left=centre[0]-box_length/2, right=centre[0]+box_length/2)

    # Rotate the axes and update
    # https://matplotlib.org/stable/gallery/mplot3d/rotate_axes3d_sgskip.html
    for cam_angle in range(0, int(180*3/rot_speed) + 1):
        # Normalize the angle to the range [-180, 180] for display
        cam_angle=cam_angle*rot_speed
        angle_norm = (cam_angle + 180) % 360 - 180

        # Cycle through a full rotation of elevation, then azimuth, roll, and all
        elev = azim = roll = 0
        if cam_angle <= 180:
            elev = angle_norm
        elif cam_angle <= 180*2:
            azim = angle_norm
        elif cam_angle <= 180*3:
            roll = angle_norm
        # else:
        #     elev = azim = roll = angle_norm

        # Update the axis view and title
        ax.view_init(elev, azim, roll)
        plt.title('Elevation: %d°, Azimuth: %d°, Roll: %d°' % (elev, azim, roll))

        plt.draw()
        plt.pause(.001)



def detect_twists(ordered_atom_lookup:OrderedAtomLookup,target_res_num:int,atom_name:str, constraints_handler:ConstraintsHandler,
                     max_solution_conformer_sep=0.1,include_CB_angles=False, verbose=False, **kwargs):
    # TODO Don't apply strong constraints, and return stats of the twists, so they can be filtered after calling this function.

    n=target_res_num # 50

    if atom_name =="N":
        angles = [
            # CA_{n-1} -- C_{n-1} -- N_n -- CA_n -- C_n 
            {"CA":n-1,"C":n-1,"N":n}, # CA_{n-1} -- C_{n-1} -- **N_n**
            {"C":n,"CA":n,"N":n,}, # #  **N_n** -- CA_n -- C_n 
        ]
    elif atom_name=="CA":
        angles = [
            {"C":n-1,"N":n,"CA":n},
            {"N":n+1,"C":n,"CA":n},
        ]
        if include_CB_angles:
            raise NotImplementedError()
    elif atom_name=="C":
        angles = [
            {"N":n,"CA":n,"C":n},
            {"CA":n+1,"N":n+1,"C":n},
        ]
        if include_CB_angles:
            raise NotImplementedError()
    else:
        print(f"Invalid atom name! {atom_name}. Please choose N, CA, or C.") # Maybe CB one day. Really need a class designed to easily go down along path of bonded atoms.



    twist_point_sets:list[list]=[]
    #protein_resnums = [r.get_id()[1] for r in struct.get_residues() if not UntangleFunctions.res_is_water(r)]
    protein_resnums = [n for n in ordered_atom_lookup.residue_nums if n not in ordered_atom_lookup.water_residue_nums]
    new_angles:list[dict[str,int]]=[]
    for angle in angles:
        if all([n in protein_resnums for n in angle.values()]):
            new_angles.append(angle) 
    angles=new_angles

    for name_resnums in angles:
        disordered_atoms={} # The atoms involved in the angle 
        for name, resnum in name_resnums.items():
            disordered_atoms[name]=[ordered_atom_lookup.better_dict[resnum][name][altloc] for altloc in ordered_atom_lookup.protein_altlocs]
        assert len(disordered_atoms)==3
        C=disordered_atoms[atom_name]
        # Figure out which atom is A and which is B
        def chain_order(dict_item):
            name,resnum = dict_item
            return resnum*4+["N","CA","CB","C"].index(name)
        atoms_in_a_row = [disordered_atoms[a[0]] for a in sorted(name_resnums.items(),key=chain_order)]
        if atoms_in_a_row[-1]!=C:
            atoms_in_a_row.reverse()
            assert atoms_in_a_row[-1]==C
            
        twist_pairs=detect_twist(atoms_in_a_row,
                                 constraints_handler=constraints_handler,
                                 verbose=verbose,
                                 **kwargs)
        twist_point_sets.append(twist_pairs) # pairs of points
        if verbose and len(twist_pairs)==1:
            print(f"candidate sep: {separation(*twist_pairs[0]):.3f} original: {separation(*C):.3f}")
            #print(np.mean(twist_pairs,axis=0))

    twist_point_combined_constraint_solutions=[]
    
    bond_needs_to_be_flipped=None
    sol = None
    if len(angles)==1: # N-terminus, or C-terminus with CB angles off
        if len(twist_point_sets[0])>0:
            sol = twist_point_sets
    else:
        seps_unflipped=[]
        seps_flipped=[]
        if len(twist_point_sets[0])==len(twist_point_sets[1])==1:
            twist_point_sets = np.array(twist_point_sets)
            # Get the deviation in the conformer coordinates for the solutions found for each geometric constraint. 
            # Since there is a bond either side of nitrogen, we need to account for the possibility that one of the bonds is "flipped" relative to the other. 
            # NOTE flipped does NOT refer to the atom site, it asks whether the bond has been flipped 
            for i in range(2):
                conformer_separation = separation(twist_point_sets[0][0][i],twist_point_sets[1][0][i])
                if verbose:
                    print(f"conformer {i} separation: {conformer_separation}")
                seps_unflipped.append(conformer_separation)
                # reverse
                conformer_separation = separation(twist_point_sets[0][0][i],twist_point_sets[1][0][1-i])
                if verbose:
                    print(f"conformer {i} separation flipped: {conformer_separation}")
                seps_flipped.append(conformer_separation)

            unflipped_sol = np.mean(twist_point_sets,axis=0)
            flipped_sol = np.zeros(shape=unflipped_sol.shape)
            for i, point in enumerate(unflipped_sol):
                flipped_sol[i]=point[::-1]

            unflipped_sol= np.mean(twist_point_sets,axis=0)
            assert twist_point_sets.shape==(2,1,2,3)
            flipped_twist_point_sets = np.zeros(shape=twist_point_sets.shape)
            flipped_twist_point_sets[0]=twist_point_sets[0]
            flipped_twist_point_sets[1]=twist_point_sets[1,:,::-1]
            flipped_sol= np.mean(flipped_twist_point_sets,axis=0)



            rmsd_unflipped = np.sqrt( np.sum(np.array(seps_unflipped)**2) )
            rmsd_flipped = np.sqrt( np.sum(np.array(seps_flipped)**2) )
            #print(rmsd_unflipped,rmsd_flipped)

            options = [(unflipped_sol,seps_unflipped,False),(flipped_sol,seps_flipped,True)]
            valid_options=[]
            for option in options:
                seps = option[1]
                for conformer_separation in seps:
                    if conformer_separation > max_solution_conformer_sep:
                        break
                else:
                    valid_options.append(option)
            
            # Choose highest sep option
            if len(valid_options)>0:
                #print(list(opt[0] for opt in valid_options))
                highest_sep=-1
                for option in valid_options:
                    option_sol=option[0]
                    sep = separation(*option_sol[0])
                    if sep > highest_sep:
                        highest_sep=sep
                        sol,_,bond_needs_to_be_flipped = option
                print(f"Separation {separation(*C)} --> {highest_sep}")

        elif len(twist_point_sets[0])>1 or len(twist_point_sets[1]) > 1:
            print("Multiple points for sets not handled")
            assert False
        else:
            pass
    if sol is not None:
        twist_point_combined_constraint_solutions.append(sol)    # XXX stupid         
    return twist_point_sets,np.array(twist_point_combined_constraint_solutions),bond_needs_to_be_flipped  


#%%
if __name__=="__main__":

    # Two-conformer implementation 

    pdb_path="/home/speno/Untangler/output/longrangetraps_TW_unrestrained2.pdb"
    #pdb_path="/home/speno/Untangler/output/longrangetraps_TW_Accepted1.pdb"

    struct = PDBParser().get_structure("struct",pdb_path)
    ordered_atom_lookup=OrderedAtomLookup(struct.get_atoms(),waters=False)
    constraints_handler = get_constraints_handler(pdb_path,geo_file_needs_generation=False)
    twists_found:list[NDArray]=[]
    twist_atom_ids:list[DisorderedTag]=[]
    bond_flips_needed:list[bool]=[]
    altlocs_involved:list[str]=[]
    debug = False
    single_res_debug=False
    print("Computing twists...")
    last_resnum=64
    first_resnum=1
    if single_res_debug:
        first_resnum=last_resnum=46
    for res_num in range(first_resnum,last_resnum+1):
        if res_num%10==1:
            print(f"residue {res_num}/{last_resnum}")
        for atom_name in "N","CA","C":
            if debug and ((res_num,atom_name) !=(46,"N")):
                continue 
            if debug:
                print(res_num)
            indiv_twist_points,solutions,bond_needs_to_be_flipped=detect_twists(ordered_atom_lookup,res_num,atom_name,constraints_handler,
                                                max_solution_conformer_sep=0.1,
                                                plot_zoomed_in=debug,
                                                plot_zoomed_out=debug,
                                                plot_rot_speed=3,
                                                take_closest=True,
                                                verbose =False
                                                #take_average=True
                                            )
            if len(solutions!=0):
                twists_found.append(solutions)
                twist_atom_ids.append(DisorderedTag(res_num,atom_name))
                bond_flips_needed.append(bond_needs_to_be_flipped)
                altlocs_involved.append("AB") # Temporary FIXME 
            if debug:
                print("====")
                print(indiv_twist_points)
                print("---")
                print(solutions)
                print("---")
                print(bond_needs_to_be_flipped)
                print("====")
    print()
    if len(twists_found)>0:
        print("Possible twists detected:")
        for _,  site_tag, tangled in zip(twists_found, twist_atom_ids, bond_flips_needed):
            print(f"Atom {site_tag}, requires flipping one bond relative to another: {tangled}")
    else:
        print("No twists detected.")
    
    solutions_chosen=np.array(twists_found)[:,0,0] # XXX
    out_handle = os.path.basename(pdb_path)[:-4]
    out_str="# resnum.name | new coords \n"
    np.set_printoptions(formatter={'float': lambda x: "{:.3f}".format(x)})
    # out_str+="\n".join([f"{site_tag} {' '.join(str(coord) for coord in coord_pair)}" \
    #                     for (coord_pair,  site_tag, tangled) in zip(solutions_chosen, twist_atom_ids, bond_flips_needed)]
    #                   )
    out_str+="\n".join([f"{site_tag} {altlocs} {' '.join(str(coord) for coord in coord_pair)} {'Y' if tangled else 'N'}" \
                        for (coord_pair,  site_tag, altlocs, tangled) in zip(solutions_chosen, twist_atom_ids,altlocs_involved, bond_flips_needed)]
                      )
    with open(f"untwist_moves_{out_handle}.txt",'w') as f:
        f.write(out_str)
    print("Done")

        
    #print(detect_twist(tuple(disordered_atoms),snap_to_ideal=False,plot=True))

#%%
if __name__=="__main__" and False:

    ideal_bond_length,ideal_angle=1.5,110
    diag_v = np.array([1.,1.,1.])
    diag_v/=norm(diag_v)
    a,b = 0.5*diag_v, 2*diag_v
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    X,Y,Z = zip(*get_arc(a,b,ideal_bond_length,ideal_angle))
    ax.scatter(X,Y,Z,s=100)
    ax.scatter(*zip(a,b),color=["red","purple"],s=100)
    
    a2,b2 = 0.5*diag_v+np.array([0.4,0.2,0]), 2*diag_v+np.array([0,0.7,-0.35])
    X,Y,Z = zip(*get_arc(a2,b2,ideal_bond_length,ideal_angle))
    ax.scatter(X,Y,Z,s=100)
    ax.scatter(*zip(a2,b2),color=["brown","orange"],s=100)
    ax.view_init(0, 0, 0)
    # Rotate the axes and update
    # https://matplotlib.org/stable/gallery/mplot3d/rotate_axes3d_sgskip.html
    rot_speed=2
    for angle in range(0, int(180*3/rot_speed) + 1):
        # Normalize the angle to the range [-180, 180] for display
        angle=angle*rot_speed
        angle_norm = (angle + 180) % 360 - 180

        # Cycle through a full rotation of elevation, then azimuth, roll, and all
        elev = azim = roll = 0
        if angle <= 180:
            elev = angle_norm
        elif angle <= 180*2:
            azim = angle_norm
        elif angle <= 180*3:
            roll = angle_norm
        # else:
        #     elev = azim = roll = angle_norm

        # Update the axis view and title
        ax.view_init(elev, azim, roll)
        plt.title('Elevation: %d°, Azimuth: %d°, Roll: %d°' % (elev, azim, roll))

        plt.draw()
        plt.pause(.001)

    



#%%

# if __name__=="__main__":
#     twists_detected=[]
#     angles_considered = get_angle_atoms()
#     for A,B,C in angles_considered:
#         twist = detect_twist()


    # from mpl_toolkits.mplot3d import axes3d

    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')

    # # Grab some example data and plot a basic wireframe.
    # X, Y, Z = axes3d.get_test_data(0.05)
    # ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)

    # # Set the axis labels
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.set_zlabel('z')

    # # Rotate the axes and update
    # for angle in range(0, 360*4 + 1):
    #     # Normalize the angle to the range [-180, 180] for display
    #     angle_norm = (angle + 180) % 360 - 180

    #     # Cycle through a full rotation of elevation, then azimuth, roll, and all
    #     elev = azim = roll = 0
    #     if angle <= 360:
    #         elev = angle_norm
    #     elif angle <= 360*2:
    #         azim = angle_norm
    #     elif angle <= 360*3:
    #         roll = angle_norm
    #     else:
    #         elev = azim = roll = angle_norm

    #     # Update the axis view and title
    #     ax.view_init(elev, azim, roll)
    #     plt.title('Elevation: %d°, Azimuth: %d°, Roll: %d°' % (elev, azim, roll))

    #     plt.draw()
    #     plt.pause(.001)