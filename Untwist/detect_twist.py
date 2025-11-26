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
    constraints_handler.load_all_constraints(pdb_file_path,ordered_atom_lookup,symmetries=None,calc_nonbonds=False)
    return constraints_handler

# Get the circle of possible values for position x that satisfy an angle between start_anchor--mid_anchor--x and distance for mid_anchor--x
def get_arc(start_anchor:tuple[float],mid_anchor:tuple[float],bond_length,atom_angle,num_points=360):
    #https://stackoverflow.com/questions/48703275/3d-truncated-cone-in-python

    # phi = cone angle
    # c = slant_height

    phi = (180-atom_angle)*np.pi/180
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
    theta = np.linspace(0, 2 * np.pi, num_points)

    points = mid_anchor+h*v + r*(np.cos(theta[...,None])*x + np.sin(theta[...,None])*y)
    return points


    
def separation_condition(atom:Atom,max_separation):
    assert type(atom)==DisorderedAtom
    A,B = [a for a in atom]
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
    bond_lengths,angles = [],[]
    for altloc, (a,b,c) in ordered_atoms.items():
        bond_length,angle=separation(b,c),get_angle(a,b,c)
        if snap_to_ideal:
            bond_length = (snap_interp*get_ideal_bond_length(constraints_handler,(b,c)) + (1-snap_interp)*bond_length)
            angle = (snap_interp*get_ideal_angle(constraints_handler,(a,b,c)) + (1-snap_interp)*angle)
            if verbose:
                print(f"Snapping to bond length {bond_length} and angle {angle}")

        arc = get_arc(a.get_coord(),b.get_coord(),bond_length,angle)
        arcs[altloc] = arc
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
        names = [atom.get_name() for atom in atoms]
        if plot_zoomed_in:
            plot_it(a_coords,b_coords,c_coords,twist_points,bond_lengths,angles,
                    names=names,
                    focus_twist=True,rot_speed=plot_rot_speed)
        if plot_zoomed_out:
            plot_it(a_coords,b_coords,c_coords,twist_points,bond_lengths,angles,
                    names=names,
                    focus_twist=False,rot_speed=plot_rot_speed)
        
    return twist_points
        

def match_midpoint(coordA,coordB,needed_midpoint,tol_frac):
    tol = tol_frac*separation(coordA,coordB)
    midpoint = (coordA+coordB)/2
    return separation(midpoint,needed_midpoint)<=tol



def plot_it(A,B,C,twist_points,bond_lengths,angles,focus_twist=True,rot_speed=1,names=[]):
    ((a,b,c),(a2,b2,c2)) = zip(A,B,C)


    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    X,Y,Z = zip(*get_arc(a,b,bond_lengths[0],angles[0]))
    ax.scatter(*zip(a,b,c),color=["red","purple","blue"],s=100)
    for i, x,y,z, name in zip(range(3),*zip(a,b,c),names):
        ax.text(x,y,z,name)

    ax.scatter(X,Y,Z,s=10,color="blue")
    
    
    X,Y,Z = zip(*get_arc(a2,b2,bond_lengths[1],angles[1]))
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



def check_nitrogen_twists(struct:Structure,target_res_num:int, constraints_handler:ConstraintsHandler,max_solution_sep=0.1, **kwargs):
    # TODO Don't apply strong constraints, and return stats of the twists, so they can be filtered after calling this function.

    n=target_res_num # 50


    angles = [
        {"CA":n-1,"C":n-1,"N":n},
        {"C":n,"CA":n,"N":n,},
    ]
    target_atom="N"
    twist_point_sets:list[list]=[]
    for atoms_resnums in angles:
        disordered_atoms={}
        for residue in struct.get_residues():
            resnum=residue.get_id()[1]
            if resnum in atoms_resnums.values():
                for atom in residue.get_atoms():
                    if atom.get_name() not in atoms_resnums:
                        continue
                    if atoms_resnums[atom.get_name()]==resnum:
                        assert atom.get_name() not in disordered_atoms
                        disordered_atoms[atom.get_name()]=atom
        assert target_atom == "N"
        
        C=disordered_atoms[target_atom]
        if target_atom == "N":
            if atoms_resnums["CA"]==atoms_resnums["C"]==atoms_resnums["N"]-1:
                A=disordered_atoms["CA"]
                B=disordered_atoms["C"]
            elif atoms_resnums["N"]==atoms_resnums["CA"]==atoms_resnums["C"]:
                # Backwards
                A=disordered_atoms["C"]
                B=disordered_atoms["CA"]

            
        twist_pairs=detect_twist((A,B,C),
                                  constraints_handler=constraints_handler,
                                **kwargs)
        twist_point_sets.append(twist_pairs) # pairs of points
        print(twist_pairs)
        if len(twist_pairs)==1:
            print(f"candidate sep: {separation(*twist_pairs[0]):.3f} original: {separation(*C):.3f}")
            #print(np.mean(twist_pairs,axis=0))

    twist_point_combined_constraint_solutions=[]
    
    seps_unflipped=[]
    seps_flipped=[]
    bond_needs_to_be_flipped=None
    if len(twist_point_sets[0])==len(twist_point_sets[1])==1:
        twist_point_sets = np.array(twist_point_sets)
        # Get the deviation in the conformer coordinates for the solutions found for each geometric constraint. 
        # Since there is a bond either side of nitrogen, we need to account for the possibility that one of the bonds is "flipped" relative to the other. 
        # NOTE flipped does NOT refer to the atom site, it asks whether the bond has been flipped 
        for i in range(2):
            conformer_separation = separation(twist_point_sets[0][0][i],twist_point_sets[1][0][i])
            print(f"conformer {i} separation: {conformer_separation}")
            seps_unflipped.append(conformer_separation)
            # reverse
            conformer_separation = separation(twist_point_sets[0][0][i],twist_point_sets[1][0][1-i])
            print(f"conformer {i} separation flipped: {conformer_separation}")
            seps_flipped.append(conformer_separation)

        unflipped_sol = np.mean(twist_point_sets,axis=0)
        flipped_sol = np.zeros(shape=unflipped_sol.shape)
        for i, point in enumerate(unflipped_sol):
            flipped_sol[i]=point[::-1]

        rmsd_unflipped = np.sqrt( np.sum(np.array(seps_unflipped)**2) )
        rmsd_flipped = np.sqrt( np.sum(np.array(seps_flipped)**2) )

        sol,seps,bond_needs_to_be_flipped = (unflipped_sol,seps_unflipped,False) if rmsd_unflipped < rmsd_flipped else (flipped_sol,seps_flipped,True)

        for conformer_separation in seps:
            if conformer_separation > max_solution_sep:
                break
        else:
            twist_point_combined_constraint_solutions.append(sol)            
    elif len(twist_point_sets[0])>1 or len(twist_point_sets[1]) > 1:
        print("Multiple points for sets not handled")
    else:
        pass
    return twist_point_sets,np.array(twist_point_combined_constraint_solutions),bond_needs_to_be_flipped  


#%%
if __name__=="__main__":

    # Tw conformer case 

    pdb_path="/home/speno/Untangler/output/longrangetraps_TW_unrestrained2_fmtd.pdb"

    struct = PDBParser().get_structure("struct",pdb_path)
    constraints_handler = get_constraints_handler(pdb_path,geo_file_needs_generation=False)
    twists_found=[]
    twist_resnums=[]
    bond_flips_needed=[]
    for res_num in range(2,64):
        print(res_num)
        indiv_twist_points,solutions,bond_needs_to_be_flipped=check_nitrogen_twists(struct,res_num,constraints_handler,
                                            max_solution_sep=0.1,
                                            #plot_zoomed_in=True,
                                            #plot_zoomed_out=True,
                                            take_closest=True,
                                            #take_average=True
                                        )
        if bond_needs_to_be_flipped is not None:
            twists_found.append(solutions)
            twist_resnums.append(res_num)
            bond_flips_needed.append(bond_needs_to_be_flipped)
        print("====")
        print(indiv_twist_points)
        print("---")
        print(solutions)
        print("---")
        print(bond_needs_to_be_flipped)
        print("====")

    print("Possible twists found:")
    for twists_found,  resnum, tangled in zip(twists_found, twist_resnums, bond_flips_needed):
        print(f"Res num {resnum}, requires flipping one bond relative to another: {tangled}")
        
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