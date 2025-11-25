#%%
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
import matplotlib.pyplot as plt


#(I) they share the same midpoint
# (II) they are farther apart
# (III) each conformer position is related via rotation of a model conformer around the vector formed by a and b.
# 
# 
# 
# 

def angle(a:Atom,b:Atom,c:Atom):
    v1 = a.get_coord() - b.get_coord()
    v2 = c.get_coord()-b.get_coord()
    def unit_vector(vector):
        return vector / np.linalg.norm(vector)
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.dot(v1_u, v2_u))*180/np.pi

def separation(v1,v2):
    return np.sqrt(np.sum((v1-v2)**2))

def separation_condition(atom:Atom,max_separation):

    def atom_sep(a:Atom,b:Atom):
        return np.sqrt(np.sum((a.get_coord()-b.get_coord())**2))
    assert type(atom)==DisorderedAtom
    A,B = [a for a in atom]
    return atom_sep(A,B)<=max_separation


def get_angle_atoms()->list[tuple[Atom,Atom,Atom]]:
    assert False
def get_ideal_angle()->float:
    assert False
def get_ideal_bond_length()->float:
    assert False


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

    points = mid_anchor+h + r*(np.cos(theta[...,None])*x + np.sin(theta[...,None])*y)
    return points


    



def detect_twist(
    atoms:tuple[Atom,Atom,Atom],
    snap_to_ideal=True  # Snap to ideal angle and bond length. Otherwise, uses current (potentially warped angle and bond length)
):
    twist_points=[]
    
    A,B,C=atoms
    # Check C is close enough. 
    max_separation=0.25 # TODO make it resolution dependent
    if not separation_condition(C,max_separation):
        return twist_points # empty list
    ordered_atoms:dict[str,list[Atom]]={}
    for i, disordered_atom in enumerate(atoms):
        for atom in disordered_atom:
            if atom.get_altloc() not in ordered_atoms:
                ordered_atoms[atom.get_altloc()]=[None,None,None]
            ordered_atoms[atom.get_altloc()][i]=atom

    midpoint = np.mean([c.get_coord() for c in C],axis=0)


    # a✓ -- b✓ -- c? 
    arcs = {altloc:[] for altloc in ordered_atoms.keys()}
    for altloc, (a,b,c) in ordered_atoms.items():
        bond_length,angle=separation(b,c),angle(a,b,c)
        if snap_to_ideal:
            bond_length,angle = get_ideal_bond_length(), get_ideal_angle()

        arc = get_arc(a.get_coord(),b.get_coord(),bond_length,angle)
        arcs[altloc].append(arc)
    
    for arcA_list,arcB_list in arcs.values():
        for arcA in arcA_list:
            for arcB in arcB_list:
                intersect_points = match_midpoint(arcA,arcB,midpoint) # May be empty list
                twist_points.extend(intersect_points)

        








#%%
if __name__=="__main__":

    ideal_bond_length,ideal_angle=1.5,110
    a,b = np.array([0.5,0.5,0.5]), np.array([2,2,2])
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    X,Y,Z = zip(*get_arc(a,b,ideal_bond_length,ideal_angle))
    ax.scatter(X,Y,Z,s=100)
    ax.scatter(*zip(a,b),color=["red","purple"],s=100)
    
    a2,b2 = np.array([0.9,1,0]), np.array([2,2.9,1.5])
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