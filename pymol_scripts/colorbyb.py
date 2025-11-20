# For use with model outputted by Measures/evaluate_tangle.py

from pymol import cmd, CmdException


def colorbyb(model):
    '''
DESCRIPTION



ARGUMENTS


EXAMPLE


    '''

    cmd.show_as('sticks', f"model {model}")
    cmd.set_bond("stick_radius",0.125,model)
    #cmd.set_bond("stick_radius",0.125,f"sidechain and {model}")
    bad_bonds = f"badBonds_{model}"
    model_copy = f"copy_{model}"
    cmd.create(model_copy, f"{model}")
    cmd.create(bad_bonds, f"{model_copy} and b > 0")
    cmd.set_bond("stick_radius",0.2,f"{bad_bonds} and b > 0")
    cmd.color('gray', model_copy)
    cmd.color('red', f"{bad_bonds} and b > 0")

    group_name=f"{model}-badbond"
    cmd.delete(group_name)
    cmd.group(group_name,f"{model_copy} or {bad_bonds}")




cmd.extend('colorbyb', colorbyb)

# vi: ts=4:sw=4:smarttab:expandtab