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
    cmd.set_bond("stick_radius",0.125,f"sidechain and {model}")
    cmd.color('gray', model)
    cmd.color('red', f"{model} and b > 0")

cmd.extend('colorbyb', colorbyb)

# vi: ts=4:sw=4:smarttab:expandtab