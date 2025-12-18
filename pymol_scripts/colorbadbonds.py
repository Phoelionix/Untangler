# For use with model outputted by Measures/evaluate_tangle.py

from pymol import cmd, CmdException


def colorbadbonds(model):
    '''
DESCRIPTION



ARGUMENTS


EXAMPLE


    '''

    def sele_exists(sele):
        sess = cmd.get_session()
        for i in sess["names"]:
            if type(i) == list:
                if sele==i[0]:
                    return 1
        return 0

    stick_factor=1.5

    #cmd.set_bond("stick_radius",0.125,f"sidechain and {model}")
    bad_bonds = f"badBonds_{model}"
    model_copy = f"copy_{model}"
    cmd.create(bad_bonds, f"{model} and b > 0")
    if not sele_exists(f"{bad_bonds}"):
        print("No bad bonds")
        return 
    cmd.create(model_copy, f"{model}")
        
    for m in [model_copy,bad_bonds]:
        cmd.show_as('sticks', f"model {m}")
    #cmd.set_bond("stick_radius",0.2*stick_factor,f"{bad_bonds} and b > 0")
    cmd.set_bond("stick_radius",0.125*stick_factor,model_copy)
    cmd.set_bond("stick_radius",0.16*stick_factor,bad_bonds)
    cmd.color('gray', model_copy)
    cmd.color('red', f"{bad_bonds} and b > 0")

    group_name=f"{model}-badbond"
    cmd.delete(group_name)
    cmd.group(group_name,f"{model_copy} or {bad_bonds}")




cmd.extend('colorbadbonds', colorbadbonds)

cmd.auto_arg[0]['colorbadbonds'] = [ cmd.object_sc, 'object', '']

# vi: ts=4:sw=4:smarttab:expandtab