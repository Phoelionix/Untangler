from pymol import cmd, CmdException
from UntangleFunctions import model_handle


def highlightconnections(model_path,tag,stick_radius=0.1,color=None,sele="b > 0"):
    def sele_exists(sele):
        sess = cmd.get_session()
        for i in sess["names"]:
            if type(i) == list:
                if sele==i[0]:
                    return 1
        return 0


    name=model_handle(model_path)
    highlighted_bonds = f"highlighted_{name}_{tag}"
    cmd.load(model_path,highlighted_bonds)
    cmd.create(highlighted_bonds, f"{highlighted_bonds} and {sele}")
    if not sele_exists(f"{highlighted_bonds}"):
        print("No highlighted bonds")
        return 

    for m in [highlighted_bonds,]:
        cmd.show_as('sticks', f"model {m}")
    cmd.set_bond("stick_radius",stick_radius,highlighted_bonds)

    if color is not None:
        cmd.color(color, f"{highlighted_bonds}")
    return highlighted_bonds