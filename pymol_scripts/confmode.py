# good default for viewing ensembles in pymol

from pymol import cmd
def confmode(water=0,white=0,mch=0):
    water=int(water); white=int(white); mch=int(mch)

    cmd.hide("everything","all")
    cmd.show("lines spheres","all")
    cmd.set("sphere_scale",0.1)
    cmd.hide("(hydro)")  # Hide hydrogen. Same as cmd.hide("everything", "hydro")  
    #cmd.hide("(hydro or solvent)")
    print(water)
    print(type(water))
    print(water=='0')
    if not water:
        cmd.hide("(resn HOH)")

    if mch:
        cmd.hide(("((byres (All))&(sc.&!(n. CA|n. N&r. PRO)))"))
        cmd.hide("(name O)")
    
    if white:
        cmd.color("white")
    else:
        for c,altloc in zip(
        ["orange","purple","green","yellow","red","blue","brown","cyan","gray","darksalmon","deepolive","hotpink","deepteal","lime","sand","ruby","marine","magenta"],
        "ABCDEFGHIJKLMNOPQR"):
            cmd.color(c,f"alt {altloc}")


    
cmd.extend('confmode', confmode)