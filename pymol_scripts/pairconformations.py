'''

Original Authors: Spencer Passmore


License:
'''

from pymol import cmd, CmdException
from time import sleep

# For two ensemble models, split into pairs to give the lowest possible RMSD

# CURRENTLY assumes 2 conformations with altlocs labelled A and B.

def pairconformations(modelA, modelB, good='ytterbium',bad='yellow',
                      colorbysep=1,sep_cutoff=0.05, doAlign=0, doPretty=1, 
                      guide=0, method='super', quiet=1):
    # try ytterbium is a moderate green
    def sele_exists(sele):
        sess = cmd.get_session()
        for i in sess["names"]:
            if type(i) == list:
                if sele==i[0]:
                    return 1
        return 0
    
    def new_name(new_model_name):
        return_name=new_model_name
        n=1
        while sele_exists(f"{return_name}"):
            return_name=f"{new_model_name}-{n}"
            n+=1
        return return_name
    def delete_name(new_model_name):
        return_name=new_model_name
        n=1
        if sele_exists(f"{return_name}"):
            cmd.delete(return_name)
        return return_name
    
    delete_name("_aln")
    delete_name("_objSelBoth")
    modelA_confs=[]
    modelB_confs=[]
    for altloc in ["A","B"]:
        modelA_confs.append(delete_name(f"{modelA}-{altloc}"))
        cmd.create(modelA_confs[-1], f"model {modelA} and altloc {altloc}")
        cmd.alter(modelA_confs[-1],'alt="X"') 
        modelB_confs.append(delete_name(f"{modelB}-{altloc}"))
        cmd.create(modelB_confs[-1], f"model {modelB} and altloc {altloc}")
        cmd.alter(modelB_confs[-1],'alt="X"') 
    
    options=[]
    options.append(((modelA_confs[0],modelB_confs[0]),(modelA_confs[1],modelB_confs[1])))
    options.append(((modelA_confs[0],modelB_confs[1]),(modelA_confs[1],modelB_confs[0])))

    lowest_distance=1e20
    best_option=None
    for option in options:
        avg_distance=0
        for pair in option:
            avg_distance+= get_dist(*pair)
        avg_distance/=len(option)
        if avg_distance< lowest_distance:
            lowest_distance=avg_distance 
            best_option=option
        print(avg_distance)
    # Group according to best option
    for i,pair in enumerate(best_option):
        # for model in pair:
        #     cmd.color({0:"orange",1:"purple"}[i],f"model {model}")
        for model in pair:
            cmd.color("grey",f"model {model}")
            cmd.set_bond("stick_radius",0.1,f"sidechain and model {model}")
            cmd.set("sphere_scale",0.15)
        if colorbysep:
            color_sep(*pair, sep_cutoff,doAlign, doPretty, guide, method, quiet=1,goodcolor=good,badcolor=bad)
        else:
            pass
            colorbyrmsd(*pair, doAlign, doPretty, guide, method, quiet=1)
        group_name=f"{modelA}-{modelB}-pair-{i}"
        cmd.delete(group_name)
        cmd.group(group_name,f"model {pair[0]} or model {pair[1]}")


    
cmd.extend('pairconformations', pairconformations)

# tab-completion of arguments
cmd.auto_arg[0]['pairconformations'] = cmd.auto_arg[0]['align']
cmd.auto_arg[1]['pairconformations'] = cmd.auto_arg[1]['align']

def color_sep(mobile, target,sep_cutoff, doAlign=0, doPretty=1, guide=0, method='super', quiet=1,goodcolor='green',badcolor='yellow'):
    from chempy import cpv

    doAlign, doPretty = int(doAlign), int(doPretty)
    guide, quiet = int(guide), int(quiet)
    aln, seleboth = '_aln', '_objSelBoth'

    try:
        align = cmd.keyword[method][0]
    except:
        print(' Error: no such method: ' + str(method))
        raise CmdException

    if guide:
        mobile = '(%s) and guide' % mobile
        target = '(%s) and guide' % target

    try:
        if doAlign:
            # superpose
            align(mobile, target)

        # get alignment without superposing
        align(mobile, target, cycles=0, transform=0, object=aln)
    except:
        print(' Error: Alignment with method %s failed' % (method))
        raise CmdException

    cmd.select(seleboth, '(%s) or (%s)' % (mobile, target))

    idx2coords = dict()
    cmd.iterate_state(-1, seleboth, 'idx2coords[model,index] = (x,y,z)', space=locals())

    if cmd.count_atoms('?' + aln, 1, 1) == 0:
        # this should ensure that "aln" will be available as selectable object
        cmd.refresh()

    b_dict = dict()
    for col in cmd.get_raw_alignment(aln):
        assert len(col) == 2
        b = cpv.distance(idx2coords[col[0]], idx2coords[col[1]])
        for idx in col:
            b_dict[idx] = b

    cmd.alter(seleboth, 'b = b_dict.get((model, index), -1)', space=locals())
    cmd.refresh()
    if doPretty:
        cmd.orient(seleboth)
        # cmd.show_as('cartoon', 'byobj ' + seleboth)
        # cmd.show('sticks', 'byobj ' + seleboth)
        cmd.color('gray',seleboth)
        cmd.color(goodcolor,  f"{seleboth} and b > -0.5")
        cmd.color(badcolor,  f"{seleboth} and b > {sep_cutoff}")

    if not quiet:
        print(" ColorByRMSD: Minimum Distance: %.2f" % (min(b_dict.values())))
        print(" ColorByRMSD: Maximum Distance: %.2f" % (max(b_dict.values())))
        print(" ColorByRMSD: Average Distance: %.2f" % (sum(b_dict.values()) / len(b_dict)))

    cmd.delete(aln)
    cmd.delete(seleboth)

def get_dist(mobile, target, doAlign=0, doPretty=1, guide=0, method='super', quiet=1):
    '''
    Taken from:
    http://pymolwiki.org/index.php/ColorByRMSD

    Original Authors: Shivender Shandilya; Jason Vertrees
    Complete rewrite by Thomas Holder

    License: BSD-2-Clause
    '''

    from chempy import cpv

    doAlign, doPretty = int(doAlign), int(doPretty)
    guide, quiet = int(guide), int(quiet)
    aln, seleboth = '_aln', '_objSelBoth'

    try:
        align = cmd.keyword[method][0]
    except:
        print(' Error: no such method: ' + str(method))
        raise CmdException

    if guide:
        mobile = '(%s) and guide' % mobile
        target = '(%s) and guide' % target

    try:
        if doAlign:
            # superpose
            align(mobile, target)

        # get alignment without superposing
        align(mobile, target, cycles=0, transform=0, object=aln)
    except:
        print(' Error: Alignment with method %s failed' % (method))
        raise CmdException

    cmd.select(seleboth, '(%s) or (%s)' % (mobile, target))

    idx2coords = dict()
    cmd.iterate_state(-1, seleboth, 'idx2coords[model,index] = (x,y,z)', space=locals())

    if cmd.count_atoms('?' + aln, 1, 1) == 0:
        # this should ensure that "aln" will be available as selectable object
        cmd.refresh()

    b_dict = dict()
    for col in cmd.get_raw_alignment(aln):
        assert len(col) == 2
        b = cpv.distance(idx2coords[col[0]], idx2coords[col[1]])
        for idx in col:
            b_dict[idx] = b

    # cmd.alter(seleboth, 'b = b_dict.get((model, index), -1)', space=locals())

    # if doPretty:
    #     cmd.orient(seleboth)
    #     cmd.show_as('cartoon', 'byobj ' + seleboth)
    #     cmd.color('gray', seleboth)
    #     cmd.spectrum('b', 'blue_red', seleboth + ' and b > -0.5')

    
    if not quiet:
        print(" ColorByRMSD: Minimum Distance: %.2f" % (min(b_dict.values())))
        print(" ColorByRMSD: Maximum Distance: %.2f" % (max(b_dict.values())))
        print(" ColorByRMSD: Average Distance: %.2f" % (sum(b_dict.values()) / len(b_dict)))

    cmd.delete(aln)
    cmd.delete(seleboth)

    mean = sum(b_dict.values()) / len(b_dict)
    return mean


##########################


'''
http://pymolwiki.org/index.php/ColorByRMSD

Original Authors: Shivender Shandilya; Jason Vertrees
Complete rewrite by Thomas Holder

License: BSD-2-Clause
'''

from pymol import cmd, CmdException


def colorbyrmsd(mobile, target, doAlign=1, doPretty=1, guide=1, method='super', quiet=1):
    '''
DESCRIPTION

    Align two structures and show the structural deviations in color to more
    easily see variable regions.

    Colors each mobile/target atom-pair by distance (the name is a bit
    misleading).

    Modifies the B-factor columns in your original structures.

ARGUMENTS

    mobile = string: atom selection for mobile atoms

    target = string: atom selection for target atoms

    doAlign = 0 or 1: Superpose selections before calculating distances
    {default: 1}

    doPretty = 0 or 1: Show nice representation and colors {default: 1}

EXAMPLE

    fetch 1ake 4ake, async=0
    remove chain B
    colorbyrmsd 1ake, 4ake
    '''
    from chempy import cpv

    doAlign, doPretty = int(doAlign), int(doPretty)
    guide, quiet = int(guide), int(quiet)
    aln, seleboth = '_aln', '_objSelBoth'

    try:
        align = cmd.keyword[method][0]
    except:
        print(' Error: no such method: ' + str(method))
        raise CmdException

    if guide:
        mobile = '(%s) and guide' % mobile
        target = '(%s) and guide' % target

    try:
        if doAlign:
            # superpose
            align(mobile, target)

        # get alignment without superposing
        align(mobile, target, cycles=0, transform=0, object=aln)
    except:
        print(' Error: Alignment with method %s failed' % (method))
        raise CmdException

    cmd.select(seleboth, '(%s) or (%s)' % (mobile, target))

    idx2coords = dict()
    cmd.iterate_state(-1, seleboth, 'idx2coords[model,index] = (x,y,z)', space=locals())

    if cmd.count_atoms('?' + aln, 1, 1) == 0:
        # this should ensure that "aln" will be available as selectable object
        cmd.refresh()

    b_dict = dict()
    for col in cmd.get_raw_alignment(aln):
        assert len(col) == 2
        b = cpv.distance(idx2coords[col[0]], idx2coords[col[1]])
        for idx in col:
            b_dict[idx] = b

    cmd.alter(seleboth, 'b = b_dict.get((model, index), -1)', space=locals())
    if doPretty:
        cmd.orient(seleboth)
        cmd.show_as('cartoon', 'byobj ' + seleboth)
        cmd.show('sticks', 'byobj ' + seleboth)
        cmd.color('gray', seleboth)
        cmd.spectrum('b', 'blue_red', seleboth + ' and b > -0.5')

    if not quiet:
        print(" ColorByRMSD: Minimum Distance: %.2f" % (min(b_dict.values())))
        print(" ColorByRMSD: Maximum Distance: %.2f" % (max(b_dict.values())))
        print(" ColorByRMSD: Average Distance: %.2f" % (sum(b_dict.values()) / len(b_dict)))

    cmd.delete(aln)
    cmd.delete(seleboth)