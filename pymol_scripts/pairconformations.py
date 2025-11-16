'''

Original Authors: Spencer Passmore


License:
'''

from pymol import cmd, CmdException


# For two ensemble models, split into pairs to give the lowest possible RMSD

# CURRENTLY assumes 2 conformations with altlocs labelled A and B.
def pairconformations(modelA, modelB, doAlign=1, doPretty=1, guide=1, method='super', quiet=1):
    modelA_confs=[]
    modelB_confs=[]
    for altloc in ["A","B"]:
        modelA_confs.append(f"modelA-{altloc}")
        cmd.create(modelA_confs[-1], f"model {modelA} and altloc {altloc}")
        modelB_confs.append(f"modelB-{altloc}")
        cmd.create(modelB_confs[-1], f"model {modelB} and altloc {altloc}")
    
    options=[]
    options.append(((modelA_confs[0],modelB_confs[0]),(modelA_confs[1],modelB_confs[1])))
    options.append(((modelA_confs[0],modelB_confs[0]),(modelA_confs[1],modelB_confs[1])))

    lowest_distance=1e20
    best_option=None
    for option in options:
        avg_distance=0
        for pair in option:
            avg_distance+= get_rmsd(*pair)
        avg_distance/=len(option)
        if avg_distance< lowest_distance:
            lowest_distance=avg_distance 
            best_option=option
    # Remove worse options
    for i,pair in enumerate(best_option):
        colorbyrmsd(*pair, doAlign, doPretty, guide, method, quiet)
        cmd.group(f"pair-{i}",f"model {pair[0]} or model {pair[1]}")


    
cmd.extend('pairconformations', pairconformations)

# tab-completion of arguments
cmd.auto_arg[0]['pairconformations'] = cmd.auto_arg[0]['align']
cmd.auto_arg[1]['pairconformations'] = cmd.auto_arg[1]['align']



def get_rmsd(mobile, target, doAlign=1, doPretty=1, guide=1, method='super', quiet=1):
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
        cmd.color('gray', seleboth)
        cmd.spectrum('b', 'blue_red', seleboth + ' and b > -0.5')

    if not quiet:
        print(" ColorByRMSD: Minimum Distance: %.2f" % (min(b_dict.values())))
        print(" ColorByRMSD: Maximum Distance: %.2f" % (max(b_dict.values())))
        print(" ColorByRMSD: Average Distance: %.2f" % (sum(b_dict.values()) / len(b_dict)))

    cmd.delete(aln)
    cmd.delete(seleboth)