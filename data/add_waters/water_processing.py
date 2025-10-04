def delete_water_copies(lines_out):
    added = []
    for i, line in enumerate(lines_out):
        resname=None
        if line.startswith("ATOM") or line.startswith("HETATM"):
            resname = line[17:20]
        if resname=="HOH":
            data=line[30:]
            if data in added:
                continue
            added.append(data) 
            #line = line[:16]+water_altloc+line[17:]
        
        lines_out[i] = line

def separate_altloc_set(lines_out):
    # move to separate altloc set.
    original_altlocs=[]
    for line_out in lines_out:
        if not ((line_out.startswith("ATOM") or line_out.startswith("HETATM"))):
            continue
        altloc = line_out[16]
        if altloc not in original_altlocs:
            original_altlocs.append(altloc)
    allowed_altlocs = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    convert_dict={}


    # TODO change to move to new altloc if vdw nonbond issues 
    for a in original_altlocs:
        for c in allowed_altlocs:
            if c not in original_altlocs and c not in convert_dict.values():
                convert_dict[a]=c
                break
    for i, line_out in enumerate(lines_out):
        if not ((line_out.startswith("ATOM") or line_out.startswith("HETATM")) and line_out[17:20]=="HOH"):
            continue
        altloc = line_out[16]
        lines_out[i] = line_out[:16]+convert_dict[altloc]+line_out[17:]
        print(lines_out[i])