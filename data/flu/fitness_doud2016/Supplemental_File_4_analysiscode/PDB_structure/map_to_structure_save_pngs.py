'''visualize site entropy of amino-acid preferences on the
1RVX pdb structure of 1934/PR8 HA.

run within pymol: `run map_to_structure.py` '''

import re

def ReadPreferences(f):
    """
    Copied from dms_tools.file_io because pymol's python search path
    might not find dms_tools for import.

    Reads the site-specific preferences written by *WritePreferences*.
    *f* is the name of an existing file or a readable file-like object.
    The return value is the tuple: *(sites, wts, pi_means, pi_95credint, h)*
    where *sites*, *wts*, *pi_means*, and *pi_95credint* will all
    have the same values used to write the file with *WritePreferences*,
    and *h* is a dictionary with *h[r]* giving the site entropy (log base
    2) for each *r* in *sites*.
    See docstring of *WritePreferences* for example usage.
    """
    charmatch = re.compile('^PI_([A-z\*\-]+)$')
    if isinstance(f, str):
        f = open(f)
        lines = f.readlines()
        f.close()
    else:
        lines = f.readlines()
    characters = []
    sites = []
    wts = {}
    pi_means = {}
    pi_95credint = {}
    h = {}
    for line in lines:
        if line.isspace():
            continue
        elif line[0] == '#' and not characters:
            entries = line[1 : ].strip().split()
            if len(entries) < 4:
                raise ValueError("Insufficient entries in header:\n%s" % line)
            if not (entries[0] in ['POSITION', 'SITE'] and entries[1][ : 2] == 'WT' and entries[2] == 'SITE_ENTROPY'):
                raise ValueError("Not the correct first three header columns:\n%s" % line)
            i = 3
            while i < len(entries) and charmatch.search(entries[i]):
                characters.append(charmatch.search(entries[i]).group(1))
                i += 1
            if i  == len(entries):
                pi_95credint = None
                linelength = len(characters) + 3
            else:
                if not len(entries) - i == len(characters):
                    raise ValueError("Header line does not have valid credible interval format:\n%s" % line)
                if not all([entries[i + j] == 'PI_%s_95' % characters[j] for j in range(len(characters))]):
                    raise ValueError("mean and credible interval character mismatch in header:\n%s" % line)
                linelength = 2 * len(characters) + 3
        elif line[0] == '#':
            continue
        elif not characters:
            raise ValueError("Found data lines before encountering a valid header")
        else:
            entries = line.strip().split()
            if len(entries) != linelength:
                raise ValueError("Line does not have expected %d entries:\n%s" % (linelength, line))
            r = entries[0]
            assert r not in sites, "Duplicate site of %s" % r
            sites.append(r)
            wts[r] = entries[1]
            assert entries[1] in characters or entries[1] == '?', "Character %s is not one of the valid ones in header. Valid possibilities: %s" % (entries[1], ', '.join(characters))
            h[r] = float(entries[2])
            pi_means[r] = dict([(x, float(entries[3 + i])) for (i, x) in enumerate(characters)])
            if pi_95credint != None:
                pi_95credint[r] = dict([(x, (float(entries[3 + len(characters) + i].split(',')[0]), float(entries[3 + len(characters) + i].split(',')[1]))) for (i, x) in enumerate(characters)])
    return (sites, wts, pi_means, pi_95credint, h)

def GetSites(sites_file):
    lines = [line for line in open(sites_file).readlines() if (not line.isspace()) and line[0] != '#']
    selectedsites = [str(int(line.split()[0])) for line in lines]
    return selectedsites

def GetEntropiesFromPreferencesFile(prefs):
    (sites, wts, pi_means, pi_95credint, h) = ReadPreferences(prefs)
    return h

def main():

    preferences_file = 'Overall-WSNHA_merged_prefs_rescaled.txt' 
    renumbered_pdb_file = '1RVX_trimer_sequentialnumbering.pdb'
    original_pdb_file = '1rvx.pdb'
    entropy = GetEntropiesFromPreferencesFile(preferences_file)
    sites = entropy.keys()
    antigenic_sites = GetSites('Caton_H1_HA_antigenic_sites.txt')
    antigeniccontact_sites = GetSites('nearby_antigenic_sites.txt')
    all_rbs_sites = GetSites('allRBS_residues.txt')
    conserved_rbs_sites = GetSites('receptor_binding_residues.txt')
    save_pngs = True

    # load file and map entropies to b-factors
    cmd.reinitialize()
    cmd.load(renumbered_pdb_file, '1rvx_trimer')
    for r in sites:
        cmd.alter("1rvx_trimer and resi %d" % int(r), "b = %g" % entropy[r])

    # copy the sialic acid receptor
    cmd.load(original_pdb_file, 'original_pdb')
    #cmd.select('all_sugars', 'original_pdb and chain A and resn GAL+NAG+SIA')
    cmd.select('receptor', 'original_pdb and chain A and resi 3021 and resn SIA')
    # make a copy of the receptor and delete the rest of the original PDB to avoid confusion with multiple numbering schemes.
    cmd.create('sialicacid', 'receptor')
    cmd.delete('original_pdb')

    # make site selections
    cmd.select('antigenic', 'resi %s' % '+'.join([r for r in antigenic_sites]))
    cmd.select('notantigenic', 'not resi %s' % '+'.join([r for r in antigenic_sites]))
    cmd.select('antigeniccontact', 'resi %s' % '+'.join([r for r in antigeniccontact_sites]))
    cmd.select('notantigenicandcontact', 'not resi %s' % '+'.join([r for r in antigeniccontact_sites]))
    cmd.select('RBS', 'resi %s' % '+'.join([r for r in all_rbs_sites]))
    cmd.select('notRBS', 'not resi %s' % '+'.join([r for r in all_rbs_sites]))

    # residues in Koel et al identified as antigenic cluster jumps in either H3 or H1. (most H3, only one is H1.)
    cmd.select('koelH3andH1clusterjumps', 'resi 202+206+169+172+171+168+158+157')
    cmd.select('notkoelH3andH1clusterjumps', 'not resi 202+206+169+172+171+168+158+157')
    # F10 epitope: residues in Sui et al colored red in figure 4c
    cmd.select('sui_F10_epitope','resi 25+45+46+47+305+306+361+362+363+364+381+384+385+388+392+395+396+399')
    cmd.select('notsui_F10_epitope','not sui_F10_epitope')
    # CR6261 epitope: residues in Ekiert et al colored darker shades in figure 3b
    cmd.select('eckiert_CR6261_epitope', 'resi 25+45+47+48+49+305+306+362+363+364+384+385+388+389+392+395+396+399')
    cmd.select('not_CR6261_epitope', 'not eckiert_CR6261_epitope')
    # residues in HA within 4.5A of FI6v3 in pdb 3ZTN, identified using 
    # the script at http://www.blopig.com/blog/2013/10/get-pdb-intermolecular-protein-contacts-and-interface-residues/
    # and converted to WSN preference numbering manually
    cmd.select('FI6v3contacts', 'resi 45+46+47+305+332+361+363+364+381+382+384+385+386+388+389+391+392+396+399+400')
    cmd.select('notI6v3contacts', 'not FI6v3contacts')
    # CR9114 epitope: residues labeled in Dreyfus et al figure 3d
    cmd.select('dreyfusCR9114epitope', 'resi 45+47+48+49+305+306+307+361+362+363+364+379+381+384+385+388+389+391+392+395+399')
    cmd.select('notCR9114epitope', 'not dreyfusCR9114epitope')
    # union of epitope footprints for these three antibodies
    cmd.select('stem_epitopes_union', 'sui_F10_epitope or eckiert_CR6261_epitope or FI6v3contacts or dreyfusCR9114epitope')
    cmd.select('not_stemepitopes', 'not stem_epitopes_union')
    # all antigenic sites (caton+contacting, koel cluster jumps, stem epitopes)
    cmd.select('all_antigen', 'antigeniccontact or koelH3andH1clusterjumps or stem_epitopes_union')
    cmd.select('not_allantigen', 'not all_antigen')

    # select the subdomains on chain B
    cmd.select('HA1fusion', 'chain B and resi 18-72+291-340')
    cmd.select('HA1ved', 'chain B and resi 73-125+279-290')
    cmd.select('HA1rbd', 'chain B and resi 126-278')
    cmd.select('HA2fusion', 'chain B and resi 344-503')

    # display settings
    cmd.hide('everything')
    cmd.show('surface', '1rvx_trimer')
    cmd.show('sticks', 'sialicacid')
    cmd.color('marine','sialicacid')
    cmd.set('bg_rgb','[1,1,1]') # white
    cmd.set('antialias','2')
    cmd.set('ray_opaque_background','off')
    cmd.deselect()
    png_width = 1600
    png_height = 1200

    # color by entropy
    cmd.spectrum('b', 'yellow_red', '1rvx_trimer')

    cmd.set_view('\
     0.232918754,    0.130249694,   -0.963816822,\
     0.972345531,   -0.052188061,    0.227946103,\
    -0.020623561,   -0.990137696,   -0.138796180,\
     0.001965282,    0.000245631, -411.707763672,\
    96.332313538,    3.910612106,   22.706035614,\
   315.693145752,  494.978698730,  -20.000000000')

    if save_pngs:
        cmd.ray(png_width,png_height)
        cmd.png('allsites_sideview.png')

    # also show with only one monomer colored by entropy
    cmd.color('gray60', 'chain B')
    cmd.color('gray30', 'chain C')

    cmd.set_view('\
     0.232918754,    0.130249694,   -0.963816822,\
     0.972345531,   -0.052188061,    0.227946103,\
    -0.020623561,   -0.990137696,   -0.138796180,\
     0.001965282,    0.000245631, -411.707763672,\
    96.332313538,    3.910612106,   22.706035614,\
   315.693145752,  494.978698730,  -20.000000000')

    if save_pngs:
        cmd.ray(png_width,png_height)
        cmd.png('allsites_monomer_sideview.png')

    for sel in ['notRBS', 'notantigenic', 'notantigenicandcontact', 'notkoelH3andH1clusterjumps', 'not_stemepitopes', 'not_allantigen']:

        file_label = sel[3:].replace('_', '')
        
        # re-spectrum whole molecule; 
        # (applying spectrum to individual selections results in different 
        # colormappings for each selection's entropy distribution - don't want that)
        cmd.spectrum('b', 'yellow_red', '1rvx_trimer')
        # color gray, by monomer, all sites not in the selection of interest
        cmd.color('gray90', '%s and chain A' % sel)
        cmd.color('gray60', '%s and chain B' % sel)
        cmd.color('gray30', '%s and chain C' % sel)
        cmd.color('marine','sialicacid')

        cmd.set_view('\
             0.083034605,    0.438198835,   -0.895136297,\
             0.995836735,   -0.000724004,    0.092045061,\
             0.039661985,   -0.898927331,   -0.436385781,\
            -0.001473397,    0.008274943, -247.293945312,\
           114.436874390,   -5.348280907,   -2.088230133,\
           141.187042236,  320.472686768,  -20.000000000')

        if save_pngs:
            cmd.ray(png_width,png_height)
            cmd.png('%s_headzoom.png'%file_label)

        if sel == 'not_stemepitopes':
            cmd.set_view('\
             0.918532729,    0.017010566,   -0.395212889,\
             0.395452172,   -0.062295828,    0.916462004,\
            -0.009029277,   -0.997950077,   -0.063960060,\
             0.012179107,   -0.004756615, -249.610641479,\
            86.267158508,  -36.820129395,   54.390720367,\
           141.187042236,  320.472686768,  -20.000000000')

            if save_pngs:
                cmd.ray(png_width,png_height)
                cmd.png('%s_stemzoom.png'%file_label)
        else:
            cmd.set_view('\
             0.232918754,    0.130249694,   -0.963816822,\
             0.972345531,   -0.052188061,    0.227946103,\
            -0.020623561,   -0.990137696,   -0.138796180,\
             0.001965282,    0.000245631, -411.707763672,\
            96.332313538,    3.910612106,   22.706035614,\
           315.693145752,  494.978698730,  -20.000000000')

            if save_pngs:
                cmd.ray(png_width,png_height)
                cmd.png('%s_sideview.png'%file_label)

main()
