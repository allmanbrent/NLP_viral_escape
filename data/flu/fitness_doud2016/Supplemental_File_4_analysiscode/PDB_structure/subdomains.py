'''According to Gamblin et al, domain numbers in WSN renumbering:

HA1:
F' fusion: D18-C72, C291-R340
Vestigial Esterase: N73-V125, S279-E290
Receptor Binding: S126-G278
HA2:
F fusion: G344-P503'''


renumbered_pdb_file = '1RVX_trimer_sequentialnumbering.pdb'
cmd.load(renumbered_pdb_file, '1rvx_trimer')

cmd.select('HA1fusion', 'chain A and resi 18-72+291-343')
cmd.select('HA1ved', 'chain A and resi 73-125+279-290')
cmd.select('HA1rbd', 'chain A and resi 126-278')
cmd.select('HA2fusion', 'chain A and resi 344-503')
cmd.set_view('\
     0.428075522,    0.056174517,   -0.902075231,\
     0.890969813,   -0.193928584,    0.410726875,\
    -0.151850626,   -0.979421794,   -0.133035392,\
     0.006424522,    0.001244783, -368.391510010,\
    74.595397949,   12.513990402,   21.364448547,\
   330.869323730,  404.743347168,  -20.000000000')

cmd.color('purple', 'HA1fusion')
cmd.color('orange', 'HA1ved')
cmd.color('marine', 'HA1rbd')
cmd.color('tv_red', 'HA2fusion')

cmd.hide('everything')
cmd.show('cartoon', '1RVX_trimer and chain A')

cmd.set('bg_rgb','[1,1,1]') # white
cmd.set('antialias','2')
cmd.set('ray_opaque_background','off')
cmd.select('None')

png_width = 1600
png_height = 1200

cmd.ray(png_width,png_height)
cmd.png('subdomains.png')