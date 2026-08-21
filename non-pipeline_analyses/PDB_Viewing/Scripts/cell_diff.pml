# load structure to map minimum scores onto
cmd.load("../../../results/mapped_scores_onto_pdb/JS-div_9MT6_mapped_pdb.pdb")

# prep protein
hide everything
show surface, chain A
show surface, chain B
show surface, chain C
show surface, chain D
show surface, chain E
show surface, chain F
show surface, chain G
show surface, chain H
show surface, chain I

# set background to white
bg_color white

# functional score colors
set_color max_color, [0.0, 0.4470588235294118, 0.6980392156862745] # default blue from polyclonal
set_color min_color, [0.9019607843137255, 0.6235294117647059, 0.000000000] # default orange from polyclonal
set_color light_gray, [0.7, 0.7, 0.7] # dark gray
set_color dark_gray, [0.5, 0.5, 0.5] # light gray

# color protein by functional score minimum
spectrum b, white max_color, maximum=0.34, minimum=0

# color two subunits gray
color light_gray, chain A
color light_gray, chain E
color light_gray, chain F
color dark_gray, chain D
color dark_gray, chain H
color dark_gray, chain I

# position protein for side view
rotate x, 270
rotate x, 270
rotate z, 60
rotate x, -45

set surface_quality, 1
set ambient, 0.5
set antialias, 4
set ray_trace_mode, 1