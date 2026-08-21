# load structure to map minimum scores onto
cmd.load("../../../results/mapped_scores_onto_pdb/C_musculinus_func_9MT6_mapped_pdb.pdb")
scores="C_musculinus_func_9MT6_mapped_pdb"


# set background to white
cmd.bg_color("white")

# functional score colors
set_color max_color, [0.0, 0.4470588235294118, 0.6980392156862745] # default blue from polyclonal
set_color min_color, [0.9019607843137255, 0.6235294117647059, 0.000000000] # default orange from polyclonal
set_color max_pH_color, [0.30196078431372547, 0.6745098039215687, 0.14901960784313725] # green
set_color min_pH_color, [0.9058823, 0.1607843137254902, 0.5411764705882353] # pink
set_color mygray, [0.9, 0.9, 0.9]
set_color light_gray, [0.7, 0.7, 0.7] # dark gray
set_color dark_gray, [0.5, 0.5, 0.5] # light gray

# get the actual min/max b-factor values present on the object
stored.b_values = []
cmd.iterate(scores, "stored.b_values.append(b)")
min_b = min(stored.b_values)
max_b = max(stored.b_values)

# base color for b == 0 (and anything untouched)
cmd.color("mygray", scores)

# negative half: min_pH_color -> mygray, anchored so 0 = mygray
cmd.spectrum("b", "min_color white", scores + " and b<0", minimum=min_b, maximum=0)

# positive half: mygray -> max_pH_color, anchored so 0 = mygray
cmd.spectrum("b", "white max_color", scores + " and b>0", minimum=0, maximum=1)

# show
hide everything
cmd.show("surface", scores + " and chain A")
cmd.show("cartoon", scores + " and chain B")
cmd.show("cartoon", scores + " and chain C")
cmd.show("surface", scores + " and chain D")
cmd.show("surface", scores + " and chain E")
cmd.show("surface", scores + " and chain F")
cmd.show("cartoon", scores + " and chain G")
cmd.show("surface", scores + " and chain H")
cmd.show("surface", scores + " and chain I")


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
# rotate z, 45
rotate x, -90
rotate y, -15


set surface_quality, 1
set ambient, 0.5
set antialias, 4
set ray_trace_mode, 1

zoom resi 269