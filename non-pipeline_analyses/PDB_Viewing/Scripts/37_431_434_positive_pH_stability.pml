# load structure to map scores onto
cmd.load("../../../results/mapped_scores_onto_pdb/positive_pH_9MT6_mapped_pdb.pdb")
scores="positive_pH_9MT6_mapped_pdb"


# set background to white
cmd.bg_color("white")

# functional score colors
set_color max_color, [0.0, 0.4470588235294118, 0.6980392156862745] # default blue from polyclonal
set_color min_color, [0.9019607843137255, 0.6235294117647059, 0.000000000] # default orange from polyclonal
set_color max_pH_color, [0.30196078431372547, 0.6745098039215687, 0.14901960784313725] # green
set_color min_pH_color, [0.9058823, 0.1607843137254902, 0.5411764705882353] # pink
set_color mygray, [0.9, 0.9, 0.9]

# get the actual min/max b-factor values present on the object
stored.b_values = []
cmd.iterate(scores, "stored.b_values.append(b)")
min_b = min(stored.b_values)
max_b = max(stored.b_values)

# positive half: mygray -> max_pH_color, anchored so 0 = mygray
cmd.spectrum("b", "mygray max_pH_color", scores, minimum=0, maximum=max_b)

hide everything
cmd.show("cartoon", scores + " and chain A")
# cmd.show("cartoon", scores + " and chain B")
cmd.show("cartoon", scores + " and chain C")
cmd.show("cartoon", scores + " and chain D")
# cmd.show("cartoon", scores + " and chain E")
cmd.show("cartoon", scores + " and chain F")
cmd.show("cartoon", scores + " and chain G")
# cmd.show("cartoon", scores + " and chain H")
cmd.show("cartoon", scores + " and chain I")

show sticks, resi 37
show sticks, resi 431
show sticks, resi 434

# position protein for side view
rotate x, 90
rotate y, 90
# rotate x, 90

zoom resi 37+431

set surface_quality, 1
set ambient, 0.5
set antialias, 4
set ray_trace_mode, 1