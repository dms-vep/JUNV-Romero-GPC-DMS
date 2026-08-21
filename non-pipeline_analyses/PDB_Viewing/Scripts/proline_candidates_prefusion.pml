# load structure
fetch 9MT6

# set background to white
cmd.bg_color("white")

# functional score colors
set_color max_color, [0.0, 0.4470588235294118, 0.6980392156862745] # default blue from polyclonal
set_color min_color, [0.9019607843137255, 0.6235294117647059, 0.000000000] # default orange from polyclonal
set_color max_pH_color, [0.30196078431372547, 0.6745098039215687, 0.14901960784313725] # green
set_color min_pH_color, [0.9058823, 0.1607843137254902, 0.5411764705882353] # pink
set_color mygray, [0.99, 0.99, 0.99]
set_color SSP, [0.667, 0.267, 0.6]
set_color GP1, [0.533, 0.8, 0.933]
set_color GP2, [0.933, 0.467, 0.2]
set_color highlight, [0.067, 0.467, 0.2]

hide everything
# show cartoon, chain A
# show cartoon, chain B
# show cartoon, chain C
show cartoon, chain D
# show cartoon, chain E
# show cartoon, chain F
# show cartoon, chain G
show cartoon, chain H
show cartoon, chain I

color SSP, chain D
color mygray, chain G
color mygray, chain A

color mygray, chain B
color mygray, chain E
color GP1, chain H

color GP2, chain I
color mygray, chain C
color mygray, chain F

show spheres, chain I and resi 317
show spheres, chain I and resi 321
show spheres, chain I and resi 322
show spheres, chain I and resi 339
show spheres, chain I and resi 352
show spheres, chain I and resi 371
show spheres, chain I and resi 384

color highlight, chain I and resi 317
color highlight, chain I and resi 321
color highlight, chain I and resi 322
color highlight, chain I and resi 339
color highlight, chain I and resi 352
color highlight, chain I and resi 371
color highlight, chain I and resi 384




# position protein for side view
rotate x, 90
rotate y, 55
rotate y, -120

set surface_quality, 1
set ambient, 0.5
set antialias, 4
set ray_trace_mode, 1