# load structure
fetch 4C53

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
show cartoon, chain A
# show cartoon, chain B
# show cartoon, chain C

color GP2, chain A
color mygray, chain B
color mygray, chain C

show spheres, chain A and resi 311 # 317
show spheres, chain A and resi 315 # 321
show spheres, chain A and resi 316 # 322
show spheres, chain A and resi 333 # 339
show spheres, chain A and resi 346 # 352
show spheres, chain A and resi 365 # 371
show spheres, chain A and resi 378 # 384

color highlight, chain A and resi 311 # 317
color highlight, chain A and resi 315 # 321
color highlight, chain A and resi 316 # 322
color highlight, chain A and resi 333 # 339
color highlight, chain A and resi 346 # 352
color highlight, chain A and resi 365 # 371
color highlight, chain A and resi 378 # 384




# position protein for side view
orient
rotate z, -90
rotate y, 45

set surface_quality, 1
set ambient, 0.5
set antialias, 4
set ray_trace_mode, 1