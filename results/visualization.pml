# PyMol script to visualize the best axis
cmd.load('../results/2moz_new.pdb')

# Create two pseudoatoms representing the best axis
cmd.pseudoatom('pt1', pos=[-0.21066172265350314, -9.338719080440962, -27.586642506581143])
cmd.pseudoatom('pt2', pos=[0.21066172265350314, 9.338719080440962, 27.586642506581143])

# Draw a dashed line between the two points
cmd.distance('axis', 'pt1', 'pt2')

# Customize the appearance of the line
cmd.set('dash_gap', '0')
cmd.set('dash_radius', '0.3')
cmd.set('dash_round_ends', '0')
cmd.set('dash_color', 'yellow', 'axis')
cmd.hide('labels', 'axis')

# Show spheres at both points representing the line endpoints
cmd.show('spheres', 'pt1')
cmd.show('spheres', 'pt2')

# Set the sphere size for both points
cmd.set('sphere_scale', '0.1', 'pt1')
cmd.set('sphere_scale', '0.1', 'pt2')
