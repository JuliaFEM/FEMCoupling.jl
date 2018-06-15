# Analytical calculation
# T=1.500
# L=210e-3
# r=16e-3
# d=2r
# G=210e9/(2*(1+0.3))
#
# Wv=pi*d^3/16
# tmax=T/Wv
# a=2*L*tmax/(G*d)
#
# println(a)

using Logging
Logging.configure(level=DEBUG)

using JuliaFEM
using JuliaFEM: add_elements!, Problem
using JuliaFEM.Preprocess
using JuliaFEM.Abaqus: create_surface_elements
using FEMBase
using FEMCoupling
using FEMCoupling: add_reference_node!, add_coupling_nodes!

# read mesh
mesh = abaqus_read_mesh("cylindrical_beam_torsion1.inp") #in this .inp file the coupling element is removed

cylinder_element = create_elements(mesh,"Body1")
# update!(cylinder_element, "geometry", mesh.nodes)
update!(cylinder_element, "youngs modulus", 210e3)
update!(cylinder_element, "poissons ratio", 0.3)
update!(cylinder_element, "density", 7.80e-9)

cylinder_problem = Problem(Elasticity,"varsinolla",3)
add_elements!(cylinder_problem, cylinder_element)
################################################################################

# Boundary conditions
bc_elements = [Element(Poi1, [j]) for j in mesh.node_sets[:Fixed_face_set]]
update!(bc_elements, "geometry", mesh.nodes)
for i=1:3
    update!(bc_elements, "fixed displacement $i", 0.0)
end

bc = Problem(Dirichlet, "fixed", 3, "displacement")
add_elements!(bc, bc_elements)

################################################################################

# Coupling

# Creating Poi1 elements for coupling nodes and numbering
# them with original node numbers.
coupling_nodes = [Element(Poi1, [j]) for j in mesh.node_sets[:Coupling_nodes_set]]
update!(coupling_nodes, "geometry", mesh.nodes)

# Creating Poi1 element for the reference node and
# applying loads for it.
reference_node = Element(Poi1, [24142])
update!(reference_node, "geometry", mesh.nodes)
update!(reference_node, "point moment 3", 1500.0)

coupling = Problem(Coupling, "cylind", 3, "displacement")
add_coupling_nodes!(coupling, coupling_nodes)
add_reference_node!(coupling, reference_node)

################################################################################

step = Analysis(Nonlinear)
add_problems!(step, [cylinder_problem, bc, coupling])
run!(step)
