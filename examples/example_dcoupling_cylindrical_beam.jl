# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

# # Cylinder in torsion using distributed coupling

# For general information about anything, see
# [this](https://en.wikipedia.org)
# wikipedia page.

# The model is a 3d cylinder, shown in picture.

# ![](example_dcoupling_cylindrical_beam/example_dcoupling_cylindrical_beam.png)

# Node sets for fixed end boundary conditions, coupling nodes and reference node
# were made to the ABAQUS input file.

using JuliaFEM
using JuliaFEM: add_elements!, Problem
using JuliaFEM.Preprocess
using JuliaFEM.Abaqus: create_surface_elements
using FEMBase
using FEMCoupling
using FEMCoupling: add_reference_node!, add_coupling_nodes!

# reading mesh from ABAQUS input file
datadir = Pkg.dir("FEMCoupling", "examples", "example_dcoupling_cylindrical_beam")
mesh = abaqus_read_mesh(joinpath(datadir, "example_dcoupling_cylindrical_beam.inp"))
println("Number of nodes in a model: ", length(mesh.nodes))

################################################################################
# Creating elements for the whole body. In the ABAQUS input file the cylinder is
# named "Body1".
cylinder_body = create_elements(mesh,"Body1")
# Updating values for the body elements.
update!(cylinder_body, "youngs modulus", 210e3)
update!(cylinder_body, "poissons ratio", 0.3)
update!(cylinder_body, "density", 7.80e-9)

# Creating an elasticity problem and adding the elements to it.
cylinder_problem = Problem(Elasticity,"varsinolla",3)
add_elements!(cylinder_problem, cylinder_body)
################################################################################

# Boundary conditions: fixed from the other end
# Creating Poi1-type elements as boundary condition elements to nodes of the
# node set Fixed_face_set.
bc_elements = [Element(Poi1, [j]) for j in mesh.node_sets[:Fixed_face_set]]
# Updating geometry for the bc_elements
update!(bc_elements, "geometry", mesh.nodes)
# Fixing all displacements for the bc elements.
for i=1:3
    update!(bc_elements, "displacement $i", 0.0)
end

# Creating a bc problem and adding the bc elements to it.
bc = Problem(Dirichlet, "fixed", 3, "displacement")
add_elements!(bc, bc_elements)

################################################################################

# Coupling

# Creating Poi1 elements to nodes in coupling nodes set.
coupling_nodes = [Element(Poi1, [j]) for j in mesh.node_sets[:Coupling_nodes_set]]
# Updating geometry for the coupling nodes.
update!(coupling_nodes, "geometry", mesh.nodes)

# Creating Poi1 element for the reference node.
reference_node_id = collect(mesh.node_sets[:ref_node_set])
reference_node = Element(Poi1, reference_node_id)
# Updating geometry for the reference node
update!(reference_node, "geometry", mesh.nodes)
#  Applying a point moment for the reference node. 
update!(reference_node, "point moment 3", 1500.0)

coupling = Problem(Coupling, "cylind", 3, "displacement")
add_coupling_nodes!(coupling, coupling_nodes)
add_reference_node!(coupling, reference_node)

################################################################################

step = Analysis(Nonlinear)
add_problems!(step, [cylinder_problem, bc, coupling])
run!(step)

################################################################################
# Tests
node_on_circle = first(mesh.node_sets[:circlenodes_set])
u = cylinder_problem("displacement", 0.0)[node_on_circle]
u_mag = norm(u)

using FEMBase.Test
@testset "displacement magnitude and rotation" begin

u_mag_expected=6.306e-4

@test isapprox(u_mag, u_mag_expected, rtol=1e-3)

radius = 16
a_expected = 6.306e-4/radius # from ABAQUS
a = u_mag/radius

@test isapprox(a,a_expected,rtol=1e-3)
end

println("reference node id = $(reference_node_id[1])")
println("node on circle id = $node_on_circle")
