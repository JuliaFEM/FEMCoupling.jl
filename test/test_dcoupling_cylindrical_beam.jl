# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE
using Revise
revise()



using JuliaFEM
using JuliaFEM: add_elements!, Problem
using JuliaFEM.Preprocess
using JuliaFEM.Abaqus: create_surface_elements
using FEMBase
using FEMCoupling
using FEMCoupling: add_reference_node!, add_coupling_nodes!


# read mesh from ABAQUS input file
datadir = Pkg.dir("FEMCoupling", "test", "test_dcoupling_cylindrical_beam")
mesh = abaqus_read_mesh(joinpath(datadir, "test_dcoupling_cylindrical_beam.inp"))
println("Number of nodes in a model: ", length(mesh.nodes))

################################################################################
# Creating body, cylindrical beam from mesh
cylinder_body = create_elements(mesh,"Body1")
# update!(cylinder_element, "geometry", mesh.nodes)
update!(cylinder_body, "youngs modulus", 210e3)
update!(cylinder_body, "poissons ratio", 0.3)
update!(cylinder_body, "density", 7.80e-9)

cylinder_problem = Problem(Elasticity,"varsinolla",3)
add_elements!(cylinder_problem, cylinder_body)
################################################################################

# Boundary conditions
bc_elements = [Element(Poi1, [j]) for j in mesh.node_sets[:Fixed_face_set]]
update!(bc_elements, "geometry", mesh.nodes)
for i=1:3
    update!(bc_elements, "displacement $i", 0.0)
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

reference_node_id = collect(mesh.node_sets[:ref_node_set])
reference_node = Element(Poi1, reference_node_id)
update!(reference_node, "geometry", mesh.nodes)
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
