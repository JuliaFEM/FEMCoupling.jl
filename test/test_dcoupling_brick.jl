# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

using JuliaFEM
using JuliaFEM: add_elements!, Problem
using JuliaFEM.Preprocess
using JuliaFEM.Abaqus: create_surface_elements
using FEMBase
using FEMCoupling
using FEMCoupling: add_reference_node!, add_coupling_nodes!

# reading mesh from ABAQUS input file

datadir = Pkg.dir("FEMCoupling", "test", "test_dcoupling_brick")
mesh = abaqus_read_mesh(joinpath(datadir, "test_dcoupling_brick.inp"))
println("Number of nodes in a model: ", length(mesh.nodes))

# Elements

brick_body = create_elements(mesh,"Body1")
update!(brick_body, "youngs modulus", 2.080e5)
update!(brick_body, "poissons ratio", 0.3)
update!(brick_body, "density", 7.80e-9)

brick_problem = Problem(Elasticity,"brick_problem",3)
add_elements!(brick_problem, brick_body)

# Boundary conditions

bc_elements = [Element(Poi1, [j]) for j in mesh.node_sets[:fixed_face]]
update!(bc_elements, "geometry", mesh.nodes)

for i=1:3
    update!(bc_elements, "displacement $i", 0.0)
end

bc = Problem(Dirichlet, "fixed", 3, "displacement")
add_elements!(bc, bc_elements)

# Distributed coupling

coupling_nodes = [Element(Poi1, [j]) for j in mesh.node_sets[:coupling_nodes]]
update!(coupling_nodes, "geometry", mesh.nodes)
reference_node_id = collect(mesh.node_sets[:ref_node])
reference_node = Element(Poi1, reference_node_id)
update!(reference_node, "geometry", mesh.nodes)
update!(reference_node, "point moment 3", -10.0e5)

coupling = Problem(Coupling, "brick", 3, "displacement")
add_coupling_nodes!(coupling, coupling_nodes)
add_reference_node!(coupling, reference_node)

# Analysis

step = Analysis(Nonlinear)
add_problems!(step, [brick_problem, bc, coupling])
run!(step)

# Results

# Comparing calculated results with ABAQUS results. Nodes in the corners should
# have the maximum displacement magnitude. The node set corner_nodes contains
# nodes which are in the corners.

corner_node1 = first(mesh.node_sets[:corner_nodes])

# Declaring displacements of a corner node at time 0.0 to variable u. The
# expected value is result from ABAQUS.

time=0.0
u = brick_problem("displacement", time)[corner_node1]
u_mag = norm(u)
u_mag_expected=norm([-3.2785241E-03,3.2785241E-03,0.0000000E+00])

# Prints
println("reference node id = $(reference_node_id[1])")
println("corner node id = $corner_node1")
println("u_mag = $u_mag")
println("u_mag_expected = $u_mag_expected")
println("relative difference[%] = $((u_mag_expected-u_mag)/u_mag_expected*100)")
# Making a testset.

using FEMBase.Test
@testset "displacement magnitude and rotation" begin

@test isapprox(u_mag, u_mag_expected, rtol=3.0e-2)
end
