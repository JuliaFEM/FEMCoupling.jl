# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

# using JuliaFEM
using FEMCoupling

# Validation of kinematic coupling: we have 1x1x1 brick which is fixed at z=0
# Reference node is 1 away from the brick in z direction and twist is applied

# This test is decoupled from JuliaFEM by assembling stiffness matrix and
# saving it as a test resource. Stiffness matrix can be read back using
# readdlm.

X = Dict(
    1 => [0.0, 0.0, 0.0],
    2 => [1.0, 0.0, 0.0],
    3 => [1.0, 1.0, 0.0],
    4 => [0.0, 1.0, 0.0],
    5 => [0.0, 0.0, 1.0],
    6 => [1.0, 0.0, 1.0],
    7 => [1.0, 1.0, 1.0],
    8 => [0.0, 1.0, 1.0],
    9 => [0.5, 0.5, 2.0])

element = Element(Hex8, [1, 2, 3, 4, 5, 6, 7, 8])
update!(element, "geometry", X)
update!(element, "youngs modulus", 2.080e5)
update!(element, "poissons ratio", 0.3)
update!(element, "density", 7.80e-9)

# problem = Problem(Elasticity, "1 hex8 element", 3)
# add_elements!(problem, [element])
# assemble!(problem, 0.0)
# K = full(problem.assembly.K)
# fn = Pkg.dir("FEMCoupling", "test", "test_dcoupling_3d", "stiffness_matrix")
# writedlm(fn, K)

bc_element = Element(Quad4, [1, 2, 3, 4])
update!(bc_element, "geometry", X)
update!(bc_element, "displacement 1", 0.0)
update!(bc_element, "displacement 2", 0.0)
update!(bc_element, "displacement 3", 0.0)

# bc = Problem(Dirichlet, "fixed", 3, "displacement")
# add_elements!(bc, [bc_element])

# Distributed coupling

coupling_node_ids = [5, 6, 7, 8]
coupling_nodes = [Element(Poi1, [j]) for j in coupling_node_ids]
reference_node_id = 9
reference_node = Element(Poi1, [reference_node_id])
update!(coupling_nodes, "geometry", X)
update!(reference_node, "geometry", X)
update!(reference_node, "point moment 3", 1.0e3)

coupling = Problem(Coupling, "brick", 3, "displacement")
add_coupling_nodes!(coupling, coupling_nodes)
add_reference_node!(coupling, reference_node)

# Analysis

# step = Analysis(Nonlinear)
# add_problems!(step, [problem, bc, coupling])
# xdmf = Xdmf("3d_coupling_results"; overwrite=true)
# add_results_writer!(step, xdmf)
# run!(step)
# close(xdmf.hdf)

K = readdlm(Pkg.dir("FEMCoupling", "test", "test_dcoupling_3d", "stiffness_matrix"))
assemble!(coupling, 0.0)
f = full(sparse(coupling.assembly.f, 24, 1))
free_dofs = (4*3+1):(8*3)
u = zeros(3*8)
u[free_dofs] = K[free_dofs, free_dofs] \ f[free_dofs]
d = Dict(j => u[3*(j-1)+1:3*(j-1)+3] for j=1:8)

# Results

using Base.Test

# Because the twisting is symmetric, we should have equal magnitude of displacement
# in all free nodes
@test isapprox(norm(d[5]), norm(d[6]))
@test isapprox(norm(d[5]), norm(d[7]))
@test isapprox(norm(d[5]), norm(d[8]))

# Comparing calculated results with ABAQUS results. The model used to get
# results is given in directory test_coupling_3d/model.inp

# for i=5:8
#     println("node $i: $(u[i])")
# end

a = 3.75e-2
@test isapprox(d[5], [ a, -a, 0.0])
@test isapprox(d[6], [ a,  a, 0.0])
@test isapprox(d[7], [-a,  a, 0.0])
@test isapprox(d[8], [-a, -a, 0.0])
