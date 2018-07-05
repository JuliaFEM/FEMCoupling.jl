# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

# using JuliaFEM
using FEMCoupling
using FEMBase.Test

# This test is decoupled from JuliaFEM by assembling stiffness matrix and
# saving it as a test resource. Stiffness matrix can be read back using
# readdlm.

# Creating 2d geometry for square shaped element and one independent node.
d = 3.0
X = Dict(1 => [0.0, 0.0],
         2 => [d, 0.0],
         3 => [d, d],
         4 => [0.0, d],
         5 => [2*d, d/2])

# Creating one plane_stress element
element1 = Element(Quad4, [1, 2, 3, 4])
update!(element1, "geometry", X)
update!(element1, "youngs modulus", 210.0e8)
update!(element1, "poissons ratio", 1/3)

# Creating Elasticity problem named "test problem" with 2 dimensions.
# problem1 = Problem(Elasticity, "test problem", 2)
# problem1.properties.formulation = :plane_stress
# add_elements!(problem1, [element1])
# assemble!(problem1, 0.0)
# K = full(problem1.assembly.K)
# fn = Pkg.dir("FEMCoupling", "test", "test_dcoupling_2d", "stiffness_matrix")
# writedlm(fn, K)

# Creating Dirichlet problem for boundary conditions with Seg2
# element which fixes nodes 4 and 1.
element2 = Element(Seg2, [4, 1])
update!(element2, "geometry", X)
update!(element2, "displacement 1", 0.0)
update!(element2, "displacement 2", 0.0)
# bc = Problem(Dirichlet, "fixed", 2, "displacement")
# add_elements!(bc, [element2])

###############################################################################
# Creating coupling problem.
coupling1 = Problem(Coupling, "test", 2, "displacement")

# Creating Poi1 elements for coupling nodes and numbering
# them with original node numbers.
coupling_node2 = Element(Poi1, [2])
coupling_node3 = Element(Poi1, [3])
update!([coupling_node2, coupling_node3], "geometry",X)

# Creating Poi1 element for the reference node and
# applying loads for it.
reference_node = Element(Poi1, [5])
update!(reference_node, "geometry", X)
update!(reference_node, "point force 1", 10.0e5)
update!(reference_node, "point force 2", 35.0e5)
update!(reference_node, "point moment 3", 80.0e5)

# Pointing out that coupling nodes and reference nodes
# are for the coupling problem. 
add_coupling_nodes!(coupling1, [coupling_node2, coupling_node3])
add_reference_node!(coupling1, reference_node)

# Making step for nonlinear analysis and adding all three problems
# to the step. Running the analysis.
# step = Analysis(Nonlinear)
# add_problems!(step, [problem1, bc, coupling1])
# run!(step)

# Testing results by comparing them against abaqus results.
# Global force vector and plane stress element's displacements
# are compared.

assemble!(coupling1, 0.0)
f = full(sparse(coupling1.assembly.f, 8, 1))
f_expected = [0.0, 0.0, 6.66667e6, 1.75e6, -5.66667e6, 1.75e6, 0.0, 0.0]
# println("f = $f")
# println("f_expected = $f_expected")
@test isapprox(f, f_expected, rtol=1.0e-3)

K = readdlm(Pkg.dir("FEMCoupling", "test", "test_dcoupling_2d", "stiffness_matrix"))
free_dofs = [3, 4, 5, 6]
u = zeros(8)
u[free_dofs] = K[free_dofs, free_dofs] \ f[free_dofs]
u2 = u[3:4]
u3 = u[5:6]
u2_expected = [0.00155379, 0.00196296]
u3_expected = [-0.00146208, 0.0019418]
@test isapprox(u2, u2_expected, rtol=1.0e-5)
@test isapprox(u3, u3_expected, rtol=1.0e-5)

# If reference node is not defined, error is thrown:

coupling2 = Problem(Coupling, "test coupling", 2, "displacement")
add_coupling_nodes!(coupling2, [coupling_node2, coupling_node3])
@test_throws Exception assemble!(coupling2, 0.0)

