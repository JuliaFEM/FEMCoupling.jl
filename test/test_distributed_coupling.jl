# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE
using Revise
revise()

using JuliaFEM
using FEMBase
using Base.Test
using FEMCoupling: add_coupling_nodes!,add_reference_node!

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
problem = Problem(Elasticity, "test problem", 2)
problem.properties.formulation = :plane_stress
add_elements!(problem, [element1])

# Creating Dirichlet problem for boundary conditions with Seg2
# element which fixes nodes 4 and 1.
bc = Problem(Dirichlet, "fixed", 2, "displacement")
element2 = Element(Seg2, [4, 1])
update!(element2, "geometry", X)
update!(element2, "displacement 1", 0.0)
update!(element2, "displacement 2", 0.0)
add_elements!(bc, [element2])

# Creating coupling problem.
coupling = Problem(Coupling, "test", 2, "displacement")

# Creating Poi1 elements for coupling nodes and numbering
# them with original node numbers.
coupling_node2 = Element(Poi1, [2])
coupling_node3 = Element(Poi1, [3])
update!([coupling_node2,coupling_node3], "geometry",X)

# Creating Poi1 element for the reference node and
# applying loads for it.
reference_node = Element(Poi1, [5])
update!(reference_node, "geometry", X)
update!(reference_node, "point force 1", 10.0e5)
update!(reference_node, "point force 2", 35.0e5)
update!(reference_node, "point moment 3", 80.0e5)

# Pointing out that coupling nodes and reference nodes
# are for the coupling problem.
add_coupling_nodes!(coupling, [coupling_node2,coupling_node3])
add_reference_node!(coupling, reference_node)

# Making step for nonlinear analysis and adding all three problems
# to the step. Running the analysis.
step = Analysis(Nonlinear)
add_problems!(step, [problem, bc, coupling])
run!(step)


# Testing results by comparing them against abaqus results.
# Global force vector and plane stress element's displacements
# are compared.
@testset "force vector and displacements" begin
f = [full(coupling.assembly.f);0;0]
f_expected = [0.0; 0.0; 6.66667e6; 1.75e6; -5.66667e6; 1.75e6; 0.0; 0.0]
@test isapprox(f,f_expected,rtol=0.001)

time = 0.0 # In this case we want to know the displacements at t=0
u = element1("displacement", time)
u1,u2,u3,u4 = u
u_expected = [-4.90255e-19, 3.6124e-19, 0.00155379, 0.00196296, -0.00146208, 0.0019418, 0.0, 0.0]
@test isapprox(u2, u_expected[3:4], rtol = 0.001)
@test isapprox(u3, u_expected[5:6], rtol = 0.001)
end
