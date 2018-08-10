# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

using Revise
using JuliaFEM
using Base.Test
using Logging
Logging.configure(level=DEBUG)
revise()
using FEMCoupling

printa(x) = show(IOContext(STDOUT, limit=true, displaysize=(1000,1000)), "text/plain", x)


# Creating 2d geometry for square shaped element and one independent node.
d = 3.0
X = Dict(4 => [0.0, 0.0],
         1 => [d, 0.0],
         2 => [d, d],
         5 => [0.0, d],
         3 => [2*d, d/2])

# Creating one plane_stress element
element1 = Element(Quad4, [4, 1, 2, 5])
update!(element1, "geometry", X)
update!(element1, "youngs modulus", 210.0e8)
update!(element1, "poissons ratio", 1/3)

# Creating Elasticity problem named "test problem" with 2 dimensions.
problem1 = Problem(Elasticity, "test problem", 2)
problem1.properties.formulation = :plane_stress
add_elements!(problem1, [element1])

# Creating Dirichlet problem for boundary conditions with Seg2
# element which fixes nodes 4 and 1.
bc = Problem(Dirichlet, "fixed", 2, "displacement")
element2 = Element(Seg2, [4, 5])
update!(element2, "geometry", X)
update!(element2, "displacement 1", 0.0)
update!(element2, "displacement 2", 0.0)
add_elements!(bc, [element2])

###############################################################################
# Creating coupling problem.
coupling1 = Problem(Coupling, "test", 2, "displacement")

# Creating Poi1 elements for coupling nodes and numbering
# them with original node numbers.
coupling_node2 = Element(Poi1, [1])
coupling_node3 = Element(Poi1, [2])
update!([coupling_node2, coupling_node3], "geometry", X)

# Creating Poi1 element for the reference node and
# applying loads for it.
reference_node = Element(Poi1, [3])
update!(reference_node, "geometry", X)
update!(reference_node, "fixed displacement 2", 0.032)

# Pointing out that coupling nodes and reference nodes
# are for the coupling problem.
add_coupling_nodes!(coupling1, [coupling_node2, coupling_node3])
add_reference_node!(coupling1, reference_node)

assemble!(coupling1, 0.0)
assemble!(problem1, 0.0)
K = full(problem1.assembly.K)[1:6, 1:6]
K += full(coupling1.assembly.K, 6, 6)
C1 = full(coupling1.assembly.C1, 6, 6)
C2 = full(coupling1.assembly.C2, 6, 6)
D = full(coupling1.assembly.D, 6, 6)
A = [K C1; C2 D]
used_rows =    [1, 2, 3, 4, 5, 6, 12]
used_columns = [1, 2, 3, 4, 5, 6, 12]
f = full(coupling1.assembly.f,6,1)
g = full(coupling1.assembly.g,6,1)
fg = [ f; g]

display([A[used_rows, used_columns] fg[used_rows]] )
# printa(full(coupling1.assembly.g))

printa([A fg])

K = A[used_rows,used_columns]
f = fg[used_rows]
u = K\f
display(u)


# Making step for nonlinear analysis and adding all three problems
# to the step. Running the analysis.
# step = Analysis(Nonlinear)
# add_problems!(step, [problem1, bc, coupling1])
# run!(step)

# Testing results by comparing them against abaqus results.
# Global force vector and plane stress element's displacements
# are compared.

# f = sparse(coupling1.assembly.f, 8, 1)
# f_expected = [0.0, 0.0, 6.66667e6, 1.75e6, -5.66667e6, 1.75e6, 0.0, 0.0]
# println("f = $f")
# println("f_expected = $f_expected")

# u = problem1("displacement", 0.0)
# println("u = $u")
