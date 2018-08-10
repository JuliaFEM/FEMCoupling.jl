# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

using JuliaFEM
using FEMCoupling

printa(x) = show(IOContext(STDOUT, limit=true, displaysize=(1000,1000)), "text/plain", x)

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

refe_node = Element(Poi1, [9])
update!(refe_node, "geomertry", X)


problem = Problem(Elasticity, "1 hex8 element", 3)
add_elements!(problem, [element,refe_node])
assemble!(problem, 0.0)
K = full(problem.assembly.K, 27,27)
# fn = Pkg.dir("FEMCoupling", "test", "test_dcoupling_3d", "stiffness_matrix")
# writedlm(fn, K)

# Distributed coupling

coupling_node_ids = [5, 6, 7, 8]
coupling_nodes = [Element(Poi1, [j]) for j in coupling_node_ids]
reference_node_id = 9
reference_node = Element(Poi1, [reference_node_id])
update!(coupling_nodes, "geometry", X)
update!(reference_node, "geometry", X)
update!(reference_node, "fixed displacement 3", 0.042)

coupling = Problem(Coupling, "brick", 3, "displacement")
add_coupling_nodes!(coupling, coupling_nodes)
add_reference_node!(coupling, reference_node)
assemble!(coupling, 0.0)

K += full(coupling.assembly.K)
K = K[13:end, 13:end]

C1 = full(coupling.assembly.C1, 27, 27)
C1 = C1[13:end, end]

C2 = full(coupling.assembly.C2, 27, 27)
C2 = C2[end:end, 13:end]

D = full(coupling.assembly.D,   27, 27)
D = D[27, 27]

f = full(coupling.assembly.f, 27, 1)
f = f[13:end, 1]

g = full(coupling.assembly.g, 27, 1)
g = g[end, 1]

KK = [K C1; C2 D]
ff = [f ; g]
u = KK\ff
# Analysis

# step = Analysis(Nonlinear)
# add_problems!(step, [problem, bc, coupling])
# xdmf = Xdmf("3d_coupling_results"; overwrite=true)
# add_results_writer!(step, xdmf)
# run!(step)
# close(xdmf.hdf)
