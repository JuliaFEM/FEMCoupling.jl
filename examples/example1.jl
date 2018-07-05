# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

using JuliaFEM
using FEMBase
using Logging
Logging.configure(level=DEBUG)

type Coupling <: BoundaryProblem
    reference_node :: Element{Poi1}
end

function Coupling()
    reference_node = Element(Poi1, Int64[])
    return Coupling(reference_node)
end

d = 3.0
X = Dict(1 => [0.0, 0.0],
         2 => [d, 0.0],
         3 => [d, d],
         4 => [0.0, d],
         5 => [3*d, 1/2*d])
element1 = Element(Quad4, [1, 2, 3, 4])
update!(element1, "geometry", X)
update!(element1, "youngs modulus", 288.0)
update!(element1, "poissons ratio", 1/3)

problem = Problem(Elasticity, "test problem", 2)
problem.properties.formulation = :plane_strain
add_elements!(problem, [element1])

bc = Problem(Dirichlet, "fixed", 2, "displacement")
element2 = Element(Seg2, [4, 1])
update!(element2, "geometry", X)
update!(element2, "displacement 1", 0.0)
update!(element2, "displacement 2", 0.0)
add_elements!(bc, [element2])

"""
    add_coupling_nodes!(problem, elements)

Add new coupling nodes into the problem. Nodes must be defined with Poi1
elements.
"""
function add_coupling_nodes!(problem, elements)
    for element in elements
        push!(problem.elements, element)
    end
    return nothing
end

"""
    add_reference_node!(problem, element)

Add new reference node into the problem. Node must be defined with a Poi1
element. If reference node is already defined, it will be replaced with new one.
"""
function add_reference_node!(problem, element)
    problem.properties.reference_node = element
    return nothing
end

function FEMBase.assemble_elements!(problem::Problem{Coupling},
                                    assembly::Assembly,
                                    elements::Vector{Element{Poi1}}, time)

    ref_node = problem.properties.reference_node
    if length(get_connectivity(ref_node)) == 0
        error("Reference node node defined. Define reference node using add_reference_node!")
    end
    ref_node_id = first(get_connectivity(ref_node))
    info("Reference node id = $ref_node_id")
    X_r = first(ref_node("geometry", time))
    info("Reference node geometry = $X_r")
    fe = zeros(2)
    for coupling_node in elements
        fill!(fe, 0.0)
        coupling_node_id = first(get_connectivity(coupling_node))
        X_n = first(coupling_node("geometry", time))
        info("Coupling node id = $coupling_node_id, geometry = $X_n")
        gdofs = get_gdofs(problem, coupling_node)
        info("gdofs = $gdofs")
        if haskey(ref_node, "point moment 3")
            M = ref_node("point moment 3")
            info("calculate point moment")
            # fe += ...
        end
        if haskey(ref_node, "point load 1")
            M = ref_node("point load 1")
            info("calculate point moment")
            # fe += ...
        end
        if haskey(ref_node, "point load 2")
            M = ref_node("point load 2")
            info("calculate point moment")
            # fe += ...
        end
        add!(assembly.f, gdofs, fe)
    end
end

coupling = Problem(Coupling, "test", 2, "displacement")
# "slave nodes"
element3 = Element(Poi1, [2])
element4 = Element(Poi1, [3])
update!([element3, element4], "geometry", X)
element5 = Element(Poi1, [5])
update!(element5, "geometry", X)
update!(element5, "point moment 3", 10.0/d)
add_coupling_nodes!(coupling, [element3, element4])
add_reference_node!(coupling, element5)

step = Analysis(Nonlinear)
#xdmf = Xdmf("example1_results"; overwrite=true)
#add_results_writer!(step, xdmf)
add_problems!(step, [problem, bc, coupling])
run!(step)
#close(xdmf.hdf)
time = 0.0
u = element1("displacement", time)
println("displacement: $u")
