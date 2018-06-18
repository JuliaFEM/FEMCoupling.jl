# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

type Coupling <: BoundaryProblem
    reference_node :: Element{Poi1}
end

export Coupling

function Coupling()
    reference_node = Element(Poi1, Int64[])
    return Coupling(reference_node)
end

function add_coupling_nodes!(problem, elements)
    for element in elements
        push!(problem.elements, element)
    end
    return nothing
end

function add_reference_node!(problem, element)
    problem.properties.reference_node = element
    return nothing
end

function FEMBase.assemble_elements!(problem::Problem{Coupling},
                                    assembly::Assembly,
                                    elements::Vector{Element{Poi1}}, time)

    function to3d(x)
        if length(x)==2
            return [x;0.0]
        end
        return x
    end
    ref_node = problem.properties.reference_node
    if length(get_connectivity(ref_node)) == 0
        error("Reference node node defined. Define reference node using add_reference_node!")
    end
    ref_node_id = first(get_connectivity(ref_node))
    # info("Reference node id = $ref_node_id")
    X_r = to3d(first(ref_node("geometry", time)))
    # info("Reference node geometry = $X_r")
    fe = zeros(3)
    weights = Dict{Int64,Float64}()
    X_n = Dict{Int64,Vector{Float64}}()
    for coupling_node in elements
        coupling_node_id = first(get_connectivity(coupling_node))
        X_n[coupling_node_id] = to3d(first(coupling_node("geometry", time)))
        weights[coupling_node_id] = 1.0
    end
    total = 0.0
    for n in keys(weights)
        total = total+weights[n]
    end
    for n in keys(weights)
        weights[n] = weights[n]/total
    end
    xbar = zeros(3)
    for n in keys(weights)
        xbar = xbar+weights[n]*X_n[n]
    end
    r=Dict{Int64,Vector{Float64}}()
    for n in keys(weights)
        r[n] = X_n[n]-xbar
    end
            rR = X_r-xbar
    T=zeros(3,3)
    for n in keys(weights)
        T = T+weights[n]*(dot(r[n],r[n])*eye(3)-(r[n]*r[n]'))
    end
    display(T)
    # T[2,2]=1
    invT = inv(T)

    for coupling_node in elements
        fill!(fe, 0.0)
        coupling_node_id = first(get_connectivity(coupling_node))
        X_n = first(coupling_node("geometry", time))
        # info("Coupling node id = $coupling_node_id, geometry = $X_n")
        gdofs = get_gdofs(problem, coupling_node)
        # info("gdofs = $gdofs")
        FR = zeros(3)
        MR = zeros(3)
        for i = 1:3
            if haskey(ref_node, "point moment $i")
                MR[i] = ref_node("point moment $i", time)
            end
            if haskey(ref_node, "point force $i")
                FR[i] = ref_node("point force $i", time)
            end
        end

        # info("FR = $FR,MR = $MR")

        # info("calculate point moment")

        display(rR)
        display(FR)
        if length(rR) == 2
            rR=[rR;0.0]
        end
            MRhat = MR + cross(rR,FR)
            n = coupling_node_id
            Fn = weights[n]*(FR+cross(invT*MRhat,r[n]))
            fe = fe+Fn[1:3]
        add!(assembly.f, gdofs, fe)
    end
end
