# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

type Coupling <: BoundaryProblem
    reference_node :: Element{Poi1}
end

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

    dimensions = problem.dimension

    function to3d(x)
        if length(x) == 2
            return [x; 0.0]
        end
        return x
    end

    ref_node = problem.properties.reference_node
    if length(get_connectivity(ref_node)) == 0
        error("Reference node node defined. Define reference node using add_reference_node!")
    end
    ref_node_id = first(get_connectivity(ref_node))
    X_r = to3d(first(ref_node("geometry", time)))
    fe = zeros(dimensions)
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
    r = Dict{Int64,Vector{Float64}}()
    for n in keys(weights)
        r[n] = X_n[n]-xbar
    end
    rR = X_r-xbar
    T = zeros(3,3)
    for n in keys(weights)
        T += weights[n]*(dot(r[n],r[n])*eye(3)-(r[n]*r[n]'))
    end
    if det(T) == 0
        indz = indmin(diag(T))
        T[indz,indz] += 1e-9
    end

    invT = inv(T)

    # handle force boundary conditions / load distribution

    for coupling_node in elements
        fill!(fe, 0.0)
        coupling_node_id = first(get_connectivity(coupling_node))
        X_n = first(coupling_node("geometry", time))
        gdofs = get_gdofs(problem, coupling_node)
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

        MRhat = MR + cross(rR,FR)
        n = coupling_node_id
        Fn = weights[n]*(FR+cross(invT*MRhat,r[n]))
        fe = fe+Fn[1:dimensions]
        add!(assembly.f, gdofs, fe)
    end

    # handle displacement boundary conditions / constraint expression

    function uRfactor(factor, rn, wn)
        un = zeros(3)
        un[factor] = 1.0
        return wn*un + cross(invT * wn*cross(rn, un), rR)
    end

    function wRfactor(factor, rn, wn)
        un = zeros(3)
        un[factor] = 1.0
        return invT*wn*cross(rn, un)
    end

    # C_all_dofs = zeros(6,3*length(r))
    # ns = []
    # for coupling_node in elements
    #     coupling_node_id = first(get_connectivity(coupling_node))
    #     push!(ns, coupling_node_id)
    # end

    rgdofs = get_gdofs(problem, ref_node)
    info("rgdofs = $rgdofs")

    for coupling_node in elements
        # couplingelementnumber = first(findin(ns,n))
        cgdofs = get_gdofs(problem, coupling_node)
        info("cgdofs = $cgdofs")

        uRn = zeros(6, 3)
        for n in keys(weights)
            for j in 1:3
                uRn[1:3, j] = uRfactor(j, r[n], weights[n])
                uRn[4:6, j] = wRfactor(j, r[n], weights[n])
            end
        end
        #C_all_dofs[1:6, couplingelementnumber*3-3)+1:(couplingelementnumber*3-3)+3]= uRn
        add!(assembly.C2, rgdofs, cgdofs, uRn[1:length(rgdofs),1:length(rgdofs)])
        # add!(assembly.C1, rgdofs, cgdofs, f)
    end

    # add -1 to diagonal of D
    for j in rgdofs
        add!(assembly.D, j, j, -1.0)
    end

    # if dimensions == 2
    #     rows=[3,4,5]
    #     sort!(rows, rev=true)
    #     for i in rows
    #         C_all_dofs=C_all_dofs[1:end .!=i, 1:end]
    #     end
    #
    #     columns=[3,6]
    #     sort!(columns, rev=true)
    #     for i in columns
    #         C_all_dofs=C_all_dofs[1:end, 1:end .!=i]
    #     end
    # end

    return nothing

end
