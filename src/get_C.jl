# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE
function get_C(refnode,slaves,dofs,ndofs,K_size)

    C_rows=size(slaves,1)*size(dofs,1)
    CC=spzeros(C_rows,K_size)
    g=spzeros(C_rows,1)
    D=spzeros(C_rows,C_rows)

    for i in 1:length(dofs),
        j in 0:length(slaves)-1
        CC[i+j*length(dofs) , refnode*ndofs-ndofs+dofs[i]]=1
        CC[i+j*length(dofs) , slaves[j+1]*ndofs-ndofs+dofs[i]]=-1
    end
    return CC,D,g
end
