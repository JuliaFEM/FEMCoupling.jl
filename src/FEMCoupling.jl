# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

""" Coupling elements for JuliaFEM. """
module FEMCoupling

using Reexport
@reexport using FEMBase

include("get_C.jl")
include("Distributed_coupling.jl")


end
