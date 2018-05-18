# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

using FEMBase
using FEMBase.Test

@testset "Test FEMCoupling.jl" begin
    @testset "test_01.jl" begin include("test_01.jl") end
end
