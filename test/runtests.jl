# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

# using FEMBase
using Base.Test

@testset "Test FEMCoupling.jl" begin
    @testset "test_dcoupling_brick.jl" begin include("test_dcoupling_brick.jl") end
    @testset "test_dcoupling_2d.jl" begin include("test_dcoupling_2d.jl") end
    @testset "test_plain_strain_kinematic_coupling.jl" begin include("test_plain_strain_kinematic_coupling.jl") end
end
