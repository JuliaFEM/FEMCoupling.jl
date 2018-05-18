# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

using FEMBase
using FEMBase.Test

using FEMCoupling

# Create coupling element

coupling = Problem(KinematicCoupling, "couple something", 2, "displacement")

# TODO ...
