# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE
using Revise
revise()

using FEMBase
using FEMBase.Test


# Create coupling element

coupling = Problem(KinematicCoupling, "couple something", 2, "displacement")

# TODO ...
