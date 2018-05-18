# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

using Documenter, FEMCoupling

makedocs(modules=[FEMCoupling],
         format = :html,
         checkdocs = :all,
         sitename = "FEMCoupling.jl",
         pages = ["index.md"]
        )
