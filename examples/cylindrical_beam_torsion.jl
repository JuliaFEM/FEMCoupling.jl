# Analytical calculation
T=1.500
L=210e-3
r=16e-3
d=2r
G=210e9/(2*(1+0.3))

Wv=pi*d^3/16
tmax=T/Wv
a=2*L*tmax/(G*d)

println(a)

using JuliaFEM
using JuliaFEM.Preprocess
using JuliaFEM.AbaqusReader.create_surface_elements
using JuliaFEM.AbaqusReader.abaqus_read_mesh

# read mesh
mesh = abaqus_read_mesh("varsinolla.inp")

# problem
cylinder_probem = Problem(Elasticity,"varsinolla",6)
cylinder_element = create_elements(mesh,"Body1")
update!(Body1, "youngs modulus", 210e3)
update!(Body1, "poissons ratio", 0.3)
update!(body1, "density", 7.80e-9)
add_elements!(cylinder_problem,cylinder_element)
