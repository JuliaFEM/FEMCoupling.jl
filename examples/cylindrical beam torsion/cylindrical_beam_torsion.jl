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

using JuliaFEM: add_elements!
using JuliaFEM.Preprocess
using JuliaFEM.Postprocess
using JuliaFEM.Abaqus: create_surface_elements
using FEMBase
# read mesh
mesh = abaqus_read_mesh("cylindrical_beam_torsion1.inp")

# problem
cylinder_problem = Problem(Elasticity,"varsinolla",6)
cylinder_element = create_elements(mesh,"Body1")
update!(cylinder_element, "youngs modulus", 210e3)
update!(cylinder_element, "poissons ratio", 0.3)
update!(cylinder_element, "density", 7.80e-9)
add_elements!(cylinder_problem,cylinder_element)
