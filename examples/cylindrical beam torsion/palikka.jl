using JuliaFEM
using JuliaFEM.Preprocess
using JuliaFEM.Postprocess
using JuliaFEM.Abaqus: create_surface_elements
using FEMBase

# read mesh
mesh = abaqus_read_mesh("palikka.inp")

# create a field problem
kalikka = Problem(Elasticity, "palikkaaa", 3)
elements = create_elements(mesh, "Solid1")
update!(elements, "youngs modulus", 207.0E3)
update!(elements, "poissons ratio", 0.30)
update!(elements, "density", 7.80E-9)
FEMBase.add_elements!(kalikka, elements)

# create boundary conditions from node sets
fixed = Problem(Dirichlet, "fixed1", 3, "displacement")
fixed_elements = create_nodal_elements(mesh, "Face_Constraint_1")
update!(fixed_elements, "displacement 1", 0.0)
update!(fixed_elements, "displacement 2", 0.0)
update!(fixed_elements, "displacement 3", 0.0)
FEMBase.add_elements!(fixed, fixed_elements)

# add the field and the boundary problems to the solver
freqs = Solver(Modal, kalikka, fixed)

# solve 6 smallest eigenvalues
freqs.properties.nev = 6
freqs.properties.which = :SM
freqs()

println("Eigenvalues: ", sqrt.(freqs.properties.eigvals) / (2*pi))
