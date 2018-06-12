using JuliaFEM
using FEMBase
using FEMBase.Test

using FEMCoupling: get_C

@testset "3D beam with point moments" begin

u_expected=zeros(6,1)
u_expected[2,1]=25
u_expected[3,1]=-25
u_expected[4,1]=416.67
u_expected[5,1]=50
u_expected[6,1]=50

# Stiffness matrix for beam X1=[0.0, 0.0, 0.0], X2=[1.0, 0.0, 0.0],
# n2=[0.0, 0.0, -1.0]
# youngs modulus= 100.0e2, shear modulus= 80.0e2
# cross-section values 2.5E-2, 1.0E-4, 0.0E-4, 1.5E-4, 3.0E-5
k_beam=
 [300.0 0.0 0.0 0.0 0.0 0.0 -300.0 0.0 0.0 0.0 0.0 0.0;
 0.0 12.0 0.0 0.0 0.0 6.0 0.0 -12.0 0.0 0.0 0.0 6.0;
 0.0 0.0 18.0 0.0 -9.0 0.0 0.0 0.0 -18.0 0.0 -9.0 0.0;
 0.0 0.0 0.0 0.24 0.0 0.0 0.0 0.0 0.0 -0.24 0.0 0.0;
 0.0 0.0 -9.0 0.0 6.0 0.0 0.0 0.0 9.0 0.0 3.0 0.0;
 0.0 6.0 0.0 0.0 0.0 4.0 0.0 -6.0 0.0 0.0 0.0 2.0;
 -300.0 0.0 0.0 0.0 0.0 0.0 300.0 0.0 0.0 0.0 0.0 0.0;
 0.0 -12.0 0.0 0.0 0.0 -6.0 0.0 12.0 0.0 0.0 0.0 -6.0;
 0.0 0.0 -18.0 0.0 9.0 0.0 0.0 0.0 18.0 0.0 9.0 0.0;
 0.0 0.0 0.0 -0.24 0.0 0.0 0.0 0.0 0.0 0.24 0.0 0.0;
 0.0 0.0 -9.0 0.0 3.0 0.0 0.0 0.0 9.0 0.0 6.0 0.0;
 0.0 6.0 0.0 0.0 0.0 2.0 0.0 -6.0 0.0 0.0 0.0 4.0]

 # Stiffness matrix for 12x12 dof beam + single 6 dof node
K=zeros(18,18)
K[1:12,1:12]=k_beam

K=K[7:end,7:end] # dofs 1-6 fixed

# Force vector: point moments to single node
f=zeros(18,1)
f[16,1]=100
f[17,1]=75
f[18,1]=50

f=f[7:end] # dofs 1-6 fixed

# get_C(refnode,slaves,dofs,ndofs,K_size)
C,D,g=get_C(3,[2],[1,2,3,4,5,6],6,18)
C=C[1:end,7:end]

K=[K C';
   C D]
f=[f;
   g]
u=lufact(K)\full(f)

@test isapprox(u[1:6],u_expected,rtol=0.0001)
end
