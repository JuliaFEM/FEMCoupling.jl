# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMCoupling.jl/blob/master/LICENSE

using FEMCoupling: get_C
using Base.Test

@testset "Plain strain kinematic Coupling" begin
# plane strain problem
# Square shaped plane strain element(8x8) + one single node(2x2)
K=zeros(10,10)
K[1:8,1:8]=
  [500 225 -350 75 -250 -225 100 -75;
   225 500 -75 100 -225 -250 75 -350;
  -350 -75 500 -225  100 75 -250 225;
   75 100 -225  500 -75 -350 225 -250;
  -250 -225 100 -75   500 225 -350 75;
  -225 -250 75  -350  225 500 -75 100;
   100 75 -250 225 -350 -75 500 -225;
  -75 -350 225 -250 75 100 -225 500]

# Removing fixed DOFs (1 and 2) from nodes 1 and 4
K=K[1:end .!=1, 1:end .!=1]
K=K[1:end .!=1, 1:end .!=1]
K=K[1:end .!=5, 1:end .!=5]
K=K[1:end .!=5, 1:end .!=5]

# Force vector
f=zeros(10,1)
f[9,1]=1200

# Removing fixed DOFs
f=f[1:end .!= 1]
f=f[1:end .!= 1]
f=f[1:end .!= 5]
f=f[1:end .!= 5]

############### Calculations by hand

# Making the C matrix by hand
A=zeros(4,10)
A[1,3]-=1;  A[1,9]=1; A[2,4]=-1; A[2,10]=1;
A[3,5]-=1;  A[3,9]=1; A[4,6]=-1; A[4,10]=1

# Eliminating fixed dofs (1,2) from nodes 1 and 4
A=A[1:end , 1:end .!=1 ]
A=A[1:end , 1:end .!=1 ]
A=A[1:end , 1:end .!=5 ]
A=A[1:end , 1:end .!=5 ]

# Renaming variables to match with get_C.jl
C_expected=A
g_expected=zeros(4,1)
D_expected=zeros(4,4)

# Assembly for solving
K_expected= [K            C_expected';
             C_expected   D_expected]

f_expected= [f;
             g_expected]

u_expected = K_expected\f_expected

############### Calculating C,D and g with get_C.jl
#                  get_C(refnode,slaves,dofs,ndofs,K_size)
K_size=size(K,1)
C,D,g= FEMCoupling.get_C(5,[2,3],[1,2],2,10)

# Removing fixed DOFs
C=C[1:end , 1:end .!=1 ]
C=C[1:end , 1:end .!=1 ]
C=C[1:end , 1:end .!=5 ]
C=C[1:end , 1:end .!=5 ]

KK=[K C';
    C D]
ff=[f;
    g]
u = lufact(KK) \ full(ff)

@test isapprox(C,C_expected,rtol=0.0001)
@test isapprox(D,D_expected,rtol=0.0001)
@test isapprox(g,g_expected,rtol=0.0001)
@test isapprox(ff,f_expected,rtol=0.0001)
@test isapprox(KK,K_expected,rtol=0.0001)
@test isapprox(u,u_expected,rtol=0.0001)
# If the last test passes, all other tests will pass too.
# Other tests are made to help tracing why the last test doesn't pass.
end
