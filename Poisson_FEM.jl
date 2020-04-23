using LinearAlgebra
using SparseArrays
using Plots

# solve Poisson's equation
#
#           -Δu = f     on Ω=[0,1]
#          u(0) = 0
#         uₓ(1) = 0


h = 0.1
x = 0:h:1
N = length(x)
f = 1

b = (f*h) * ones(N)

# build stiffness matrix Aᵢⱼ = ∫ϕ′ᵢ ⋅ ϕ′ⱼ
Id = 1:N
Jd = copy(Id)
Vd = (2/h)*ones(N)

Il = 2:N
Jl = 1:N-1
Vl = (-1/h)*ones(N-1)

Iu = 1:N-1
Ju = 2:N
Vu = (-1/h)*ones(N-1)

A = sparse(Id, Jd, Vd, N, N) + sparse(Il, Jl, Vl, N, N) + sparse(Iu, Ju, Vu, N, N)

# basis function ϕₙ is only defined for x < xₙ
A[N,N] = 0.5A[N,N]

# basis function ϕ₁ is only defined for x₁ < x
A[1,1] = 0.5A[1,1]

# use a Dirichlet extension to impose non-homogeneous boundary conditions

u∂ = zeros(N)
u∂[1] = 0
u∂[N] = 0

A∂ = A * u∂

b̃ = b - A∂

# omit basis functions on boundary elements
A[1,1] = 1
A[1,2] = 0
A[2,1] = 0


# Neumann(natural) conditions are built in to the variational form
# imposing them is just a matter of NOT imposing Dirichlet(essential)
# conditions at a boundary point.

# imposing Dirichlet conditon at x=1

A[N,N] = 1
A[N,N-1] = 0
A[N-1,N] = 0


# set RHS data
b̃[1] = 0
b̃[N] = 0

u = A\b̃

plot(u)
