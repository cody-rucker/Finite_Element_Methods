using LinearAlgebra
using SparseArrays

# solve Poisson's equation
#
#           -Δu = f     on Ω=[0,1]
#          u(0) = 0
#         uₓ(1) = 0

N = 10
h = 0.1
x = 0:h:1

f = 1

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

# Dirichlet condition u(0) = 0 gives ϕ₁ = 0
#A[1,1] = 0
#A[1,2] = 0
#A[2,1] = 0

b = (f*h) * ones(N)

u = A\b

@show u
