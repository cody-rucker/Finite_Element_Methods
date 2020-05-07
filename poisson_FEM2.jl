using LinearAlgebra
using Plots
#using Printf
include("convergence_rates.jl")

function trapezoid_quad(A, h)
    N = size(A)[2]
    M = size(A)[1]
    I = 0.0
    w = 1.0*ones(size(A))
    w[1,1] = 0.25
    w[M,1] = 0.25
    w[1,N] = 0.25
    w[M,N] = 0.25

    w[1, 2:N-1] = 0.5*ones(N-2,1)
    w[M, 2:N-1] = 0.5*ones(N-2,1)
    w[2:M-1, 1] = 0.5*ones(N-2,1)
    w[2:M-1, N] = 0.5*ones(N-2,1)

    for i = 1:N
        for j = 1:M
            I += w[i,j] * A[i,j] * h^2
        end
    end
    return I
end



function f(x,y)
    return 2*(pi^2)*sin(pi*x)*sin(pi*y)
end

function poisson_fem(h)
    #h = 1/32           # grid spacing
    M = Int(1/h)         # number of nodes in each direction
    N = (M-1)^2          # number of interior nodes
    A = zeros(N, N)
    V = -1*ones(M+1, M+1)
    b = zeros(N, 1)

    # grid for plotting later on
    x = 0:h:1
    y = 0:h:1

    for i = 1:M-1
        for j = 1:M-1
            V[j+1, i+1] = (i-1)*(M-1) + j
        end
    end

    lower_nodes = zeros(M^2, 3)
    upper_nodes = zeros(M^2, 3)
    lower_barycenters = zeros(M^2, 2)
    upper_barycenters = zeros(M^2, 2)

    for i = 1:M
        for j = 1:M
            lower_nodes[(i-1)*M+j,:] = [V[j,i], V[j, i+1], V[j+1, i+1]]
            upper_nodes[(i-1)*M+j,:] = [V[j,i], V[j+1, i+1], V[j+1, i]]
            lower_barycenters[(i-1)*M+j, :] = [2*h/3 + (i-1)*h, h/3+(j-1)*h]
            upper_barycenters[(i-1)*M+j, :] = [h/3+(i-1)*h, 2*h/3+(j-1)*h]
        end
    end

    Ke_lower = [ 0.5 -0.5 0;
                -0.5 1 -0.5;
                0 -0.5 0.5]
    Ke_upper = [0.5 0 -0.5;
                0 0.5 -0.5;
               -0.5 -0.5 1]

    s = M^2

    function boundary_exile(node_triplet)
        indices = [1, 2, 3]
        index = [-1]

        for j = 1:3
            if node_triplet[j] > -1
                append!(index, indices[j])
            end
        end
        sidx = Int(size(index)[1])
        index = index[2:sidx]
        nodes = node_triplet[index[:]]
        number_of_nodes = Int(size(nodes)[1])

        return nodes, number_of_nodes, index

    end

    #nodes, num_nodes, idx = boundary_exile([1 2 3])

    for i = 1:s
        nodes, num_nodes, idx = boundary_exile(lower_nodes[i,:])
        for m = 1:num_nodes
            for p = 1:num_nodes
                A[Int(nodes[m]), Int(nodes[p])] += Ke_lower[idx[m], idx[p]]
            end

            b[Int(nodes[m])] += h^2/(6)  * f(lower_barycenters[i,1], lower_barycenters[i,2])

            # I think here we can compute portion of b from
            # lower triangles. b[Int(nodes[m])] = f(lower_barycenters[m,:])*ϕ(lower_bartcenters[m,:])|K|
        end
        nodes, num_nodes, idx = boundary_exile(upper_nodes[i,:])
        for m = 1:num_nodes
            for p = 1:num_nodes
                A[Int(nodes[m]), Int(nodes[p])] += Ke_upper[idx[m], idx[p]]
            end

            b[Int(nodes[m])] += h^2/(6)  * f(upper_barycenters[i,1], upper_barycenters[i,2])

        end
    end

    u = A\b

    uexact(x,y) = sin(pi*x)*sin(pi*y)

    u = reshape(u, (M-1,M-1))
    Ã = zeros(M+1,M+1)
    Ã[2:M, 2:M] = u[:,:]

    U = zeros(size(Ã))
    for i = 1:length(x)
        for j = 1:length(y)
            U[i,j] = uexact(x[i], y[j])
        end
    end

    error = trapezoid_quad(Ã - U, h)
    return error

    #=
    pyplot()
    p1 = surface(x, y, uexact, title="exact")
    p2 = surface(x, y, Ã, title="approximation")

    plot(p1,p2, layout=(1,2))
    =#
end

convergence_rates(poisson_fem, 0.25, 5)
