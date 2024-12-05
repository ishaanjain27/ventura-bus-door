using Plots
using LinearAlgebra

n = 1000  # Number of nodes 
L = 5.0  # Beam length
EI = 70e6 * 5e-6  # Flexural rigidity

x = range(0.0, L, length=n)
dx = L / (n-1)

function assemble_matrix(n, dx)
    # Fourth order derivative discretization
    A = zeros(n, n)
    for i in 3:n-2
        A[i, i-2] = 1.0
        A[i, i-1] = -4.0
        A[i, i] = 6.0
        A[i, i+1] = -4.0
        A[i, i+2] = 1.0
    end
    return A / dx^4
end

function apply_boundary_conditions_simply_supported(A, rhs, n, dx)
    # w(0) = 0, w(L) = 0
    A[1, :] .= 0.0
    A[1, 1] = 1.0
    A[end, :] .= 0.0
    A[end, end] = 1.0

    # w''(0) = 0, w''(L) = 0
    A[2, :] .= 0.0
    A[2, 1:3] = [1.0, -2.0, 1.0] / dx^2
    A[end-1, :] .= 0.0
    A[end-1, end-2:end] = [1.0, -2.0, 1.0] / dx^2

    rhs[1] = 0.0
    rhs[2] = 0.0
    rhs[end-1] = 0.0
    rhs[end] = 0.0

    return A, rhs
end

function apply_boundary_conditions_clamped(A, rhs, n, dx)
    # w(0) = 0, w'(0) = 0

    A[1, :] .= 0.0
    A[1, 1] = 1.0  

    A[2, :] .= 0.0
    A[2, 1] = -1.0 / dx  
    A[2, 2] = 1.0 / dx
    # w''(L) = 0, w'''(L) = 0

    A[end, :] .= 0.0
    A[end, end-2:end] = [-1.0, 2.0, -1.0] / dx^2 
    A[end-1, :] .= 0.0
    A[end-1, end-3:end] = [1.0, -3.0, 3.0, -1.0] / dx^3  

    rhs[1] = 0.0  
    rhs[2] = 0.0 
    rhs[end-1] = 0.0  
    rhs[end] = 0.0  

    return A, rhs
end

function uniform_load(n, L, EI, x)
    rhs = fill(200.0, n)  
    return -rhs / EI
end

function linear_load(n, L, EI, x)
    q0 = 0.0  
    slope = 40.0  
    rhs = q0 .+ slope .* x  
    return -collect(rhs) / EI
end

A0 = assemble_matrix(n, dx)
rhs = uniform_load(n, L, EI, x)
A, rhs = apply_boundary_conditions_clamped(A0, rhs, n, dx)
w = A \ rhs

max_deflection = minimum(w)
max_deflection_location = x[argmin(w)]


println("Maximum Deflection: ", max_deflection, " at x = ", max_deflection_location)

plot(x, w, label="Deflection w",xlabel="x (m)", ylabel="Deflection (mm)",title="Beam Deflection",legend=:top)

eigvals, eigvecs = eigen(A)

plot_eigenshapes = plot(x, eigvecs[:, 1], label="Mode 1", legend=:topright)
for i in 2:5
    plot!(plot_eigenshapes, x, eigvecs[:, i], label="Mode $(i)")
end

title!("Modes of the Beam")
xlabel!("x (m)")
ylabel!("Amplitude")

          