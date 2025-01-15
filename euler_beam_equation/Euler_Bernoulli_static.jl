using Plots
using LinearAlgebra

n = 100 # Number of nodes 
L = 5  # Beam length
EI = 70e6 * 5e-6  # Flexural rigidity
spring_stiffness_left = 100
spring_stiffness_right = 1000

x = range(0.0, L, length=n)
dx = L / (n-1)

function assemble_matrix(n, dx)
    # Fourth order derivative discretization
    A = zeros(n, n)
    for i in 3:n-2
        A[i, i-2:i+2] = [1.0, -4.0, 6.0, -4.0, 1.0]
    end
    A[1, 1:3] = [6.0, -4.0, 1.0]
    A[end, end-2: end] = [1.0, -4.0, 6.0]
    A[2, 1:4] = [-4.0, 6.0, -4.0, 1.0]
    A[end-1, end-3:end] = [1.0, -4.0, 6.0, -4.0]
    return A 
end

# Apply boundary conditions (cl = clmaped, ss = simply supported, fr = free end, spring = with a spring at the end)
function apply_boundary_conditions(left, right, A, load, n, dx)
    k1, k2, k3 = 6, -4, 1

    if left == "cl"
        A[1, 2] += k2; A[1, 3] += k3; A[2, 2] += k3
        A[1, :] .= 0; A[1, 1] = 1; load[1] = 0
    elseif left == "fr"
        A[1, 1] += k2 + 4k3; A[1, 2] += k2 - 4k3; A[1, 3] += -k2 + k3
        A[2, 1] += k3; A[2, 2] += k3; A[2, 3] += -k3
    elseif left == "ss"
        A[1, 2] += -k2; A[1, 3] += -k3; A[2, 2] += -k3
        A[1, :] .= 0; A[1, 1] = 1; load[1] = 0
    elseif left == "spring"
        A[1, 1] += k2 + 4k3 - k2*dx^3*spring_stiffness_left/EI + 2*k3*dx^3*spring_stiffness_left/EI; A[1, 2] += k2 - 4k3; A[1, 3] += -k2 + k3
        A[2, 1] += k3 - spring_stiffness_left*dx^3*k3/EI; A[2, 2] += k3; A[2, 3] += -k3
    end

    if right == "cl"
        A[end-1, end-1] += k3; A[end, end-2] += k3; A[end, end-1] += k2
        A[end, :] .= 0; A[end, end] = 1; load[end] = 0
    elseif right == "fr"
        A[end, end] += k2 + 4k3; A[end, end-1] += k2 - 4k3; A[end, end-2] += -k2 + k3
        A[end-1, end] += k3; A[end-1, end-1] += k3; A[end-1, end-2] += -k3
    elseif right == "ss"
        A[end-1, end-1] += -k3; A[end, end-2] += -k3; A[end, end-1] += -k2
        A[end, :] .= 0; A[end, end] = 1; load[end] = 0
    elseif right == "spring"
        A[end, end] += k2 + 4k3 - k2*dx^3*spring_stiffness_right/EI + 2*k3*dx^3*spring_stiffness_right/EI ; A[end, end-1] += k2 - 4k3; A[end, end-2] += -k2 + k3
        A[end-1, end] += k3 - spring_stiffness_right*dx^3*k3/EI; A[end-1, end-1] += k3; A[end-1, end-2] += -k3
    end
    return A, load
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

rhs = uniform_load(n, L, EI, x) #right hand side vector
A, rhs = apply_boundary_conditions("ss", "ss", A0, rhs, n, dx) 
A = A/dx^4

####   MODAL ANALYSIS


eigvalues, eigvectors = eigen(A)

# Normalize eigenvector signs
for i in 1:size(eigvectors, 2)
    idx = findfirst(abs.(eigvectors[:, i]) .> 1e-8)  
    if eigvectors[idx, i] < 0
        eigvectors[:, i] .= -eigvectors[:, i]  
    end
end
# Compute frequencies and sort eigenvalues and eigenvectors
frequencies = sqrt.(abs.(eigvalues))  
sorted_indices = sortperm(frequencies)  
eigvalues_sorted = eigvalues[sorted_indices]
eigvectors_sorted = eigvectors[:, sorted_indices]
frequencies_sorted = frequencies[sorted_indices]

plot(x, eigvectors_sorted[:, 1], label="Mode 1, Frequency=$(round(frequencies_sorted[1], digits=2))", legend=:bottomleft)
for i in 2:6
    plot!(x, eigvectors_sorted[:, i], label="Mode $(i), Frequency=$(round(frequencies_sorted[i], digits=2))")
end
title!("Eigenvectors with Lowest Frequencies")
xlabel!("x (m)")
ylabel!("Amplitude")

####   STATIC DEFLECTION

w = A \ rhs
max_deflection = minimum(w)
max_deflection_location = x[argmin(w)]


println("Maximum Deflection: ", max_deflection, " at x = ", max_deflection_location)

plot(x, w, label="Deflection w",xlabel="x (m)", ylabel="Deflection (mm)",title="Beam Deflection",legend=:top)

