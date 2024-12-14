using Plots

n = 20000  # Number of nodes 
L = 5.0  # Beam length
EI = 70e6 * 5e-6  # Flexural rigidity

x = range(0.0, L, length=n)
dx = L / (n-1)

function assemble_matrix(n, dx)
    A = zeros(n, n)
    for i in 2:n-1
        A[i, i-1] = 1.0
        A[i, i] = -2.0
        A[i, i+1] = 1.0
    end
    return A / dx^2
end

function apply_boundary_conditions_simply_supported(A, rhs, n, dx)
    # Simply supported at x = 0 and x = L (w(0) = 0, w(L) = 0)
    A[1, :] .= 0.0
    A[1, 1] = 1.0  
    A[end, :] .= 0.0
    A[end, end] = 1.0  
    rhs[1] = 0.0  
    rhs[end] = 0.0  
    return A, rhs
end



function uniform_load(n, L, EI, x)
    rhs = fill(200.0, n)  
    return -rhs / EI
end
 

A = assemble_matrix(n, dx)

rhs = uniform_load(n, L, EI, x)

A, rhs = apply_boundary_conditions_simply_supported(A, rhs, n, dx)
A2 = A * A

w = A2 \ rhs

max_deflection = minimum(w)
max_deflection_location = x[argmin(w)]


println("Maximum Deflection: ", max_deflection, " at x = ", max_deflection_location)

plot(x, w, label="Deflection w",
     xlabel="x (m)", ylabel="Deflection (mm)",
     title="Beam Deflection",
     legend=:top)

eigvalues, eigvectors = eigen(A2)

plot_eigenshapes = plot(x, eigvectors[:, 1], label="Mode 1", legend=:topright)
for i in 2:5
    plot!(plot_eigenshapes, x, eigvectors[:, i], label="Mode $(i)")
end

title!("Modes of the Beam")
xlabel!("x (m)")
ylabel!("Amplitude")
