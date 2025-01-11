using Plots
using DifferentialEquations
using LinearAlgebra

# Parameters
n = 100  # Number of nodes
L = 5.0  # Beam length
EI = 70e6 * 5e-6  # Flexural rigidity (kPa * m^4)
mu = 1.0  # Mass density
c = 0.5  

x = range(0.0, L, length=n)
dx = L / (n-1)

# Assemble matrix for fourth-order derivative discretization
function assemble_matrix(n, dx)
    A = zeros(n, n)
    for i in 3:n-2
        A[i, i-2:i+2] = [1.0, -4.0, 6.0, -4.0, 1.0]
    end
    A[1, 1:3] = [6.0, -4.0, 1.0]
    A[end, end-2:end] = [1.0, -4.0, 6.0]
    A[2, 1:4] = [-4.0, 6.0, -4.0, 1.0]
    A[end-1, end-3:end] = [1.0, -4.0, 6.0, -4.0]
    return A
end

# Apply boundary conditions
function apply_boundary_conditions(left, right, A, load, n, dx)
    k1, k2, k3 = 6, -4, 1

    if left == "cl"
        A[1, 2] += k2; A[1, 3] += k3; A[2, 2] += k3
        A[1, :] .= 0; A[1, 1] = 1
    elseif left == "fr"
        A[1, 1] += k2 + 4k3; A[1, 2] += k2 - 4k3; A[1, 3] += -k2 + k3
        A[2, 1] += k3; A[2, 2] += k3; A[2, 3] += -k3
    elseif left == "ss"
        A[1, 2] += -k2; A[1, 3] += -k3; A[2, 2] += -k3
        A[1, :] .= 0; A[1, 1] = 1
    end

    if right == "cl"
        A[end-1, end-1] += k3; A[end, end-2] += k3; A[end, end-1] += k2
        A[end, :] .= 0; A[end, end] = 1
    elseif right == "fr"
        A[end, end] += k2 + 4k3; A[end, end-1] += k2 - 4k3; A[end, end-2] += -k2 + k3
        A[end-1, end] += k3; A[end-1, end-1] += k3; A[end-1, end-2] += -k3
    elseif right == "ss"
        A[end-1, end-1] += -k3; A[end, end-2] += -k3; A[end, end-1] += -k2
        A[end, :] .= 0; A[end, end] = 1
    end

    load[2] = 0.0; load[end-1] = 0.0
    return A, load
end

# Uniform load
function uniform_load(n)
    load = fill(200.0, n)  # Uniform load (N/m)
    return -load
end

# Beam solver
function beam_solver(A, q, mu, c, n, x, EI)
    function beam_ode!(du, u, p, t)
        EI, mu, q, A, c = p
        w = u[1:n]        # Displacement
        v = u[n+1:2n]     # Velocity
        du[1:n] .= v      # w' = v
        du[n+1:2n] .= (-EI * A * w + q - c*v) / mu  # v' = -(EI * A * w + q) / mu
    end

    # Initial conditions
    u0 = zeros(2n)
    tspan = (0.0, 10.0)
    params = (EI, mu, q, A, c)

    prob = ODEProblem(beam_ode!, u0, tspan, params)
    sol = solve(prob, ImplicitEuler(), saveat=0.1)

    return sol
end

# Main program
A0 = assemble_matrix(n, dx)
q = uniform_load(n)
A, q = apply_boundary_conditions("ss", "ss", A0, q, n, dx)
A = A / dx^4

sol = beam_solver(A, q, mu, c, n, x, EI)

# Animation of beam deflection
w = sol[1:n, :]
times = sol.t

anim = @animate for i in 1:length(times)
    plot(x, w[:, i], xlabel="x (m)", ylabel="w (mm)", label="Deflection",
         title="Beam Deflection at t = $(round(times[i], digits=2)) s", legend=false, ylims = (-100, 100), xlims = (0, 6))
end

# Save as GIF
gif(anim, "beam_deflection.gif", fps=30)
