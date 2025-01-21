using Plots
using DifferentialEquations
using LinearAlgebra
using SparseArrays

# Parameters
n = 50  # Number of nodes
L = 5.0  # Beam length
EI = 70e9 * 5e-6  # Flexural rigidity (GPa * m^4)
mu = 27.0  # Mass density
c = 0  #damping coefficient (critically damping c = 6750)
spring_stiffness_left = 100
spring_stiffness_right = 100

x = range(0.0, L, length=n)
dx = L / (n-1)

# Assemble matrix for fourth-order derivative discretization
function assemble_matrix(n, dx)
    e = ones(n)
    A = spdiagm(-2 => 1*e[3:end], -1 => -4*e[2:end], 0 => 6*e, 1 => -4*e[2:end], 2 => 1*e[3:end])
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

# Uniform load
function uniform_load(n)
    load = fill(200.0, n)  # Uniform load (N/m)
    return -load
end

# Beam solver
function beam_solver(A, q, mu, c, n, x, EI)
    function beam_ode!(du, u, p, t)
        EI, mu, q, A, c = p
        w = u[1:n]        
        v = u[n+1:2n]     
        du[1:n] .= v      # w' = v
        du[n+1:2n] .= (-EI * A * w + q - c*v) / mu  # v' = -(EI * A * w + q) / mu
    end

    u0 = zeros(2n) # Initial conditions
    tspan = (0.0, 2.0)
    params = (EI, mu, q, A, c)

    prob = ODEProblem(beam_ode!, u0, tspan, params)
    sol = solve(prob, Tsit5())

    return sol
end

# Main program

@time begin
    A0 = assemble_matrix(n, dx)
    q = uniform_load(n)
    A, q = apply_boundary_conditions("cl", "fr", A0, q, n, dx)
    A = A / dx^4

    sol = beam_solver(A, q, mu, c, n, x, EI)
end
# Animation of beam deflection
w = sol[1:n, :] * 10^3
times = sol.t
frame_count = 500  # Desired number of frames
step_size = max(1, div(length(times), frame_count))

anim = @animate for i in 1:step_size:length(times)
    plot(x, w[:, i], xlabel="x (m)", ylabel="w (mm)", label="Deflection",
         title="Beam Deflection at t = $(round(times[i], digits=2)) s", legend=false, ylims = (-100, 100), xlims = (0, 6))
end
gif(anim, "beam_deflection.gif", fps=30)
