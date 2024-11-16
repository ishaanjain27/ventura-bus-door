using DifferentialEquations
using Plots
using SparseArrays
using LinearAlgebra

# constants
n = 10
m = 1.0 * ones(n)
c = 0.2 * ones(n+1)
k = 5.0 * ones(n+1)
F = 0.0 * ones(n)

# define the differential equation
function n_masses!(ddu, du, u, p, t)
    m, c, k, F = p
    n = length(u)

    for i in 1:n
        if i == 1
            # left boundary
            ddu[i] = (F[i] - c[i]*(du[i] - 0) - c[i+1]*(du[i] - du[i+1]) - k[i]*(u[i] - 0) - k[i+1]*(u[i] - u[i+1])) / m[i]
        elseif i == n
            # right boundary
            ddu[i] = (F[i] - c[i]*(du[i] - du[i-1]) - c[i+1]*(du[i] - 0) - k[i]*(u[i] - u[i-1]) - k[i+1]*(u[i] - 0)) / m[i]
        else
            # middle masses
            ddu[i] = (F[i] - c[i]*(du[i] - du[i-1]) - c[i+1]*(du[i] - du[i+1]) - k[i]*(u[i] - u[i-1]) - k[i+1]*(u[i] - u[i+1])) / m[i]
        end
    end
end

# initial conditions
u0 = ones(n)
du0 = zeros(n)
tspan = (0.0, 20.0)
p = (m, c, k, F)

# solve the system
prob = SecondOrderODEProblem(n_masses!, du0, u0, tspan, p)
sol = solve(prob)

# plot the solution
plot(sol, vars=n+1:2*n)
