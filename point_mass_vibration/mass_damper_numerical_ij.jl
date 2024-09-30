using LinearAlgebra
using SparseArrays
using Plots

# define domain and mesh
t_start = 0.0
t_end = 50.0
n_intervals = 30000
h = (t_end - t_start) / n_intervals
t_domain = range(t_start, stop=t_end - h, length=n_intervals)

# sefine constants in SI units
m = 1.0
c = 0.2
k = 5.0

# define system matrix coefficients
R = m / (h^2) + c / h + k
S = -2 * m / (h^2) - c / h
T = m / (h^2)

# create system matrix
A = spdiagm(-2 => fill(T, n_intervals-2),
            -1 => fill(S, n_intervals-1),
             0 => fill(R, n_intervals))

# define vector b
b = zeros(n_intervals)
b[1] = -1 * (T + S)
b[2] = -1 * T

# solve the system
x = A \ b

# plot the solution
plot(t_domain, x, label="Displacement", xlabel="Time (s)", ylabel="Displacement (m)", grid=true)
