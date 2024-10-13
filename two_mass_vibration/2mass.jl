using LinearAlgebra
using Plots
using DifferentialEquations

# Parameters for ODE and Numerical solutions
k1 = 5.0      
m1 = 1.0       
c1 = 0.2  
k2 = 10.0
m2 = 2.0
c2 = 0.4     
x10 = 2.0 
x20 = 1.0  
x1dot0 = 0.0 
x2dot0 = 0.0 
dt = 0.01  
tfinal = 20.0
time = collect(0:dt:tfinal)

f(t) = 0

# ODE SOLVER
function equation!(du, u, p, t)
    x1, x2, x3, x4 = u  
    m1, m2, k1, k2, c1, c2 = p
    du[1] = x3  
    du[2] = x4  
    du[3] = (-k1*x1 - c1*x3 + k2*(x2 - x1) + c2*(x4 - x3)) / m1  
    du[4] = (-k2*(x2 - x1) - c2*(x4 - x3)) / m2  
end

u0 = [x10, x20, x1dot0, x2dot0]  
param = [m1, m2, k1, k2, c1, c2]

ode = ODEProblem(equation!, u0, (0.0, tfinal), param)
ode_sol = solve(ode, saveat=time)

ode_time = ode_sol.t
x1_position_ode = ode_sol[1, :]  
x2_position_ode = ode_sol[2, :]

# NUMERICAL SOLUTION (Euler's Method)
A = [0.0 1.0 0.0 0.0;
     -(k1+k2)/m1 -(c1 +c2)/m1 k2/m1 c2/m1;
     0.0 0.0 0.0 1.0;
     k2/m2 c2/m2 -k2/m2 -c2/m2]

# Initial conditions for mass 1 and mass 2
u = [x10;x1dot0;x20;x2dot0]
disp1 = [x10]
disp2 = [x20]
vel1 = [x1dot0]
vel2 = [x2dot0]

for t in time[2:end] 
    global u
    uprime = A*u
    u_new = u + dt * uprime  
    push!(disp1, u_new[1]) 
    push!(disp2, u_new[3]) 
    push!(vel1, u_new[2])   
    push!(vel2, u_new[4])   

    u = u_new
end

# Plot results - comparison between ODE and numerical for mass 1
plot1 = plot(ode_time, x1_position_ode, label="Mass 1 (ODE)", color=:blue, lw=2, 
    title="Mass 1 Position", xlabel="Time (s)", ylabel="Position", legend=:topright)
plot!(plot1, time, disp1, label="Mass 1 (Numerical)", color=:red, lw=2)

# Plot results - comparison between ODE and numerical for mass 2
plot2 = plot(ode_time, x2_position_ode, label="Mass 2 (ODE)", color=:blue, lw=2,
    title="Mass 2 Position", xlabel="Time (s)", ylabel="Position", legend=:topright)
plot!(plot2, time, disp2, label="Mass 2 (Numerical)", color=:red, lw=2)

# Combine the two plots
plot(plot1, plot2, layout = (2, 1), size=(800, 600))
