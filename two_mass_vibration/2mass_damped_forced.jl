using Plots
using DifferentialEquations
using SymPy

# Parameters for ODE, Numerical, and Analytical solutions
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

# External force function
f(t) = 0

# ----------------------------------------------
# ODE SOLVER with DifferentialEquations.jl
# ----------------------------------------------
function equation!(du, u, p, t)
    x1, x2, x3, x4 = u  
    m1, m2, k1, k2, c1, c2 = p
    du[1] = x3  
    du[2] = x4  
    du[3] = (-k1*x1 - c1*x3 + k2*(x2 - x1) + c2*(x4 - x3)) / m1  + f(t) / m1  # force on mass 1
    du[4] = (-k2*(x2 - x1) - c2*(x4 - x3)) / m2  + f(t) / m2      # force on mass 2
end

u0 = [x10, x20, x1dot0, x2dot0]  
param = [m1, m2, k1, k2, c1, c2]

ode = ODEProblem(equation!, u0, (0.0, tfinal), param)
ode_sol = solve(ode, Trapezoid(), dt=dt, saveat=time)

ode_time = ode_sol.t
x1_position_ode = ode_sol[1, :]  
x2_position_ode = ode_sol[2, :]

# ----------------------------------------------
# NUMERICAL SOLUTION (Euler's Method)
# ----------------------------------------------

A = [0.0 1.0 0.0 0.0;
     -(k1+k2)/m1 -(c1 +c2)/m1 k2/m1 c2/m1;
     0.0 0.0 0.0 1.0;
     k2/m2 c2/m2 -k2/m2 -c2/m2]

# Initial conditions for mass 1 and mass 2
u = [x10; x1dot0; x20; x2dot0]
x1_position_numerical = [x10]
x2_position_numerical= [x20]
vel1 = [x1dot0]
vel2 = [x2dot0]

for (i, t) in enumerate(time[2:end]) 
    global u
    force = [0; f(t)/m1; 0; f(t)/m2]  # External force vector

    # Euler's method update
    uprime = A * u + force
    u_new = u + dt * uprime  

    # Store results
    push!(x1_position_numerical, u_new[1]) 
    push!(x2_position_numerical, u_new[3]) 
    push!(vel1, u_new[2])   
    push!(vel2, u_new[4])   

    u = u_new
end

# ----------------------------------------------
# ANALYTICAL SOLUTION (Using SymPy.jl)
# ----------------------------------------------
t_sym = symbols("t")
x1_sym = SymFunction("x1")(t_sym)  
x2_sym = SymFunction("x2")(t_sym)  

eq1 = Eq(m1*diff(x1_sym, t_sym, t_sym), -k1*x1_sym - c1*diff(x1_sym, t_sym) + k2*(x2_sym - x1_sym) + c2*(diff(x2_sym, t_sym) - diff(x1_sym, t_sym)) + f(t_sym))
eq2 = Eq(m2*diff(x2_sym, t_sym, t_sym), -k2*(x2_sym - x1_sym) - c2*(diff(x2_sym, t_sym) - diff(x1_sym, t_sym)) + f(t_sym))

analytical_sol = dsolve([eq1, eq2])

x1_general = analytical_sol[1].rhs
x2_general = analytical_sol[2].rhs

ics = (
    Eq(x1_general.subs(t_sym, 0), x10),       
    Eq(diff(x1_general, t_sym).subs(t_sym, 0), x1dot0), 
    Eq(x2_general.subs(t_sym, 0), x20),        
    Eq(diff(x2_general, t_sym).subs(t_sym, 0), x2dot0)   
)

constants = solve(ics)

x1 = x1_general.subs(constants)
x2 = x2_general.subs(constants)

x1_position_analytical = [x1(t_val) for t_val in time]
x2_position_analytical = [x2(t_val) for t_val in time]

# ----------------------------------------------
# PLOTTING: Compare ODE, Numerical, and Analytical
# ----------------------------------------------

# Plot results - Mass 1 comparison
plot1 = plot(ode_time, x1_position_ode, label="Mass 1 (ODE)", color=:blue, lw=2, 
    title="Mass 1 Position", xlabel="Time (s)", ylabel="Position", legend=:topright)
plot!(plot1, time, x1_position_numerical, label="Mass 1 (Numerical)", color=:red, lw=2)
plot!(plot1, time, x1_position_analytical, label="Mass 1 (Analytical)", color=:green, lw=2)

# Plot results - Mass 2 comparison
plot2 = plot(ode_time, x2_position_ode, label="Mass 2 (ODE)", color=:blue, lw=2,
    title="Mass 2 Position", xlabel="Time (s)", ylabel="Position", legend=:topright)
plot!(plot2, time, x2_position_numerical, label="Mass 2 (Numerical)", color=:red, lw=2)
plot!(plot2, time, x2_position_analytical, label="Mass 2 (Analytical)", color=:green, lw=2)

#RMSE
rmse1_numvsana = sqrt.((x1_position_numerical .- x1_position_analytical).^2)
rmse1_odevsana = sqrt.((x1_position_ode .- x1_position_analytical).^2)

p1 = plot(time, rmse1_numvsana, label="Numerical vs Analytical", color=:red, lw=2,
           title="RMSE between Numerical and Analytical Solutions for X1",
           xlabel="Time (s)", ylabel="RMSE",
           legend=:topright, grid=true)
plot!(p1, time, rmse1_odevsana, label="ODE vs Analytical", color=:blue, lw=2)

rmse2_numvsana = sqrt.((x2_position_numerical .- x2_position_analytical).^2)
rmse2_odevsana = sqrt.((x2_position_ode .- x2_position_analytical).^2)

p2 = plot(time, rmse2_numvsana, label="Numerical vs Analytical", color=:red, lw=2,
           title="RMSE between Numerical and Analytical Solutions for X2",
           xlabel="Time (s)", ylabel="RMSE",
           legend=:topright, grid=true)
plot!(p2, time, rmse2_odevsana, label="ODE vs Analytical", color=:blue, lw=2)


plot(plot1, p1, plot2, p2, layout=(2, 2), size=(1400, 1000))
