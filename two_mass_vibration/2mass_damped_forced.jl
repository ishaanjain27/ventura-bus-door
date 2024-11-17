using DifferentialEquations
using SymPy
using Statistics
using Printf  # For formatted output
using Plots

elapsed_time = @elapsed begin

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

    f(t) = 0.0

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

    # Measure computation time using @elapsed
    println("Solving ODE...")
    ode_sol = solve(ode, abstol = 1e-10, reltol = 1e-10)

    ode_time = ode_sol.t
    x1_position_ode = ode_sol[1, :]  
    x2_position_ode = ode_sol[2, :]

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

    x1_position_analytical = [x1(t_val) for t_val in ode_time]
    x2_position_analytical = [x2(t_val) for t_val in ode_time]

    # ----------------------------------------------
    # PLOTTING: Compare ODE, Numerical, and Analytical
    # ----------------------------------------------

    # Plot results - Mass 1 comparison
    plot1 = plot(ode_time, x1_position_ode, label="Mass 1 (Numerical)", color=:blue, lw=2, 
        title="Mass 1 Position", xlabel="Time (s)", ylabel="Displacement (m)", legend=:topright)
    plot!(plot1, ode_time, x1_position_analytical, label="Mass 1 (Analytical)", color=:green, lw=2)

    # Plot results - Mass 2 comparison
    plot2 = plot(ode_time, x2_position_ode, label="Mass 2 (Numerical)", color=:blue, lw=2,
        title="Mass 2 Position", xlabel="Time (s)", ylabel="Displacement (m)", legend=:topright)
    plot!(plot2, ode_time, x2_position_analytical, label="Mass 2 (Analytical)", color=:green, lw=2)

    # ----------------------------------------------
    # ERRORS
    # ----------------------------------------------

    # Absolute errors
    abs_error_x1 = abs.(x1_position_ode - x1_position_analytical)
    abs_error_x2 = abs.(x2_position_ode - x2_position_analytical)

    # Relative errors
    rel_error_x1 = abs.(x1_position_ode - x1_position_analytical) ./ abs.(x1_position_analytical)
    rel_error_x2 = abs.(x2_position_ode - x2_position_analytical) ./ abs.(x2_position_analytical)

    # Plot absolute error
    plot_abs_error = plot(ode_time, (abs_error_x1), label="Mass 1 Absolute Error", color=:red, lw=2,
        title="Absolute Error", xlabel="Time (s)", ylabel="Error", legend=:bottomright)
    plot!(plot_abs_error, ode_time, (abs_error_x2), label="Mass 2 Absolute Error", color=:orange, lw=2)

    # Plot relative error
    plot_rel_error = plot(ode_time, log10.(rel_error_x1), label="Mass 1 Relative Error", color=:red, lw=2,
        title="Relative Error", xlabel="Time (s)", ylabel="Log10(Error)", legend=:bottomright)
    plot!(plot_rel_error, ode_time, log10.(rel_error_x2), label="Mass 2 Relative Error", color=:orange, lw=2)

    # Calculate the time steps
    time_steps = diff(ode_sol.t)

    # Use the midpoints of the time intervals for plotting
    time_midpoints = ode_sol.t[1:end-1] .+ diff(ode_sol.t) / 2

    # Plot the time steps vs. time midpoints
    plot_time_steps = plot(time_midpoints, time_steps, label="Time Step", color=:purple, lw=2, xlabel="Time (s)", ylabel="Time Step (s)", legend=:topright)

    
end

@printf("Computation time for numerical solution: %.3f seconds\n", elapsed_time)

plot(plot_abs_error)
