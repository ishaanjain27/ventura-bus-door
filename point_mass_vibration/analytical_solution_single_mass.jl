using SymPy
using Plots

# Define the time variable
t_sym = symbols("t")
x_sym = SymFunction("x")(t_sym)

# Define parameters directly (no symbols here)
m = 1       # mass
c = 0.2     # damping coefficient
k = 5       # spring constant
F0 = 1      # Force amplitude
t0 = 10     # Time for the center of the Gaussian pulse
sigma = 0.01 # Small width for Gaussian pulse

@syms DiracDelta

# Define the Gaussian force function (symbolically using a small sigma)
f(t) = F0 * DiracDelta(t-t0)

# Define the differential equation (m * x''(t) + c * x'(t) + k * x(t) = f(t))
eq = Eq(m * diff(x_sym, t_sym, 2) + c * diff(x_sym, t_sym) + k * x_sym, f(t_sym))

# Solve the differential equation symbolically
x_solution = dsolve(eq)
x_general = rhs(x_solution)

# Define initial conditions
x0 = 1      # Initial position
x_dot0 = 0  # Initial velocity

# Apply initial conditions to solve for constants
ics = (
    Eq(x_general.subs(t_sym, 0), x0),
    Eq(diff(x_general, t_sym).subs(t_sym, 0), x_dot0)
)

# Solve for the constants (C1, C2)
constants = solve(ics)
final_solution = x_general.subs(constants)

# Print the analytical solution
println("Analytical solution: ", final_solution)

# Define the time domain (discrete points for numerical evaluation)
t_domain = range(0, stop=50, length=500)

# Evaluate the symbolic expression for displacement at each time value
x_vals = [float(final_solution.subs(t_sym, t_val)) for t_val in t_domain]

# Plot the displacement function
plot(t_domain, real(x_vals), label="Impulsive Force F(t) = 2 δ(t-10)", xlabel="Time (s)", ylabel="Displacement (m)", legend=:topright, grid=true, linewidth=2)
