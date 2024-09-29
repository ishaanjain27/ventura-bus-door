import sympy as sm
import numpy as np
import matplotlib.pyplot as plt

# variables
t = sm.Symbol('t')
x = sm.Function('x')(t)

# constants
m = sm.Symbol('m')
c = sm.Symbol('c')
k = sm.Symbol('k')

# mass damper de in 1 dimension
mass_damper = sm.Eq(m * x.diff(t, 2) + c * x.diff(t) + k * x, 0)

# solve the de
x_solution = sm.dsolve(mass_damper, x)
x_general = x_solution.rhs
x_prime = x_general.diff(t)

# boundary conditions
boundary_eq1 = sm.Eq(x_general.subs(t, 0), 1)
boundary_eq2 = sm.Eq(x_prime.subs(t, 0), 0)

# solve for constants
C1, C2 = sm.symbols('C1 C2')
solution = sm.solve([boundary_eq1, boundary_eq2], (C1, C2))
final_solution = x_general.subs(solution)

# print the analytical solution
print(final_solution)

# numerical values of constants in si units
m_val = 1
c_val = 0.2
k_val = 5

# substitute constants
numerical_solution = final_solution.subs({m: m_val, c: c_val, k: k_val})

# plotting
t_domain = np.linspace(0, 50, 500)
x_numeric = sm.lambdify(t, numerical_solution, 'numpy')
x_vals = x_numeric(t_domain)

# Plot the function
plt.plot(t_domain, x_vals)
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('mass damper system with initial dismlacement')
plt.grid(True)
plt.show()
