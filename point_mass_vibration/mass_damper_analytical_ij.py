import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# variables
t = sp.Symbol('t')
x = sp.Function('x')(t)

# constants
m = sp.Symbol('m')
c = sp.Symbol('c')
k = sp.Symbol('k')

# mass damper de in 1 dimension
mass_damper = sp.Eq(m * x.diff(t, 2) + c * x.diff(t) + k * x, 0)

# solve the de
x_solution = sp.dsolve(mass_damper, x)
x_general = x_solution.rhs
x_prime = x_general.diff(t)

# boundary conditions
boundary_eq1 = sp.Eq(x_general.subs(t, 0), 1)
boundary_eq2 = sp.Eq(x_prime.subs(t, 0), 0)

# solve for constants
C1, C2 = sp.symbols('C1 C2')
solution = sp.solve([boundary_eq1, boundary_eq2], (C1, C2))
final_solution = x_general.subs(solution)

print(final_solution)

# plotting
t_domain = np.linspace(0, 10, 50)
x_numeric = sp.lambdify(t, final_solution, 'numpy')


# print the solution
