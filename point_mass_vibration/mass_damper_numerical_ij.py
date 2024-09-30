import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# first doing it in python and then converting to julia

# define domain and mesh
t_start = 0
t_end = 5
n_intervals = 100000
h = (t_end-t_start) / n_intervals
t_domain = np.linspace(t_start, t_end-h, n_intervals)

# define constants in si units
m = 1
c = 0.2
k = 5

# define system matrix
R = m/(h*h) + c/h + k
S = -2*m/(h*h) - c/h
T = m/(h*h)

A = sp.sparse.diags([T, S, R], [-2, -1, 0], shape=(n_intervals, n_intervals)).todense()
b = np.zeros(n_intervals)
b[0] = -1*(T+S)
b[1] = -1*T

# solve
x = np.linalg.solve(A, b)

# plot
plt.figure()
plt.plot(t_domain, x)
plt.grid(True)
plt.show()