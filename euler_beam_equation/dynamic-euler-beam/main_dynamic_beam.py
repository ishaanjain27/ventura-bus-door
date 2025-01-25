import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from forcing_functions import *
from modal_analysis import *
from matrix_builder_v3 import *

#Equation being modeled: https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory 
#(and we also added a damping term to the RHS: -c * dw/dt)

##########################
#    GLOBAL VARIABLES    #
##########################

# Mesh Parameters
N = 300
beam_length = 10 
dx = beam_length/N

# Material Properties (vector for young's modulus * moment of inertia at each node)
EI_vec = 20*np.ones(N+1)

MEW = 1  # mass per unit length
c = 0.5    # damping coefficient

# Spring constants
spring_stiffness_left = 4
spring_stiffness_right = 10

# Time marching
TStart = 0 
TEnd = 100
dt = 0.0001

# Force Magnitude
F = 0.001

##########################
#   INITIAL CONDITIONS   # 
##########################

u0 = np.zeros((N+1)) 
uprime = np.zeros((N+1))

##########################
#     SYSTEM MATRIX A    # 
##########################

# Make raw 4th derivative matrix
A, k1, k2, k3 = A_matrix_builder(N, dx)

# Apply boundary conditions (HERE YOU SPECIFY THE BOUNDARY CONDITIONS- cl, fr, ss, spring)
A = BC_left(A, dx, EI_vec,"spring", k1, k2, k3, spring_stiffness_left)
A_withBCs = BC_right(A, dx, EI_vec, "spring", k1, k2, k3, spring_stiffness_right)
A_withBCs = A_withBCs * EI_vec[:, np.newaxis]

modal_analysis(A_withBCs, N)



##########################
#    FE TIME MARCHING    # 
##########################

t = TStart
u = u0.copy() + dt * uprime
uminus1 = u0.copy()

solution_matrix = [u0.copy()]
t_lst = [t]

while t <= TEnd:
    # Calculate velocity using central difference
    velocity = (u - uminus1) / dt
    
    # Time-marching update (Forward Euler) with damping term (CHANGE q_RHS_uniform to desired forcing function!)
    uplus1 = (-dt**2 / MEW) * (A_withBCs.dot(u) - q_RHS_uniform(t, F, N, dx) + c * velocity) + 2*u - uminus1

    solution_matrix.append(u)
    t_lst.append(t)

    # Update time-stepping variables
    uminus1 = u
    u = uplus1
    t += dt

##########################
#       GRAPHING         # 
##########################

# Animation
fig, ax = plt.subplots()
x = np.linspace(0, beam_length, N + 1)
line, = ax.plot([], [], lw=2)

# Set plot limits
ax.set_xlim(0, beam_length)
ax.set_ylim(-0.2, 0.2)
ax.set_xlabel('Beam Length (X)')
ax.set_ylabel('Displacement (u)')

def init():
    line.set_data([], [])
    return line,

def update(frame):
    line.set_data(x, solution_matrix[frame])
    ax.set_title(f"Vibrating Beam Over Time {t_lst[frame]:.3f}")
    return line,

# Create the animation
frames_to_show = range(0, len(t_lst), 2000)
ani = FuncAnimation(fig, update, frames=frames_to_show, init_func=init, interval=1)

plt.show()
