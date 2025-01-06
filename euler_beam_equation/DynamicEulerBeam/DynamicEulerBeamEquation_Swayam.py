import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from MatrixBuilder import *
from ForcingFunction import *
from ModalAnalysis import modal_analysis



##########################
#    GLOBAL VARIABLES    #
##########################

# Mesh Parameters
N = 300
beam_length = 10 
dx = beam_length/N

E_MODULUS = 200*10**3
INERTIA = 0.000001650388 
EI = 20
MEW  = 8.3 # mass per unit length

# Time marching
TStart = 0 
TEnd = 10
dt = 0.000001

# Force Magnitude
F = 300


##########################
#   INITIAL CONDITIONS   # 
##########################

u0 = np.zeros((N+1)) 
uprime = np.zeros((N+1))


##########################
#     SYSTEM MATRIX A    # 
##########################

A_withoutBCs = A_matrix_builder(N,E_MODULUS,INERTIA,dx)
A_withBCs = BC_clamped_clamped(A_withoutBCs)
print(A_withBCs)

A_withBCs = A_withBCs*((EI)/dx**4) 
#A_withBCs = sparse.csr_matrix(A_withBCs)

##########################
#    FE TIME MARCHING    # 
##########################

t = TStart
u = u0.copy() + dt * uprime
uminus1 = u0.copy()

solution_matrix = [u0.copy()]  # Start with initial condition
t_lst = [t]

while t <= TEnd:
    # Time-marching update (Forward Euler)
    uplus1 = (-dt**2 / MEW) * (A_withBCs.dot(u) - q_RHS_uniform(t, F, N, dx)) + 2*u - uminus1

    solution_matrix.append(u)
    t_lst.append(t)

    # Update time-stepping variables
    uminus1 = u
    u = uplus1
    t += dt



##########################
#       GRAPHING         # 
##########################

#modal_analysis(A_withBCs,N)

# Animation
fig, ax = plt.subplots()
x = np.linspace(0, beam_length, N + 1)
line, = ax.plot([], [], lw=2)  # Line plot

# Set plot limits
ax.set_xlim(0, beam_length)
ax.set_ylim(min(solution_matrix[-1]),max(solution_matrix[-1]))
ax.set_xlabel('Beam Length (X)')
ax.set_ylabel('Displacement (u)')


def init():
    """Initialize the plot."""
    line.set_data([], [])
    return line,

def update(frame):
    """Update the plot for each frame."""
    line.set_data(x, solution_matrix[frame])
    ax.set_title(f"Vibrating Beam Over Time {t_lst[frame]}")
    return line,

# Create the animation
frames_to_show = range(0, len(t_lst), 6000)  # Show every 1000th frame
ani = FuncAnimation(fig, update, frames=frames_to_show, init_func=init, interval=1)

plt.show()
