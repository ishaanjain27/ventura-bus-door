import numpy as np
import matplotlib.pyplot as plt

#Grid Parameters
N = 2000 #number of subdivisions
n = N-1 #number of interior nodes
L = 5.0
dx = L/N

EI = EI = 24

#Creating A matrix (second derivate) (nxn) TODO: Make it sparse
A = np.zeros((n,n))
np.fill_diagonal(A[:-1, 1:], -4)
np.fill_diagonal(A[1:, :-1], -4)
np.fill_diagonal(A[:-2, 2:], 1)
np.fill_diagonal(A[2:, :-2], 1)
B = 6 * np.ones((n))
np.fill_diagonal(A,B)

#Expand size of the matrix of nxn to (n+2 x n+2)
A = np.pad(A, pad_width=1, mode='constant', constant_values=0)

#Applying clamped clamped BC's
def clamped_clamped_bc(A):
    A[0,0] = 1
    A[-1,-1] = 1
    A[1,1] = 7
    A[-2,-2] = 7
    A[1,0] = -4 #ask abt this one?
    print(A)
    return A
#A = clamped_clamped_bc(A)

#Apply Canteliver beam BC's
def clamped_free_bc(A):
    A[0,0] = 1
    A[1,1] = 7
    A[1,0] = -4
    
    #free end aspects 
    #2nd last row using i = 4 and 2nd derivative BC
    A[-3,-1] = 1
    A[-2,-1] = -2
    A[-2,-2] = 5
    A[-2,-3] = -4
    A[-2,-4] = 1
    
    #last row using 3rd derivative BC
    A[-1,-1]  = 1
    A[-1,-2] = -3
    A[-1,-3] = 3
    A[-1,-4] = -1
    return A
A = clamped_free_bc(A)
print(A)

#RHS vector size n + 2
f = np.ones((n))
f = np.concatenate(([0], f, [0]))
#print(f)

#Solving
u_vector = (dx**4/EI) * np.dot(np.linalg.inv(A),f)
#print(u_vector)

#Graphing
plt.plot(np.arange(0,N+1,1),u_vector)


# Modal Analysis
eigenvalues, eigenvectors = np.linalg.eig(A)

# Normalize eigenvector signs
eigenvectors = eigenvectors.T  # Transpose for easier iteration over rows
for i in range(len(eigenvectors)):
    # Check the first nonzero element's sign and normalize
    first_nonzero = np.argmax(np.abs(eigenvectors[i]) > 1e-8)  # Find the first significant non-zero
    if eigenvectors[i, first_nonzero] < 0:
        eigenvectors[i] *= -1  # Flip the sign of the entire eigenvector
eigenvectors = eigenvectors.T  # Transpose back to original shape

# Sort eigenvalues and eigenvectors by frequency (sqrt(eigenvalue))
frequencies = np.sqrt(np.abs(eigenvalues))  # Frequencies are proportional to sqrt(eigenvalues)
sorted_indices = np.argsort(frequencies)    # Indices for sorted frequencies

# Sort eigenvalues and eigenvectors
eigenvalues_sorted = eigenvalues[sorted_indices]
eigenvectors_sorted = eigenvectors[:, sorted_indices]

# Plot the eigenvectors with the lowest frequencies
plt.figure(figsize=(10, 6))
for i in range(4):  # Plot the first three eigenvectors with the lowest frequencies
    plt.plot(np.arange(0, N+1, 1), eigenvectors_sorted[:, i], label=f"Eigenvector {i+1}, Frequency={frequencies[sorted_indices[i]]:.2f}")
    
plt.xlabel("Node Index")
plt.ylabel("Amplitude")
plt.title("Eigenvectors with Lowest Frequencies")
plt.legend()
plt.grid()
plt.show()