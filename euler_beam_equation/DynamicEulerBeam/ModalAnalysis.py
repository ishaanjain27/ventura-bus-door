import numpy as np
import matplotlib.pyplot as plt


##########################
#     MODAL ANALYSIS     #
##########################

def modal_analysis(A,N):
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
    for i in range(6):  # Plot the first x eigenvectors with the lowest frequencies
        plt.plot(np.arange(0, N+1, 1), eigenvectors_sorted[:, i], label=f"Eigenvector {i+1}, Frequency={frequencies[sorted_indices[i]]:.2f}")
        print(frequencies[sorted_indices[i]])
        
        
    plt.xlabel("Node Index")
    plt.ylabel("Amplitude")
    plt.title("Eigenvectors with Lowest Frequencies")
    plt.legend()
    plt.grid()
    plt.show()

