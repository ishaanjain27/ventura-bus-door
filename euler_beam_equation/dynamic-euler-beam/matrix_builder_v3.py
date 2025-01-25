import numpy as np

# reference: https://fenbildergi.aku.edu.tr/wp-content/uploads/2016/12/035601-693-710.pdf


# Builds "raw" N+1 size 4th derivative matrix
def A_matrix_builder(N, dx):
    """
    k1 = (EI*6)/dx**4 
    k2 = (EI*-4)/dx**4
    k3 = (EI*1)/dx**4
    """

    k1 = (6)/dx**4 
    k2 = (-4)/dx**4
    k3 = (1)/dx**4

    A_withoutBCs = np.zeros((N+1, N+1))
    
    np.fill_diagonal(A_withoutBCs[:-1, 1:], k2)
    np.fill_diagonal(A_withoutBCs[1:, :-1], k2)
    np.fill_diagonal(A_withoutBCs[:-2, 2:], k3)
    np.fill_diagonal(A_withoutBCs[2:, :-2], k3)
    B = k1 * np.ones((N+1))
    np.fill_diagonal(A_withoutBCs, B)    
    
    return A_withoutBCs, k1, k2, k3

def BC_left(A, dx, EI_vec, text, k1, k2, k3, spring_stiffness_left=None):
    if text == "cl":
        A[0,1] += k2
        A[0,2] += k3
        A[1,1] += k3
        A[0,:] = 0
        #A[:,0] = 0
        A[0,0] = 1
    elif text == "fr":
        A[0,0] += k2 + 4*k3
        A[0,1] += k2 - 4*k3
        A[0,2] += -k2 + k3
        A[1,0] += k3
        A[1,1] += k3
        A[1,2] += -k3
    elif text == "ss":
        A[0,1] += -k2
        A[0,2] += -k3
        A[1,1] += -k3
        A[0,:] = 0
        A[:,0] = 0
        A[0,0] = 1
    elif text == "spring":
        A[0,0] += k2 + 4*k3 - k2*dx**3*spring_stiffness_left/EI_vec[0] + 2*k3*dx**3*spring_stiffness_left/EI_vec[0]
        A[0,1] += k2 - 4*k3
        A[0,2] += -k2 + k3
        A[1,0] += k3 - spring_stiffness_left*dx**3*k3/EI_vec[1]
        A[1,1] += k3
        A[1,2] += -k3
    return A

def BC_right(A, dx, EI_vec, text, k1, k2, k3, spring_stiffness_right=None):
    if text == "cl":
        A[-2,-2] += k3
        A[-1,-3] += k3
        A[-1,-2] += k2
        A[-1,:] = 0
        #A[:,-1] = 0
        A[-1,-1] = 1
    elif text == "fr":
        A[-1,-1] += k2 + 4*k3
        A[-1,-2] += k2 - 4*k3
        A[-1,-3] += -k2 + k3
        A[-2,-1] += k3
        A[-2,-2] += k3
        A[-2,-3] += -k3
    elif text == "ss":
        A[-2,-2] += -k3
        A[-1,-3] += -k3
        A[-1,-2] += -k2
        A[-1,:] = 0
        A[:,-1] = 0
        A[-1,-1] = 1
    elif text == "spring":
        A[-1,-1] += k2 + 4*k3 - k2*dx**3*spring_stiffness_right/EI_vec[-1] + 2*k3*dx**3*spring_stiffness_right/EI_vec[-1]
        A[-1,-2] += k2 - 4*k3
        A[-1,-3] += -k2 + k3
        A[-2,-1] += k3 - spring_stiffness_right*dx**3*k3/EI_vec[-2]
        A[-2,-2] += k3
        A[-2,-3] += -k3
    return A