import numpy as np

# reference: https://fenbildergi.aku.edu.tr/wp-content/uploads/2016/12/035601-693-710.pdf

# Builds "raw" N+1 size 4th derivative matrix
def A_matrix_builder(N,EI,dx): 

    k1 = (EI*6)/dx**4 
    k2 = (EI*-4)/dx**4
    k3 = (EI*1)/dx**4

    A_withoutBCs = np.zeros((N+1,N+1))
    
    np.fill_diagonal(A_withoutBCs[:-1, 1:], k2)
    np.fill_diagonal(A_withoutBCs[1:, :-1], k2)
    np.fill_diagonal(A_withoutBCs[:-2, 2:], k3)
    np.fill_diagonal(A_withoutBCs[2:, :-2], k3)
    B = k1 * np.ones((N+1))
    np.fill_diagonal(A_withoutBCs,B)    
    
    return  A_withoutBCs, k1, k2, k3


# Left side BC
def BC_left(A,text,k1,k2,k3):
    if text == "cl":
        A[0,1] += k2
        A[0,2] += k3
        A[1,1] += k3

        A[0,:] = 0
        A[:,0] = 0
        A[0,0] = 1

    if text == "fr":
        A[0,0] += k2 + 4*k3
        A[0,1] += k2 - 4*k3
        A[0,2] += -k2 + k3
        A[1,0] += k3
        A[1,1] += k3
        A[1,2] += -k3

    if text == "ss":
        A[0,1] += -k2
        A[0,2] += -k3
        A[1,1] += -k3
        
        A[0,:] = 0
        A[:,0] = 0
        A[0,0] = 1
        
    return A

# Right side BC
def BC_right(A,text,k1,k2,k3):
    if text == "cl":
        A[-2,-2] += k3
        A[-1,-3] += k3
        A[-1,-2] += k2

        A[-1,:] = 0
        A[:,-1] = 0
        A[-1,-1] = 1

    if text == "fr":
        A[-1,-1] += k2 + 4*k3
        A[-1,-2] += k2 - 4*k3
        A[-1,-3] += -k2 + k3
        A[-2,-1] += k3
        A[-2,-2] += k3
        A[-2,-3] += -k3

    if text == "ss":
        A[-2,-2] += -k3
        A[-1,-3] += -k3
        A[-1,-2] += -k2
        A[-1,:] = 0
        A[:,-1] = 0
        A[-1,-1] = 1
    
    return A
