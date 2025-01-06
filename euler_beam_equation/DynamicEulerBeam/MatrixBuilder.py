import numpy as np

# Builds "raw" N+1 size 4th derivative matrix
def A_matrix_builder(N,E_MODULUS,INERTIA,dx): 

    A_withoutBCs = np.zeros((N+1,N+1))
    
    np.fill_diagonal(A_withoutBCs[:-1, 1:], -4)
    np.fill_diagonal(A_withoutBCs[1:, :-1], -4)
    np.fill_diagonal(A_withoutBCs[:-2, 2:], 1)
    np.fill_diagonal(A_withoutBCs[2:, :-2], 1)
    B = 6 * np.ones((N+1))
    np.fill_diagonal(A_withoutBCs,B)    
    print(A_withoutBCs)
    return  A_withoutBCs


##########################
#  BOUNDARY CONDITIONS   #
##########################

def BC_clamped_free(A_withoutBCs):

    # Applying clamp to left side
    A_withoutBCs[0,:] = 0
    A_withoutBCs[1,:] = 0

    A_withoutBCs[0,0] = 1   # u(x=0) = 0

    A_withoutBCs[1,0] = -1
    A_withoutBCs[1,1] = 1   # u'(x=0) = 0

    # Applying free BC to right side
    A_withoutBCs[-1,:] = 0
    A_withoutBCs[-2,:] = 0

    #A_withoutBCs[-2,-1] = 1 
    #A_withoutBCs[-2,-2] = -2   # this part applies u''(x=N+1) = 0 using central dif
    #A_withoutBCs[-2,-3] = 1

    A_withoutBCs[-2, -4] = 1
    A_withoutBCs[-2, -3] = -4
    A_withoutBCs[-2, -2] = 5
    A_withoutBCs[-2, -1] = -2

    A_withoutBCs[-1,-1] = 1
    A_withoutBCs[-1,-2] = -3    # this part applies u'''(x=N+1) = 0 using backward dif
    A_withoutBCs[-1,-3] = 3
    A_withoutBCs[-1,-4] = -1

    A_withBCs = A_withoutBCs

    return A_withBCs 


def BC_clamped_clamped(A_withoutBCs):

    A_withoutBCs[0,:] = 0
    A_withoutBCs[1,:] = 0  # this part wipes two and bottom 2 rows to 0
    A_withoutBCs[-1,:] = 0 
    A_withoutBCs[-2,:] = 0
    
    A_withoutBCs[0,0] = 1  # this part applies u(x=0,x=N+1) = 0
    A_withoutBCs[-1,-1] = 1

    A_withoutBCs[1,0] = -1
    A_withoutBCs[1,1] = 1  # this part applies u'(x=0,x=N+1) = 0
    A_withoutBCs[-2,-1] = 1
    A_withoutBCs[-2,-2] = -1

    return A_withoutBCs 

