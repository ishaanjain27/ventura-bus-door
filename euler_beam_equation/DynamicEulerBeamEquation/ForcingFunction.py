import numpy as np 

def q_RHS_uniform(t, F, N, dx):

    q = (F/dx) * np.ones((N-1))
    
    return np.concatenate(([0], q, [0]))  # Ensure boundary values are 0


def q_RHS_diracdelta(t, F, N, dx):

    q = np.zeros(N-1)
    q[int(np.floor(N/2))] = F/dx *np.sin(5*t)
    
    return np.concatenate(([0], q, [0]))  # Ensure boundary values are 0


def q_RHS_sine(t, F, N):

    q = F*np.sin(4*np.pi*np.linspace(0,2*np.pi,N-1) *t)
    
    return np.concatenate(([0], q, [0]))  # Ensure boundary values are 0

