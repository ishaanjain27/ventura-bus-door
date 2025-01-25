import numpy as np 

# Define the q(x,t) function for the RHS of the PDE

def q_RHS_uniform(t, F, N, dx):

    q = F * np.ones((N-1))
    #print(f"q = {q}")

    return np.concatenate(([0], q, [0]))  # Ensure boundary values are 0 (might be unnecessary idk)


def q_RHS_diracdelta(t, F, N, dx):

    q = np.zeros(N-1)
    q[int(np.floor(N/2))] = F/dx *np.sin(2*t)
    #print(f"q = {q}")

    return np.concatenate(([0], q, [0]))  # Ensure boundary values are 0 (might be unnecessary idk)


def q_RHS_sine(t, F, N):

    q = F*np.sin(4*np.pi*np.linspace(0,2*np.pi,N-1) *t)
    #print(f"q = {q}")


    return np.concatenate(([0], q, [0]))  # Ensure boundary values are 0 (might be unnecessary idk)

