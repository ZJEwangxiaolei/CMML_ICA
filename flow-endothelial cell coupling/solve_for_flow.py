import numpy as np

def solve_for_flow(G, Pin, Pout, H=None):
    Nn = 41
    Nseg = 40
    G[G == 0] = 1e-25
    
    P = np.zeros(Nn)
    Q = np.zeros(Nseg)
    C = np.zeros((Nn, Nn))
    B = np.zeros(Nn)
    
    connections = [
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5),  # vessel1
        (5, 6), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 12), (12, 13), (13, 14), (14, 15),  # vessel2
        (15, 16), (16, 17), (17, 18), (18, 19), (19, 20),  # vessel3
        (5, 21), (21, 22), (22, 23), (23, 24), (24, 25),  # vessel4
        (25, 26), (26, 27), (27, 28), (28, 29), (29, 30), (30, 31), (31, 32), (32, 33), (33, 34), (34, 35),  # vessel5
        (35, 36), (36, 37), (37, 38), (38, 39), (39, 40)  # vessel6
    ]
    #(25, 35),(20, 26),
    C[0, 0] = G[0]
    B[0] = G[0] * Pin
    for seg, (i, j) in enumerate(connections):
        C[i, i] += G[seg]
        C[j, j] += G[seg]
        C[i, j] -= G[seg]
        C[j, i] -= G[seg]
    C[20, 20] = G[19]
    B[20] = G[19] * Pout
    C[40, 40] = G[39]
    B[40] = G[39] * Pout
    
    P = np.linalg.solve(C, B)
    for seg, (i, j) in enumerate(connections):
        Q[seg] = G[seg] * (P[j] - P[i])
    
    if H is not None:
        tau = H * Q
        return P, Q, tau
    return P, Q