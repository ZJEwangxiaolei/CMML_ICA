import numpy as np
from make_segments import make_segments
from solve_for_flow import solve_for_flow
from realign_polarity import realign_polarity
from cell_migration import cell_migration,apply_branch_rule
from plot_network import plot_network

# Set random seed for reproducibility
np.random.seed(123456789)

# Input parameters
Nt = 20  # Number of time steps
Pin = 4 * 98  # Inlet pressure (Pa)
Pout = 1 * 98  # Outlet pressure (Pa)

mu = 3.5e-3  # Dynamic viscosity of blood (Pa-s)
Nn = 40  # Number of nodes
Nseg = 40  # Number of segments
num_cell = 10  # Initial number of cells per segment
cell_size = 5e-6  # Size of each cell (m)

branch_rule = 1  # Branching rule
branch_alpha = 0.5  # For BR5

# Polarization re-alignment weights
w2 = 1  # Flow component weight
w3 = 0.00  # Neighbor re-alignment weight
w4 = 0.00  # Random re-alignment weight
w1 = 1 - w2 - w3 - w4  # Persistence component

def compute_conductance(Nseg, Ncell, cell_size, mu, L):
    D = np.zeros(Nseg)
    G = np.zeros(Nseg)
    for seg in range(Nseg):
        D[seg] = max(Ncell[seg] * cell_size / np.pi, 1e-6)
        G[seg] = (np.pi * D[seg]**4) / (128 * mu * L[seg]) # Hagen-Poiseuille conductance
    return D, G

# Initialize segment cell structures
def initialize_segments(Nseg, num_cell):
    seg_cells = [{} for _ in range(Nseg)]
    for seg in range(Nseg):
        seg_cells[seg]['num'] = int(num_cell)
        seg_cells[seg]['polarity'] = [np.random.randn(2) for _ in range(num_cell)]
        for v in seg_cells[seg]['polarity']:
            v /= np.linalg.norm(v)
        seg_cells[seg]['migration'] = [0] * num_cell
    return seg_cells


L = np.ones(Nseg) * 10e-6  # Segment length (m)
vessels = make_segments(L)
Ncell = np.ones(Nseg) * num_cell
seg_cells = initialize_segments(Nseg, num_cell)

# Compute initial segment conductance
D, G = compute_conductance(Nseg, Ncell, cell_size, mu, L)
P, Q = solve_for_flow(G, Pin, Pout)
tau = (8 * mu * Q) / (np.pi * D**3)  # Initial shear stress (tau = 4 * mu * |Q| / (pi * D^3))


print("Initial network state:")
# plot_network(vessels, D, P, Q, seg_cells, tau)

# Time stepping for migration process
for t in range(Nt):
    print(f'Time step {t+1}/{Nt}')
    migrate = np.zeros(Nseg)
    new_seg_cells = [dict(seg_cells[i]) for i in range(Nseg)]  
    
    for seg in range(Nseg):
        seg_cells, new_seg_cells = realign_polarity(seg, Q, seg_cells, new_seg_cells, w1, w2, w3, w4)
        seg_cells, new_seg_cells = cell_migration(seg, seg_cells, new_seg_cells, migrate, Q, branch_rule, branch_alpha, tau, cell_size)

        # print(f"Segment {seg}: cells before migration = {seg_cells[seg]['num']}")
        # print(f"Segment {seg}: cells after migration = {new_seg_cells[seg]['num']}")
    
    # 
    #for seg in range(Nseg):
    #    assert len(new_seg_cells[seg]['polarity']) == new_seg_cells[seg]['num'], f"Seg {seg}: polarity length mismatch"
    #    assert len(new_seg_cells[seg]['migration']) == new_seg_cells[seg]['num'], f"Seg {seg}: migration length mismatch"
    
    seg_cells = new_seg_cells
    Ncell = np.array([seg_cells[seg]['num'] for seg in range(Nseg)])
    # print(Ncell)
    D, G = compute_conductance(Nseg, Ncell, cell_size, mu, L)
    P, Q = solve_for_flow(G, Pin, Pout)
    tau = (8 * mu * Q) / (np.pi * D**3)
    
    if (t + 1) % 2 == 0:
    #print(D,Q)
        plot_network(vessels, D, P, Q, seg_cells, tau)
    if t == 20:
        break

print("Simulation completed.")