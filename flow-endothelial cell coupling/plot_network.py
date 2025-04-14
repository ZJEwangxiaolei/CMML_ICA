import numpy as np
import matplotlib.pyplot as plt

def plot_network(segments, D, P, Q, seg_cells, tau=None):  # Add tau=None
    """Plot the vessel network along with pressure, flow, and cell polarity vectors."""
    
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.title('Pressure, Flow, Diameter of Network')
    
    seg_idx = 0
    for vessel in segments:  # 遍历每个血管段
        for i in range(len(vessel) - 1):
            if Q[seg_idx] > 0:
                color = "red"
            else:
                color = "blue"
            plt.plot([vessel[i, 0], vessel[i + 1, 0]], 
                     [vessel[i, 1], vessel[i + 1, 1]], 
                     color=color, linewidth=D[seg_idx] * 1e6 / 2)
            seg_idx += 1
    
    plt.grid()
    
    # Polarity distribution plot
    plt.subplot(1, 2, 2)
    plt.title('Distribution of Cell Polarity')
    plt.axis([-1, 1, -1, 1])
    plt.grid()
    
    for seg in range(len(seg_cells)):
        for cell in range(seg_cells[seg]['num']):
            polarity = seg_cells[seg]['polarity'][cell]
            plt.plot([0, polarity[0]], [0, polarity[1]], 'b-')
    
    plt.show()
