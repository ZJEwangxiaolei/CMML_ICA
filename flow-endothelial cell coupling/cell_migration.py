import numpy as np


def cell_migration(seg, seg_cells, new_seg_cells, migrate, Q, branch_rule, branch_alpha=None, tau=None, cell_size=5e-6):
    """Handle cellular migration in the agent-based model."""
    mchance = 1   # Assume full migration probability for now

    if seg_cells[seg]['num'] != 0:
        # Check if the segment contains cells
        cells_to_remove = []
        for cell in range(seg_cells[seg]['num']):
            mcell = np.random.rand()
            if mcell <= mchance:
                polar_vect = seg_cells[seg]['polarity'][cell]
                migrate_vect = cell_size * polar_vect
                
                # Migration logic of non-bifurcated segments
                if seg <= 4 or (21 <= seg <= 24):
                    # Vessel1 (0–4), Vessel4 (21–24): Vertical migration (up/down)
                    if migrate_vect[1] >= cell_size / 2:
                        # up
                        new_seg_cells[seg+1]['num'] += 1
                        new_seg_cells[seg+1]['polarity'].append(polar_vect)
                        new_seg_cells[seg+1]['migration'].append(0)
                        cells_to_remove.append(cell)
                    elif migrate_vect[1] <= -cell_size / 2:
                        # down
                        new_seg_cells[seg-1]['num'] += 1
                        new_seg_cells[seg-1]['polarity'].append(polar_vect)
                        new_seg_cells[seg-1]['migration'].append(0)
                        cells_to_remove.append(cell)
                        migrate[seg] += 1
                # Bifurcation points (seg == 5, 25)
                elif seg in [5, 25]:
                    target_seg = apply_branch_rule(seg, seg_cells, Q, tau, branch_rule, branch_alpha)
                    new_seg_cells[target_seg]['num'] += 1
                    new_seg_cells[target_seg]['polarity'].append(polar_vect)
                    new_seg_cells[target_seg]['migration'].append(0)
                    cells_to_remove.append(cell)
                elif 35 <= seg <= 38:
                    if migrate_vect[1] <= -cell_size / 2:  
                        new_seg_cells[seg+1]['num'] += 1
                        new_seg_cells[seg+1]['polarity'].append(polar_vect)
                        new_seg_cells[seg+1]['migration'].append(0)
                        cells_to_remove.append(cell)
                    elif migrate_vect[0] <= -cell_size / 2:  
                        new_seg_cells[seg-1]['num'] += 1
                        new_seg_cells[seg-1]['polarity'].append(polar_vect)
                        new_seg_cells[seg-1]['migration'].append(0)
                        cells_to_remove.append(cell)
                # new cells were added to segment 1 when migrating 40/41
                elif seg == 40 or seg ==39:
                    if migrate_vect[1] <= -cell_size / 2:
                        cells_to_remove.append(cell)  
                        
                        new_seg_cells[1]['num'] += 1
                        new_seg_cells[1]['polarity'].append(np.array([0, 1]))  
                        new_seg_cells[1]['migration'].append(0)
            # Migration of horizontal segments
                else:  
                    if migrate_vect[0] >= cell_size / 2:  
                        new_seg_cells[seg+1]['num'] += 1
                        new_seg_cells[seg+1]['polarity'].append(polar_vect)
                        new_seg_cells[seg+1]['migration'].append(0)
                        cells_to_remove.append(cell)
                    elif migrate_vect[0] <= -cell_size / 2: 
                        new_seg_cells[seg-1]['num'] += 1
                        new_seg_cells[seg-1]['polarity'].append(polar_vect)
                        new_seg_cells[seg-1]['migration'].append(0)
                        cells_to_remove.append(cell)
            
        # Remove cells after iteration to avoid index issues
        for cell in sorted(cells_to_remove, reverse=True):  
            new_seg_cells[seg]['num'] -= 1
            new_seg_cells[seg]['polarity'].pop(cell)
            new_seg_cells[seg]['migration'].pop(cell)
    
    return seg_cells, new_seg_cells

def apply_branch_rule(seg, seg_cells, Q, tau, branch_rule, branch_alpha=None):
    """Apply branching rule at bifurcation."""
    # Define bifurcation branches
    if seg == 5:
        proximal, distal = (6, 21)
    elif seg == 25:
         proximal, distal = (34, 35)
    else:
        proximal, distal = (26, 36)

    if branch_rule == 1:  # BR1: Higher shear stress
        return proximal if tau[proximal] > tau[distal] else distal
    elif branch_rule == 2: 
        return distal 
    elif branch_rule == 3: # BR3: Random (50% probability)
        return proximal if np.random.rand() < 0.5 else distal
    elif branch_rule == 4: # BR4: Biased random (70% proximal)
        return proximal if np.random.rand() < 0.7 else distal
    elif branch_rule == 5: # BR5: Weighted tau and cell count
        if tau is None or branch_alpha is None:
            raise ValueError("BR5 requires tau and branch_alpha")
        tau_prox, tau_dist = tau[proximal], tau[distal]
        n_prox, n_dist = seg_cells[proximal]['num'], seg_cells[distal]['num']
        P_prox = branch_alpha * (tau_prox / (tau_prox + tau_dist)) + \
                 (1 - branch_alpha) * (n_prox / (n_prox + n_dist))
        return proximal if np.random.rand() < P_prox else distal
    else:
        return proximal