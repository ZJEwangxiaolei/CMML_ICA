# Agent-Based Model of Endothelial Cell Migration

This project contains a Python implementation of an Agent-Based Model (ABM) (original code written in MATLAB, see LICENSE file in the parent directory), which simulates endothelial cell (EC) migration in response to blood flow within a bifurcating vessel network as below:
####################################################################
---------- distal branch
|        |
|        |
---------- proximal branch
|        |
|        |
inlet    outlet
####################################################################
This model follows an iterative approach based-on time steps where cells (the "agents") realign their polarity, migrate, and update the blood flow dynamics accordingly.

## How the Model Works

1. **Initialization**: 
   - The vessel network is created with a defined number of segments.
   - Each segment is assigned an initial number of endothelial cells.
   - The cells have an initial polarity vector that dictates their migration direction.

2. **Flow Computation**:
   - The pressure difference across the network drives blood flow.
   - Vessel wall Shear stress values are computed for each segment.

3. **Cell Migration**:
   - Cells move based on their polarity, flow influence, and random walk behavior.
   - The migration of cells is guided by different rules (e.g., persistence, flow alignment, or random choice).

4. **Realignment of Polarity**:
   - The cells adjust their polarity after each step to reflect flow changes and random re-orientation.

5. **Updating the Network**:
   - New network conductance values are computed based on the updated cell distribution.
   - The blood flow is recomputed, and the process is repeated for the defined number of time steps.

## Description of Each Python File


### 1. **abm_ec_simulation_v2.py** (Updated Simulation with Branching Rules)
   - A modified version of the main simulation with additional **branching rules**.
   - Allows for different cell migration behaviors depending on blood flow and segment population.

### 2. **solve_for_flow.py** (Solves Flow in the Vessel Network)
   - Computes blood flow and pressure distribution in the vessel network.
   - Uses **matrix operations** to solve for nodal pressures and segment flow rates.

### 3. **cell_migration.py** (Handles Endothelial Cell Migration)
   - Determines how cells move within the network.
   - Migration decisions are based on:
     - Flow direction (e.g., against flow, or along flow)
     - Random movement (i.e., stochastic motion)
     - Persistence of previous movement at last time step

### 4. **realign_polarity.py** (Adjusts Cell Polarity Vectors)
   - Commands the cells to align their polarity based on **flow direction** and **random influences**.
   - Ensures directional migration patterns.

### 5. **plot_network.py** (Visualizes the Vessel Network)
   - Creates plots to display the:
     - Vessel structure.
     - Blood pressure and flow distribution.
     - Endothelial cell polarity vectors.

### 6. **make_segments.py** (Creates the Vessel Network Segments)
   - Defines the structure of the bifurcating vessel network.
   - Assigns positions to vessel segments.

## Running the Model

1. Ensure you have **Python 3+** installed.
2. Install dependencies using:
   ```sh
   pip install numpy matplotlib scipy
   ```


For an alternative version with enhanced **branching behavior**, run:
   ```sh
   python abm_ec_simulation_v2.py
   ```
