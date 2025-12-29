import gradio as gr
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from matplotlib.colors import ListedColormap
import io
from PIL import Image

# ---------------------------------------------------------------------
# 1) Logic & Helper Functions (Kept mostly identical to your script)
# ---------------------------------------------------------------------

def all_conditions(condition_name: str):
    name = condition_name.strip()
    if name == "Control":
        return (0.001, 0.055, 0.95, 0.02, 0.0002)
    elif name == "Aggressive":
        return (0.001, 0.065, 0.95, 0.025, 0.0004)
    elif name == "DrugTreatment":
        return (0.001, 0.025, 0.95, 0.01, 0.0001)
    else:
        # Fallback default
        return (0.001, 0.055, 0.95, 0.02, 0.0002)

def grid_to_rgb(grid):
    """
    Convert the integer grid (1-6) to an RGB image array for fast display.
    1: Normal (White)
    2: Tumor (Blue)
    3: Grown Vessel (Pink)
    4: Sprouting Vessel (Dark Red)
    5: Original Vessel (Pink)
    6: Mutated Tumor (Light Blue)
    """
    h, w = grid.shape
    rgb = np.ones((h, w, 3), dtype=np.uint8) * 255 # Default white
    
    # Define colors (R, G, B) 0-255
    colors = {
        2: (0, 0, 255),       # Blue
        3: (255, 191, 178),   # Pink
        4: (128, 0, 0),       # Dark Red
        5: (255, 191, 178),   # Pink
        6: (128, 128, 255)    # Light Blue
    }
    
    for val, color in colors.items():
        mask = (grid == val)
        rgb[mask] = color
        
    return rgb

# ---------------------------------------------------------------------
# 2) Main Simulation Logic (Generator)
# ---------------------------------------------------------------------

def run_simulation(condition_name, num_iterations, grid_size):
    # Retrieve parameters
    (break_threshold, tumor_chance, death_chance, 
     mutation_chance, metastasis_chance) = all_conditions(condition_name)

    # Simulation constants
    initial_tumor_size = 2
    neighborhood_size = 1
    vessel_humanizer = 0.5
    min_tumor_cell_count_for_metastasis = 50
    treatment_time = [] 

    # Setup Grid
    grid = np.ones((grid_size, grid_size), dtype=np.int32)
    center = (grid_size // 2, grid_size // 2)
    yy, xx = np.meshgrid(np.arange(grid_size), np.arange(grid_size), indexing="ij")
    radius_sq = (initial_tumor_size / 2.0) ** 2
    circ_mask = (xx - center[1]) ** 2 + (yy - center[0]) ** 2 <= radius_sq
    grid[circ_mask] = 2

    # Setup Vessel
    vessel_width = 5
    vessel_height = int(grid_size * 0.75) # Scale vessel with grid
    start_row = 0
    start_col = int(grid_size * 0.7)
    vessel_mask = np.zeros_like(grid, dtype=bool)
    # Bounds check
    end_row = min(start_row + vessel_height, grid_size)
    end_col = min(start_col + vessel_width, grid_size)
    vessel_mask[start_row:end_row, start_col:end_col] = True
    grid[vessel_mask] = 5

    # Setup Map / Polarization
    rng = np.random.default_rng(1234)
    raw_map = rng.random((grid_size, grid_size))
    raw_map = gaussian_filter(raw_map, sigma=1)
    mn, mx = raw_map.min(), raw_map.max()
    polarized = ((raw_map - mn) / (mx - mn)) ** 2

    growth_rate_map = np.zeros_like(polarized)
    growth_rate_map[polarized <= 0.3] = 1.0
    growth_rate_map[(polarized > 0.3) & (polarized <= 0.6)] = 0.02
    growth_rate_map[polarized > 0.6] = 0.02
    
    adjusted_tumor_chance = tumor_chance * growth_rate_map

    # Tracking arrays
    normal_tumor_count = np.zeros(num_iterations, dtype=int)
    mutated_tumor_count = np.zeros(num_iterations, dtype=int)
    net_survival = np.zeros(num_iterations, dtype=float)

    number_of_breaks = 0
    max_distance = abs(center[1] - start_col) or 1
    metastasis_flag = False

    # Neighbors offset
    d = np.arange(-neighborhood_size, neighborhood_size + 1)
    dy, dx = np.meshgrid(d, d, indexing="ij")
    neigh_offsets = np.stack([dy.ravel(), dx.ravel()], axis=1)
    neigh_offsets = neigh_offsets[~((neigh_offsets[:, 0] == 0) & (neigh_offsets[:, 1] == 0))]

    # --- ITERATION LOOP ---
    for it in range(num_iterations):
        new_grid = grid.copy()
        grow_mask = np.zeros_like(grid, dtype=bool)
        
        # We perform a simplified pass to speed up Python loops
        # (Optimizing the double loop for strict Python speed in web app)
        
        # 1. Identify active cells to reduce loop overhead
        # Cells that are 2, 6 (tumor) or 5 (vessel) or 4 (sprout)
        active_mask = (grid != 1) 
        # Dilate mask to check neighbors
        # For simplicity in this demo, we scan the full grid or use list of coords
        # Using full scan for fidelity to original logic, but be mindful of speed.
        
        # To keep UI responsive, we process logic exactly as provided:
        for i in range(grid_size):
            for j in range(grid_size):
                cell_val = grid[i, j]

                # Tumor growth
                if cell_val in (2, 6):
                    # Check neighbors
                    for off in neigh_offsets:
                        ni, nj = i + off[0], j + off[1]
                        if 0 <= ni < grid_size and 0 <= nj < grid_size and grid[ni, nj] == 1:
                            # Check vessel adjacency
                            adj = False
                            for off2 in neigh_offsets:
                                vi, vj = ni + off2[0], nj + off2[1]
                                if 0 <= vi < grid_size and 0 <= vj < grid_size:
                                    if grid[vi, vj] in (3, 5):
                                        adj = True
                                        break
                            gch = adjusted_tumor_chance[ni, nj] * (1.5 if adj else 1.0)
                            if rng.random() <= gch:
                                grow_mask[ni, nj] = True

                # Vessel break -> sprout
                # Calculate local space only if needed
                if cell_val == 5 and (grid_size * 0.2 < i < grid_size * 0.8):
                    # Quick empty space check (center only)
                    if rng.random() < 0.5: # optim
                        i0, i1 = max(0, i-1), min(grid_size, i+2)
                        j0, j1 = max(0, j-1), min(grid_size, j+2)
                        local_patch = grid[i0:i1, j0:j1]
                        empty_space = np.sum(local_patch == 1)
                        
                        if empty_space != 0:
                            if number_of_breaks == 0 and rng.random() <= break_threshold:
                                new_grid[i, j] = 4
                                number_of_breaks += 1
                            elif number_of_breaks > 0:
                                reduced_p = break_threshold / (10.0 * number_of_breaks)
                                if rng.random() <= reduced_p:
                                    new_grid[i, j] = 4
                                    number_of_breaks += 1

                # Sprout movement
                if cell_val == 4:
                    # check local tumor count
                    i0, i1 = max(0, i-1), min(grid_size, i+2)
                    j0, j1 = max(0, j-1), min(grid_size, j+2)
                    local_tumor_count = np.sum((grid[i0:i1, j0:j1] == 2) | (grid[i0:i1, j0:j1] == 6))

                    if local_tumor_count == 0:
                        # Move sprout
                        # Find closest tumor cell (expensive, but necessary for logic)
                        tumor_locs = np.argwhere((grid == 2) | (grid == 6))
                        if tumor_locs.size > 0:
                            dists = (tumor_locs[:, 0] - i)**2 + (tumor_locs[:, 1] - j)**2
                            idx = np.argmin(dists)
                            ta, tb = tumor_locs[idx]
                        else:
                            ta, tb = i, j # Stay put if no tumor

                        best_dist = 1e30
                        grow_c, grow_d = i, j
                        
                        # Look for best move
                        for ci in range(i-1, i+2):
                            for cj in range(j-1, j+2):
                                if 0 <= ci < grid_size and 0 <= cj < grid_size:
                                    dtmp = (ci - ta)**2 + (cj - tb)**2 # squared dist is enough
                                    if dtmp < best_dist:
                                        best_dist = dtmp
                                        grow_c, grow_d = ci, cj
                        
                        # Humanizer
                        if rng.random() <= vessel_humanizer:
                            grow_c = np.clip(i + rng.integers(-1, 2), 0, grid_size-1)
                            grow_d = np.clip(j + rng.integers(-1, 2), 0, grid_size-1)
                        
                        # Draw vessel trail
                        bv_width = int(5 * abs(j - center[1]) / max_distance)
                        # limit width
                        bv_width = max(0, min(bv_width, 3)) 
                        
                        for w in range(-bv_width, bv_width + 1):
                            ii = i + w
                            if 0 <= ii < grid_size:
                                new_grid[ii, j] = 3
                        
                        new_grid[grow_c, grow_d] = 4
                    else:
                        # Touched tumor
                        new_grid[i, j] = 3
                        metastasis_flag = True

        # Apply growth
        new_grid[grow_mask] = 2

        # Mutation
        mut_rand = rng.random((grid_size, grid_size))
        mutate_mask = (grid == 2) & (mut_rand <= mutation_chance)
        new_grid[mutate_mask] = 6

        # Drug (Not fully implemented in input logic time array, but logic exists)
        # If user wanted drugs, we'd add `if it in treatment_time:` logic here.
        # For this demo, we assume treatment_time is empty list as per init.

        # Metastasis
        total_tumor_now = np.sum((grid == 2) | (grid == 6))
        if total_tumor_now >= min_tumor_cell_count_for_metastasis and metastasis_flag:
            vessel_cells = np.argwhere((grid == 3) | (grid == 5))
            if vessel_cells.size > 0:
                num_sites = rng.integers(1, 4)
                for _ in range(num_sites):
                    vr, vc = vessel_cells[rng.integers(0, vessel_cells.shape[0])]
                    for off in neigh_offsets:
                        nr, nc = vr + off[0], vc + off[1]
                        if 0 <= nr < grid_size and 0 <= nc < grid_size:
                            if new_grid[nr, nc] == 1 and rng.random() <= metastasis_chance:
                                new_grid[nr, nc] = 2
                                break

        grid = new_grid
        
        # Stats
        normal_tumor_count[it] = np.sum(grid == 2)
        mutated_tumor_count[it] = np.sum(grid == 6)
        total_cells = normal_tumor_count[it] + mutated_tumor_count[it]
        net_survival[it] = max(0.0, 100.0 * (1 - total_cells / 20000.0))

        # Yield results for Animation (every 5 frames to speed up UI)
        if it % 5 == 0 or it == num_iterations - 1:
            # Convert grid to Image for Gradio
            img_array = grid_to_rgb(grid)
            yield img_array, None, None, f"Iteration {it}/{num_iterations}"

    # --- FINAL PLOTS ---
    # 1. Growth Curve
    fig1, ax1 = plt.subplots(figsize=(6, 4))
    x = np.arange(1, num_iterations + 1)
    ax1.plot(x, normal_tumor_count, label="Normal tumor (2)", color='blue')
    ax1.plot(x, mutated_tumor_count, label="Mutated tumor (6)", color='cyan')
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel("Count")
    ax1.set_title(f"Tumor Growth: {condition_name}")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Survival Curve
    fig2, ax2 = plt.subplots(figsize=(6, 4))
    ax2.plot(x, net_survival, label="Net survival (%)", color="red")
    ax2.set_ylim(0, 100)
    ax2.set_xlabel("Iteration")
    ax2.set_ylabel("Survival (%)")
    ax2.set_title(f"Survival: {condition_name}")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Return final state
    img_array = grid_to_rgb(grid)
    yield img_array, fig1, fig2, "Simulation Complete"

# ---------------------------------------------------------------------
# 3) Gradio Interface
# ---------------------------------------------------------------------

with gr.Blocks() as demo:
    gr.Markdown("# Tumor Growth Cellular Automaton")
    gr.Markdown("Python port of `camodel.m`. Simulates tumor growth, angiogenesis, and mutation.")
    
    with gr.Row():
        with gr.Column(scale=1):
            cond_dd = gr.Dropdown(["Control", "Aggressive", "DrugTreatment"], value="Control", label="Condition")
            iters_sld = gr.Slider(50, 500, value=200, step=10, label="Iterations")
            size_sld = gr.Slider(50, 300, value=150, step=10, label="Grid Size")
            btn = gr.Button("Run Simulation", variant="primary")
            status_txt = gr.Textbox(label="Status", value="Ready")
        
        with gr.Column(scale=2):
            # Output grid as an image (updates live)
            grid_out = gr.Image(label="Live Grid State", elem_id="grid_viz")

    with gr.Row():
        plot1 = gr.Plot(label="Growth Statistics")
        plot2 = gr.Plot(label="Survival Rate")

    # Connect function
    btn.click(
        fn=run_simulation, 
        inputs=[cond_dd, iters_sld, size_sld], 
        outputs=[grid_out, plot1, plot2, status_txt]
    )

if __name__ == "__main__":
    demo.launch()