# Tumor Growth Cellular Automaton

A Python-based cellular automaton simulation for modeling tumor growth, angiogenesis, and mutation dynamics. Features an interactive web interface built with Gradio for real-time visualization.

## Overview

This project is a Python port of `camodel.m`, implementing a 2D cellular automaton that simulates:

- **Tumor proliferation** — Cancer cells divide and spread based on local microenvironment conditions
- **Angiogenesis** — Blood vessels sprout and grow toward the tumor mass
- **Mutation** — Normal tumor cells can mutate into more aggressive variants
- **Metastasis** — Tumor cells can spread via the vascular system

## Features

- 🔬 **Three simulation conditions**: Control, Aggressive, and Drug Treatment scenarios
- 📊 **Real-time visualization**: Watch the tumor evolve iteration by iteration
- 📈 **Statistical tracking**: Growth curves and survival rate plots
- ⚡ **Configurable parameters**: Adjust grid size and iteration count
- 🌐 **Web interface**: Easy-to-use Gradio UI

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/tumor-ca-simulation.git
cd tumor-ca-simulation

# Install dependencies
pip install gradio numpy matplotlib scipy pillow
```

## Usage

```bash
python app.py
```

Then open your browser to the local URL displayed (typically `http://127.0.0.1:7860`).

### Interface Controls

| Parameter | Description | Range |
|-----------|-------------|-------|
| **Condition** | Simulation scenario | Control, Aggressive, DrugTreatment |
| **Iterations** | Number of simulation steps | 50–500 |
| **Grid Size** | Simulation grid dimensions | 50–300 |

## Simulation Conditions

| Condition | Tumor Growth | Death Rate | Mutation | Metastasis | Description |
|-----------|--------------|------------|----------|------------|-------------|
| **Control** | Moderate | High | Low | Low | Baseline tumor behavior |
| **Aggressive** | High | High | Medium | High | Fast-growing, highly metastatic |
| **DrugTreatment** | Low | High | Very Low | Very Low | Simulates therapeutic intervention |

## Cell Types & Visualization

| Color | Cell Type | Description |
|-------|-----------|-------------|
| ⬜ White | Normal tissue | Healthy cells |
| 🟦 Blue | Tumor | Primary tumor cells |
| 🩵 Light Blue | Mutated tumor | Mutated variant cells |
| 🩷 Pink | Blood vessel | Original and grown vasculature |
| 🟥 Dark Red | Sprouting vessel | Actively growing vessel tip |

## How It Works

1. **Initialization**: A small tumor seed is placed at the grid center; a blood vessel runs along one edge.

2. **Tumor Growth**: Each iteration, tumor cells probabilistically convert neighboring normal cells. Growth rate is modulated by a spatially-varying microenvironment map.

3. **Angiogenesis**: Vessel cells can "break" and sprout toward the tumor. Sprouts navigate using a combination of chemotaxis (toward tumor) and random movement.

4. **Mutation**: Tumor cells have a small probability of mutating into a more aggressive phenotype each iteration.

5. **Metastasis**: Once the tumor reaches a critical size and contacts vasculature, cells can seed new tumors along blood vessels.

## Output

The simulation produces:

- **Live grid visualization** — Updated every 5 iterations
- **Growth statistics plot** — Tracks normal and mutated tumor cell counts
- **Survival curve** — Estimates net survival based on tumor burden

## Project Structure

```
tumor-ca-simulation/
├── app.py          # Main application with simulation logic and Gradio UI
└── README.md       # This file
```

## Requirements

- Python 3.8+
- NumPy
- Matplotlib
- SciPy
- Pillow
- Gradio

## License

MIT License — see [LICENSE](LICENSE) for details.

## Acknowledgments

Based on the MATLAB cellular automaton model `camodel.m` for tumor growth simulation.
