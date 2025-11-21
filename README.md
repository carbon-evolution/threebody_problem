# Relativistic Three-Body Problem

A Python simulation of the Three-Body Problem that goes beyond Newtonian gravity to include **General Relativistic (GR)** effects using the Einstein-Infeld-Hoffmann (EIH) approximation.

## Features

- **Relativistic Physics**: Uses the **Einstein-Infeld-Hoffmann (EIH)** approximation (1st Post-Newtonian order) to simulate:
  - Time dilation
  - Spatial curvature
  - Nonlinear gravity effects
  - Orbital precession (like Mercury's perihelion advance)

- **Animated 3D Visualization**:
  - Real-time orbital motion animation (1-minute loops)
  - True 3D trajectories in all three spatial dimensions
  - Persistent trajectory trails showing complete orbital paths
  - Color-coded bodies with matching start markers

- **Interactive Controls**:
  - **Scroll Wheel**: Zoom in/out
  - **Left Click + Drag**: Rotate the 3D view
  - **Spacebar or Button**: Play/Pause animation
  - **Dynamic Axis Scaling**: Automatically expands view as bodies move apart

- **N-Body Support**: Generalized to simulate any number of bodies

## Physics

This simulation solves the N-body equations of motion derived from the **Einstein Field Equations**:

$$ G_{\mu\nu} = 8\pi G T_{\mu\nu} $$

Instead of simple Newtonian forces ($F = GmM/r^2$), it calculates the relativistic acceleration:

$$ \mathbf{a} = \mathbf{a}_{Newtonian} + \mathbf{a}_{1PN} $$

For a detailed explanation of the math and physics, see [RELATIVITY.md](RELATIVITY.md).

## Usage

1. Install dependencies:
   ```bash
   pip install numpy matplotlib scipy
   ```

2. Run the simulation:
   ```bash
   python threebody_problem.py
   ```

## Scenarios

The simulation runs two scenarios side-by-side:

### 1. Solar System (Sun + 3 Planets)
- A hierarchical system with a massive central "Sun" (mass 10.0) and three orbiting planets
- Demonstrates stable orbits with **Schwarzschild precession** (slow rotation of elliptical orbits)
- All bodies exhibit true 3D motion including the central star

### 2. Chaotic 3-Body (No Sun)
- Three bodies of comparable mass (1.0, 1.2, 0.9) interacting in free space
- Demonstrates the classic "Three-Body Problem" chaos
- Highly sensitive to initial conditions with unpredictable long-term behavior
- Non-linear relativistic effects are most pronounced

## Controls

- **Scroll Wheel**: Zoom in/out
- **Left Click + Drag**: Rotate the 3D view
- **Spacebar**: Toggle Play/Pause
- **Play/Pause Button**: Click to toggle animation

## Technical Details

- **Simulation Duration**: 375 time units (1-minute animation)
- **Integration Method**: RK45 adaptive step-size solver
- **Precision**: `rtol=1e-6`, `atol=1e-9` for accurate GR effects
- **Animation**: 3000 frames at 20ms intervals with continuous looping
- **Dynamic Scaling**: Starts zoomed in, automatically expands as needed

## License

MIT License