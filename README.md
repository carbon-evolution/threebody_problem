# Relativistic Three-Body Problem

A Python simulation of the Three-Body Problem that goes beyond Newtonian gravity to include **General Relativistic (GR)** effects.

## Features

- **Relativistic Physics**: Uses the **Einstein-Infeld-Hoffmann (EIH)** approximation (1st Post-Newtonian order) to simulate time dilation, spatial curvature, and nonlinear gravity.
- **Interactive 3D Plotting**: Zoom and pan around the trajectories using the mouse scroll wheel.
- **Precession**: Observe orbital precession (like Mercury's perihelion advance) in hierarchical systems.

## Physics

This simulation solves the N-body equations of motion derived from the **Einstein Field Equations**:

$$ G_{\mu\nu} = 8\pi G T_{\mu\nu} $$

Instead of simple Newtonian forces ($F = GmM/r^2$), it calculates the relativistic acceleration:

$$ \mathbf{a} = \mathbf{a}_{Newtonian} + \mathbf{a}_{1PN} $$

For a detailed explanation of the math and physics, see [RELATIVITY.md](RELATIVITY.md).

## Usage

Run the simulation with default parameters:

```
python programs/threebody_problem/threebody_problem.py
```

The simulation will display two side-by-side plots:
- Left: Equal masses with balanced initial velocities
- Right: Random masses with random initial velocities

## How It Works

The simulation uses the Runge-Kutta 4th-order method (through SciPy's `solve_ivp`) to numerically solve the coupled differential equations of motion for the three bodies. The gravitational force between any two bodies is calculated using Newton's law of universal gravitation.

## Physics Background

In the three-body problem, we have three masses exerting gravitational forces on each other. The equations of motion are derived from Newton's laws and the gravitational force law:

F = G * m₁ * m₂ / r²

Where:
- F is the force between the masses
- G is the gravitational constant
- m₁ and m₂ are the masses
- r is the distance between the centers of the masses

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Feel free to open an issue or submit a pull request. 