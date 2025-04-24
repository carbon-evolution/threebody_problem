# Three-Body Problem Simulation

A Python-based simulation of the classical three-body gravitational problem in 3D space, using numerical integration methods.

## Overview

This project simulates the motion of three bodies under mutual gravitational attraction. It's a classic problem in physics that has no general analytical solution and exhibits chaotic behavior.

## Features

- 3D visualization of three-body trajectories
- Configurable masses and initial velocities
- Multiple simulation scenarios:
  - Equal masses with balanced velocities
  - Random masses with random velocities
- Interactive 3D plots for trajectory analysis

## Requirements

- Python 3.6+
- NumPy
- Matplotlib
- SciPy

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/threebody-simulation.git
   cd threebody-simulation
   ```

2. Install the required dependencies:
   ```
   pip install -r requirements.txt
   ```

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

![image](https://github.com/user-attachments/assets/7188ec97-de42-4b51-ad3f-d8d2884aa6b4)

![image](https://github.com/user-attachments/assets/05b7281e-2813-4f4c-91a8-083c2fdb4bd2)

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
