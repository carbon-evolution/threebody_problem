import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import random
from mpl_toolkits.mplot3d import Axes3D

# Constants
G = 1  # Gravitational constant (simplified)

def equations(t, y, m1, m2, m3):
    # Unpack positions and velocities (now including z coordinates)
    x1, y1, z1, vx1, vy1, vz1, x2, y2, z2, vx2, vy2, vz2, x3, y3, z3, vx3, vy3, vz3 = y
    
    # 3D distance calculations
    r12 = ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)**1.5
    r13 = ((x3 - x1)**2 + (y3 - y1)**2 + (z3 - z1)**2)**1.5
    r23 = ((x3 - x2)**2 + (y3 - y2)**2 + (z3 - z2)**2)**1.5

    # 3D Gravitational forces
    ax1 = G * m2 * (x2 - x1) / r12 + G * m3 * (x3 - x1) / r13
    ay1 = G * m2 * (y2 - y1) / r12 + G * m3 * (y3 - y1) / r13
    az1 = G * m2 * (z2 - z1) / r12 + G * m3 * (z3 - z1) / r13

    ax2 = G * m1 * (x1 - x2) / r12 + G * m3 * (x3 - x2) / r23
    ay2 = G * m1 * (y1 - y2) / r12 + G * m3 * (y3 - y2) / r23
    az2 = G * m1 * (z1 - z2) / r12 + G * m3 * (z3 - z2) / r23

    ax3 = G * m1 * (x1 - x3) / r13 + G * m2 * (x2 - x3) / r23
    ay3 = G * m1 * (y1 - y3) / r13 + G * m2 * (y2 - y3) / r23
    az3 = G * m1 * (z1 - z3) / r13 + G * m2 * (z2 - z3) / r23

    return [vx1, vy1, vz1, ax1, ay1, az1, 
            vx2, vy2, vz2, ax2, ay2, az2,
            vx3, vy3, vz3, ax3, ay3, az3]

def simulate_system(m1, m2, m3, initial_velocities, title, ax):
    # Initial positions and velocities
    vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3 = initial_velocities
    
    initial_conditions = [
        -1.0, 0.0, 0.0,    # x1, y1, z1
        vx1, vy1, vz1,     # velocities for body 1
        1.0, 0.0, 0.0,     # x2, y2, z2
        vx2, vy2, vz2,     # velocities for body 2
        0.0, 0.5, 0.0,     # x3, y3, z3
        vx3, vy3, vz3      # velocities for body 3
    ]

    print(f"\nSystem parameters for {title}:")
    print(f"Masses: m1={m1:.2f}, m2={m2:.2f}, m3={m3:.2f}")
    print(f"Body 1: vx={vx1:.2f}, vy={vy1:.2f}, vz={vz1:.2f}")
    print(f"Body 2: vx={vx2:.2f}, vy={vy2:.2f}, vz={vz2:.2f}")
    print(f"Body 3: vx={vx3:.2f}, vy={vy3:.2f}, vz={vz3:.2f}")

    # Reduced simulation parameters
    t_span = (0, 50)
    t_eval = np.linspace(0, 50, 5000)
    
    solution = solve_ivp(
        lambda t, y: equations(t, y, m1, m2, m3),
        t_span,
        initial_conditions,
        t_eval=t_eval,
        method='RK45'
    )

    # Extract positions
    x1, y1, z1 = solution.y[0], solution.y[1], solution.y[2]
    x2, y2, z2 = solution.y[6], solution.y[7], solution.y[8]
    x3, y3, z3 = solution.y[12], solution.y[13], solution.y[14]

    # Plot trajectories
    ax.plot(x1, y1, z1, label=f"Body 1 (m={m1:.2f})", linewidth=0.8)
    ax.plot(x2, y2, z2, label=f"Body 2 (m={m2:.2f})", linewidth=0.8)
    ax.plot(x3, y3, z3, label=f"Body 3 (m={m3:.2f})", linewidth=0.8)
    ax.scatter([-1, 1, 0], [0, 0, 0], [0, 0, 0], color='red', s=100, label='Initial positions')
    
    ax.set_xlabel("X Position")
    ax.set_ylabel("Y Position")
    ax.set_zlabel("Z Position")
    ax.set_title(title)
    ax.set_xlim(-15, 15)
    ax.set_ylim(-15, 15)
    ax.set_zlim(-15, 15)
    ax.view_init(elev=30, azim=45)
    ax.grid(True, linestyle='--', alpha=0.3)

# Create figure
fig = plt.figure(figsize=(20, 8))

# Equal masses and balanced velocities for interesting interactions
ax1 = fig.add_subplot(121, projection='3d')
equal_velocities = [0.0, 1.0, 0.0,    # v1: moving up
                   0.0, -1.0, 0.0,    # v2: moving down
                   1.0, 0.0, 0.0]     # v3: moving right
simulate_system(1.0, 1.0, 1.0, equal_velocities, 
               "Equal Masses & Balanced Velocities", ax1)

# Random masses and random velocities
ax2 = fig.add_subplot(122, projection='3d')
random_velocities = [random.uniform(-1.5, 1.5) for _ in range(9)]
m1 = random.uniform(0.5, 2.0)
m2 = random.uniform(0.5, 2.0)
m3 = random.uniform(0.5, 2.0)
simulate_system(m1, m2, m3, random_velocities,
               "Random Masses & Random Velocities", ax2)

plt.tight_layout()
plt.show()
