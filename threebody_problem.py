import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import random

# Constants
G = 1  # Gravitational constant (simplified)

# Random masses between 0.5 and 2.0
m1 = random.uniform(0.5, 2.0)
m2 = random.uniform(0.5, 2.0)
m3 = random.uniform(0.5, 2.0)

print(f"Random masses: m1={m1:.2f}, m2={m2:.2f}, m3={m3:.2f}")

# Initial positions and velocities
initial_conditions = [
    -1.0, 0.0, 0.0, 1.0,  # x1, y1, vx1, vy1
     1.0, 0.0, 0.0, -1.0, # x2, y2, vx2, vy2
     0.0, 0.5, -1.0, 0.0  # x3, y3, vx3, vy3
]

def detect_repetition(positions, tolerance=0.1):
    """Check if the system has started repeating its pattern"""
    if len(positions) < 1000:  # Need enough points to detect pattern
        return False
    
    # Check every 500 points instead of 100 to reduce computation
    step = 500
    window = 1000  # Only look at last 1000 points for comparison
    
    # Get current position
    current_pos = positions[-1]
    
    # Compare with previous positions in chunks
    for i in range(len(positions)-window, len(positions)-step, step):
        prev_pos = positions[i]
        if np.all(np.abs(current_pos - prev_pos) < tolerance):
            return True
    
    return False

def equations(t, y):
    # Unpack positions and velocities
    x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3 = y
    
    # Distance calculations
    r12 = ((x2 - x1)**2 + (y2 - y1)**2)**1.5
    r13 = ((x3 - x1)**2 + (y3 - y1)**2)**1.5
    r23 = ((x3 - x2)**2 + (y3 - y2)**2)**1.5

    # Gravitational forces
    ax1 = G * m2 * (x2 - x1) / r12 + G * m3 * (x3 - x1) / r13
    ay1 = G * m2 * (y2 - y1) / r12 + G * m3 * (y3 - y1) / r13

    ax2 = G * m1 * (x1 - x2) / r12 + G * m3 * (x3 - x2) / r23
    ay2 = G * m1 * (y1 - y2) / r12 + G * m3 * (y3 - y2) / r23

    ax3 = G * m1 * (x1 - x3) / r13 + G * m2 * (x2 - x3) / r23
    ay3 = G * m1 * (y1 - y3) / r13 + G * m2 * (y2 - y3) / r23

    return [vx1, vy1, ax1, ay1, vx2, vy2, ax2, ay2, vx3, vy3, ax3, ay3]

# Modified simulation parameters
t_span = (0, 200)  # Reduced maximum time
t_eval = np.linspace(0, 200, 20000)  # Reduced number of points
solution = solve_ivp(equations, t_span, initial_conditions, t_eval=t_eval)

# Store positions for repetition detection
positions = np.vstack([solution.y[0:2], solution.y[4:6], solution.y[8:10]]).T

# Find where pattern starts repeating
repetition_index = len(positions)
for i in range(1000, len(positions)):
    if detect_repetition(positions[:i]):
        repetition_index = i
        break

# Trim the solution to just before repetition
x1, y1 = solution.y[0][:repetition_index], solution.y[1][:repetition_index]
x2, y2 = solution.y[4][:repetition_index], solution.y[5][:repetition_index]
x3, y3 = solution.y[8][:repetition_index], solution.y[9][:repetition_index]

# Enhanced plotting
plt.figure(figsize=(12, 12))  # Increased from 8x8
plt.plot(x1, y1, label="Body 1", linewidth=0.8)
plt.plot(x2, y2, label="Body 2", linewidth=0.8)
plt.plot(x3, y3, label="Body 3", linewidth=0.8)
plt.scatter([-1, 1, 0], [0, 0, 0.5], color='red', s=100, label='Initial positions')  # Keeping your original initial positions
plt.title("Three-Body Problem Simulation (Extended Time and Space)", pad=20)
plt.xlabel("X Position")
plt.ylabel("Y Position")
plt.xlim(-30, 30)  # Extended view
plt.ylim(-30, 30)  # Extended view
plt.grid(True, linestyle='--', alpha=0.3)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()
