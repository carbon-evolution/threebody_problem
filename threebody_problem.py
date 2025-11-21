import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import random
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button

# =============================================================================
# RELATIVISTIC PHYSICS CONSTANTS
# =============================================================================
G = 1.0
c = 50.0 
c2 = c * c 

def get_acceleration(positions, velocities, masses):
    """
    Calculates acceleration for all bodies using the Einstein-Infeld-Hoffmann (EIH)
    approximation (1st Post-Newtonian order).
    """
    n = len(masses)
    accelerations = np.zeros_like(positions)
    
    for a in range(n):
        acc_newton = np.zeros(3)
        acc_eih = np.zeros(3)
        
        pos_a = positions[a]
        vel_a = velocities[a]
        v_a_sq = np.dot(vel_a, vel_a)
        
        for b in range(n):
            if a == b:
                continue
            
            pos_b = positions[b]
            vel_b = velocities[b]
            mass_b = masses[b]
            
            r_vec = pos_b - pos_a
            r = np.linalg.norm(r_vec)
            n_vec = r_vec / r 
            
            # 1. NEWTONIAN GRAVITY
            term_newton = (G * mass_b / (r**2)) * n_vec
            acc_newton += term_newton
            
            # 2. RELATIVISTIC CORRECTIONS (1PN)
            v_b_sq = np.dot(vel_b, vel_b)
            v_a_dot_v_b = np.dot(vel_a, vel_b)
            n_dot_v_a = np.dot(n_vec, vel_a)
            n_dot_v_b = np.dot(n_vec, vel_b)
            
            term_A = (1/c2) * (G * mass_b / (r**2)) * n_vec * (
                4 * (G * mass_b / r) +
                v_a_sq +
                2 * v_b_sq -
                4 * v_a_dot_v_b -
                1.5 * (n_dot_v_b**2)
            )
            
            vec_correction = (1/c2) * (G * mass_b / (r**2)) * n_dot_v_a * (4 * vel_a - 3 * vel_b)
            
            acc_eih += term_A + vec_correction

        accelerations[a] = acc_newton + acc_eih
        
    return accelerations

def equations(t, y, masses):
    n = len(masses)
    pos = np.zeros((n, 3))
    vel = np.zeros((n, 3))
    
    for i in range(n):
        idx = i * 6
        pos[i] = y[idx : idx+3]
        vel[i] = y[idx+3 : idx+6]
    
    acc = get_acceleration(pos, vel, masses)
    
    dydt = np.zeros(6 * n)
    for i in range(n):
        idx = i * 6
        dydt[idx : idx+3] = vel[i]
        dydt[idx+3 : idx+6] = acc[i]
    
    return dydt

class SimulationAnimation:
    def __init__(self, fig, ax, solution, masses, title):
        self.fig = fig
        self.ax = ax
        self.solution = solution
        self.masses = masses
        self.n = len(masses)
        self.title = title
        self.is_paused = False
        
        # Data setup
        self.lines = []
        self.points = []
        self.history_len = 200 # How much trail to show
        
        # Initialize plots
        for i in range(self.n):
            # Line for trajectory trail
            line, = ax.plot([], [], [], linewidth=1.0, label=f"Body {i+1} (m={masses[i]:.2f})")
            self.lines.append(line)
            # Point for current position
            point, = ax.plot([], [], [], marker='o', markersize=6, color=line.get_color())
            self.points.append(point)
            
            # Start marker
            idx = i * 6
            ax.scatter(solution.y[idx][0], solution.y[idx+1][0], solution.y[idx+2][0], 
                      color=line.get_color(), s=50, marker='o', alpha=0.5)

        # Set limits
        self.set_limits()
        
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title(title)
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.3)

    def set_limits(self):
        # Start zoomed in on initial positions
        all_x = []
        all_y = []
        all_z = []
        for i in range(self.n):
            idx = i * 6
            # Use only first few points for initial zoom
            all_x.append(self.solution.y[idx][0])
            all_y.append(self.solution.y[idx+1][0])
            all_z.append(self.solution.y[idx+2][0])
            
        # Add some padding around initial positions
        margin = 1.5
        max_range = max(max(all_x)-min(all_x), max(all_y)-min(all_y), max(all_z)-min(all_z)) + margin
        if max_range < 2.0:
            max_range = 2.0
            
        mid_x = (max(all_x)+min(all_x)) * 0.5
        mid_y = (max(all_y)+min(all_y)) * 0.5
        mid_z = (max(all_z)+min(all_z)) * 0.5
        
        self.ax.set_xlim(mid_x - max_range, mid_x + max_range)
        self.ax.set_ylim(mid_y - max_range, mid_y + max_range)
        self.ax.set_zlim(mid_z - max_range, mid_z + max_range)

    def update(self, frame):
        if self.is_paused:
            return self.lines + self.points
            
        # Speed up: skip frames
        idx_frame = frame * 5 
        if idx_frame == 0:
            idx_frame = 1
        if idx_frame >= len(self.solution.t):
            idx_frame = len(self.solution.t) - 1
            
        for i in range(self.n):
            idx = i * 6
            x = self.solution.y[idx][:idx_frame]
            y = self.solution.y[idx+1][:idx_frame]
            z = self.solution.y[idx+2][:idx_frame]
            
            # Update trail (show full history)
            self.lines[i].set_data(x, y)
            self.lines[i].set_3d_properties(z)
            
            # Update current point
            self.points[i].set_data([x[-1]], [y[-1]])
            self.points[i].set_3d_properties([z[-1]])
        
        # Dynamic axis expansion: check if any body is near boundaries
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        zlim = self.ax.get_zlim()
        
        # Get current positions
        needs_expansion = False
        for i in range(self.n):
            idx = i * 6
            if idx_frame < len(self.solution.t):
                curr_x = self.solution.y[idx][idx_frame]
                curr_y = self.solution.y[idx+1][idx_frame]
                curr_z = self.solution.y[idx+2][idx_frame]
                
                # Check if close to boundaries (within 10% margin)
                x_margin = (xlim[1] - xlim[0]) * 0.1
                y_margin = (ylim[1] - ylim[0]) * 0.1
                z_margin = (zlim[1] - zlim[0]) * 0.1
                
                if (curr_x < xlim[0] + x_margin or curr_x > xlim[1] - x_margin or
                    curr_y < ylim[0] + y_margin or curr_y > ylim[1] - y_margin or
                    curr_z < zlim[0] + z_margin or curr_z > zlim[1] - z_margin):
                    needs_expansion = True
                    break
        
        if needs_expansion:
            # Expand limits by 20%
            x_center = (xlim[0] + xlim[1]) / 2
            y_center = (ylim[0] + ylim[1]) / 2
            z_center = (zlim[0] + zlim[1]) / 2
            
            x_range = (xlim[1] - xlim[0]) * 1.2
            y_range = (ylim[1] - ylim[0]) * 1.2
            z_range = (zlim[1] - zlim[0]) * 1.2
            
            self.ax.set_xlim(x_center - x_range/2, x_center + x_range/2)
            self.ax.set_ylim(y_center - y_range/2, y_center + y_range/2)
            self.ax.set_zlim(z_center - z_range/2, z_center + z_range/2)
            
        return self.lines + self.points

    def toggle_pause(self, event=None):
        self.is_paused = not self.is_paused

def run_simulation(masses, initial_velocities, initial_pos, title, ax, fig):
    n = len(masses)
    initial_conditions = []
    for i in range(n):
        initial_conditions.extend(initial_pos[i])
        initial_conditions.extend(initial_velocities[i])
    
    print(f"\nSystem parameters for {title}:")
    for i in range(n):
        print(f"Body {i+1}: m={masses[i]:.2f}")
    
    t_span = (0, 375)  # Extended for 1-minute animation
    t_eval = np.linspace(0, 375, 15000)
    
    solution = solve_ivp(
        lambda t, y: equations(t, y, masses),
        t_span,
        initial_conditions,
        t_eval=t_eval,
        method='RK45',
        rtol=1e-6,
        atol=1e-9
    )
    
    anim = SimulationAnimation(fig, ax, solution, masses, title)
    return anim

def on_scroll(event):
    """
    Zoom in/out on scroll event for 3D axes.
    """
    if event.inaxes is None:
        return

    ax = event.inaxes
    if not isinstance(ax, Axes3D):
        return

    # Get current limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    zlim = ax.get_zlim()

    # Determine scale factor
    base_scale = 1.1
    if event.button == 'up':
        # Zoom in
        scale_factor = 1 / base_scale
    elif event.button == 'down':
        # Zoom out
        scale_factor = base_scale
    else:
        return

    # Calculate new limits
    x_center = (xlim[0] + xlim[1]) / 2
    y_center = (ylim[0] + ylim[1]) / 2
    z_center = (zlim[0] + zlim[1]) / 2

    x_range = (xlim[1] - xlim[0]) * scale_factor
    y_range = (ylim[1] - ylim[0]) * scale_factor
    z_range = (zlim[1] - zlim[0]) * scale_factor

    ax.set_xlim([x_center - x_range / 2, x_center + x_range / 2])
    ax.set_ylim([y_center - y_range / 2, y_center + y_range / 2])
    ax.set_zlim([z_center - z_range / 2, z_center + z_range / 2])

    event.canvas.draw_idle()

def on_key(event):
    if event.key == ' ':
        anim1.toggle_pause()
        anim2.toggle_pause()

# Create figure
fig = plt.figure(figsize=(16, 8))
fig.canvas.mpl_connect('key_press_event', on_key)
fig.canvas.mpl_connect('scroll_event', on_scroll)

# =============================================================================
# SCENARIO 1: Solar System-like
# =============================================================================
ax1 = fig.add_subplot(121, projection='3d')
m_sun = 10.0
m_p1 = 0.05
m_p2 = 0.1
m_p3 = 0.2
masses_1 = [m_sun, m_p1, m_p2, m_p3]
pos_sun = [0.0, 0.0, 0.5]  # Sun at Z=0.5 for 3D motion
pos_p1 = [1.2, 0.0, 0.0]
pos_p2 = [2.0, 0.0, 0.5]
pos_p3 = [3.5, 0.0, -0.3]
v_p1 = np.sqrt(G * m_sun / 1.2)
v_p2 = np.sqrt(G * m_sun / 2.0)
v_p3 = np.sqrt(G * m_sun / 3.5)
vel_sun = [0.0, 0.0, 0.1]  # Sun has slight Z-motion
vel_p1 = [0.0, v_p1, 0.3]
vel_p2 = [0.0, v_p2, -0.5]
vel_p3 = [0.0, v_p3, 0.2]

anim1 = run_simulation(masses_1, [vel_sun, vel_p1, vel_p2, vel_p3], 
                      [pos_sun, pos_p1, pos_p2, pos_p3],
                      "Solar System (Sun + 3 Planets)", ax1, fig)

# =============================================================================
# SCENARIO 2: Chaotic 3-Body
# =============================================================================
ax2 = fig.add_subplot(122, projection='3d')
m1 = 1.0
m2 = 1.2
m3 = 0.9
masses_2 = [m1, m2, m3]
pos1 = [1.0, 0.0, 0.0]
pos2 = [-0.5, 0.866, 0.0]
pos3 = [-0.5, -0.866, 0.0]

# Velocities (Rotational setup to keep them bound initially)
# Tangential velocities with significant Z-components for 3D motion
vel1 = [0.0, 0.8, 0.4]
vel2 = [-0.7, -0.4, -0.3]
vel3 = [0.7, -0.4, 0.5]

anim2 = run_simulation(masses_2, [vel1, vel2, vel3], 
                      [pos1, pos2, pos3],
                      "Chaotic 3-Body (No Sun)", ax2, fig)

# Add Play/Pause Button
ax_play = plt.axes([0.45, 0.05, 0.1, 0.075])
btn_play = Button(ax_play, 'Play/Pause')
def on_click(event):
    anim1.toggle_pause()
    anim2.toggle_pause()
btn_play.on_clicked(on_click)

# Create animations (1 minute = 3000 frames at 20ms/frame)
ani1 = FuncAnimation(fig, anim1.update, frames=3000, interval=20, blit=False, repeat=True)
ani2 = FuncAnimation(fig, anim2.update, frames=3000, interval=20, blit=False, repeat=True)

plt.tight_layout()
plt.show()
