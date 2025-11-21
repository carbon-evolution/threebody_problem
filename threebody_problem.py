import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import random
from mpl_toolkits.mplot3d import Axes3D

# =============================================================================
# RELATIVISTIC PHYSICS CONSTANTS
# =============================================================================
# G: Gravitational Constant (set to 1 for normalized units)
G = 1.0

# c: Speed of Light
# In this simulation, 'c' controls the strength of relativistic effects.
# - High 'c' (e.g., 1000.0) -> Newtonian Gravity (v/c is small)
# - Low 'c' (e.g., 10.0)   -> Strong Relativistic Effects (v/c is significant)
# We choose a moderately low value to make effects like precession visible.
c = 50.0 
c2 = c * c # c squared

def get_acceleration(positions, velocities, masses):
    """
    Calculates acceleration for all bodies using the Einstein-Infeld-Hoffmann (EIH)
    approximation (1st Post-Newtonian order).
    
    This approximation is derived from the Einstein Field Equations:
    G_mu_nu = 8*pi*G * T_mu_nu
    
    It adds relativistic correction terms to the standard Newtonian gravity.
    """
    n = len(masses)
    accelerations = np.zeros_like(positions)
    
    for a in range(n):
        # Newtonian term accumulator
        acc_newton = np.zeros(3)
        # Relativistic correction term accumulator
        acc_eih = np.zeros(3)
        
        pos_a = positions[a]
        vel_a = velocities[a]
        v_a_sq = np.dot(vel_a, vel_a) # v_a^2
        
        for b in range(n):
            if a == b:
                continue
            
            pos_b = positions[b]
            vel_b = velocities[b]
            mass_b = masses[b]
            
            # Relative position vector r_ab = x_a - x_b
            # Note: Standard formula often uses x_b - x_a for force direction.
            # Here we want force ON a FROM b.
            # Vector pointing from a to b: r_vec = pos_b - pos_a
            r_vec = pos_b - pos_a
            r = np.linalg.norm(r_vec)
            n_vec = r_vec / r # Unit vector pointing to b
            
            # -----------------------------------------------------------------
            # 1. NEWTONIAN GRAVITY (The "0th Order" Approximation)
            # -----------------------------------------------------------------
            # F = G * m_a * m_b / r^2
            # a = F / m_a = G * m_b / r^2
            term_newton = (G * mass_b / (r**2)) * n_vec
            acc_newton += term_newton
            
            # -----------------------------------------------------------------
            # 2. RELATIVISTIC CORRECTIONS (The "1st Post-Newtonian" Order)
            # -----------------------------------------------------------------
            # These terms come from expanding the metric tensor g_mu_nu.
            # They account for:
            # - Time dilation (clocks running slower near masses)
            # - Spatial curvature (space being "stretched")
            # - Nonlinearity (gravity creating gravity)
            
            # Dot products needed for EIH
            v_b_sq = np.dot(vel_b, vel_b)       # v_b^2
            v_a_dot_v_b = np.dot(vel_a, vel_b)  # v_a . v_b
            n_dot_v_a = np.dot(n_vec, vel_a)    # n . v_a
            n_dot_v_b = np.dot(n_vec, vel_b)    # n . v_b
            
            # EIH Equation Terms (simplified standard form):
            # A_scalar * n_vec + B_scalar * (v_a - v_b)
            
            # Term A: Scalar corrections along the line of sight (n_vec)
            # - 4 * G * m_b / r : Gravitational potential contribution
            # - v_a^2 : Kinetic energy contribution (mass-energy equivalence)
            # - 1.5 * (n . v_b)^2 : Retardation effect
            term_A = (1/c2) * (G * mass_b / (r**2)) * n_vec * (
                4 * (G * mass_b / r) +  # Self-gravity potential
                v_a_sq +                # Time dilation / Kinematics
                2 * v_b_sq -            # Source motion
                4 * v_a_dot_v_b -       # Frame dragging coupling
                1.5 * (n_dot_v_b**2)    # Retardation
            )
            
            # Term B: Vector corrections along the velocity difference
            # This represents "gravito-magnetic" effects (like magnetic force but for gravity)
            term_B = (1/c2) * (G * mass_b / (r**2)) * (pos_b - pos_a) * (
                # Actually, the vector part is usually along (v_a - v_b) or similar.
                # Let's use the standard EIH vector form:
                # 4 * (n . v_a) * v_a  - 3 * (n . v_a) * v_b ...
                # A common simplified form for simulation is:
                # (n . v_a) * (4*v_a - 3*v_b)
                0 # Placeholder, let's implement the specific vector term below
            )
            
            # Correct Vector Term Implementation:
            # The force has a component along the velocity vector difference.
            # F_vector = (G m_a m_b / r^2) * (1/c^2) * (n . v_a) * (4 v_a - 3 v_b)
            vec_correction = (1/c2) * (G * mass_b / (r**2)) * n_dot_v_a * (4 * vel_a - 3 * vel_b)
            
            acc_eih += term_A + vec_correction

        accelerations[a] = acc_newton + acc_eih
        
    return accelerations

def equations(t, y, m1, m2, m3):
    """
    Differential equations for the system.
    y contains [x1,y1,z1, vx1,vy1,vz1, x2... , x3...]
    """
    # Reshape y into positions and velocities
    # y is a 1D array of size 18 (3 bodies * 6 coordinates)
    
    # Extract positions (indices 0-2, 6-8, 12-14)
    pos = np.array([y[0:3], y[6:9], y[12:15]])
    # Extract velocities (indices 3-5, 9-11, 15-17)
    vel = np.array([y[3:6], y[9:12], y[15:18]])
    
    masses = [m1, m2, m3]
    
    # Calculate acceleration using Relativistic EIH approximation
    acc = get_acceleration(pos, vel, masses)
    
    # Flatten output
    # [vx1, vy1, vz1, ax1, ay1, az1, ...]
    dydt = np.zeros(18)
    dydt[0:3] = vel[0]
    dydt[3:6] = acc[0]
    dydt[6:9] = vel[1]
    dydt[9:12] = acc[1]
    dydt[12:15] = vel[2]
    dydt[15:18] = acc[2]
    
    return dydt

def simulate_system(m1, m2, m3, initial_velocities, title, ax):
    # Initial positions (same as before)
    initial_pos = [
        [-1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.5, 0.0]
    ]
    
    vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3 = initial_velocities
    
    initial_conditions = [
        initial_pos[0][0], initial_pos[0][1], initial_pos[0][2], vx1, vy1, vz1,
        initial_pos[1][0], initial_pos[1][1], initial_pos[1][2], vx2, vy2, vz2,
        initial_pos[2][0], initial_pos[2][1], initial_pos[2][2], vx3, vy3, vz3
    ]

    print(f"\nSystem parameters for {title}:")
    print(f"Masses: m1={m1:.2f}, m2={m2:.2f}, m3={m3:.2f}")
    print(f"Relativistic Limit: c={c:.1f} (Lower c = Stronger GR effects)")
    
    # Simulation parameters
    t_span = (0, 50)
    t_eval = np.linspace(0, 50, 5000)
    
    solution = solve_ivp(
        lambda t, y: equations(t, y, m1, m2, m3),
        t_span,
        initial_conditions,
        t_eval=t_eval,
        method='RK45',
        rtol=1e-6, # Higher precision needed for GR effects
        atol=1e-9
    )

    # Extract positions for plotting
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
    ax.set_xlim(-2, 2) # Zoomed in a bit more by default
    ax.set_ylim(-2, 2)
    ax.set_zlim(-2, 2)
    ax.view_init(elev=30, azim=45)
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.legend()

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

# Create figure
fig = plt.figure(figsize=(16, 8))
fig.canvas.mpl_connect('scroll_event', on_scroll)

# Scenario 1: Equal Masses (Chaotic + Relativistic)
ax1 = fig.add_subplot(121, projection='3d')
equal_velocities = [0.0, 0.8, 0.0,    
                   0.0, -0.8, 0.0,    
                   0.8, 0.0, 0.0]     
simulate_system(1.0, 1.0, 1.0, equal_velocities, 
               "Relativistic 3-Body (Equal Masses)", ax1)

# Scenario 2: Hierarchical (Star + 2 Planets) - Good for seeing precession
ax2 = fig.add_subplot(122, projection='3d')
# Heavy center, lighter orbiters
m1_h = 10.0 # "Star"
m2_h = 0.1  # "Planet 1"
m3_h = 0.1  # "Planet 2"
# Velocities for stable-ish orbits
hierarchical_velocities = [
    0.0, 0.0, 0.0,      # Star stationary-ish
    0.0, 2.5, 0.0,      # Planet 1 fast orbit
    0.0, -3.0, 0.5      # Planet 2 fast orbit inclined
]
simulate_system(m1_h, m2_h, m3_h, hierarchical_velocities,
               "Relativistic Hierarchical System", ax2)

plt.tight_layout()
plt.show()
