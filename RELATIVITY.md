# Relativistic Physics in N-Body Simulations

This document explains the physics behind the relativistic corrections implemented in this simulation. We move beyond Newtonian gravity to a **Post-Newtonian (PN)** approximation derived directly from the **Einstein Field Equations (EFE)**.

## 1. The Einstein Field Equation (EFE)

The core of General Relativity is the Einstein Field Equation:

$$ G_{\mu\nu} = 8\pi G T_{\mu\nu} $$

Where:
- **$G_{\mu\nu}$ (Einstein Tensor)**: Describes the curvature of spacetime. It is constructed from the **Metric Tensor $g_{\mu\nu}$** and its derivatives.
- **$T_{\mu\nu}$ (Stress-Energy Tensor)**: Describes the distribution of matter and energy (mass, momentum, pressure).
- **$8\pi G$**: Coupling constant (in units where $c=1$).

In simple terms: **Mass tells spacetime how to curve, and curved spacetime tells mass how to move.**

## 2. The Metric Tensor $g_{\mu\nu}$

Instead of a simple "force," gravity is the geometry of 4D spacetime. This geometry is defined by the metric tensor $g_{\mu\nu}$, which measures distances in spacetime:

$$ ds^2 = g_{\mu\nu} dx^\mu dx^\nu $$

In the **Weak-Field Approximation** (where gravity is not extremely strong, like near a black hole horizon), we can write the metric as a small perturbation $h_{\mu\nu}$ on top of flat Minkowski space $\eta_{\mu\nu}$:

$$ g_{\mu\nu} \approx \eta_{\mu\nu} + h_{\mu\nu} $$

- $\eta_{\mu\nu} = \text{diag}(-1, 1, 1, 1)$ (Flat spacetime)
- $h_{\mu\nu} \ll 1$ (Small gravitational perturbation)

## 3. From EFE to Motion (Geodesics)

Objects in GR follow **geodesics**—the straightest possible paths in curved spacetime. The equation of motion is:

$$ \frac{d^2 x^\mu}{d\tau^2} + \Gamma^\mu_{\alpha\beta} \frac{dx^\alpha}{d\tau} \frac{dx^\beta}{d\tau} = 0 $$

Here, $\Gamma^\mu_{\alpha\beta}$ are the **Christoffel symbols**, which represent the "gravitational force" field and are calculated from the derivatives of the metric $g_{\mu\nu}$.

## 4. The Einstein-Infeld-Hoffmann (EIH) Approximation

Solving the full EFE for N moving bodies is extremely difficult (nonlinear and coupled). However, for solar-system-like scales, we can expand the equations in powers of $(v/c)^2$. This is called the **Post-Newtonian (PN)** expansion.

The **1PN (First Post-Newtonian)** approximation gives us the **Einstein-Infeld-Hoffmann (EIH)** equations of motion. This adds relativistic correction terms to the standard Newtonian acceleration.

### The EIH Equation of Motion

The acceleration $\mathbf{a}_a$ of body $a$ due to other bodies $b$ is given by:

$$ \mathbf{a}_a = \mathbf{a}_{Newtonian} + \mathbf{a}_{1PN} $$

$$ \mathbf{a}_a = \sum_{b \neq a} \frac{Gm_b}{r_{ab}^2} \mathbf{n}_{ab} + \frac{1}{c^2} \sum_{b \neq a} \frac{Gm_b}{r_{ab}^2} \left( \mathbf{A}_{ab} + \mathbf{B}_{ab} \right) $$

Where the relativistic correction terms are:

**Term A (Scalar corrections along the radial vector $\mathbf{n}_{ab}$):**
$$ \mathbf{A}_{ab} = \mathbf{n}_{ab} \left[ v_a^2 (1 + \dots) - v_b^2 (1 + \dots) - \frac{3}{2}(\mathbf{n}_{ab} \cdot \mathbf{v}_b)^2 + \dots \right] $$
*(Includes time dilation and length contraction effects)*

**Term B (Vector corrections along velocity vectors):**
$$ \mathbf{B}_{ab} = (\mathbf{n}_{ab} \cdot \mathbf{v}_a) (4 \mathbf{v}_a - 3 \mathbf{v}_b) $$
*(Includes "gravito-magnetic" frame-dragging-like effects)*

*(Note: The full EIH equation is quite long; the code implements the standard form used in celestial mechanics.)*

## 5. Physical Interpretation of Corrections

1.  **Perihelion Precession**: Orbits are no longer perfectly closed ellipses. They precess (rotate) over time. This explains the famous **Precession of Mercury**.
2.  **Time Dilation**: Clocks run slower in strong gravity and at high speeds. This affects how "fast" a body appears to accelerate.
3.  **Non-Linearity**: Gravity interacts with itself. The gravitational field has energy, and that energy creates more gravity.

## 6. Limitations

- **Weak Field**: This model breaks down near black hole event horizons or neutron star mergers.
- **Low Velocity**: Valid for $v \ll c$.
- **No Gravitational Waves**: This 1PN approximation does not include the energy loss due to gravitational radiation (which appears at 2.5PN).

## 7. Why This Code?

This code calculates the EIH acceleration at every time step. It allows us to simulate relativistic phenomena like orbital precession without the immense computational cost of full Numerical Relativity.
