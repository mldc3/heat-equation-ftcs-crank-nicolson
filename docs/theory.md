# Theory Background

This document develops the physical and numerical theory behind the one-dimensional heat-equation simulations in this repository. The goal is to explain the mathematical model, the physical interpretation, and the numerical methods in a way that is clearer and more complete than a short classroom derivation.

The practice studies a parabolic partial differential equation: the heat equation. It compares explicit and implicit finite-difference time integration, constant and spatially variable diffusion, stability restrictions, conservative discretization, numerical error, and computational cost.

---

## 1. Physical meaning of diffusion

Diffusion is the process by which a quantity spreads from regions of high concentration or high temperature toward regions of lower concentration or temperature. In this project the evolving scalar field is the temperature $T(x,t)$, but the same mathematical structure appears in heat conduction, particle diffusion, chemical concentration, probability spreading and smoothing processes.

The essential physical feature of diffusion is that gradients are reduced over time. If the initial temperature has a sharp peak, the peak gradually becomes lower and wider. Heat does not move as a wave with a fixed propagation speed. Instead, the temperature field relaxes by smoothing spatial variations.

This behaviour is very different from hyperbolic equations such as the wave equation or the advection equation. Hyperbolic equations transport information with finite-speed characteristics. Parabolic equations such as the heat equation are dissipative: they remove sharp structures and damp high-frequency spatial components.

---

## 2. The one-dimensional heat equation

For constant diffusion coefficient, the one-dimensional heat equation is:

$$
\frac{\partial T}{\partial t} =
D
\frac{\partial^2 T}{\partial x^2}.
$$

Here $T(x,t)$ is the temperature field, $x$ is the spatial coordinate, $t$ is time and $D$ is the diffusion coefficient.

The parameter $D$ determines how quickly heat spreads. If $D$ is large, temperature gradients are smoothed quickly. If $D$ is small, the profile evolves more slowly.

The right-hand side contains the second spatial derivative. This is important physically. The second derivative measures curvature. A point where the temperature profile is curved downward, such as the top of a peak, tends to decrease. A point where the profile is curved upward tends to increase. This is how the equation redistributes temperature and smooths the profile.

The heat equation is parabolic. This classification matters because parabolic equations are strongly connected to stability restrictions in explicit numerical schemes. In particular, the stable timestep for a simple explicit method scales as $\Delta x^2$, not as $\Delta x$. This makes explicit heat-equation simulations increasingly expensive when the grid is refined.

---

## 3. Conservation law and heat flux

The heat equation can be understood from two physical principles: local conservation and Fourier's law.

A local conservation law states that the amount of heat inside a small interval can only change because heat flows through the boundaries of that interval. In one dimension this can be written as:

$$
\frac{\partial T}{\partial t} +
\frac{\partial J}{\partial x} = 0,
$$

where $J(x,t)$ is the heat flux.

The flux $J$ measures the amount of heat crossing a point per unit time. Fourier's law states that heat flows from hot regions to cold regions:

$$
J=-D\frac{\partial T}{\partial x}.
$$

The minus sign is physically important. If temperature increases with $x$, then $\partial T/\partial x>0$, so the heat flux is negative, meaning heat flows toward smaller $x$. In other words, heat flows opposite to the temperature gradient.

Substituting Fourier's law into the conservation equation gives:

$$
\frac{\partial T}{\partial t} =
\frac{\partial}{\partial x}
\left(
D\frac{\partial T}{\partial x}
\right).
$$

If $D$ is constant, then $D$ can be taken outside the derivative:

$$
\frac{\partial T}{\partial t} =
D
\frac{\partial^2 T}{\partial x^2}.
$$

This derivation is crucial because it shows that the physically fundamental form of diffusion is a flux-divergence equation. That becomes especially important when the diffusion coefficient varies in space.

---

## 4. Constant diffusion and variable diffusion

When $D$ is constant, the heat equation is:

$$
\frac{\partial T}{\partial t} =
D
\frac{\partial^2 T}{\partial x^2}.
$$

This is the simplest model. Every point in the material conducts heat equally. The same temperature gradient produces the same flux everywhere.

However, real materials may not be homogeneous. The thermal diffusivity may depend on position. In that case, the diffusion coefficient is written as $D(x)$, and the correct equation is:

$$
\frac{\partial T}{\partial t} =
\frac{\partial}{\partial x}
\left(
D(x)
\frac{\partial T}{\partial x}
\right).
$$

This is not the same as simply writing:

$$
\frac{\partial T}{\partial t} =
D(x)
\frac{\partial^2 T}{\partial x^2}.
$$

To see the difference, expand the derivative:

$$
\begin{aligned}
\frac{\partial}{\partial x}
\left(
D(x)
\frac{\partial T}{\partial x}
\right)
&=
D(x)
\frac{\partial^2 T}{\partial x^2} +
\frac{dD}{dx}
\frac{\partial T}{\partial x}.
\end{aligned}
$$

The extra term involving $dD/dx$ is missing if one simply multiplies the Laplacian by $D(x)$. This is why the practice distinguishes between a first non-conservative variable-diffusion attempt and the corrected conservative discretization.

The conservative form is physically preferable because it describes fluxes through interfaces between neighbouring cells. When diffusivity changes in space, heat flux depends not only on the local temperature curvature but also on how the material property changes from one region to another.

---

## 5. Initial condition

The simulations use a localized Gaussian initial condition:

$$
T(x,0)=100e^{-20x^2}.
$$

This profile is useful because it is smooth, has a clear maximum at the centre, decays toward the boundaries and makes the diffusion process visually obvious.

At $t=0$, most of the heat is concentrated near $x=0$. As time evolves, the central peak decreases and the temperature spreads outward. This is the expected physical behaviour of the heat equation.

The Gaussian is also numerically convenient. Unlike a discontinuous initial condition, it does not introduce artificial jumps at the beginning of the simulation. This makes it easier to focus on diffusion, stability and the behaviour of the time integration schemes.

---

## 6. Boundary conditions

The simulations use homogeneous Dirichlet boundary conditions:

$$
T(-1,t)=0,
\qquad
T(1,t)=0.
$$

Physically, these boundaries can be interpreted as fixed-temperature reservoirs held at zero temperature. Heat reaching the boundaries can be absorbed by the reservoirs.

Numerically, boundary conditions must be imposed carefully. If the boundary values are not fixed after every update, they can drift due to matrix operations or roundoff errors. This would change the physical problem being solved.

Dirichlet boundaries also affect the long-time behaviour. Because the boundaries are held at zero, the temperature profile eventually decays toward zero throughout the domain. Heat is not conserved inside the finite domain; it can leave through the boundaries.

This is different from insulated Neumann boundaries, where the derivative would be set to zero and heat would remain in the system. Therefore, boundary conditions are not just numerical details. They define the physical model.

---

## 7. Spatial discretization

The continuous spatial interval is replaced by a finite grid:

$$
x_j=x_{\min}+j\Delta x.
$$

The temperature at each grid point and time level is written as:

$$
T_j^n \approx T(x_j,t^n).
$$

The second derivative is approximated using the centered finite-difference stencil:

$$
\frac{\partial^2 T}{\partial x^2}
\bigg|_{x_j}
\approx \frac{
T_{j+1}^n -
2T_j^n
+
T_{j-1}^n
}{\Delta x^2}.
$$

This formula is second-order accurate in space. It uses only the nearest neighbours of each point. That locality is the reason the discrete Laplacian matrix is tridiagonal.

In matrix form, the second derivative operator can be written as:

$$
L=
\frac{1}{\Delta x^2}
\begin{pmatrix}
-2 & 1 & 0 & \cdots & 0 \\
1 & -2 & 1 & \cdots & 0 \\
0 & 1 & -2 & \cdots & 0 \\
\cdots & \cdots & \cdots & \cdots & \cdots \\
0 & 0 & 0 & 1 & -2
\end{pmatrix}.
$$

The matrix is sparse because almost all entries are zero. Sparse matrices are useful for computational physics because they allow large systems to be stored and manipulated efficiently.

---

## 8. FTCS explicit method

FTCS stands for Forward Time, Centered Space. The time derivative is approximated with a forward difference:

$$
\frac{\partial T}{\partial t}
\approx
\frac{T_j^{n+1}-T_j^n}{\Delta t}.
$$

The spatial derivative is approximated with the centered second-derivative stencil:

$$
\frac{\partial^2 T}{\partial x^2}
\approx
\frac{
T_{j+1}^n -
2T_j^n
+
T_{j-1}^n
}{\Delta x^2}.
$$

Substituting into the heat equation gives:

$$
\frac{T_j^{n+1}-T_j^n}{\Delta t} =
D
\frac{
T_{j+1}^n -
2T_j^n
+
T_{j-1}^n
}{\Delta x^2}.
$$

Solving for the new temperature gives:

$$
\begin{aligned}
T_j^{n+1} &=
T_j^n +
r
\left(
T_{j+1}^n -
2T_j^n
+
T_{j-1}^n
\right),
\\
r
&=
\frac{D\Delta t}{\Delta x^2}.
\end{aligned}
$$

The method is explicit because all quantities on the right-hand side are known at time level $n$. This makes each timestep cheap. No linear system must be solved.

The disadvantage is that explicit heat-equation schemes are conditionally stable. If $\Delta t$ is too large, the numerical solution develops nonphysical oscillations and can blow up.

---

## 9. FTCS stability condition

For the one-dimensional heat equation, FTCS is stable only if:

$$
r\leq\frac{1}{2}.
$$

Since:

$$
r=\frac{D\Delta t}{\Delta x^2},
$$

the timestep must satisfy:

$$
\Delta t
\leq
\frac{\Delta x^2}{2D}.
$$

This condition has an important computational consequence. If the spatial grid is refined by making $\Delta x$ smaller, the timestep must decrease as $\Delta x^2$. Therefore, doubling the spatial resolution can require roughly four times more timesteps to reach the same final time.

For variable diffusion, the most restrictive part of the domain is where $D(x)$ is largest. A safe condition is:

$$
\Delta t
\leq
\frac{\Delta x^2}{2\max(D(x))}.
$$

This explains why stronger diffusion makes explicit simulations more expensive: increasing $D$ forces a smaller stable timestep.

---

## 10. Von Neumann stability analysis

Von Neumann analysis studies the evolution of Fourier modes under a numerical scheme. A discrete perturbation can be written as:

$$
T_j^n = Q^n e^{ikj\Delta x},
$$

where $Q$ is the amplification factor.

A method is stable if all modes satisfy:

$$
|Q|\leq1.
$$

For FTCS applied to the heat equation, substitution gives:

$$
Q = 1 -
4r
\sin^2
\left(
\frac{k\Delta x}{2}
\right).
$$

The strongest restriction comes from the highest spatial frequencies, for which the sine term can reach 1. Requiring $|Q|\leq1$ gives:

$$
r\leq\frac{1}{2}.
$$

This analysis explains why numerical instability often appears as point-to-point oscillations. High-frequency modes are the first to become unstable.

For Crank-Nicolson, the amplification factor has the form:

$$
Q = \frac{ 1 -
2r
\sin^2
\left(
\frac{k\Delta x}{2}
\right)
}{
1
+
2r
\sin^2
\left(
\frac{k\Delta x}{2}
\right)
}.
$$

For positive $r$, the denominator is always larger in magnitude than the numerator, so $|Q|\leq1$. This is the origin of the unconditional stability of Crank-Nicolson for the linear heat equation.

---

## 11. Crank-Nicolson method

Crank-Nicolson is obtained by averaging the diffusion operator between the old and new time levels. Instead of evaluating the spatial derivative only at time $n$, it evaluates the average of the derivative at times $n$ and $n+1$.

For constant diffusion:

$$
\frac{T_j^{n+1}-T_j^n}{\Delta t} =
\frac{D}{2}
\left[
\left(
\frac{\partial^2 T}{\partial x^2}
\right)_j^n
+
\left(
\frac{\partial^2 T}{\partial x^2}
\right)_j^{n+1}
\right].
$$

In matrix form:

$$
\left(
I-\frac{D\Delta t}{2}L
\right)
\mathbf{T}^{n+1} =
\left(
I+\frac{D\Delta t}{2}L
\right)
\mathbf{T}^{n}.
$$

This method is implicit because $\mathbf{T}^{n+1}$ appears inside a matrix equation. Each timestep requires solving a linear system.

The advantages are second-order accuracy in time, unconditional stability for the linear heat equation, better behaviour for larger timesteps and suitability for longer simulations. The disadvantage is computational cost per timestep because solving a linear system is more expensive than applying an explicit update.

---

## 12. Explicit versus implicit trade-off

FTCS and Crank-Nicolson represent two different numerical philosophies.

FTCS is simple, direct and cheap per timestep. However, it is limited by a strict stability condition. It may require many small timesteps.

Crank-Nicolson is more complex and expensive per timestep. However, it can use larger timesteps without becoming unstable. This makes it powerful for stiff diffusion problems or long-time simulations.

The best method depends on the goal. If the grid is small and the final time is short, FTCS may be efficient. If the grid is fine, the diffusion coefficient is large, or the final time is long, Crank-Nicolson can be more practical despite the cost of solving linear systems.

This is why the practice includes both runtime and error comparisons.

---

## 13. Conservative discretization for variable diffusion

For spatially variable diffusion, the physically correct equation is:

$$
\frac{\partial T}{\partial t} =
\frac{\partial}{\partial x}
\left(
D(x)
\frac{\partial T}{\partial x}
\right).
$$

The conservative finite-difference approach discretizes heat fluxes at cell interfaces. Define:

$$
D_{j+1/2} =
\frac{D_j+D_{j+1}}{2}.
$$

The flux through the interface between $j$ and $j+1$ is approximated by:

$$
J_{j+1/2} = -
D_{j+1/2}
\frac{T_{j+1}-T_j}{\Delta x}.
$$

The divergence of the flux is then:

$$
-\frac{J_{j+1/2}-J_{j-1/2}}{\Delta x}.
$$

Substituting the flux approximation gives the conservative operator:

$$
\begin{aligned}
\left(
\nabla\cdot D\nabla T
\right)_j
&\approx
\frac{1}{\Delta x^2}
\left[
D_{j+1/2}
\left(
T_{j+1}-T_j
\right) -
D_{j-1/2}
\left(
T_j-T_{j-1}
\right)
\right].
\end{aligned}
$$

This form is better than simply using $D_jL$ because it respects the flux balance between neighbouring cells. It is also more physically meaningful at material interfaces or in regions where diffusivity changes rapidly.

---

## 14. Error analysis

Numerical methods are not judged only by whether they remain stable. A stable method can still be inaccurate.

The project compares numerical solutions obtained with different timesteps. A highly resolved solution, usually obtained with a small timestep, can be used as a reference. The error of coarser simulations is then measured relative to that reference.

The maximum error is:

$$
E_{\max}
=
\max_j
\left|
T_j^{\mathrm{num}} -
T_j^{\mathrm{ref}}
\right|.
$$

The mean error is:

$$
E_{\mathrm{mean}} =
\frac{1}{N}
\sum_j
\left|
T_j^{\mathrm{num}} -
T_j^{\mathrm{ref}}
\right|.
$$

The maximum error identifies the worst pointwise discrepancy. The mean error measures the global average discrepancy. Both are useful because a method may have a localized large error while remaining accurate on average, or it may have moderate error spread across the whole domain.

---

## 15. Runtime analysis

Runtime analysis connects numerical methods to practical computation.

FTCS has low cost per timestep:

$$
\mathbf{T}^{n+1} =
\left(
I+D\Delta t L
\right)
\mathbf{T}^{n}.
$$

This is essentially a matrix-vector update or a local stencil update.

Crank-Nicolson requires solving:

$$
M_{\mathrm{left}}
\mathbf{T}^{n+1} =
\mathbf{b}.
$$

Solving this system is more expensive per timestep. However, Crank-Nicolson may need fewer timesteps because it is not limited by the explicit FTCS stability condition.

Therefore, the runtime comparison is not trivial. The faster method depends on grid size, final time, timestep choice, implementation, and whether sparse matrix methods are used efficiently.

---

## 16. Diffusion-amplitude sweeps

The practice studies variable diffusion of the form:

$$
D(x)=\alpha
\left(
1+e^{-x^2}
\right).
$$

The parameter $\alpha$ controls the overall strength of diffusion. Increasing $\alpha$ increases $D(x)$ everywhere, especially near the centre where $e^{-x^2}$ is largest.

Physically, larger $\alpha$ means faster smoothing. The Gaussian peak should flatten more rapidly, and the temperature should spread more strongly toward the boundaries.

Numerically, increasing $\alpha$ also tightens the FTCS stability condition:

$$
\Delta t
\leq
\frac{\Delta x^2}{2\max(D(x))}.
$$

Thus the alpha sweep connects physics and computation: stronger diffusion produces faster physical relaxation but requires smaller explicit timesteps.

---

## 17. Main theoretical conclusions

The key conclusions are:

1. The heat equation is parabolic and dissipative.
2. Diffusion smooths gradients and damps sharp spatial structures.
3. FTCS is simple and explicit but conditionally stable.
4. The FTCS stability limit scales as $\Delta x^2$.
5. Crank-Nicolson is implicit, more expensive per timestep, but unconditionally stable for the linear heat equation.
6. Stability does not automatically imply accuracy.
7. Variable diffusion must be discretized in conservative flux form for physical consistency.
8. Sparse matrices are natural because finite-difference operators are local.
9. Runtime comparisons are essential because different schemes have different per-step cost and timestep restrictions.
10. The numerical method must be chosen according to both the physics of the PDE and the computational cost.
