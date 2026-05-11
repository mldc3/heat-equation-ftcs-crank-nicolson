# Theory Background

This document explains the physical and numerical theory behind the one-dimensional heat-equation simulations in this repository.

---

## 1. Physical meaning of the heat equation

The heat equation describes the time evolution of a scalar field that spreads through space. In this project the scalar field is the temperature $T(x,t)$, but the same equation also appears in diffusion of particles, concentration fields, probability densities and relaxation processes.

In one spatial dimension, the constant-diffusion heat equation is:

$$
\frac{\partial T}{\partial t}
=
D
\frac{\partial^2 T}{\partial x^2}.
$$

Here $D$ is the diffusion coefficient. A large value of $D$ means that temperature gradients are smoothed quickly. A small value of $D$ means that the initial profile relaxes more slowly.

The equation is parabolic. This means that it is time-dependent and dissipative. Unlike a wave equation, it does not transport a disturbance at a fixed velocity. Instead, it smooths gradients. Sharp peaks become wider and lower as time evolves.

---

## 2. Conservation law and Fourier law

The heat equation can be derived from a local conservation law. If heat is neither created nor destroyed inside the domain, the change of temperature in a small region must be explained by the net heat flux entering or leaving that region.

The local conservation form is:

$$
\frac{\partial T}{\partial t}
+
\frac{\partial J}{\partial x}
=
0.
$$

The quantity $J$ is the heat flux. Fourier's law states that heat flows from hot regions to cold regions:

$$
J=-D\frac{\partial T}{\partial x}.
$$

Substituting Fourier's law into the conservation equation gives:

$$
\frac{\partial T}{\partial t}
=
\frac{\partial}{\partial x}
\left(
D\frac{\partial T}{\partial x}
\right).
$$

If $D$ is constant, this reduces to:

$$
\frac{\partial T}{\partial t}
=
D
\frac{\partial^2 T}{\partial x^2}.
$$

This derivation is important because it shows that the physically correct equation is naturally written in conservative flux form. This becomes essential when $D$ depends on position.

---

## 3. Constant diffusion versus variable diffusion

For constant $D$, the diffusion coefficient can be taken outside the spatial derivative:

$$
\frac{\partial}{\partial x}
\left(
D\frac{\partial T}{\partial x}
\right)
=
D
\frac{\partial^2 T}{\partial x^2}.
$$

For spatially variable diffusion $D(x)$, this simplification is not valid. The correct equation is:

$$
\frac{\partial T}{\partial t}
=
\frac{\partial}{\partial x}
\left(
D(x)\frac{\partial T}{\partial x}
\right).
$$

Expanding the derivative gives:

$$
\frac{\partial T}{\partial t}
=
D(x)
\frac{\partial^2 T}{\partial x^2}
+
\frac{dD}{dx}
\frac{\partial T}{\partial x}.
$$

Therefore, replacing the constant coefficient $D$ by a diagonal matrix $D(x)$ multiplying the usual Laplacian is only an approximation. It ignores the additional gradient term. This is why the repository distinguishes between a first non-conservative variable-diffusion attempt and a corrected conservative discretization.

---

## 4. Initial and boundary conditions

The simulations use a one-dimensional spatial domain:

$$
x\in[-1,1].
$$

The initial temperature profile is a localized Gaussian peak:

$$
T(x,0)=100e^{-20x^2}.
$$

This profile is high near the centre and small near the boundaries. It is a good test because it allows the simulation to show the expected diffusion process clearly: the central peak decreases while heat spreads toward neighbouring points.

The main boundary conditions are homogeneous Dirichlet conditions:

$$
T(-1,t)=0,
\qquad
T(1,t)=0.
$$

These conditions represent boundaries held at zero temperature. Physically, the boundaries act as thermal reservoirs that absorb heat. Numerically, these values must be imposed at every time step; otherwise, the boundary values can drift due to roundoff or matrix operations.

---

## 5. Spatial discretization

The continuous domain is replaced by a uniform grid:

$$
x_j=x_{\min}+j\Delta x.
$$

The temperature field becomes a vector:

$$
\mathbf{T}^n=
\left(
T_0^n,
T_1^n,
\ldots,
T_{N-1}^n
\right)^T.
$$

The second derivative is approximated by the centred finite-difference stencil:

$$
\frac{\partial^2 T}{\partial x^2}
\bigg|_{x_j}
\approx
\frac{
T_{j+1}^n
-
2T_j^n
+
T_{j-1}^n
}{\Delta x^2}.
$$

This stencil is second order in space. It is local: each grid point is coupled only to its two nearest neighbours. In matrix form, this produces a tridiagonal discrete Laplacian.

The discrete Laplacian has the structure:

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

This matrix representation is useful because the whole time evolution can be written compactly using linear algebra.

---

## 6. FTCS explicit method

FTCS means Forward Time, Centered Space. The time derivative is approximated with a forward difference, while the spatial derivative uses the centred second-derivative stencil.

The method is:

$$
\frac{T_j^{n+1}-T_j^n}{\Delta t}
=
D
\frac{
T_{j+1}^n
-
2T_j^n
+
T_{j-1}^n
}{\Delta x^2}.
$$

Solving for $T_j^{n+1}$ gives:

$$
\begin{aligned}
T_j^{n+1}
&=
T_j^n
+
r
\left(
T_{j+1}^n
-
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

FTCS is explicit because the new value is computed only from known values at time level $n$. This makes the method simple and fast per step.

The disadvantage is stability. For the one-dimensional heat equation, FTCS is stable only if:

$$
r\leq\frac{1}{2}.
$$

Equivalently:

$$
\Delta t
\leq
\frac{\Delta x^2}{2D}.
$$

This condition becomes very restrictive when the grid is refined, because the allowed timestep scales as $\Delta x^2$.

---

## 7. Crank-Nicolson method

Crank-Nicolson averages the diffusion operator between the old time level $n$ and the new time level $n+1$. For constant diffusion, the scheme can be written as:

$$
\frac{T_j^{n+1}-T_j^n}{\Delta t}
=
\frac{D}{2}
\left[
\left(
\frac{\partial^2 T}{\partial x^2}
\right)_j^{n}
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
\mathbf{T}^{n+1}
=
\left(
I+\frac{D\Delta t}{2}L
\right)
\mathbf{T}^{n}.
$$

The method is implicit because $\mathbf{T}^{n+1}$ appears inside a linear system. Therefore, each time step requires solving a matrix equation.

Crank-Nicolson is more expensive per time step than FTCS, but it is much more stable. It is second order in time and is unconditionally stable for the linear heat equation. In practice, this means that the timestep can be chosen mainly from accuracy requirements rather than from a strict stability limit.

---

## 8. Von Neumann stability analysis

Von Neumann stability analysis studies how each Fourier mode evolves under a numerical scheme. A small numerical error can be decomposed into modes:

$$
T_j^n
=
Q^n e^{ikj\Delta x}.
$$

The factor $Q$ is the amplification factor. Stability requires:

$$
|Q|\leq 1.
$$

For FTCS applied to the heat equation, substituting a Fourier mode gives:

$$
Q
=
1
-
4r
\sin^2
\left(
\frac{k\Delta x}{2}
\right).
$$

The worst case occurs for the highest-frequency mode represented by the grid. Enforcing $|Q|\leq1$ for all modes gives:

$$
r\leq\frac{1}{2}.
$$

This explains why FTCS becomes unstable when $\Delta t$ is too large. The instability usually appears as point-to-point oscillations that grow instead of being smoothed.

For Crank-Nicolson, the amplification factor is:

$$
Q
=
\frac{
1
-
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

This satisfies $|Q|\leq1$ for all positive $r$, explaining its unconditional stability.

---

## 9. Conservative discretization for variable diffusion

For variable diffusion, the physically correct form is:

$$
\frac{\partial T}{\partial t}
=
\frac{\partial}{\partial x}
\left(
D(x)\frac{\partial T}{\partial x}
\right).
$$

A conservative finite-difference discretization evaluates diffusion coefficients at cell interfaces:

$$
D_{j+1/2}
=
\frac{D_j+D_{j+1}}{2}.
$$

The discrete operator becomes:

$$
\frac{1}{\Delta x^2}
\left[
D_{j+1/2}
\left(
T_{j+1}-T_j
\right)
-
D_{j-1/2}
\left(
T_j-T_{j-1}
\right)
\right].
$$

This form is better than simply multiplying the standard Laplacian by $D_j$, because it describes heat fluxes between neighbouring cells. It respects the conservation-law origin of the heat equation.

---

## 10. Numerical interpretation

This project compares several important numerical ideas:

- explicit versus implicit time integration,
- stability versus accuracy,
- constant versus spatially variable diffusion,
- non-conservative versus conservative discretization,
- timestep restriction versus computational cost,
- matrix-based implementation using sparse operators,
- error behaviour as $\Delta t$ changes.

The main physical expectation is that temperature profiles should smooth over time. The main numerical challenge is to reproduce that smoothing without generating nonphysical oscillations, excessive artificial diffusion or unstable growth.
