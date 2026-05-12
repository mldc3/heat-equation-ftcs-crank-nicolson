# Numerical Method

This document explains how the heat-equation solvers are implemented in the code.

---

## 1. Computational domain

The code solves a one-dimensional diffusion problem on the interval:

$$
x\in[-1,1].
$$

The domain is discretized with a uniform grid:

$$
x_j=x_{\min}+j\Delta x.
$$

The initial condition is:

$$
T(x,0)=100e^{-20x^2}.
$$

After creating the initial profile, the boundary values are explicitly set to zero:

$$
T_0^n=0, \qquad T_{N-1}^n=0.
$$

These are homogeneous Dirichlet boundary conditions.

---

## 2. Discrete Laplacian construction

The second derivative is approximated by:

$$
\frac{\partial^2 T}{\partial x^2} \bigg|_{x_j} \approx \frac{ T_{j+1} - 2T_j + T_{j-1} }{\Delta x^2}.
$$

In the code, this operator is represented as a sparse tridiagonal matrix. The diagonal contains $-2$, and the upper and lower diagonals contain $1$.

Using sparse matrices is important because most entries of the derivative matrix are zero. Sparse storage avoids wasting memory and makes matrix-vector products more efficient.

---

## 3. FTCS implementation for constant diffusion

For constant diffusion, the FTCS update is:

$$
\begin{aligned}
T_j^{n+1} &= T_j^n +
\frac{D\Delta t}{\Delta x^2}
\left( T_{j+1}^n - 2T_j^n
+
T_{j-1}^n
\right).
\end{aligned}
$$

The code writes this in matrix form as:

$$
\mathbf{T}^{n+1} = \mathbf{T}^{n} +
D\Delta t L\mathbf{T}^{n}.
$$

Equivalently:

$$
\mathbf{T}^{n+1} = \left( I+D\Delta t L
\right)
\mathbf{T}^{n}.
$$

The boundary values are reset after each update. This prevents the numerical update from modifying the fixed-temperature boundaries.

---

## 4. Stability timestep for FTCS

The FTCS timestep is chosen from the stability condition:

$$
\frac{D\Delta t}{\Delta x^2}\leq\frac{1}{2}.
$$

The code uses a safety factor:

$$
\Delta t =
0.98
\frac{\Delta x^2}{2D}.
$$

When $D$ depends on position, the maximum value of the diffusion coefficient controls the stability limit:

$$
\Delta t =
0.98
\frac{\Delta x^2}{2\max(D(x))}.
$$

This is necessary because the region with largest diffusion is the most restrictive part of the domain.

---

## 5. First variable-diffusion attempt

One implementation uses a diagonal matrix containing $D(x)$ and multiplies it by the standard Laplacian. This corresponds roughly to:

$$
D(x)\frac{\partial^2 T}{\partial x^2}.
$$

This is computationally simple, but it is not the full conservative diffusion operator. The correct continuous equation is:

$$
\frac{\partial}{\partial x}
\left(
D(x)\frac{\partial T}{\partial x}
\right).
$$

The difference matters because:

$$
\frac{\partial}{\partial x} \left(
D(x)\frac{\partial T}{\partial x}
\right) = D(x) \frac{\partial^2 T}{\partial x^2} + \frac{dD}{dx} \frac{\partial T}{\partial x}.
$$

The first attempt is therefore useful as a diagnostic comparison, but the corrected conservative discretization is physically preferable.

---

## 6. Conservative variable-diffusion operator

The conservative discretization computes diffusion at interfaces:

$$
D_{j+1/2} = \frac{D_j+D_{j+1}}{2}.
$$

The discrete operator is:

$$
\begin{aligned} \left(
\nabla\cdot D\nabla T
\right)_j
&\approx
\frac{1}{\Delta x^2}
\left[
D_{j+1/2}
\left(
T_{j+1}-T_j \right) -
D_{j-1/2}
\left(
T_j-T_{j-1}
\right)
\right].
\end{aligned}
$$

In matrix form, this produces a tridiagonal operator with coefficients depending on neighbouring interface diffusivities.

This is the version that best matches the physical conservation-law derivation of the heat equation.

---

## 7. Crank-Nicolson implementation

The Crank-Nicolson scheme is written as:

$$
\left(
I-\frac{\Delta t}{2}A
\right)
\mathbf{T}^{n+1} =
\left(
I+\frac{\Delta t}{2}A
\right)
\mathbf{T}^{n},
$$

where $A$ is the discrete diffusion operator. For constant diffusion:

$$
A=DL.
$$

For variable diffusion, $A$ is the conservative variable-diffusion operator.

The code constructs two sparse matrices:

$$
M_{\mathrm{left}} =
I-\frac{\Delta t}{2}A,
$$

and

$$
M_{\mathrm{right}} =
I+\frac{\Delta t}{2}A.
$$

At each timestep, it computes:

$$
\mathbf{b} =
M_{\mathrm{right}}
\mathbf{T}^{n},
$$

and solves:

$$
M_{\mathrm{left}}
\mathbf{T}^{n+1} =
\mathbf{b}.
$$

This is done using sparse linear algebra.

---

## 8. Runtime comparison

FTCS has a cheap timestep, but the timestep must be small for stability. Crank-Nicolson has a more expensive timestep because it solves a linear system, but it can use larger timesteps without becoming unstable.

The comparison therefore studies a real numerical trade-off:

- FTCS: cheap per step, many steps required.
- Crank-Nicolson: expensive per step, fewer steps can be used.

---

## 9. Error versus timestep

The repository studies the dependence of the numerical error on $\Delta t$.

The error measures include maximum error and mean error:

$$
E_{\max} =
\max_j
\left|
T_j^{\mathrm{num}}
-
T_j^{\mathrm{ref}}
\right|.
$$

$$
E_{\mathrm{mean}} =
\frac{1}{N}
\sum_j
\left|
T_j^{\mathrm{num}} -
T_j^{\mathrm{ref}}
\right|.
$$

This allows the project to compare not only visual behaviour but also quantitative accuracy.

---

## 10. Alpha parameter sweeps

The coefficient $D(x)$ is controlled through a parameter $\alpha$:

$$
D(x)=\alpha
\left(
1+e^{-x^2}
\right).
$$

Increasing $\alpha$ increases the overall strength of diffusion. Physically, the temperature profile should spread faster and decay more strongly.

Numerically, increasing $\alpha$ also makes the FTCS stability restriction more severe, because the stable timestep decreases as $\max(D(x))$ increases.

---

## 11. Summary of implementation skills

The code demonstrates finite-difference discretization, sparse tridiagonal matrix construction, explicit FTCS time stepping, implicit Crank-Nicolson time stepping, Dirichlet boundary conditions, constant diffusion, variable diffusion, conservative flux discretization, sparse linear-system solving, timestep stability control, numerical error measurement, runtime benchmarking and generation of publication-style plots.
