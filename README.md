# Heat Equation Diffusion Solver: FTCS and Crank-Nicolson Methods

## Scientific motivation
Diffusion and heat-transport models are foundational in computational physics. This repository is structured as a professional portfolio project for studying one-dimensional parabolic partial differential equations with finite-difference methods.

## Heat equation and diffusion model
The project studies the one-dimensional diffusion equation:

$$
\frac{\partial u}{\partial t} = \frac{\partial}{\partial x}\left(D(x)\frac{\partial u}{\partial x}\right),
$$

where $u(x,t)$ is the transported scalar field and $D(x)$ is the diffusion coefficient.

## Constant diffusion and spatially variable diffusion $D(x)$
Two model classes are planned:
- Constant diffusion, $D(x)=\alpha$.
- Spatially variable diffusion, $D(x)$ varying across the domain.

## FTCS explicit method
The Forward-Time Centered-Space (FTCS) scheme will be used as an explicit baseline for the transient diffusion problem.

## Crank-Nicolson method
The Crank-Nicolson method will be used as an implicit second-order approach in time for improved stability and accuracy.

## Stability and Von Neumann analysis
A dedicated stability study will compare method behavior under timestep variation, including Von Neumann analysis for the constant-diffusion case.

## Numerical error versus timestep
A systematic error analysis will quantify temporal discretization error trends versus $\Delta t$.

## Runtime comparison
Runtime benchmarks will compare explicit and implicit approaches under matched simulation settings.

## Repository structure
- `src/` — source code (to be added).
- `docs/` — theory, methods, and summary notes.
- `figures/constant_diffusion/` — figures for constant-$D$ experiments.
- `figures/variable_diffusion/` — figures for variable-$D(x)$ experiments.
- `figures/stability/` — stability-analysis plots.
- `figures/error_analysis/` — error-versus-timestep plots.
- `figures/runtime/` — runtime comparison plots.
- `figures/alpha_sweeps/` — parameter-sweep plots.
- `raw_upload/` — uploaded raw project assets and outputs.

## Skills demonstrated
- Mathematical modelling of diffusion PDEs.
- Finite-difference discretization design.
- Explicit and implicit time integration methods.
- Stability, error, and runtime analysis workflows.
- Scientific project organization and reproducible documentation.

## Project status
This is the initial scaffold only. Python implementation files and generated figures will be added after code and simulation outputs are uploaded.

## Author
**María Lourdes Domínguez Cacho**  
Final-semester Physics student, University of Alicante  
GitHub: [mldc3](https://github.com/mldc3)
