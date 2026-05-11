# Figure Manifest

This file defines how uploaded figures must be organized, renamed and used in the documentation. Use the closest matching uploaded filename if accents, spacing or capitalization differ. Do not invent missing figures.

For `README.md`, use paths like:

```markdown
![Caption](figures/constant_diffusion/ftcs_constant_diffusion.png)
```

For files inside `docs/`, use paths like:

```markdown
![Caption](../figures/constant_diffusion/ftcs_constant_diffusion.png)
```

## Constant diffusion

### FTCS constant diffusion

Original filename:

```text
difusion FTCS D constante.png
```

Target path:

```text
figures/constant_diffusion/ftcs_constant_diffusion.png
```

Use in `README.md` and `docs/results_summary.md`.

### Crank-Nicolson constant diffusion

Original filename:

```text
difusion D no variable Crank Nicolson.png
```

Target path:

```text
figures/constant_diffusion/crank_nicolson_constant_diffusion.png
```

Use in `README.md` and `docs/results_summary.md`.

## Variable diffusion

### Non-conservative FTCS variable-diffusion attempt

Original filename:

```text
difusion 1D con D(x) variable mal hecho FCTS.png
```

Target path:

```text
figures/variable_diffusion/ftcs_variable_diffusion_nonconservative_attempt.png
```

Use in `docs/results_summary.md`.

### Conservative FTCS variable diffusion

Original filename:

```text
difusion 1D FTCS con derivada bienhecha.png
```

Target path:

```text
figures/variable_diffusion/ftcs_variable_diffusion_conservative.png
```

Use in `README.md` and `docs/results_summary.md`.

### Crank-Nicolson variable diffusion

Original filename:

```text
difucion 1D con D(x) variable CN.png
```

Target path:

```text
figures/variable_diffusion/crank_nicolson_variable_diffusion.png
```

Use in `README.md` and `docs/results_summary.md`.

### Additional Crank-Nicolson variable diffusion profile

Original filename:

```text
difusion 1D con crank nicolson con d(x) variable.png
```

Target path:

```text
figures/variable_diffusion/crank_nicolson_variable_diffusion_profile.png
```

Use in `docs/results_summary.md`.

## Stability

Original filename:

```text
convergenciametodoCrankNicolson.png
```

Target path:

```text
figures/stability/crank_nicolson_convergence.png
```

Use in `README.md` and `docs/results_summary.md`.

## Runtime

Original filename:

```text
comparacion tiempo de simulacion FTCS vs Crank nicolson.png
```

Target path:

```text
figures/runtime/ftcs_vs_crank_nicolson_runtime.png
```

Use in `README.md` and `docs/results_summary.md`.

## Error analysis

Original filename:

```text
error numerico ftcs con d(x) varibale al cambiar dt.png
```

Target path:

```text
figures/error_analysis/ftcs_error_vs_dt_variable_diffusion.png
```

Use in `README.md` and `docs/results_summary.md`.

Original filename:

```text
error numerico cn con D(x) variable al variar dt.png
```

Target path:

```text
figures/error_analysis/crank_nicolson_error_vs_dt_variable_diffusion.png
```

Use in `README.md` and `docs/results_summary.md`.

## Alpha sweeps

```text
solucion final para t=0.1s alpha =0.05.png -> figures/alpha_sweeps/final_profile_t01_alpha005.png
solucion final t = 0.1 s alpha 0.1.png -> figures/alpha_sweeps/final_profile_t01_alpha01.png
solucion 1s alpoha 0.05.png -> figures/alpha_sweeps/final_profile_t1_alpha005.png
solucion 1s alpha 0.1.png -> figures/alpha_sweeps/final_profile_t1_alpha01.png
solucion final t = 1s alpha = 0.2.png -> figures/alpha_sweeps/final_profile_t1_alpha02.png
solucion final t = 1s alpha 0.5.png -> figures/alpha_sweeps/final_profile_t1_alpha05.png
```

Use alpha sweep figures in `docs/results_summary.md`.

## Duplicates and files to avoid

Do not include:

```text
ifusion FTCS D constante.png
difusion 1D con crank nicolson con d(x) variable.png.png
comparacion de tiempo simulacion FTCS vs crank nicolson.png
.DS_Store
__MACOSX/
._*
*.aux
*.log
*.out
*.toc
*.fls
*.fdb_latexmk
*.synctex.gz
```
