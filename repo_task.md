# Repository Build Task

Rebuild this repository as a clean computational-physics portfolio project for a numerical modelling practice on the one-dimensional heat equation.

The project studies the 1D heat/diffusion equation, explicit FTCS time stepping, Crank-Nicolson time stepping, constant diffusion, spatially variable diffusion $D(x)$, conservative and non-conservative discretizations, homogeneous Dirichlet boundary conditions, Von Neumann stability, numerical error, runtime comparison and diffusion-amplitude sweeps.

## Rules

- Do not invent results, figures or plots.
- Do not modify the numerical logic of the Python code.
- Do not change physical parameters.
- Do not summarize the source documentation files.
- Copy the source files exactly into the target documentation files.
- Use GitHub-compatible Markdown math only: inline `$...$` and display `$$ ... $$` with delimiters on separate lines.
- Multi-line equations must use `\begin{aligned}...\end{aligned}`.
- Do not use `\( ... \)` or `\[ ... \]`.
- Do not leave broken image links.
- If an image listed in the manifest is missing, remove only that image block.
- Do not include `.DS_Store`, `__MACOSX`, `._*`, LaTeX auxiliary files or duplicate images.

## Target structure

```text
heat-equation-ftcs-crank-nicolson/
├── README.md
├── requirements.txt
├── .gitignore
├── src/
│   ├── heat_equation_ftcs_crank_nicolson.py
│   └── ftcs_matrix_helper.py
├── docs/
│   ├── theory.md
│   ├── numerical_method.md
│   ├── results_summary.md
│   └── sources_and_notes.md
├── figures/
│   ├── constant_diffusion/
│   ├── variable_diffusion/
│   ├── stability/
│   ├── error_analysis/
│   ├── runtime/
│   └── alpha_sweeps/
└── raw_upload/
```

## Copy these source files exactly

- `raw_upload/README_source.md` -> `README.md`
- `raw_upload/theory_source.md` -> `docs/theory.md`
- `raw_upload/numerical_method_source.md` -> `docs/numerical_method.md`
- `raw_upload/results_source.md` -> `docs/results_summary.md`
- `raw_upload/sources_and_notes_source.md` -> `docs/sources_and_notes.md`

## Code organization

Move the main Python file to:

```text
src/heat_equation_ftcs_crank_nicolson.py
```

Possible original filename:

```text
ultimaversionpractica4modelizacion.py
```

If `ftcs_diffusion.py` exists, move it to:

```text
src/ftcs_matrix_helper.py
```

If no helper file exists, do not create one.

## requirements.txt

Create `requirements.txt` with exactly:

```text
numpy>=1.24
scipy>=1.10
matplotlib>=3.7
```

## .gitignore

Create `.gitignore` with exactly:

```gitignore
# Python
__pycache__/
*.py[cod]
*.pyo
*.pyd
.venv/
venv/
env/

# macOS
.DS_Store
__MACOSX/
._*

# Jupyter
.ipynb_checkpoints/

# LaTeX auxiliary files
*.aux
*.log
*.out
*.toc
*.fls
*.fdb_latexmk
*.synctex.gz

# Generated numerical outputs
output/
outputs/
tmp/
*.npy
*.npz
```

## Final verification

Verify that all docs render correctly, all image links point to existing files or have been removed, all math uses GitHub-compatible syntax, the Python numerical logic is unchanged, and the final repository looks like a professional computational-physics portfolio project. Make one clean commit.
