# Figure code

Python scripts that lay out the published figures from the analysis outputs. The MATLAB scripts in the repository
root compute every number and also save their own versions of the plots (`out/figures/`). These scripts redraw the
same data in the journal layout: Helvetica Neue Light, panel letters, and vector artwork for the schematic panels.

Run them from this folder. They need Python 3 with numpy, scipy, matplotlib, openpyxl, scikit-learn and Pillow.
`fig1b_orbitals.py` also needs PySCF, and `fig1b_render.py` needs PyVista. Placing the artwork uses Swift and Core
Graphics, so it runs on macOS only.

| Script | Figure | Inputs |
|---|---|---|
| `fig1_composite.py .. <out.pdf>` | Fig. 1 | `out/fig1c_values.csv`, `out/metrics_distribution.csv`, `out/statistics_table4_*.csv`, `fig1b_data/`, `artwork/Fig1a.pdf` |
| `fig1b_orbitals.py`, `fig1b_render.py` | Fig. 1b orbitals | glycine geometry in `quantum_chemistry/geometries/` |
| `main_figures.py .. <out dir>` | Figs. 2–4, Suppl. Fig. S3 | `out/` CSVs, `out/simulations/` (Fig. 4b), `artwork/Fig2a.pdf`, `Fig4a.pdf`, `Fig4c.pdf` |
| `sim_figures.py .. <out prefix> --fig4b <path>` | Suppl. Figs. S4–S6, Fig. 4b | `out/simulations/` from `simulations.m` |
| `fig_s1.py <out.pdf>` | Suppl. Fig. S7 | database, `quantum_chemistry/amino_acids.csv`, `out_with_S1/features_per_sample.csv` |
| `fig_s1_lumos.py .. <out.pdf> --field S1_TDDFT_eV` | Suppl. Fig. S8 | `out/simulations/S1_TDDFT_eV/` from `simulations.m` with `descriptor = 'S1'` |

## Notes

- **Simulation outputs are not included.** Each `sim.mat` is about 500 MB. Run `simulations.m` first, about 30 min per metric.
- **Font:** the Helvetica Neue font files are not redistributed. On macOS, `hnfont.py` extracts them from the system font (`/System/Library/Fonts/HelveticaNeue.ttc`) on first use.
- **SVG output:** each `.svg` written next to a PDF has live text for editing. The artwork is left out, and a comment in the file gives where to place it.
