# Life/Nonlife Discrimination via Orbital Energy Gaps
Code and data used to reproduce the results in "Distinguishing life from non-life via molecular frontier orbital energy gaps".

# Citation
Ramírez-Colón JL, Ni Z, Carr CE. Distinguishing life from non-life via molecular frontier orbital energy gaps. Submitted.

# Compatibility
Tested with MATLAB R2025b on macOS. Requires the Statistics and Machine Learning, the Predictive Maintenance, and the System Identification toolboxes. Third-party code is included in `code/3rdparty/` with its license.

The figure scripts in `figure_code/` need Python 3 and macOS (see `figure_code/README.md`).

# Installation

```bash
git clone https://github.com/jlramirezcolon/hlg-life-detection
cd hlg-life-detection
```

Open MATLAB in the repository folder. Every script uses paths relative to it.

# Data

**`data/Amino.Acid.Database.Release.2026-11-21.xlsx`** holds all amino acid (AA) samples and molecular descriptors used in the study.
- **MATLAB:** version of 'Database' read by the abundance-weighted analysis (`abundance_weighted_analysis.m`).
- **MATLAB_Distributions:** presence/absence data read by the unweighted analysis (`molecular_descriptor_distribution.m`).
- **Database:** AA profiles (columns) compiled from the literature. Samples are color-coded: pink (abiotic), orange (biotic) and blue (abiotic simulations).
- **Properties:** molecular descriptors of every AA in 'Database'. The two electronic descriptors are:
  - `HLG_wB97XD_TZVP_pH_7_eV`: HOMO-LUMO gap, ωB97X-D3/def2-TZVP/COSMO(water), Boltzmann-weighted over conformers;
  - `S1_TDDFT_eV`: optical gap, TDDFT ωB97X/def2-TZVPD/CPCM(water), highest-weight conformer.
  - `Dipole_moment_Debye`: dipole moment (Debye) from the same single points as the HLG, Boltzmann-weighted mean of the magnitude over conformers.
- **Properties Distributions:** 'Properties' in the layout read by `molecular_descriptor_distribution.m`. It also covers the seven AAs that occur only in the presence/absence data.
- **Statistics Tables:** summary tables cited in the paper (Supplemental Dataset). They are generated from the CSV files in `out/`.
- **Extraterrestrial Samples Info / Terrestrial Samples Info:** for each sample, its ID, category, classification, amount extracted, analytical method and similar details.
- **Legends:** abbreviations and color codes used in 'Database' and 'MATLAB_Distributions'.
- **pH Analysis:** dominant protonation state of each AA at pH 7.
- **References:** literature source of each sample.

**`quantum_chemistry/`** holds the conformer geometries, Boltzmann weights, and per-conformer HOMO, LUMO, HLG and S1 values behind the electronic descriptors (HLG, S1 and dipole moment) (see `quantum_chemistry/README.md`).

**`settings/`** holds the LUMOS simulation settings: `sim_settings_HLG.xlsx`, and `sim_settings_S1.xlsx` with the S1 bin edges.

# Usage

Run the analyses in this order:

1. **Descriptor distributions (unweighted)**
   ```matlab
   export_distribution_results
   ```
   Compares the distributions of 10 molecular descriptors between biotic and abiotic samples (relative entropy, and a Wilcoxon rank-sum test with each AA counted once per subcategory).

2. **Abundance-weighted analysis and classifiers**
   ```matlab
   export_weighted_results
   export_rankings
   ```
   Computes the abundance-weighted mean, variance (wvar) and Gini coefficient of each descriptor (30 features), ranks them for class separation, and trains the classifiers. Two switches near the top of the script control it:
   - `include_S1 = true` adds the optical gap (S1) as an 11th descriptor (33 features).
   - `ml_predictors` picks the classifier inputs: `'HLG'` (HLG wvar), `'S1'` (S1 wvar, needs `include_S1 = true`) or `'all'`. Together these give the three columns of Suppl. Table S5.

3. **LUMOS simulations**
   ```matlab
   simulations
   ```
   Monte Carlo simulations of the Bayesian framework:
   - P(E|A): probability of the evidence given an abiotic sample;
   - P(E|B): probability of the evidence given a biotic sample;
   - P(B|E): posterior probability that the sample is biotic.

   Set `descriptor = 'HLG'` or `'S1'` near the top. Each run takes about 30 minutes per metric and writes `out/simulations/`. Those outputs are about 500 MB per file and are not included in the repository.

4. **Simulation plots**
   ```matlab
   simulation_plots
   ```
   Plots the simulation results. Set the same `descriptor` as in step 3. The script reads the 1,000,000-replicate outputs (`1000000_class_freq_uniform_wor_excl`). If you change the number of replicates in the settings file, change `suffix_dir` in `simulation_plots.m` to match. Exporting the vector versions of Suppl. Figs. S4–S6 can take about an hour.

# Outputs

`export_distribution_results` runs `molecular_descriptor_distribution.m` and `export_weighted_results` runs `abundance_weighted_analysis.m`; each then saves the numbers below to `out/`. `export_rankings` adds the chi-squared and MRMR scores. Running the two analysis scripts on their own gives the same figures without the CSV files.

`out/` holds the results of steps 1 and 2 (HLG-only run):

| File | Contents |
|---|---|
| `metrics_distribution.csv` | Class statistics and relative entropy (nats) per descriptor (Fig. 1c–d) |
| `statistics_table4_*.csv`, `subcategory_ranges_*.csv`, `ranksum_p_values.csv`, `fig1c_values.csv` | Subcategory statistics and pairwise rank-sum tests |
| `features_per_sample.csv` | The 30 abundance-weighted features of each sample |
| `metrics_weighted.csv`, `rankings_chi2_mrmr.csv` | Relative entropy, AUC, chi-squared and MRMR scores of each feature (Fig. 3a, Suppl. Fig. S3) |
| `ml/`, `ml_all_features/`, `ml_S1/` | Classifier results with HLG wvar, all features and S1 wvar as predictors (Suppl. Table S5) |
| `figures/` | MATLAB versions of the figures (Fig1e and the extra histograms are not in the paper) |

`out_with_S1/` holds the features from the run with `include_S1 = true` (Suppl. Fig. S7b).

## Figures

The published layouts are drawn by the Python scripts in `figure_code/` (see its README):
- **Fig. 1:** b–d from `molecular_descriptor_distribution` outputs and `quantum_chemistry/`.
- **Figs. 2–3 and Suppl. Figs. S1–S3:** from `abundance_weighted_analysis` outputs.
- **Fig. 4b and Suppl. Figs. S4–S6:** from the HLG simulations.
- **Suppl. Fig. S7:** S1 against the HLG.
- **Suppl. Fig. S8:** from the S1 simulations.

# Helper Functions
- `addFeature.m`: computes abundance-weighted statistical features from molecular properties.
- `PlotFeature.m`: histogram of one feature.
- `PlotFeatures2d.m`: two features against each other.
- `LookupProperties.m`: retrieves property values for molecules in the data table.
- `sampleWithoutConflicts.m`: draws AAs without replacement, keeping D- and L- entries mutually exclusive.

# Licensing

This repository is dual-licensed by content type.

- **Code** (`*.m`, `code/`, `figure_code/*.py`, `figure_code/*.swift`): GNU Affero General Public License v3.0 only. See [`LICENSE.md`](LICENSE.md). Third-party components keep their own terms: `code/3rdparty/gini/` is BSD-2-Clause (© 2010 Yvan Lengwiler).
- **Data** (`data/`, `settings/`, `quantum_chemistry/`, `out/`, `out_with_S1/`): Creative Commons Attribution 4.0 International (CC BY 4.0). See [`LICENSE-DATA.md`](LICENSE-DATA.md).

This algorithm is a product of the Planetary eXploration Lab (PXL) at Georgia Tech. We are committed to Open Science and have released this implementation under AGPLv3 so that the life-detection community benefits from transparent, peer-reviewed methods.

Interested in collaboration? If you are interested in integrating this algorithm into flight hardware, proprietary software suites, or NASA mission proposals where AGPLv3 compliance may be complex, please reach out. We would love to discuss science team partnerships, validation support, or alternative licensing. Visit [pxl.earth](https://www.pxl.earth/).
