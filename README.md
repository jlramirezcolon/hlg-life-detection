# Life/Nonlife Discrimination via Orbital Energy Gaps
Code and data used to reproduce the results outlined in "Distinguishing life from non-life via molecular frontier orbital energy gaps". 

# Citation
Ramírez-Colón JL, Ni Z, Carr CE. Distinguishing life from non-life via molecular frontier orbital energy gaps. In preparation.

# Compatibility
Tested with MATLAB R2024b on OS X 13.6.5 and 14.8.3. Requires the Statistics and Machine Learning, the Predictive Maintenance, and the System Identification toolboxes. All third-party code is included in the `3rd_party/` folder with appropriate licenses.

# Installation

## Get Scripts
```bash
# Clone the repository
git clone https://github.com/jlramirezcolon/hlg-master.git
cd hlg-master
```

Open MATLAB and add the repository to your path:
```matlab
addpath(genpath('path/to/hlg-master'))
```
## Data
**Amino.Acid.Database.Release.2026-07-20.xlsx**: Contains all amino acid samples and calculated molecular descriptors used in this study. 
- MATLAB: Version of 'Database' that is MATLAB compatible and used for the abundance-weighted analysis (abundance_weighted_analysis.m).
- MATLAB_Distributions: Database used for the unweighted analysis (molecular_descriptor_distribution.m).
- Database: Shows all amino acids profiles (columns) compiled from the literature. Samples are color-coded: pink (abiotic), orange (biotic), and blue (abiotic simulations).
- Properties: Molecular properties calculated for all identified amino acids in 'Database'.
- Properties Distributions: Version of 'Properties' compatible with the script in molecular_descriptor_distribution.m.
- Extraterrestrial Samples Info: Contains information relevant to each extraterrestrial sample (e.g., ID, category, classification,	amount extracted (g), analytical method, etc.). 
- Terrestrial Samples Info: Contains information relevant to each terrestrial sample (e.g., category, location,  amount extracted (g), analytical method, etc.). 
- Legends: Key for all the abbreviations and colors identifications used in 'MATLAB', 'MATLAB_Distirbutions', and 'Database'.
- References: Contains all the literature references for each of the samples within the database.

**sim_settings_2026-08-11.xlsx**: Configuration parameters for Bayesian simulation framework. 

## Usage

### Analysis Workflow

Run analyses in the following order to reproduce paper results:

1. **Molecular descriptor distribution analysis**
```matlab
molecular_descriptor_distribution
```
   Analyzes the distribution of 10 molecular descriptors across biotic and abiotic environments.

2. **Abundance-weighted classification**
```matlab
abundance_weighted_analysis
```
   Computes weighted molecular descriptor statistics, evaluates class separation metrics, and performs machine learning classification to identify the most discriminative features.

3. **Bayesian framework simulations**
```matlab
simulations
```
   Calculates Bayesian probabilities:
   - P(E|A): Probability of evidence given abiotic sample
   - P(E|B): Probability of evidence given biotic sample  
   - P(B|E): Posterior probability that sample is biotic

4. **Generate simulation plots**
```matlab
simulation_plots
```
   Produces publication figures from Bayesian simulation results.

### Reproducing Paper Figures

PDF versions of all figures saved from MATLAB figures can be found in /out/figures.

molecular_descriptor_distribution.m:
- Figure 1c. Box plots showing amino acid's HOMO-LUMO gap distributions by class
- Figure 1d. Bar chart of relative entropy values for different molecular descriptors
- Figure 1e. Statistical significance matrix showing pairwise comparisons between categories

abundance_weighted_analysis.m: 
- Figure 2b. Abundance-weighted molecular descriptors improve classification of biotic and abiotic amino acid samples
- Figure 3.  Performance of abundance-weighted molecular descriptors in biotic–abiotic classification
- Figure S1. Abundance of Amino Acids Across Biotic and Abiotic Classes
- Figure S2. Amino Acid Abundance Patterns by Biotic and Abiotic Subcategories
- Figure S3. Evaluation of Alternative Class Separation Methods Using Abundance-Weighted Amino Acid Descriptors to Differentiate Biotic and Abiotic Samples

simulations_plots.m: 
- Figure 4b. Heatmap of P (B|E) as function of the number of amino acids measured
- Figure S4. Probability of evidence given abiotic and biotic datasets
- Figure S5. Confidence of biogenicity at prior of 0.5 across metrics
- Figure S6. Effect of prior on evidence distribution and posterior confidence in biogenicity

### Helper Functions
- `addFeature.m`: Computes abundance-weighted statistical features from molecular properties
- `PlotFeature.m`: Creates histogram visualizations of features
- `PlotFeatures2d.m`: Creates 2d visualization of two features
- `LookupProperties.m`: Retrieves property values for molecules in data table

# License
This software is distributed under the GNU Affero General Public License v3 (AGPLv3). This algorithm is a product of the Planetary eXploration Lab (PXL) at Georgia Tech. We are committed to Open Science and have released this implementation under AGPLv3 to ensure the life-detection community benefits from transparent, peer-reviewed methods.

Interested in Collaboration? If you are interested in integrating this algorithm into flight hardware, proprietary software suites, or NASA mission proposals where AGPLv3 compliance may be complex, please reach out. We would love to discuss science team partnerships, validation support, or alternative licensing. Visit [pxl.earth](https://www.pxl.earth/).
