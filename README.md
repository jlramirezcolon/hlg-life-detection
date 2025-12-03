# LUMOS (Life Unveiled via Molecular Orbital Signatures)
Code and data used to reproduce the results outlined in "Distinguishing life from non-life via molecular frontier orbital energy gaps". 

# Citation
Ramírez-Colón JL, Ni Z, Carr CE. Distinguishing life from non-life via molecular frontier orbital energy gaps. In preparation.

# Compatibility
Tested with MATLAB R2024b on OS X 13.6.5. Requires the Statistics and Machine Learning, the Predictive Maintenance, and the System Identification toolboxes. All third-party code is included in the `3rd_party/` folder with appropriate licenses.

# Installation
## Get Scripts
Download: <https://github.com/jlramirezcolon/LUMOS/archive/master.zip> or use command line ```git clone git@github.com:jlramirezcolon/LUMOS.git```.

Unzip to preferred location, here denoted ```/LUMOS-master```.

Open MATLAB and add the repository to your path:
```matlab
addpath(genpath('path/to/LUMOS-master'))
```
## Data
**Amino.Acid.Database.v1.3.Release.2025-07-19.xlsx**: Contains all amino acid samples and calculated molecular descriptors used in this study. 
- MATLAB: Version of 'Database' that is MATLAB compatible and used for the abundance-weighted analysis (abundance_weighted_analysis.m).
- MATLAB_Distributions: Database used for the unweighted analysis (molecular_descriptor_distribution.m).
- Database: Shows all amino acids profiles (columns) compiled from the literature. Samples are color-coded: pink (abiotic), orange (biotic), and blue (abiotic simulations).
- Properties: Molecular properties calculated for all identified amino acids in 'Database'.
- Properties_Distributions: Version of 'Properties' compatible with the script in molecular_descriptor_distribution.m.
- Extraterrestrial_Samples_Info: Contains information relevant to each extraterrestrial sample (e.g., ID, category, classification,	amount extracted (g), analytical method, etc.). 
- Terrestrial_Samples_Info: Contains information relevant to each terrestrial sample (e.g., category, location,	amount extracted (g), analytical method, etc.). 
- Abiotic_Simulations: Contains information relevant to each simulated sample (e.g., conditions). 
- Legends: Key for all the abbreviations and colors identifications used in 'MATLAB', 'MATLAB_Distirbutions', and 'Database'.
- References: Contains all the literature references for each of the samples within the database.

## Run Analysis
In MATLAB, go to your ```/zerogseq-master``` path, and run the main script: ```zerogseq```. This will perform the same analysis as in the publication (see citation, above), running each script in turn, including running some twice, once each on the "Flight" or "Ground" datasets. See each script for details and instructions.

The analysis may take up to 5-6 hours.

The results of running this analysis in MATLAB include a series of EPS and/or PDF figures, replicating those in the paper, and various tab-delimited files. All times are elapsed time, and for reference, the start time is: 2017-11-17 18:28:51 UTC.

# License
Distributed under an MIT license. See LICENSE for details.
