# Quantum-chemical data

Inputs and per-conformer results behind the electronic descriptors in the database
(`data/Amino.Acid.Database.Release.2026-11-21.xlsx`, sheets *Properties* and *Properties Distributions*):

- `HLG_wB97XD_TZVP_pH_7_eV`: HOMO-LUMO gap, Boltzmann-weighted over conformers.
- `S1_TDDFT_eV`: optical gap (lowest singlet vertical excitation energy), highest-weight conformer.
- `Dipole_moment_Debye`: dipole moment, Boltzmann-weighted over conformers.

All calculations were run through [Rowan](https://rowansci.com). Psi4 1.11 was used for ground-state orbitals and GPU4PySCF for TDDFT. Full methods are in the paper (Methods, "Quantum-chemical calculations").

## Files

| File | Contents |
|---|---|
| `geometries/` | 226 conformer geometries (`<amino acid>__<state>__c<conformer>.xyz`), AIMNet2/ALPB(water)-optimized |
| `conformers.csv` | One row per conformer: Boltzmann weight, SMILES, geometry file, HOMO, LUMO, HLG and dipole moment at ωB97X-D3/def2-TZVP/COSMO(water), and HLG at ωB97X-3c/vDZP/COSMO(water) |
| `amino_acids.csv` | One row per amino acid (62): database names, protonation state, number of conformers, Boltzmann HLG with both methods, Boltzmann dipole moment, S1 (ωB97X/def2-TZVPD/CPCM(water)) with its conformer and oscillator strength, and S1 at ωB97X-3c/vDZP |
| `s1_conformers.csv` | S1 on every conformer of the eight amino acids used for the conformer-sensitivity test (Suppl. Table S3) |

## Conventions

- **Protonation state:** the dominant microstate at pH 7 (Rowan Starling workflow). Most amino acids are zwitterions; those with an ionizable side chain carry its charge (database sheet *pH Analysis*).
- **Conformers:** ETKDG generation, AIMNet2 (ωB97M-D3) optimization with ALPB(water), and Boltzmann weights at 298.15 K. Only conformers with weight above 0.01 were kept, so the stored weights of an amino acid sum to between 0.96 and 1. The HLG is the weighted mean using the weights renormalized over the retained conformers.
- **S1:** linear-response TDDFT without the Tamm-Dancoff approximation, six roots. Solvation is CPCM(water) with the engine's default non-equilibrium treatment (optical dielectric constant 1.78). The functional is ωB97X with D3(0) dispersion, because GPU4PySCF does not implement ωB97X-D3; the dispersion correction does not change excitation energies.
- **Dipole moment:** magnitude of the total dipole from Psi4, about the centre of mass (the reference point matters for the charged amino acids), converted from atomic units with 1 e·a₀ = 2.5417465 D. The database value is the Boltzmann-weighted mean of the magnitudes, with the same weights as the HLG.
- **Enantiomers:** D-, L- and DL- forms share the values of one representative structure; `database_names` lists every database entry that uses each row.
- **The seven extra amino acids:** these occur only in the presence/absence dataset (`in_abundance_database = False`) and were not part of the ωB97X-3c comparison.
- **Energy units:** eV. Orbital energies were converted from hartree with 1 Eh = 27.211386 eV.
