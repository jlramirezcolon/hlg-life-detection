"""Fig. 1b: frontier orbitals of zwitterionic glycine (the pH-7 state used throughout the paper).

The paper's HLG for glycine is a Psi4 wB97X-D3/def2-TZVP/COSMO(water) single point on its one AIMNet2 conformer
(ensemble_tzvp; HOMO -0.345319, LUMO 0.105802 Eh). PySCF has no wB97X-D3, so the orbitals for the picture are
recomputed locally on the same geometry with wB97X-D's exchange-correlation part (libxc; the dispersion term does not
change orbitals) and def2-TZVP/COSMO(water). The script prints both sets of energies so the match can be checked;
panel b labels use the Psi4 values (the database values).

Usage: ~/orbitals-venv/bin/python fig1b_orbitals.py   -> fig1b_data/glycine.{homo,lumo}.cube, energies.json
"""
import json
from pathlib import Path

from pyscf import dft, gto
from pyscf.tools import cubegen

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent   # repository root
OUT = HERE / 'fig1b_data'
EH = 27.211386
XYZ = ROOT / 'quantum_chemistry' / 'geometries' / 'glycine__zwitterionic__c1.xyz'
PSI4 = {'homo_Eh': -0.345319, 'lumo_Eh': 0.105802}   # ensemble_tzvp/results.csv, glycine

lines = XYZ.read_text().strip().splitlines()[2:]
mol = gto.M(atom='\n'.join(lines), charge=0, spin=0, basis='def2-tzvp', unit='Angstrom', verbose=0)
mf = dft.RKS(mol, xc='hyb_gga_xc_wb97x_d').density_fit().PCM()
mf.with_solvent.method = 'COSMO'
mf.with_solvent.eps = 78.3553
mf.kernel()
assert mf.converged
h = mol.nelectron // 2 - 1
OUT.mkdir(exist_ok=True)
for name, idx in (('homo', h), ('lumo', h + 1)):
    cubegen.orbital(mol, str(OUT / f'glycine.{name}.cube'), mf.mo_coeff[:, idx], resolution=0.12, margin=4.0)
res = {'pyscf_homo_eV': mf.mo_energy[h] * EH, 'pyscf_lumo_eV': mf.mo_energy[h + 1] * EH,
       'psi4_homo_eV': PSI4['homo_Eh'] * EH, 'psi4_lumo_eV': PSI4['lumo_Eh'] * EH}
res['pyscf_gap_eV'] = res['pyscf_lumo_eV'] - res['pyscf_homo_eV']
res['psi4_gap_eV'] = res['psi4_lumo_eV'] - res['psi4_homo_eV']
(OUT / 'energies.json').write_text(json.dumps(res, indent=1))
print(json.dumps({k: round(v, 3) for k, v in res.items()}, indent=1))
