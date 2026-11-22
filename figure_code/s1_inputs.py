"""Inputs for fig_s1.py, read from this repository: per-sample amino-acid abundances (database, MATLAB sheet) and
per-amino-acid HLG and S1 values (quantum_chemistry/amino_acids.csv).

  ens_name   amino-acid key (quantum_chemistry name) for each abundance row of the MATLAB sheet
  eh, d3, c3 the 55 abundance-database amino acids and their Boltzmann HLG at wB97X-D3/def2-TZVP and wB97X-3c
  s1, s1_3c  S1 at wB97X/def2-TZVPD and at wB97X-3c, by amino acid
  auc_with   wvar AUC (biotic vs abiotic) for a list of values aligned with ens_name
"""
import csv
from pathlib import Path

import numpy as np
import openpyxl
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
DB = ROOT / 'data' / 'Amino.Acid.Database.Release.2026-11-21.xlsx'
QC = list(csv.DictReader(open(ROOT / 'quantum_chemistry' / 'amino_acids.csv', encoding='utf-8')))

wb = openpyxl.load_workbook(DB, read_only=True, data_only=True)
rows = list(wb['MATLAB'].iter_rows(values_only=True))
S0 = 8                                                     # first sample column
lab = np.array([float(x) if x not in (None, 'NaN') else np.nan for x in rows[1][S0:]])
keep = ~np.isnan(lab)
cls = lab[keep].astype(bool)
body = [r for r in rows[3:] if r[5] is not None]


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


A = np.nan_to_num(np.array([[num(v) for v in r[S0:]] for r in body])[:, keep])
sym2name = {r[1]: r[0] for i, r in enumerate(wb['Properties'].iter_rows(values_only=True)) if i >= 3 and r[0]}
name2key = {n: q['amino_acid'] for q in QC for n in q['database_names'].split('; ')}
ens_name = [name2key.get(sym2name.get(r[5])) for r in body]

eh = [q for q in QC if q['in_abundance_database'] == 'True']
d3 = np.array([float(q['hlg_wB97XD3_TZVP_boltzmann_eV']) for q in eh])
c3 = np.array([float(q['hlg_wB97X3c_boltzmann_eV']) for q in eh])
s1 = {q['amino_acid']: float(q['s1_wB97X_TZVPD_eV']) for q in QC}
s1_3c = {q['amino_acid']: float(q['s1_wB97X3c_eV']) for q in QC if q['s1_wB97X3c_eV']}


def auc_with(vals):
    v = np.array(vals, dtype=float)
    m = ~np.isnan(v)
    Am, vm = A[m], v[m][:, None]
    t = Am.sum(0)
    mu = (Am * vm).sum(0) / t
    return roc_auc_score(cls, (Am * (vm - mu) ** 2).sum(0) / t)
