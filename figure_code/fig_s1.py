"""Figure: the optical gap (S1) against the HLG. Helvetica Neue Light, panel letters.

  a  S1 vs HLG for the 55 amino acids of the abundance database (aromatic and sulfur-containing amino acids marked)
  b  ROC curves of the abundance-weighted variance (wvar) for the HLG and S1 (MATLAB features, out_with_S1/)
  c  wvar AUC with two functionals on identical conformers (main level vs wB97X-3c), HLG vs S1
     main level: wB97X-D3 for the HLG (Psi4); wB97X with D3(0) for S1 (GPU4PySCF has no wB97X-D3)
  d  epsilon-aminocaproic acid relative to glycine for each descriptor and functional

Values are the database values (3 decimals), so AUCs match the MATLAB pipeline. Also writes fig_s1_stats.csv.
Inputs: s1_inputs.py (database + quantum_chemistry/amino_acids.csv), out_with_S1/features_per_sample.csv.
Usage: python3 fig_s1.py <out.pdf>
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent   # repository root
sys.path.insert(0, str(HERE))
from hnfont import use_helvetica_neue_light  # noqa: E402

use_helvetica_neue_light()
FS = 8
plt.rcParams.update({'font.size': FS, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': True,
                     'axes.spines.right': True, 'legend.frameon': True})   # closed boxes
# tones of the main figures: Fig. 4's LUMOS blue and criteria gold, Fig. 3a's wine red, and a pastel olive
HLG_C, S1_C = '#7e9bb1', '#d45f69'      # blue = HLG, wine = S1
AROM_C, SULF_C = '#cba95c', '#b5bd6c'   # panel a: gold = aromatic, olive = thiol/thioether
GREY, INK = '#9aa0a6', '#222222'
AROM = {'phenylalanine': 'Phe', 'tyrosine': 'Tyr', 'tryptophan': 'Trp', 'histidine': 'His'}
SULF = {'methionine': 'Met', 'cysteine': 'Cys'}
OTHER_LABELS = {'taurine': 'Tau', 'aminocaproic-acid': 'ε-aca', 'glycine': 'Gly'}


def r3(d):
    return {k: round(float(v), 3) for k, v in d.items()}


def main(out):
    import s1_inputs as g
    auc, names = g.auc_with, g.ens_name
    hlg = r3({r['amino_acid']: a for r, a in zip(g.eh, g.d3)})
    hlg3c = r3({r['amino_acid']: b for r, b in zip(g.eh, g.c3)})
    s1 = r3(g.s1)
    s13c = r3(g.s1_3c)
    A = lambda d: auc([d.get(n, np.nan) for n in names])
    aucs = {('HLG', 'ωB97X-D3'): A(hlg), ('HLG', 'ωB97X-3c'): A(hlg3c), ('S1', 'ωB97X-D3'): A(s1),
            ('S1', 'ωB97X-3c'): A(s13c)}

    F = list(csv.DictReader(open(ROOT / 'out_with_S1' / 'features_per_sample.csv')))
    y = np.array([f['Class'] == '1' for f in F])
    feats = {'HLG': np.array([float(f['HLG_wB97XD_TZVP_pH_7_eV_wvar']) for f in F]),
             'S1': np.array([float(f['S1_TDDFT_eV_wvar']) for f in F])}

    fig = plt.figure(figsize=(7.0, 5.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    ax = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(4)]

    # a: S1 vs HLG
    aa = sorted(set(hlg) & set(s1))
    x = np.array([hlg[a] for a in aa]); yv = np.array([s1[a] for a in aa])
    r = np.corrcoef(x, yv)[0, 1]
    for a, xi, yi in zip(aa, x, yv):
        c = AROM_C if a in AROM else (SULF_C if a in SULF else GREY)
        ax[0].scatter(xi, yi, s=30 if c != GREY else 20, color=c, edgecolor='white', linewidth=0.4, zorder=3)
        lab = AROM.get(a) or SULF.get(a) or OTHER_LABELS.get(a)
        if lab == 'ε-aca':
            ax[0].annotate(lab, (xi, yi), xytext=(-34, -22), textcoords='offset points', fontsize=FS - 1, color=INK,
                           arrowprops=dict(arrowstyle='-', lw=0.5, color=INK))
            ax[0].scatter(xi, yi, s=30, facecolor='none', edgecolor=INK, linewidth=0.6, zorder=4)
        elif lab == 'Tau':   # next to the right frame: label on the left
            ax[0].annotate(lab, (xi, yi), xytext=(-4, -2), textcoords='offset points', ha='right',
                           fontsize=FS - 1, color=INK)
        elif lab:
            ax[0].annotate(lab, (xi, yi), xytext=(4, -2 if a != 'glycine' else 3), textcoords='offset points',
                           fontsize=FS - 1, color=INK)
    ax[0].scatter([], [], s=30, color=AROM_C, label='Aromatic'); ax[0].scatter([], [], s=30, color=SULF_C, label='Thiol/thioether')
    ax[0].scatter([], [], s=20, color=GREY, label='Other')
    ax[0].legend(loc='upper left', fontsize=FS - 1, handlelength=1.0, handletextpad=0.35, borderpad=0.3, labelspacing=0.25, borderaxespad=0.4, frameon=True, edgecolor='#cccccc', framealpha=0.95)
    ax[0].set_xlabel('HLG, ωB97X-D3/def2-TZVP (eV)'); ax[0].set_ylabel('Optical gap S1, TDDFT (eV)')
    ax[0].text(0.97, 0.05, f'r = {r:.2f}, n = {len(aa)}', transform=ax[0].transAxes, ha='right', fontsize=FS - 0.5)

    # b: ROC curves
    for k, c in (('HLG', HLG_C), ('S1', S1_C)):
        fpr, tpr, _ = roc_curve(y, feats[k])
        ax[1].plot(fpr, tpr, color=c, lw=1.4, label=f'{k} wvar (AUC = {roc_auc_score(y, feats[k]):.3f})',
                   zorder=5, clip_on=False)   # drawn over the y-axis line
    ax[1].plot([0, 1], [0, 1], color=GREY, lw=0.6, ls='--')
    ax[1].set_xlim(0, 1); ax[1].set_ylim(0, 1.01)
    ax[1].set_xlabel('False positive rate'); ax[1].set_ylabel('True positive rate')
    ax[1].legend(loc='lower right', fontsize=FS - 1, handlelength=1.0, handletextpad=0.35, borderpad=0.3, labelspacing=0.25, borderaxespad=0.4, frameon=True, edgecolor='#cccccc', framealpha=0.95)

    # c: functional robustness
    xs = np.arange(2); w = 0.36
    for j, (k, c) in enumerate((('HLG', HLG_C), ('S1', S1_C))):
        vals = [aucs[(k, f)] for f in ('ωB97X-D3', 'ωB97X-3c')]
        bars = ax[2].bar(xs + (j - 0.5) * w, vals, w * 0.92, color=c, label=k)
        for b, v in zip(bars, vals):
            ax[2].text(b.get_x() + b.get_width() / 2, v + 0.008, f'{v:.3f}', ha='center', fontsize=FS - 1)
    ax[2].set_xticks(xs); ax[2].set_xticklabels(['ωB97X(-D3)', 'ωB97X-3c']); ax[2].set_ylim(0.5, 1.1); ax[2].set_yticks(np.arange(0.5, 1.01, 0.1))
    ax[2].set_ylabel('wvar AUC (same conformers)'); ax[2].legend(loc='upper center', ncol=2, fontsize=FS - 1, columnspacing=0.8, handlelength=1.0, handletextpad=0.35, borderpad=0.3, labelspacing=0.25, borderaxespad=0.4, frameon=True, edgecolor='#cccccc', framealpha=0.95)

    # d: epsilon-aminocaproic acid relative to glycine
    e, gly = 'aminocaproic-acid', 'glycine'
    vals = {('HLG', 'ωB97X-D3'): hlg[e] - hlg[gly], ('HLG', 'ωB97X-3c'): hlg3c[e] - hlg3c[gly],
            ('S1', 'ωB97X-D3'): s1[e] - s1[gly], ('S1', 'ωB97X-3c'): s13c[e] - s13c[gly]}
    for j, (k, c) in enumerate((('HLG', HLG_C), ('S1', S1_C))):
        v = [vals[(k, f)] for f in ('ωB97X-D3', 'ωB97X-3c')]
        bars = ax[3].bar(xs + (j - 0.5) * w, v, w * 0.92, color=c, label=k)
        for b, vv in zip(bars, v):
            ax[3].text(b.get_x() + b.get_width() / 2, vv - 0.03, f'{vv:+.2f}', ha='center', va='top', fontsize=FS - 1)
    ax[3].axhline(0, color=INK, lw=0.6)
    ax[3].set_xticks(xs); ax[3].set_xticklabels(['ωB97X(-D3)', 'ωB97X-3c']); ax[3].set_ylim(-1.05, 0.1)
    ax[3].set_ylabel('ε-aminocaproic acid − glycine (eV)')
    ax[3].set_ylim(-1.05, 0.22)   # room above zero for the legend
    ax[3].set_yticks(np.arange(-1.0, 0.01, 0.2))
    ax[3].legend(loc='upper center', ncol=2, columnspacing=0.8, fontsize=FS - 1, handlelength=1.0, handletextpad=0.35, borderpad=0.3, labelspacing=0.25, borderaxespad=0.4, frameon=True, edgecolor='#cccccc', framealpha=0.95)   # one row on top, like c

    for a_, l in zip(ax, 'abcd'):
        a_.text(-0.16, 1.03, l, transform=a_.transAxes, fontsize=FS + 3, fontweight='bold', va='bottom')
    fig.savefig(out)
    print('wrote', out)
    with open(Path(out).with_name('fig_s1_stats.csv'), 'w', newline='') as f:
        wr = csv.writer(f); wr.writerow(['quantity', 'value'])
        wr.writerow(['pearson_r_S1_vs_HLG', f'{r:.3f}'])
        for (k, fn), v in aucs.items():
            wr.writerow([f'AUC_wvar_{k}_{fn}', f'{v:.3f}'])
        for (k, fn), v in vals.items():
            wr.writerow([f'eaca_minus_gly_{k}_{fn}', f'{v:.3f}'])
    print(dict((f'{k} {fn}', round(v, 3)) for (k, fn), v in aucs.items()), 'r', round(r, 3),
          {f'{k} {fn}': round(v, 2) for (k, fn), v in vals.items()})


if __name__ == '__main__':
    main(sys.argv[1])
