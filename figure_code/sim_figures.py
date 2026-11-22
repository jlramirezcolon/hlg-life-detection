"""LUMOS simulation figures (Figs. S4-S6 and Fig. 4b) redrawn from simulations.m output (sim.mat), in Helvetica
Neue Light with panel letters. A matplotlib port of simulation_plots.m: same data, same panels, same contour
levels and colour map; differences are styling only, plus
  * a fixed 0-1 colour scale (the MATLAB version auto-scaled empty panels to -1..1),
  * y-axis units by metric: the Gini coefficient is dimensionless and the weighted variance is in eV^2.

Usage: python3 sim_figures.py <run folder> <output prefix> [--field <field>] [--fig4b <path>] [--s6-range EV POST]
  e.g. HLG: python3 sim_figures.py .. ../out/SI. --fig4b ../out/Fig4b.pdf
  run folder: the repository root; it reads out/simulations/<field>/<metric>/<sub>/sim.mat (from simulations.m)
  writes <prefix>PofE.pdf, <prefix>EqualPrior(Metrics).pdf, <prefix>EffectPrior.pdf
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hnfont import use_helvetica_neue_light  # noqa: E402

use_helvetica_neue_light()
FS = 8.5
plt.rcParams.update({'font.size': FS, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'xtick.major.size': 2.5, 'ytick.major.size': 2.5})
LEVELS = [0.01, 0.1, 0.5, 0.9, 0.99]
# MATLAB drew the surface at 70% opacity on white; pre-blend the jet map with white instead (no seams between cells)
_jet = plt.get_cmap('jet')(np.linspace(0, 1, 256))
_jet[:, :3] = 0.7 * _jet[:, :3] + 0.3
HEAT_ALPHA = 0.95   # extra heatmap transparency over white (5% lighter)
GRID_ALPHA = 0.15   # strength of the grid seen through the heatmap
_jet[:, :3] = HEAT_ALPHA * _jet[:, :3] + (1 - HEAT_ALPHA)
CMAP = matplotlib.colors.ListedColormap(_jet, name='jet70')
CMAP.set_bad('white')
FIELD = 'HLG_wB97XD_TZVP_pH_7_eV'
SUB = '1000000_class_freq_uniform_wor_excl'


def load(root, metric, field=FIELD):
    m = sio.loadmat(str(Path(root) / 'out' / 'simulations' / field / metric / SUB / 'sim.mat'),
                    squeeze_me=True, struct_as_record=False)
    return m['op'], m['sim']


def label(op):
    name = str(op.field_label).split(' (')[0]          # 'HLG (eV)' -> 'HLG'
    metric = str(op.metric)
    unit = {'gini': '', 'wmean': ' (eV)', 'wvar': ' (eV$^2$)'}[metric]
    return f"{str(op.metric_label).capitalize()} of {name}{unit}"


def grid(op):
    x = np.atleast_1d(op.N_AA).astype(float)
    e = np.asarray(op.metric_edges, float)
    return x, (e[:-1] + e[1:]) / 2


def fit_top(op, cdf, ylim, level=0.99, margin=1.15):
    """Raise the top of ylim, if needed, so the `level` iso-probability contour fits inside with a margin."""
    _, y = grid(op)
    z = np.asarray(cdf, float)
    tops = [y[np.flatnonzero(row >= level)[0]] for row in z if np.any(row >= level)]
    need = ylim[0] + margin * (max(tops) - ylim[0]) if tops else ylim[1]   # margin measured from the axis bottom
    return (ylim[0], max(ylim[1], need))


def heat(ax, op, z, ylim, contours=True, title=None, ylab=None, levels=LEVELS, interior_labels=False, cax=None):
    x, y = grid(op)
    z = np.asarray(z, float)                             # (n_AA, bins)
    keep = (y >= ylim[0] - (y[1] - y[0])) & (y <= ylim[1] + (y[1] - y[0]))
    pm = ax.pcolormesh(x, y[keep], np.ma.masked_invalid(z[:, keep].T), cmap=CMAP, vmin=0, vmax=1, shading='nearest',
                       rasterized=True, linewidth=0, antialiased=False, zorder=1)
    if contours and np.nanmax(z[:, keep]) > 0:
        cs = ax.contour(x, y[keep], z[:, keep].T, levels=levels, colors='k', linewidths=0.5, linestyles='dashed')
        manual = None
        if interior_labels:
            # one label per level, away from the panel edges and staggered across n so neighbours don't stack
            lo, hi = ylim[0] + 0.01 * (ylim[1] - ylim[0]), ylim[1] - 0.04 * (ylim[1] - ylim[0])
            xs = np.linspace(x.min() + 0.25 * np.ptp(x), x.max() - 0.15 * np.ptp(x), len(cs.levels))
            manual = []
            for segs, xt in zip(cs.allsegs, xs):
                pts = np.vstack([sg for sg in segs if len(sg)]) if any(len(sg) for sg in segs) else np.empty((0, 2))
                pts = pts[(pts[:, 1] > lo) & (pts[:, 1] < hi)] if len(pts) else pts
                if len(pts):
                    manual.append(tuple(pts[np.argmin(np.abs(pts[:, 0] - xt))]))
        ax.clabel(cs, fmt=lambda v: f'{v:g}', fontsize=FS - 1.5, inline_spacing=2, manual=manual or False)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(*ylim)
    # grid as seen through a heatmap at HEAT_ALPHA opacity on white: the colours are pre-lightened (CMAP) and the
    # grid is drawn on top at the strength that would show through (exact for a white background, works in every PDF viewer)
    ax.grid(True, color='#4d4d4d', lw=0.4, alpha=GRID_ALPHA, zorder=1.5)
    ax.set_xlabel('Number of amino acids')
    ax.set_ylabel(ylab or label(op))
    if title:
        ax.set_title(title, fontsize=FS + 0.5)
    cb = plt.colorbar(pm, cax=cax) if cax is not None else plt.colorbar(pm, ax=ax, pad=0.03, fraction=0.06)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(width=0.5, size=2)
    return ax


def letters(axes):
    for ax, l in zip(axes, 'abcdefgh'):
        ax.text(-0.27, 1.06, l, transform=ax.transAxes, fontsize=FS + 3, fontweight='bold', va='bottom', ha='left')


def save(fig, path):
    fig.savefig(path, dpi=300)
    plt.close(fig)
    print('wrote', path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root'); ap.add_argument('prefix'); ap.add_argument('--fig4b')
    ap.add_argument('--field', default=FIELD, help='database field folder under out/simulations '
                    '(HLG_wB97XD_TZVP_pH_7_eV or S1_TDDFT_eV)')
    ap.add_argument('--s6-range', type=float, nargs=2, help='y range of Fig. S6 (P(E), posterior); default 0.15 0.15')
    a = ap.parse_args()
    runs = {m: load(a.root, m, a.field) for m in ('gini', 'wmean', 'wvar')}

    # Fig. S4: P(E|A) and P(E|B)
    fig, axs = plt.subplots(3, 2, figsize=(7.0, 8.4), constrained_layout=True)
    for r, m in enumerate(('gini', 'wmean', 'wvar')):
        op, sim = runs[m]
        heat(axs[r, 0], op, sim.P_of_E_given_A.cdf, fit_top(op, sim.P_of_E_given_A.cdf, op.metric_range_abiotic),
             title='P(E|A) (CDF)' if r == 0 else None)
        heat(axs[r, 1], op, sim.P_of_E_given_B.cdf, op.metric_range_biotic, title='P(E|B) (CDF)' if r == 0 else None)
    letters(axs.flat)
    save(fig, f'{a.prefix}PofE.pdf')

    # Fig. S5: P(E) and P(B|E) at prior 0.5
    fig, axs = plt.subplots(3, 2, figsize=(7.0, 8.4), constrained_layout=True)
    for r, m in enumerate(('gini', 'wmean', 'wvar')):
        op, sim = runs[m]
        k = int(np.flatnonzero(np.atleast_1d(sim.P_of_B) == 0.5)[0])
        pe, pbe = np.atleast_1d(sim.P_of_E)[k], np.atleast_1d(sim.P_of_B_given_E)[k]
        lev = [0.5, 0.9, 0.99] if m == 'wvar' else LEVELS   # S5e: 0.01, 0.1 and 0.5 nearly coincide; keep 0.5
        heat(axs[r, 0], op, pe.cdf, op.metric_range_evidence, title='P(E) (CDF), P(B) = 0.5' if r == 0 else None,
             levels=lev, interior_labels=True)
        heat(axs[r, 1], op, pbe.pdf, op.metric_range_posterior, contours=False,
             title='Posterior P(B|E), P(B) = 0.5' if r == 0 else None)
    letters(axs.flat)
    save(fig, f'{a.prefix}EqualPrior.pdf')

    # Fig. S6: effect of the prior (wvar)
    op, sim = runs['wvar']
    ev, post = a.s6_range or (0.15, 0.15)
    fig, axs = plt.subplots(3, 2, figsize=(7.0, 8.4), constrained_layout=True)
    for r, pb in enumerate((0.1, 0.01, 0.001)):
        k = int(np.flatnonzero(np.isclose(np.atleast_1d(sim.P_of_B), pb))[0])
        pe, pbe = np.atleast_1d(sim.P_of_E)[k], np.atleast_1d(sim.P_of_B_given_E)[k]
        # same y-range as the posterior beside it; at P(B) = 0.1 the 0.99 contour lies near 0.8-1.0 eV2, off this scale
        heat(axs[r, 0], op, pe.cdf, (0, ev), title='P(E) (CDF)' if r == 0 else None, interior_labels=True,
             ylab=f'P(B) = {pb:g}\n{label(op)}')
        heat(axs[r, 1], op, pbe.pdf, (0, post), contours=False, title='Posterior P(B|E)' if r == 0 else None)
    letters(axs.flat)
    save(fig, f'{a.prefix}EffectPrior.pdf')

    # Fig. 4b: posterior at prior 0.01 (a single panel of Fig. 4, so no letter)
    if a.fig4b:
        k = int(np.flatnonzero(np.isclose(np.atleast_1d(sim.P_of_B), 0.01))[0])
        fig, ax = plt.subplots(figsize=(3.4, 3.0), constrained_layout=True)
        heat(ax, op, np.atleast_1d(sim.P_of_B_given_E)[k].pdf, (0, post), contours=False)
        save(fig, a.fig4b)


if __name__ == '__main__':
    main()
