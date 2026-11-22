"""Figures 2-4 and Supplementary Fig. S3 (ClassificationMethods) from the MATLAB outputs (HLG-only run: 10 descriptors,
30 features). Schematic panels are vector artwork (artwork/); data panels are drawn in Helvetica Neue Light with
bold panel letters. Needs macOS (the Swift helpers place the artwork with Core Graphics).

  Figure 2: a from the vector artwork PANEL_2A; b redrawn (Gini, weighted mean, weighted variance of the ensemble HLG).
  Figure 3: a, b, c redrawn (30 features; HLG vs MA index wvar; ROC curves).
  Figure 4: a and c from the vector artwork (PANEL_4A, PANEL_4C); b redrawn (posterior at prior 0.01, n = 3-35); in c, the "Robust" text now says 91%.
  Suppl. ClassificationMethods: a-d redrawn (rank heatmap, AUC, chi-square, MRMR).

Usage: python main_figures.py <repository root> <output dir>
Needs out/rankings_chi2_mrmr.csv and, for Fig. 4b, out/simulations/ from simulations.m (HLG).
"""
import csv
import math
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sim_figures as SF  # noqa: E402  (sets the Helvetica Neue Light font)

LN2 = math.log(2)
PANEL_4A = Path(__file__).resolve().parent / 'artwork' / 'Fig4a.pdf'
PANEL_4C = Path(__file__).resolve().parent / 'artwork' / 'Fig4c.pdf'
PANEL_2A = Path(__file__).resolve().parent / 'artwork' / 'Fig2a.pdf'
DESC = {'HLG_wB97XD_TZVP_pH_7_eV': 'HOMO-LUMO Gap', 'MA_index': 'MA index', 'TopoPSA_A2': 'TPSA',
        'Gibbs_Energy_kJ_mol': 'Gibbs free energy', 'Heat_Formation_kJ_mol': 'Heat of formation',
        'Molar_mass_g_mol': 'Molar mass', 'Carbon_number': 'Carbon number', 'LogP': 'LogP', 'LogS': 'LogS',
        'Dipole_moment_Debye': 'Dipole moment'}
SHORT = dict(DESC, HLG_wB97XD_TZVP_pH_7_eV='HLG')
MCOL = {'wvar': '#d45f69', 'gini': '#f0686b', 'wmean': '#f8bfb3'}
MNAME = {'wvar': 'Weighted variance', 'gini': 'Gini coefficient', 'wmean': 'Weighted mean'}
ABIO, BIO = '#f5a142', '#b9dca6'
ABIO_E, BIO_E = '#c97a1e', '#6fa45a'


def split(f):
    d, m = f.rsplit('_', 1)
    return d, m


def style(fs):
    plt.rcParams.update({'font.size': fs, 'axes.linewidth': 0.8, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
                         'xtick.major.size': 4, 'ytick.major.size': 4, 'legend.frameon': True,
                         'legend.edgecolor': '#999999', 'legend.fancybox': False})


class Page:
    """A transparent page the size of a composite; axes and text placed in pixels of a 1600-px-wide render."""

    def __init__(self, w, h, ref=1600):
        self.w, self.h, self.px = w, h, w / ref
        self.fig = plt.figure(figsize=(w / 72, h / 72))
        self.fig.patch.set_alpha(0)

    def ax(self, x0, y0, x1, y1, **kw):
        p = self.px
        return self.fig.add_axes([x0 * p / self.w, 1 - y1 * p / self.h, (x1 - x0) * p / self.w, (y1 - y0) * p / self.h], **kw)

    def letter(self, l, x, y, fs):
        self.fig.text(x * self.px / self.w, 1 - y * self.px / self.h, l, fontsize=fs, fontweight='bold',
                      va='center', ha='left')

    def rect_pt(self, x0, y0, x1, y1):
        """pixel box -> 'x,y,w,h' in PDF points (origin bottom-left) for merge_keep.swift"""
        p = self.px
        return f'{x0 * p:.2f},{self.h - y1 * p:.2f},{(x1 - x0) * p:.2f},{(y1 - y0) * p:.2f}'

    def fit(self, art, x0, y0, x1, y1):
        """Queue vector artwork scaled into a pixel box (aspect kept, top-left aligned). Returns (scale, x_pt, top_pt)
        so overlays can be positioned in the artwork's own coordinates."""
        r = subprocess.run(['swift', str(HERE / 'pdf_size.swift'), str(art)], capture_output=True, text=True, check=True)
        aw, ah = map(float, r.stdout.split())
        p = self.px
        sc = min((x1 - x0) * p / aw, (y1 - y0) * p / ah)
        self.placed = getattr(self, 'placed', []) + [(art, x0 * p, y0 * p, ah, sc)]
        return sc, x0 * p, y0 * p

    def save(self, out, original=None, keep=(), place=None, svg=False):
        """place = (artwork.pdf, x0, y0, x1, y1), or artwork queued with fit(): vector artwork instead of clipping the
        original. svg=True also writes out.svg with live text (artwork left out; a comment gives its position)."""
        if place:
            self.fit(*place)
        ov = Path(out).with_suffix('.overlay.pdf')
        self.fig.savefig(ov, transparent=True)
        notes = []
        if getattr(self, 'placed', None):
            specs = [f'{art}@{x:.2f},{self.h - top - ah * sc:.2f},{sc:.4f}' for art, x, top, ah, sc in self.placed]
            subprocess.run(['swift', str(HERE / 'place_pdf.swift'), str(ov), str(out), *specs], check=True)
            notes = [f'<!-- Place "{Path(art).name}" with its top-left corner at x = {x:.1f} pt, y = {top:.1f} pt '
                     f'from the top-left of the artboard, scaled to {sc * 100:.1f}% -->' for art, x, top, ah, sc in self.placed]
        else:
            subprocess.run(['swift', str(HERE / 'merge_keep.swift'), str(original) if original else '-', str(ov),
                            str(out), *keep], check=True)
        if svg:
            with plt.rc_context({'svg.fonttype': 'none'}):
                self.fig.savefig(Path(out).with_suffix('.svg'), format='svg', transparent=True)
            if notes:
                t = Path(out).with_suffix('.svg').read_text()
                k = t.index('>', t.index('<svg')) + 1
                Path(out).with_suffix('.svg').write_text(t[:k] + ''.join(notes) + t[k:])
        plt.close(self.fig)
        ov.unlink()
        print('wrote', out)


def load(repo):
    rows = list(csv.DictReader(open(repo / 'out' / 'features_per_sample.csv')))
    cls = np.array([int(float(r['Class'])) for r in rows])
    feats = [k for k in rows[0] if k not in ('Samples', 'Class', 'Class_Label', '_sum')]
    X = {k: np.array([float(r[k]) for r in rows]) for k in feats}
    met = {r['feature']: r for r in csv.DictReader(open(repo / 'out' / 'metrics_weighted.csv'))}
    rk = {r['feature']: r for r in csv.DictReader(open(repo / 'out' / 'rankings_chi2_mrmr.csv'))}
    return cls, feats, X, met, rk


def logbins(v, n=40):
    v = v[v > 0]
    return np.logspace(np.log10(v.min()), np.log10(v.max()), n)


# ---------------------------------------------------------------- Figure 2
def figure2(repo, orig, out):
    cls, feats, X, met, rk = load(repo)
    style(16)
    P = Page(1061.13, 505.918)
    cols = [(70, 505), (620, 1052), (1167, 1597)]
    spec = [('gini', 'Gini coefficient (HLG)', True), ('wmean', 'Weighted mean (HLG), eV', False),
            ('wvar', 'Weighted variance (HLG), eV$^2$', True)]
    for (x0, x1), (m, xl, lg) in zip(cols, spec):
        v = X[f'HLG_wB97XD_TZVP_pH_7_eV_{m}']
        a, b = v[cls == 0], v[cls == 1]
        bins = logbins(v) if lg else np.linspace(v.min(), v.max(), 40)
        lo, hi = bins[0], bins[-1]
        h = P.ax(x0, 422, x1, 685)
        h.hist(a, bins=bins, color=ABIO, alpha=0.8, edgecolor=ABIO_E, lw=0.4, label='Abiotic')
        h.hist(b, bins=bins, color=BIO, alpha=0.8, edgecolor=BIO_E, lw=0.4, label='Biotic')
        h.set_xscale('log' if lg else 'linear')
        h.set_xlim(lo, hi)
        h.set_xlabel(xl)
        h.set_ylabel('Number of Samples')
        h.legend(loc='best', fontsize=13, handlelength=0.9, handleheight=0.9, handletextpad=0.4, borderpad=0.3,
                 labelspacing=0.2, borderaxespad=0.3, framealpha=0.9, edgecolor='#cccccc')
        h.grid(color='#eeeeee', lw=0.5)
        h.set_axisbelow(True)   # grid behind the bars and boxes
        s = P.ax(x0, 300, x1, 370, sharex=h)
        bp = s.boxplot([a, b], vert=False, widths=0.55, patch_artist=True, whis=1.5,
                       flierprops=dict(marker='.', markersize=3))
        for i, (fc, ec) in enumerate(((ABIO, ABIO_E), (BIO, BIO_E))):
            bp['boxes'][i].set(facecolor=fc, edgecolor='#333333', lw=0.8)
            bp['medians'][i].set(color='#333333', lw=0.8)
            bp['fliers'][i].set(markeredgecolor=ec, markerfacecolor=ec)
            for k in (2 * i, 2 * i + 1):
                bp['whiskers'][k].set(color='#333333', lw=0.8, ls=(0, (4, 2)))
                bp['caps'][k].set(color='#333333', lw=0.8)
        s.set_yticks([]); s.tick_params(labelbottom=False)
        s.set_ylim(0.4, 2.6)
        s.grid(axis='x', color='#eeeeee', lw=0.5)
        s.set_axisbelow(True)   # grid behind the bars and boxes
        s.tick_params(which='both', length=0)
        for vals, y, col, va in ((b, 2.95, BIO_E, 'bottom'), (a, -0.0, ABIO_E, 'top')):
            s.plot([vals.min(), vals.max()], [y, y], color='#333333', lw=0.9, clip_on=False)
            for xx in (vals.min(), vals.max()):
                s.plot([xx, xx], [y - 0.12, y + 0.12], color='#333333', lw=0.9, clip_on=False)
            mid = math.sqrt(vals.min() * vals.max()) if lg and vals.min() > 0 else (vals.min() + vals.max()) / 2
            s.text(mid, y + (0.2 if va == 'bottom' else -0.2), f'{vals.max() - vals.min():.3f}', color=col,
                   ha='center', va=va, fontsize=16, clip_on=False)
    P.letter('b', 0, 270, 25)
    # the panel-a artwork carries its own letter 'a', so it is placed flush left and no letter is drawn here
    P.save(out, place=(PANEL_2A, 0, 10, 1600, 240), svg=True)


# ---------------------------------------------------------------- Figure 3
def roc(score, y):
    o = np.argsort(-score, kind='mergesort')
    s, yy = score[o], y[o]
    tp = np.cumsum(yy); fp = np.cumsum(1 - yy)
    keep = np.r_[np.diff(s) != 0, True]
    tpr = np.r_[0, tp[keep] / yy.sum()]; fpr = np.r_[0, fp[keep] / (1 - yy).sum()]
    return fpr, tpr


def figure3(repo, out):
    cls, feats, X, met, rk = load(repo)
    style(16)
    P = Page(1085.5, 935.326)
    # a: relative entropy of all 30 features
    re = {f: float(met[f]['relent_nats']) / LN2 for f in met}
    order = sorted(re, key=re.get)
    a = P.ax(290, 18, 790, 1300)   # wider, closer to the panel letter
    y = np.arange(len(order))
    a.barh(y, [re[f] for f in order], color=[MCOL[split(f)[1]] for f in order], edgecolor='#555555', lw=0.4,
           height=0.75)
    a.set_yticks(y); a.set_yticklabels([DESC[split(f)[0]] for f in order])
    a.set_xscale('log'); a.set_ylim(-0.6, len(order) - 0.4)
    a.set_xlabel('Relative entropy (bits)'); a.set_ylabel('Relative entropy rank', labelpad=12)
    a.tick_params(axis='y', length=0)
    a.legend([plt.Rectangle((0, 0), 1, 1, fc=MCOL[m], ec='#555555', lw=0.4) for m in ('wvar', 'gini', 'wmean')],
             [MNAME[m] for m in ('wvar', 'gini', 'wmean')], loc='lower right', fontsize=13, handlelength=1.0, handleheight=0.9, handletextpad=0.35, borderpad=0.3, labelspacing=0.25, borderaxespad=0.4)
    # b: HLG wvar vs MA index wvar, marginal histograms
    hx, hy = X['HLG_wB97XD_TZVP_pH_7_eV_wvar'], X['MA_index_wvar']
    b = P.ax(955, 180, 1450, 668)
    for c, col, ec, mk, lab in ((0, ABIO, ABIO_E, 'o', 'Abiotic'), (1, BIO, BIO_E, 's', 'Biotic')):
        b.scatter(hx[cls == c], hy[cls == c], s=85 if c == 0 else 50, marker=mk, color=col, edgecolor='black', lw=0.3,
                  alpha=0.9, label=lab, zorder=2 + c)   # biotic on top, as in the marginal histograms
    b.set_xscale('log'); b.set_yscale('log')
    b.set_xlabel('Weighted variance (HLG, eV$^2$)'); b.set_ylabel('Weighted variance (MA index)')
    b.legend(loc='upper left', fontsize=14, markerscale=0.9, handlelength=1.0, handleheight=0.9, handletextpad=0.35, borderpad=0.3, labelspacing=0.25, borderaxespad=0.4)
    xlo, xhi = hx[hx > 0].min() / 1.3, hx.max() * 1.3
    ylo, yhi = hy[hy > 0].min() / 1.3, hy.max() * 1.3
    bx, by = np.logspace(np.log10(xlo), np.log10(xhi), 40), np.logspace(np.log10(ylo), np.log10(yhi), 40)
    b.set_xlim(xlo, xhi); b.set_ylim(ylo, yhi)
    t = P.ax(955, 32, 1450, 172, sharex=b)          # 140 px thick, like the side histogram
    r = P.ax(1460, 180, 1600, 668, sharey=b)
    for c, col, ec in ((0, ABIO, ABIO_E), (1, BIO, BIO_E)):
        nx, _ = np.histogram(hx[cls == c], bins=bx)
        ny, _ = np.histogram(hy[cls == c], bins=by)
        kx, ky = nx > 0, ny > 0                       # empty bins are not drawn (no zero-height outlines)
        t.bar(bx[:-1][kx], nx[kx], width=np.diff(bx)[kx], align='edge', color=col, alpha=0.9, edgecolor='black', lw=0.25)
        r.barh(by[:-1][ky], ny[ky], height=np.diff(by)[ky], align='edge', color=col, alpha=0.9, edgecolor='black', lw=0.25)
    for ax, base in ((t, 'bottom'), (r, 'left')):
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.spines[base].set_visible(True); ax.spines[base].set_color('black'); ax.spines[base].set_linewidth(0.8)
        ax.tick_params(which='both', length=0, labelbottom=False, labelleft=False)   # shared axes: hide, don't clear
    cmax = 1.05 * max(t.get_ylim()[1], r.get_xlim()[1])   # one count scale for both marginals
    t.set_ylim(0, cmax); r.set_xlim(0, cmax)
    # c: ROC curves of the wvar of four descriptors
    c = P.ax(955, 810, 1450, 1300)
    styles = [('HLG_wB97XD_TZVP_pH_7_eV', 'HOMO-LUMO Gap', '#1aa3d9', (0, (4, 1.5))),
              ('MA_index', 'MA index', '#1b8ad0', (0, (2, 1))), ('Molar_mass_g_mol', 'Molar mass', '#2160a8', (0, (6, 2, 1, 2))),
              ('Carbon_number', 'Carbon number', '#1d3b70', (0, (1, 1)))]
    for d, lab, col, ls in styles:
        f = f'{d}_wvar'
        sgn = 1 if float(met[f]['direction']) > 0 else -1
        fpr, tpr = roc(sgn * X[f], cls)
        c.plot(fpr, tpr, color=col, lw=1.6, ls=ls, label=lab,
               drawstyle='steps-post')
    c.plot([0, 0, 1], [0, 1, 1], color='#bcd9ea', lw=1.6, label='Perfect classifier', zorder=10)   # above curves and frame
    c.plot([0, 1], [0, 1], color='#333333', lw=1.2, ls=(0, (1, 1.5)), label='Random classifier')
    c.set_xlim(0, 1); c.set_ylim(0, 1)
    for ln in c.get_lines():
        ln.set_clip_on(False)
    c.set_xlabel('False positive rate'); c.set_ylabel('True positive rate')
    c.grid(color='#eeeeee', lw=0.5)
    c.legend(loc='lower right', fontsize=13, handlelength=1.8, handletextpad=0.4, borderpad=0.3, labelspacing=0.25, borderaxespad=0.4)
    for l, (x, yy) in {'a': (0, 20), 'b': (875, 20), 'c': (875, 785)}.items():
        P.letter(l, x, yy, 25)
    P.save(out, svg=True)


# ---------------------------------------------------------------- Figure 4
def figure4(repo_sim, orig, out):
    style(16)
    P = Page(1228.26, 563.405)
    op, sim = SF.load(repo_sim, 'wvar')
    k = int(np.flatnonzero(np.isclose(np.atleast_1d(sim.P_of_B), 0.01))[0])
    z = np.asarray(np.atleast_1d(sim.P_of_B_given_E)[k].pdf, float)
    x, yc = SF.grid(op)
    keep = yc <= 0.15 + (yc[1] - yc[0])
    ax = P.ax(407, 118, 985, 612)
    pm = ax.pcolormesh(x, yc[keep], np.ma.masked_invalid(z[:, keep].T), cmap=SF.CMAP, vmin=0, vmax=1,
                       shading='nearest', rasterized=True, linewidth=0, antialiased=False)
    ax.set_xlim(x.min(), x.max()); ax.set_ylim(0, 0.15)
    ax.set_xlabel('Number of amino acids'); ax.set_ylabel('Weighted variance of HLG (eV$^2$)')
    cax = P.ax(1012, 118, 1040, 612)
    cb = plt.colorbar(pm, cax=cax)
    cb.outline.set_linewidth(0.6)
    cb.set_label('Confidence of biogenicity, P(B) = 0.01', rotation=270, va='bottom', labelpad=8)   # reads top to bottom
    # panels a and c: vector artwork (PANEL_4A, PANEL_4C)
    P.fit(PANEL_4A, 40, 40, 275, 725)
    sc, cx, ctop = P.fit(PANEL_4C, 1150, 40, 1600, 725)
    # panel c: the "Robust" box states the classifier accuracy (91%); drawn over the artwork's box, in its coordinates
    # (box x 77.6-315.0 pt, y 366.1-420.1 pt from the artwork's top; lines centred at 377.1, 393.9, 410.4 pt)
    pt = lambda ax_, ay_: ((cx + sc * ax_) / P.px, (ctop + sc * ay_) / P.px)
    (bx0, by0), (bx1, by1) = pt(77.6, 366.1), pt(315.0, 420.1)
    box = P.ax(bx0, by0, bx1, by1)
    box.set_xlim(0, 1); box.set_ylim(0, 1); box.axis('off')
    box.add_patch(plt.Rectangle((0, 0), 1, 1, color=(245 / 255, 238 / 255, 222 / 255), transform=box.transAxes))
    for line, yl in zip(('Classifies biotic and abiotic origins', 'with 91% accuracy using varied', 'analytical techniques'),
                        (377.1, 393.9, 410.4)):
        box.text((79.7 - 77.6) / (315.0 - 77.6), 1 - (yl - 366.1) / (420.1 - 366.1), line, fontsize=14.6 * sc,
                 color='#1a1a1a', va='center', ha='left', transform=box.transAxes)
    P.letter('a', 0, 40, 25)
    P.letter('b', 303, 40, 25)
    P.letter('c', 1112, 40, 25)
    P.save(out, svg=True)


# ---------------------------------------------------------------- Suppl. ClassificationMethods
def classification(repo, out):
    cls, feats, X, met, rk = load(repo)
    style(14)
    P = Page(996.669, 1028)            # shorter page so the figure plus caption fits on one page at text width
    f30 = list(met)
    scores = {'Relative Entropy': {f: float(met[f]['relent_nats']) for f in f30},
              'AUC': {f: float(met[f]['auc']) for f in f30},
              'Chi-Square': {f: float(rk[f]['chi2_score']) for f in f30},
              'MRMR': {f: float(rk[f]['mrmr_score']) for f in f30}}
    ranks = {m: {f: i + 1 for i, f in enumerate(sorted(f30, key=lambda f: -s[f]))} for m, s in scores.items()}
    order = sorted(f30, key=lambda f: np.mean([ranks[m][f] for m in ranks]))
    lab = lambda f: f'{split(f)[1]} ({SHORT[split(f)[0]]})'
    # a: rank heatmap
    # all four panels share one size: 460 x 690 px (row 1: y 30-720, row 2: y 880-1570); 690/30 rows keeps 12.5-pt labels apart
    a = P.ax(272, 30, 732, 720)        # left column moved toward its letters
    R = np.array([[ranks[m][f] for m in ranks] for f in order])
    im = a.imshow(R, cmap='bone_r', vmin=1, vmax=30, aspect='auto')
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            a.text(j, i, R[i, j], ha='center', va='center', fontsize=11.5, color='white' if R[i, j] > 15 else 'black')
    a.set_yticks(range(len(order))); a.set_yticklabels([lab(f) for f in order], fontsize=12.5)
    a.set_xticks(range(4)); a.set_xticklabels(list(ranks), rotation=35, ha='right')
    a.tick_params(length=0)
    for xv in (0.5, 1.5, 2.5):                       # thin lines between the method columns
        a.axvline(xv, color='black', lw=0.4)
    for yv in np.arange(0.5, len(order) - 1, 1.0):   # and between the feature rows
        a.axhline(yv, color='black', lw=0.4)
    a.set_ylabel('Features', labelpad=4)
    cb = plt.colorbar(im, cax=P.ax(747, 30, 767, 720)); cb.ax.invert_yaxis()
    cb.outline.set_linewidth(0.6)
    cb.set_label('Rank', rotation=270, va='bottom', labelpad=4)
    # the bar is inverted (rank 1, light, at the top); axes-fraction positions are not, so 1 = top
    cb.ax.text(0.5, 1.008, 'Better', transform=cb.ax.transAxes, ha='center', va='bottom', fontsize=12)
    cb.ax.text(0.5, -0.03, 'Worse', transform=cb.ax.transAxes, ha='center', va='top', fontsize=12)

    def bars(ax, m, xl, title):
        s = scores[m]
        o = sorted(f30, key=lambda f: s[f])
        yy = np.arange(len(o))
        ax.barh(yy, [s[f] for f in o], color=[MCOL[split(f)[1]] for f in o], edgecolor='#555555', lw=0.4, height=0.75)
        ax.set_yticks(yy); ax.set_yticklabels([SHORT[split(f)[0]] for f in o], fontsize=12.5)
        ax.set_ylim(-0.6, len(o) - 0.4); ax.tick_params(axis='y', length=0)
        ax.set_xlabel(xl)
        ax.set_ylabel(f'{m if m != "Chi-Square" else "Chi-squared"} rank')
        ax.legend([plt.Rectangle((0, 0), 1, 1, fc=MCOL[k], ec='#555555', lw=0.4) for k in ('wvar', 'gini', 'wmean')],
                  ['wvar', 'gini', 'wmean'], loc='lower right', fontsize=12, handlelength=1.0, handleheight=0.9,
                  handletextpad=0.35, borderpad=0.3, labelspacing=0.2, borderaxespad=0.4)

    bars(P.ax(1100, 30, 1560, 720), 'AUC', 'Area under the ROC curve', 'AUC method')
    bars(P.ax(272, 880, 732, 1570), 'Chi-Square', 'Predictor importance score', 'Chi-squared method')
    bars(P.ax(1100, 880, 1560, 1570), 'MRMR', 'Predictor importance score', 'MRMR method')
    for l, (x, yy) in {'a': (0, 15), 'b': (870, 15), 'c': (0, 862), 'd': (870, 862)}.items():
        P.letter(l, x, yy, 22)
    P.save(out, svg=True)


if __name__ == '__main__':
    repo, outdir = (Path(p) for p in sys.argv[1:3])
    figure2(repo, None, outdir / 'Figure2.Features.pdf')
    figure3(repo, outdir / 'Figure3.Separation.pdf')
    figure4(repo, None, outdir / 'Figure4.Framework.pdf')
    classification(repo, outdir / 'SI.ClassificationMethods.pdf')
