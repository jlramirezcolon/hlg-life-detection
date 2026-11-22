"""Figure 1: uses the vector panel-a artwork (PANEL_A) and redraws
b (frontier orbitals of zwitterionic glycine; fig1b_orbitals.py + fig1b_render.py), c (ensemble HLG by subcategory)
and d (unweighted relative entropy, 10 descriptors) in the composite's layout.
Panel e (pairwise rank-sum matrix) is dropped, so d takes the whole right column.

Step 1 (this script) draws c and d on a transparent page the size of the original composite.
Step 2 (place_pdf.swift) places the vector panel-a artwork (PANEL_A), then this overlay, into one PDF.

Usage: python fig1_composite.py <repository root> <output.pdf>
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
from hnfont import use_helvetica_neue_light  # noqa: E402

use_helvetica_neue_light()
W, H = 1586.57, 692.24          # original composite page, points
PX = W / 1800                   # the layout below is measured on a 1800-px-wide render of the original
PANEL_A = Path(__file__).resolve().parent / 'artwork' / 'Fig1a.pdf'   # panel a artwork (vector)
HOMO_EV, LUMO_EV = -0.345319 * 27.211386, 0.105802 * 27.211386   # glycine, Psi4 (ensemble_tzvp/results.csv)
FS = 19
plt.rcParams.update({'font.size': FS, 'axes.linewidth': 0.8, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8})

GROUPS = [('Biotic', '#1aab5a', [('Hydrothermal', '#5fd3a0'), ('Organism', '#3cbf86'), ('Terrestrial biome', '#22a56c')]),
          ('Abiotic (Simulated)', '#2a87c8', [('Asteroid Parent Body', '#8ec5ef'), ('Hadean Earth', '#73afe6'),
                                               ('Titan', '#5895dc'), ('Noachian Mars', '#3f7dcc'), ('Europa', '#2563b6')]),
          ('Abiotic', '#f39325', [('Meteorite', '#f8ac5c'), ('Asteroid', '#f19544'), ('Lunar', '#e47f30')])]
LABELS = {'HLG_wB97XD_TZVP_pH_7_eV': 'HOMO-LUMO Gap', 'Carbon_number': 'Carbon number',
          'Gibbs_Energy_kJ_mol': 'Gibbs free energy', 'MA_index': 'MA index', 'Heat_Formation_kJ_mol': 'Heat of Formation',
          'Molar_mass': 'Molar mass', 'LogP': 'LogP', 'Dipole_moment_Debye': 'Dipole moment', 'TopoPSA_A2': 'TPSA',
          'LogS': 'LogS'}


def ax_px(fig, x0, y0, x1, y1):
    """Axes from pixel box (x0, y0 top-left; x1, y1 bottom-right) of the 1800-px layout."""
    return fig.add_axes([x0 * PX / W, 1 - y1 * PX / H, (x1 - x0) * PX / W, (y1 - y0) * PX / H])


def write_svg(out, x_pt, top_pt, scale):
    """out.svg: panels b-d with live text (Helvetica Neue Light) for editing in Illustrator. Panel a is left out
    (SVG import does not reproduce its masked gradient); place PANEL_A at the position recorded in the SVG comment."""
    tmp = out.with_suffix('.svg.tmp')
    m = tmp.read_text()
    k = m.index('>', m.index('<svg')) + 1
    note = (f'<!-- Panel a: place "{PANEL_A.name}" with its top-left corner at x = {x_pt:.1f} pt, y = {top_pt:.1f} pt '
            f'from the top-left of the artboard, scaled to {scale * 100:.1f}% -->')
    out.with_suffix('.svg').write_text(m[:k] + note + m[k:])
    tmp.unlink()


def main(repo, out):
    repo = Path(repo)
    vals = {}
    for r in csv.DictReader(open(repo / 'out' / 'fig1c_values.csv')):
        vals.setdefault(r['subcategory'], []).append(float(r['hlg_eV']))
    samples = {r['Category']: int(r['Samples'])
               for r in csv.DictReader(open(repo / 'out' / 'statistics_table4_HLG_wB97XD_TZVP_pH_7_eV.csv'))}
    re_bits = {r['variable']: float(r['mean']) / math.log(2)
               for r in csv.DictReader(open(repo / 'out' / 'metrics_distribution.csv')) if r['class'].startswith('RE')}

    fig = plt.figure(figsize=(W / 72, H / 72))
    fig.patch.set_alpha(0)

    # ---- c: three stacked box-plot panels, same vertical layout as the original
    boxes = [(37, 218), (231, 535), (543, 727)]
    xlim = (9.0, 12.6)
    drawn = []
    for (gname, gcol, subs), (top, bot) in zip(GROUPS, boxes):
        ax = ax_px(fig, 662, top, 1268, bot)
        data = [vals[s] for s, _ in subs][::-1]
        cols = [c for _, c in subs][::-1]
        bp = ax.boxplot(data, vert=False, widths=0.55, patch_artist=True, whis=1.5,
                        flierprops=dict(marker='o', markersize=5, markerfacecolor='none', markeredgewidth=0.9))
        for i, c in enumerate(cols):
            bp['boxes'][i].set(facecolor=matplotlib.colors.to_rgba(c, 0.75), edgecolor=c, linewidth=1.3)
            bp['medians'][i].set(color=matplotlib.colors.to_rgba('k', 0.0))
            ax.plot([np.median(data[i])] * 2, [i + 1 - 0.275, i + 1 + 0.275], color=c, lw=1.6, zorder=3,
                    solid_capstyle='butt')
            for k in (2 * i, 2 * i + 1):
                bp['whiskers'][k].set(color=c, linewidth=1.3)
                bp['caps'][k].set(color=c, linewidth=1.3)
            bp['fliers'][i].set(markeredgecolor=c)
        ax.set_yticks(range(1, len(subs) + 1))
        ax.set_yticklabels([s for s, _ in subs][::-1])
        ax.set_xlim(*xlim)
        ax.set_ylim(0.4, len(subs) + 0.6)
        ax.tick_params(axis='y', length=0, pad=8)
        ax.grid(axis='x', color='#e6e6e6', lw=0.6)
        ax.set_axisbelow(True)
        if gname != 'Abiotic':
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('HOMO-LUMO Gap (eV)')
        drawn.append((ax, gname, gcol, subs, top, bot))

    # group bars sit just left of the longest subcategory name (measured), so no gap opens between them
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    FigureCanvasAgg(fig)
    rend = fig.canvas.get_renderer()
    left_disp = min(t.get_window_extent(rend).x0 for ax, *_ in drawn for t in ax.get_yticklabels())
    bar_r = left_disp / fig.dpi * 72 / PX - 12
    for ax, gname, gcol, subs, top, bot in drawn:
        n = sum(samples[s] for s, _ in subs)
        # coloured group bar with the group name, and n beside it (as in the original)
        bx = ax_px(fig, bar_r - 28, top + 2, bar_r, bot - 2)
        bx.add_patch(matplotlib.patches.FancyBboxPatch((0.08, 0.01), 0.84, 0.98, boxstyle='round,pad=0,rounding_size=0.3',
                                                       transform=bx.transAxes, color=gcol, mutation_aspect=0.05))
        bx.text(0.5, 0.5, gname, rotation=90, ha='center', va='center', color='white', fontsize=FS + 1,
                transform=bx.transAxes)
        bx.text(-0.75, 0.5, f'n = {n}', rotation=90, ha='center', va='center', fontsize=FS, transform=bx.transAxes)
        bx.axis('off')

    # ---- d: relative entropy, now the whole right column (panel e removed)
    order = sorted(re_bits, key=re_bits.get)
    ax = ax_px(fig, 1545, 37, 1790, 727)   # same top and bottom as panel c
    y = np.arange(len(order))
    ax.barh(y, [re_bits[k] for k in order], color='#f9c22e', edgecolor='#3a3a3a', linewidth=0.6, height=0.72)
    ax.set_yticks(y)
    ax.set_yticklabels([LABELS[k] for k in order])
    ax.set_xscale('log')
    ax.set_xlim(1e-2, 5)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.tick_params(axis='y', length=0, pad=6)
    ax.set_xlabel('Relative entropy (bits)')
    ax.set_ylabel('Molecular descriptor', labelpad=10)
    # boxed, like the panel c axes

    # ---- b: glycine frontier orbitals (zwitterion), energies from the database (Psi4)
    from PIL import Image
    eb = ax_px(fig, 118, 372, 362, 725)
    eb.set_ylim(-16.2, 9.5); eb.set_xlim(0, 1)
    for sname in ('top', 'right', 'bottom'):
        eb.spines[sname].set_visible(False)
    eb.set_xticks([])
    eb.set_yticks([5, 0, -5, -10, -15]); eb.set_yticklabels(['5.0', '0.0', '-5.0', '-10.0', '-15.0'])
    eb.set_ylabel('Energy (eV)')
    eb.plot([0], [9.5], marker='^', color='k', markersize=6, clip_on=False)
    lev = '#f5b400'
    for e, lab, va, dy in ((LUMO_EV, 'LUMO', 'bottom', 0.6), (HOMO_EV, 'HOMO', 'top', -0.6)):
        eb.plot([0.04, 0.30], [e, e], color=lev, lw=2.4, solid_capstyle='butt')
        eb.text(0.17, e + dy, lab, ha='center', va=va)
    eb.plot([0.17, 0.17], [HOMO_EV + 0.4, LUMO_EV - 0.4], color=lev, lw=1.6)
    eb.text(0.20, (HOMO_EV + LUMO_EV) / 2, f'{LUMO_EV - HOMO_EV:.2f} eV', ha='left', va='center')
    for orb, e in (('lumo', LUMO_EV), ('homo', HOMO_EV)):
        im = Image.open(HERE / 'fig1b_data' / f'glycine.{orb}.png')
        im = im.crop(im.getbbox())
        wpx = 142; hpx = wpx * im.size[1] / im.size[0]
        ycen = 434 + (5 - e) * (708 - 434) / 20      # same energy-to-pixel map as the axis above
        ia = ax_px(fig, 218, ycen - hpx / 2, 218 + wpx, ycen + hpx / 2)
        ia.imshow(im); ia.axis('off')

    for letter, (x, y0) in {'a': (0, 20), 'b': (8, 368), 'c': (372, 20), 'd': (1300, 20)}.items():
        fig.text(x * PX / W, 1 - y0 * PX / H, letter, fontsize=FS + 9, fontweight='bold', va='center', ha='left')

    overlay = Path(out).with_suffix('.overlay.pdf')
    fig.savefig(overlay, transparent=True)
    # editable SVG for Illustrator: live text (not outlines), panel a inserted from its own SVG (Inkscape)
    with plt.rc_context({'svg.fonttype': 'none'}):
        fig.savefig(Path(out).with_suffix('.svg.tmp'), format='svg', transparent=True)
    plt.close(fig)
    # panel a: the vector artwork (PANEL_A), scaled into the a box (x 18-355 px, y 15-345 px of the layout)
    aw, ah = 296.308, 289.862
    sc = min((355 - 18) * PX / aw, (345 - 15) * PX / ah)
    place = f'{PANEL_A}@{18 * PX:.2f},{H - 345 * PX:.2f},{sc:.4f}'
    subprocess.run(['swift', str(HERE / 'place_pdf.swift'), str(overlay), out, place], check=True)
    write_svg(Path(out), 18 * PX, 345 * PX - ah * sc, sc)
    overlay.unlink()
    print('wrote', out)
    print('relative entropy (bits):', {LABELS[k]: round(v, 2) for k, v in sorted(re_bits.items(), key=lambda kv: -kv[1])})


if __name__ == '__main__':
    main(*sys.argv[1:3])
