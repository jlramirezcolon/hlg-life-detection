"""Supplementary Fig. S8: the optical gap (S1) in the LUMOS framework, weighted variance only.
Combines the wvar rows of the former S1 P(E|A)/P(E|B) and prior-0.5 figures into one row of three panels:
  a  P(E|A) (CDF)     b  P(E|B) (CDF)     c  posterior P(B|E) at P(B) = 0.5
Same data, colour map and contours as sim_figures.py.

Usage: python fig_s1_lumos.py <run folder> <output.pdf> [--field FIELD]
  e.g. python fig_s1_lumos.py .. ../out/SI.S1LUMOS.pdf --field S1_TDDFT_eV
  (needs simulations.m run with descriptor = 'S1', which uses the S1 bin edges in settings/sim_settings_S1.xlsx)
"""
import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sim_figures import FS, FIELD, heat, letters, load, save  # noqa: E402

YLAB = 'Weighted variance of S1 (eV$^2$)'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root'); ap.add_argument('out')
    ap.add_argument('--field', default=FIELD)
    a = ap.parse_args()
    op, sim = load(a.root, 'wvar', a.field)
    k = int(np.flatnonzero(np.atleast_1d(sim.P_of_B) == 0.5)[0])
    pbe = np.atleast_1d(sim.P_of_B_given_E)[k]

    # three square panels, each with a colour bar of the same height (positions in inches)
    W, H, side, bottom, cbw = 7.2, 2.3, 1.3, 0.55, 0.07
    fig = plt.figure(figsize=(W, H))
    lefts = (0.68, 3.08, 5.48)   # each block: y labels, square, colour bar and its ticks (2.4 in)
    axs = [fig.add_axes([l / W, bottom / H, side / W, side / H]) for l in lefts]
    caxs = [fig.add_axes([(l + side + 0.07) / W, bottom / H, cbw / W, side / H]) for l in lefts]
    heat(axs[0], op, sim.P_of_E_given_A.cdf, op.metric_range_abiotic, title='P(E|A) (CDF)', ylab=YLAB,
         interior_labels=True, cax=caxs[0])
    heat(axs[1], op, sim.P_of_E_given_B.cdf, op.metric_range_biotic, title='P(E|B) (CDF)', ylab=YLAB,
         interior_labels=True, cax=caxs[1], levels=[0.1, 0.5, 0.9, 0.99])   # 0.01 runs into 0.1 at the bottom
    heat(axs[2], op, pbe.pdf, op.metric_range_posterior, contours=False, title='P(B|E), P(B) = 0.5', ylab=YLAB,
         cax=caxs[2])
    for l, lab in zip(lefts, 'abc'):   # letters at the left edge of each block, level with the titles
        fig.text((l - 0.66) / W, (bottom + side + 0.13) / H, lab, fontsize=FS + 3, fontweight='bold', ha='left', va='bottom')
    save(fig, a.out)


if __name__ == '__main__':
    main()
