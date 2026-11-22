"""Render the Fig. 1b glycine HOMO/LUMO isosurfaces (cubes from fig1b_orbitals.py) as transparent PNGs.

Same ball-and-stick style as the chapter-2 orbital figures (orbitals_S4-7/render.py), but with the blue and
yellow lobes that the Fig. 1 caption describes (blue = positive phase, yellow = negative). HOMO and LUMO are
drawn from the same viewpoint (molecule on its principal axes).
Usage: ~/orbitals-venv/bin/python fig1b_render.py [--size PX] [--out DIR]   (default 2000 px into fig1b_data/)
"""
from pathlib import Path

import argparse

import numpy as np
import pyvista as pv

HERE = Path(__file__).resolve().parent
BOHR = 0.529177
COLORS = {'H': '#ffffff', 'C': '#8a8a8a', 'N': '#2f55c8', 'O': '#d42a2a'}     # as orbitals_S4-7/render.py
RADII = {'H': 0.22, 'C': 0.36, 'N': 0.36, 'O': 0.36}
COVALENT = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66}


def principal_frame(coords):
    xyz = np.asarray(coords); centre = xyz.mean(axis=0)
    _, _, vt = np.linalg.svd(xyz - centre)
    if np.linalg.det(vt) < 0:
        vt[2] *= -1
    return centre, vt


def read_cube(path):
    lines = Path(path).read_text().splitlines()
    n_atoms, *origin = lines[2].split()
    n_atoms = abs(int(n_atoms)); origin = np.array(origin[:3], float)
    axes, counts = [], []
    for k in range(3):
        n, *v = lines[3 + k].split()
        counts.append(int(n)); axes.append(np.array(v[:3], float))
    vals = np.array(' '.join(lines[6 + n_atoms:]).split(), float).reshape(counts)
    i, j, k = np.meshgrid(*(np.arange(n) for n in counts), indexing='ij')
    pts = origin + i[..., None] * axes[0] + j[..., None] * axes[1] + k[..., None] * axes[2]
    grid = pv.StructuredGrid(*(pts[..., c] * BOHR for c in range(3)))
    grid['psi'] = vals.ravel(order='F')
    return {'grid': grid}

pv.OFF_SCREEN = True
D = HERE / 'fig1b_data'
ISO = 0.05
LOBES = (('#4f8fd1', 1), ('#e8c44a', -1))     # blue positive, yellow negative
ap = argparse.ArgumentParser(); ap.add_argument('--size', type=int, default=2000); ap.add_argument('--out', default=None)
ARGS = ap.parse_args()
SIZE = ARGS.size
OUTDIR = Path(ARGS.out) if ARGS.out else D
OUTDIR.mkdir(parents=True, exist_ok=True)

lines = (HERE.parents[1] / 'ensemble_tzvp' / 'geometries' / 'glycine__zwitterionic__c1.xyz').read_text().strip().splitlines()[2:]
sym = [l.split()[0] for l in lines]
xyz0 = np.array([[float(v) for v in l.split()[1:4]] for l in lines])
centre, rot = principal_frame(xyz0)
xyz = (xyz0 - centre) @ rot.T
for orb in ('homo', 'lumo'):
    grid = read_cube(D / f'glycine.{orb}.cube')['grid']
    grid.points = (grid.points - centre) @ rot.T
    p = pv.Plotter(off_screen=True, window_size=(SIZE, SIZE))
    p.set_background('white')
    p.enable_anti_aliasing('ssaa')
    for s, x in zip(sym, xyz):
        p.add_mesh(pv.Sphere(radius=RADII[s], center=x, theta_resolution=48, phi_resolution=48),
                   color=COLORS[s], smooth_shading=True, ambient=0.35, diffuse=0.7, specular=0.25)
    for i in range(len(sym)):
        for j in range(i + 1, len(sym)):
            if np.linalg.norm(xyz[i] - xyz[j]) < 1.15 * (COVALENT[sym[i]] + COVALENT[sym[j]]):
                p.add_mesh(pv.Tube(pointa=xyz[i], pointb=xyz[j], radius=0.09, n_sides=28),
                           color='#9a9a9a', smooth_shading=True, ambient=0.3)
    for col, sign in LOBES:
        surf = grid.contour([sign * ISO], scalars='psi')
        if surf.n_points:
            p.add_mesh(surf, color=col, opacity=0.55, smooth_shading=True, ambient=0.55, diffuse=0.5, specular=0.1)
    p.camera_position = 'xy'
    p.camera.zoom(0.95)
    p.screenshot(str(OUTDIR / f'glycine.{orb}.png'), transparent_background=True)
    p.close()
    print('wrote', OUTDIR / f'glycine.{orb}.png')
