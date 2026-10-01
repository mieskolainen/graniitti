#!/usr/bin/env python3
# Plot self-describing Sudakov and Shuvaev interpolation arrays
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from core.io.serialize import decode_doubles, load_json_file

SCRIPT_DIR = Path(__file__).resolve().parent
LABELS = {'q2': r'$Q^2$ (GeV$^2$)', 'mu': r'$\mu$ (GeV)', 'x': r'$x$'}


# Read grid coordinates and undo only the transformations recorded by the generator
def load_array(path):
    meta = load_json_file(path)
    if meta.get('type') != 'IArray2D' or meta.get('version') != 1:
        raise ValueError(f'{path}: unsupported interpolation cache')
    data = decode_doubles(meta['data'], 'data').reshape(meta['shape'])
    if data.ndim != 3 or data.shape[-1] != 4:
        raise ValueError(f'{path}: expected four columns per grid cell')
    data = data.reshape(-1, 4)
    axes = [np.unique(data[:, i]) for i in range(2)]
    grid = data.reshape(len(axes[0]), len(axes[1]), 4)
    if not (np.allclose(grid[:, :, 0], axes[0][:, None])
            and np.allclose(grid[:, :, 1], axes[1][None, :])):
        raise ValueError(f'{path}: coordinates do not form an ordered rectangular grid')
    if meta['axes'] not in (['q2', 'mu'], ['q2', 'x']):
        raise ValueError(f'{path}: unsupported physical axes {meta["axes"]}')
    if len(meta['log']) != 2 or any(type(flag) is not bool for flag in meta['log']):
        raise ValueError(f'{path}: invalid coordinate transformations')
    axes = [np.exp(axis) if log else axis for axis, log in zip(axes, meta['log'], strict=True)]
    if any(not np.all(np.isfinite(axis)) for axis in axes):
        raise ValueError(f'{path}: nonfinite physical coordinates')
    values = grid[:, :, 2].T
    if meta.get('quantity') == 'sudakov_radiator':
        log_range = np.maximum(0.0, np.log(axes[1][:, None] ** 2 / axes[0][None, :]))
        values = np.exp(-values * log_range ** 2)
    return meta, axes, values


# Label physical axes and select their scales from the stored grid transformations
def label_axes(ax, meta):
    ax.set_xlabel(LABELS[meta['axes'][0]])
    ax.set_ylabel(LABELS[meta['axes'][1]])
    ax.set_xscale('log' if meta['log'][0] else 'linear')
    ax.set_yscale('log' if meta['log'][1] else 'linear')
    ax.set_title(rf'$\sqrt{{s}} = {meta["sqrts"]:g}$ GeV, {meta["pdf"]}')


# Plot the full physical domain and slices selected across the stored grid
def plot_array(path, fig_dir, slices):
    meta, axes, values = load_array(path)
    output = fig_dir / path.name
    output.mkdir(parents=True, exist_ok=True)
    symbol = 'T' if meta['axes'][1] == 'mu' else 'H'
    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(*axes, values, shading='auto')
    label_axes(ax, meta)
    fig.colorbar(mesh, ax=ax, label=rf'${symbol}$')
    fig.savefig(output / 'map.pdf', bbox_inches='tight')
    plt.close(fig)
    for fixed in range(2):
        varied = 1 - fixed
        fig, ax = plt.subplots()
        for index in np.unique(np.linspace(0, len(axes[fixed]) - 1, slices, dtype=int)):
            curve = values[:, index] if fixed == 0 else values[index, :]
            ax.plot(axes[varied], curve, label=f'{meta["axes"][fixed]} = {axes[fixed][index]:.4g}')
        ax.set_xlabel(LABELS[meta['axes'][varied]])
        ax.set_ylabel(rf'${symbol}$')
        ax.set_xscale('log' if meta['log'][varied] else 'linear')
        if np.all(values > 0):
            ax.set_yscale('log')
        ax.legend()
        ax.set_title(rf'$\sqrt{{s}} = {meta["sqrts"]:g}$ GeV, {meta["pdf"]}')
        fig.savefig(output / f'{meta["axes"][varied]}.pdf', bbox_inches='tight')
        plt.close(fig)
    print(f'{path}: {values.shape[1]} x {values.shape[0]} grid -> {output}')


# Plot explicitly selected data files without assuming an energy, PDF or grid
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('arrays', nargs='+', type=Path, help='Sudakov or Shuvaev JSON cache files')
    parser.add_argument('--fig-dir', type=Path, default=SCRIPT_DIR / 'figs')
    parser.add_argument('--slices', type=int, default=7, help='Number of curves per coordinate')
    args = parser.parse_args()
    if args.slices < 1:
        parser.error('--slices must be positive')
    for path in args.arrays:
        plot_array(path, args.fig_dir, args.slices)


if __name__ == '__main__':
    main()
