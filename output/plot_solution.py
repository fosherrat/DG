#!/usr/bin/env python3
"""Plot scalar or compressible DG conservative coefficients."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


GAMMA = 1.4


def read_metadata(path):
    metadata = {}
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            break
        if line.startswith("# nvar =") or line.startswith("# ndof ="):
            key, value = line[1:].split("=", 1)
            metadata[key.strip()] = int(value)
    return metadata


def main():
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input", nargs="?", default=script_dir / "solution_coefficients.dat"
    )
    parser.add_argument("-o", "--output", default=script_dir / "solution.png")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    input_path = Path(args.input)
    metadata = read_metadata(input_path)
    ndof = metadata["ndof"]
    nvar = metadata["nvar"]
    data = np.loadtxt(input_path, comments="#")

    if nvar == 1:
        fig, axes = plt.subplots(1, 1, figsize=(8, 3))
        axes = (axes,)
        labels = (r"$\rho$",)
    elif nvar == 3:
        fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 7))
        labels = (r"$\rho$", r"$u$", r"$p$")
    else:
        raise ValueError(f"unsupported number of variables: {nvar}")

    for element_id in np.unique(data[:, 0].astype(int)):
        rows = data[data[:, 0] == element_id]
        rows = rows[np.argsort(rows[:, 1])]
        if len(rows) != ndof:
            raise ValueError(f"element {element_id} has {len(rows)} coefficients, expected {ndof}")

        x_center, jacobian = rows[0, 2:4]
        coefficients = rows[:, 4 : 4 + nvar]
        xi = np.linspace(-1.0, 1.0, 100)
        basis = np.polynomial.legendre.legvander(xi, ndof - 1)
        solution = basis @ coefficients

        x = x_center + jacobian * xi
        rho = solution[:, 0]
        axes[0].plot(x, rho, "k-")
        if nvar == 3:
            velocity = solution[:, 1] / rho
            pressure = (GAMMA - 1.0) * (solution[:, 2] - 0.5 * rho * velocity**2)
            axes[1].plot(x, velocity, "k-")
            axes[2].plot(x, pressure, "k-")

    for axis, label in zip(axes, labels):
        axis.set_ylabel(label)
        axis.grid(True, alpha=0.3)
    axes[-1].set_xlabel("x")

    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"saved {args.output}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
