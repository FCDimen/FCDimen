#! /usr/bin/env python

from argparse import ArgumentParser
from fcdimen.functions.reader import read_data
from fcdimen.functions.fcmat import force_matrix
import os
import numpy as np
from fcdimen.functions.welcome import welcome_logo
from fcdimen.functions.score import calc_score, scanner


def _print_pair(label, supercell, forceconstants_zero_diagonal, value):
    """Print the first atom pair (row-major order) whose force constant equals value."""
    print(label)
    rows, cols = np.nonzero(forceconstants_zero_diagonal == value)
    if len(rows) != 0:
        i, j = rows[0], cols[0]
        print("{}, {} = {:.2f} eV/A^2".format(supercell.symbols[i], supercell.symbols[j],
                                              forceconstants_zero_diagonal[i][j]))


def run():
    parser = ArgumentParser(description="Example: fcdimen -p PATH -i phonon.yaml")
    parser.set_defaults(
        input_file="phonon.yaml",
        path=os.getcwd(),
        threshold=None
    )
    parser.add_argument("-i", "--input_file", help="Phonopy yaml file name, default= phonon.yaml")
    parser.add_argument("-p", "--path", help="Path to phonopy yaml file with FORCE_SET file,"
                                             " default= Current directory")
    parser.add_argument("-t", "--threshold", help="Choose a specific threshold")
    args = parser.parse_args()
    if args.path is not None:
        os.chdir(args.path)

    welcome_logo()

    phonon, unitcell, forceconstants, supercell = read_data(filename=args.input_file)
    natom, forceconstants_zero_diagonal, forceconstants_reshaped = force_matrix(forceconstants)
    column_max = forceconstants_zero_diagonal.max(axis=0)
    maximumforce = column_max.max()
    # Generating thresholds for scanning
    pl0 = np.round(np.arange(0, 50, 0.1), 2).tolist()
    pl1 = np.arange(50, maximumforce, 1).tolist()
    pl = pl0 + pl1

    _print_pair("MaxFC between atoms: ", supercell, forceconstants_zero_diagonal, maximumforce)
    _print_pair("MinFC between atoms: ", supercell, forceconstants_zero_diagonal, column_max.min())

    if args.threshold is None:
        ths = [x for x in pl if x <= maximumforce]
        fcdimlist, dimensionality1, dimensionality2 = scanner(natom, forceconstants_zero_diagonal, supercell, ths)
        if dimensionality1 is not None:
            print("Dimensionality with Class. I: {}".format(dimensionality1))

        if dimensionality2 is not None:
            print("Dimensionality with Class. II: {}".format(dimensionality2))

        if fcdimlist is not None:
            scores = calc_score(fcdimlist, maximumforce)
            print("Dimensionality with Class. III: {}".format(max(scores, key=scores.get)))
            print("Scores: {}".format(scores))

    else:
        assert float(args.threshold) >= 0.0, "ERROR: negative Negative threshold"
        ths = [float(args.threshold)]
        fcdimlist, dimensionality1, dimensionality2 = scanner(natom, forceconstants_zero_diagonal, supercell, ths)
        print("Dimensionality for t = {} eV/A^2 is: {}".format(ths[0], dimensionality1))


if __name__ == "__main__":
    run()
