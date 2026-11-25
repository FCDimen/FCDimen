import numpy as np
import networkx as nx
import pytest
from ase import Atoms

from fcdimen.functions.cluster import doubled_fc, mic_distance_vectors
from fcdimen.functions.score import scanner
from fcdimen.functions.dimensionality import calc_dimension


def reference_doubled_fc(natom, forceconstants, supercell):
    """Original (slow) loop implementation, kept as reference."""
    ions = supercell.get_positions()
    basisvec = supercell.cell
    distances = supercell.get_all_distances(mic=True, vector=True)
    ions_doubled = np.block([[ions], [ions + basisvec[0, :]]])
    ions_doubled = np.block([[ions_doubled], [ions_doubled + basisvec[1, :]]])
    ions_doubled = np.block([[ions_doubled], [ions_doubled + basisvec[2, :]]])
    doubled = Atoms(positions=ions_doubled, cell=np.array(basisvec) * 2, pbc=True)
    distances_new = doubled.get_all_distances(mic=True, vector=True)
    a = np.array(basisvec)[0]
    fd = np.zeros([8 * natom, 8 * natom])
    for k in range(8 * natom):
        nk = k % natom
        for m in range(8 * natom):
            diff = np.sum((distances_new[k, m, :] - distances[nk, :, :]) ** 2, axis=1)
            j = np.nonzero(diff == diff.min())[0][0]
            ok = abs(np.sum(distances[nk, j, :] * a) / np.linalg.norm(a)) <= 0.48 * np.linalg.norm(a)
            if diff.min() < 1e-8 and ok:
                fd[k, m] = forceconstants[nk, j]
    return fd


def reference_dimension(fc, fd, threshold):
    g = nx.from_numpy_array(fc >= threshold)
    gd = nx.from_numpy_array(fd >= threshold)
    return calc_dimension(list(nx.connected_components(g)), list(nx.connected_components(gd)))


def random_case(seed):
    rng = np.random.default_rng(seed)
    if seed % 2:
        cell = np.diag(rng.uniform(4, 9, 3)) + rng.uniform(-2, 2, (3, 3)) * (1 - np.eye(3))
        frac = rng.random((rng.integers(3, 12), 3))
    else:  # lattice-like supercell: many pairs exactly half a cell apart
        reps = rng.integers(1, 4, 3)
        base = rng.random((2, 3))
        frac = np.array([(b + [x, y, z]) / reps for x in range(reps[0]) for y in range(reps[1])
                         for z in range(reps[2]) for b in base])
        cell = np.diag(rng.uniform(3, 6, 3) * reps)
    atoms = Atoms('H%d' % len(frac), scaled_positions=frac, cell=cell, pbc=True)
    n = len(atoms)
    d = atoms.get_all_distances(mic=True)
    fc = 10 * np.exp(-d) * rng.uniform(0.5, 1.5, (n, n))
    fc = (fc + fc.T) / 2
    fc[rng.random((n, n)) < 0.3] = 0
    np.fill_diagonal(fc, 0)
    return n, fc, atoms


@pytest.mark.parametrize("seed", range(12))
def test_doubled_fc_matches_reference(seed):
    n, fc, atoms = random_case(seed)
    reference = reference_doubled_fc(n, fc, atoms)
    assert np.allclose(doubled_fc(n, fc, atoms).toarray(), reference, rtol=0, atol=1e-12)

    thresholds = sorted(set(np.round(np.linspace(0, fc.max(), 15), 2).tolist()))
    dims, _, _ = scanner(n, fc, atoms, thresholds)
    for t in thresholds:
        assert dims[str(t)] == reference_dimension(fc, reference, t)


def test_mic_vectors_match_ase():
    _, _, atoms = random_case(4)
    atoms = atoms.repeat((3, 3, 3))
    assert np.array_equal(mic_distance_vectors(atoms.positions, np.array(atoms.cell)),
                          atoms.get_all_distances(mic=True, vector=True))
