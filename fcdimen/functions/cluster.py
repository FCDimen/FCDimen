import numpy as np
from scipy.sparse import coo_matrix, issparse
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree, breadth_first_order
from ase.geometry.geometry import general_find_mic, get_distances

# Upper bound on the number of vectors handled at once in the minimum image
# search (memory ~ 27 * 3 * 8 bytes per vector).
_MIC_CHUNK = 200000
_PBC = np.array([True, True, True])


# ----------------------------------------------------------------------------
# Geometry of the initial supercell
# ----------------------------------------------------------------------------
def _mic_upper_pairs(positions, cell):
    """Yield blocks (ii, jj, D_ij) of minimum-image vectors for all pairs j > i.

    Same convention and bitwise the same vectors as ASE
    ``Atoms.get_all_distances(mic=True, vector=True)``, which evaluates
    mic(p_j - p_i) for j > i (D[j, i] = -D[i, j]). Small cells use ASE in a
    single batch; large cells are processed in blocks of pairs (every block
    has more than one vector, which keeps the rounding identical to a single
    batch) so memory stays bounded.
    """
    positions = np.asarray(positions, dtype=float)
    cell = np.asarray(cell, dtype=float)
    natom = len(positions)
    if natom < 2:
        return
    if natom * natom <= 2 * _MIC_CHUNK:
        D, _ = get_distances(positions, cell=cell, pbc=_PBC)
        ii, jj = np.triu_indices(natom, k=1)
        yield ii, jj, D[ii, jj]
        return
    # row ranges holding about _MIC_CHUNK pairs each
    pairs_per_row = natom - 1 - np.arange(natom)
    cumulative = np.cumsum(pairs_per_row)
    edges = [0]
    while edges[-1] < natom - 1:
        start = edges[-1]
        done = cumulative[start - 1] if start > 0 else 0
        stop = int(np.searchsorted(cumulative, done + _MIC_CHUNK, side="right"))
        stop = min(max(stop, start + 1), natom - 1)
        if cumulative[natom - 2] - cumulative[stop - 1] < 2:
            stop = natom - 1  # do not leave a block with a single vector
        edges.append(stop)
    for start, stop in zip(edges[:-1], edges[1:]):
        rows = np.arange(start, stop)
        ii, jj = np.nonzero(np.arange(natom)[None, :] > rows[:, None])
        ii = rows[ii]
        vmin, _ = general_find_mic(positions[jj] - positions[ii], cell, pbc=_PBC)
        yield ii, jj, vmin


def mic_distance_vectors(positions, cell):
    """Minimum-image distance vectors between all atoms of a periodic cell.

    Equivalent to ASE ``Atoms.get_all_distances(mic=True, vector=True)``,
    but evaluated in chunks so memory stays bounded for large supercells.

    Parameters:
    positions: ndarray (N, 3)
     Cartesian positions
    cell: ndarray (3, 3)
     Lattice vectors (rows)

    Returns:
    D : ndarray (N, N, 3)
     D[i, j] = minimum-image vector from atom i to atom j
    """
    D = np.zeros((len(positions), len(positions), 3))
    for ii, jj, v in _mic_upper_pairs(positions, cell):
        D[ii, jj] = v
        D[jj, ii] = -v
    return D


def check_dist(distances, basisvec):
    """Check that interatomic vectors are not close to the supercell boundary.

    A vector fails the check when its projection on the first lattice
    vector is larger than 48% of that lattice vector length.
    NOTE: only the first lattice vector is tested. This reproduces the
    behaviour of the original loop, which returned on its first iteration.

    Parameters:
    distances: ndarray (..., 3)
     interatomic distance vector(s)
    basisvec: ndarray (3, 3)
     lattice vectors of the supercell

    Returns:
    Boolean or ndarray of booleans
     True if the vector is far enough from the boundary
    """
    basisvec = np.asarray(basisvec, dtype=float)
    a = basisvec[0]
    norm_a = np.linalg.norm(a)
    projection = np.abs(np.asarray(distances, dtype=float) @ a) / norm_a
    ok = ~(projection > 0.48 * norm_a)
    if np.ndim(ok) == 0:
        return bool(ok)
    return ok


def doubled_cell_data(natom, forceconstants, supercell):
    """Describe the doubled (2x2x2) supercell without building it.

    The doubled supercell is periodic in the initial supercell. For every pair
    (i, j) of the initial supercell, exactly one image of j in the doubled cell
    has, seen from i, the same minimum-image vector d_ij as in the initial cell.
    Its sub-cell index (0..7) follows from the lattice translation
    p_i + d_ij - p_j. Only pairs passing check_dist are kept, as in the
    original implementation.

    Parameters:
    natom: integer
     Number of atoms in the initial supercell
    forceconstants: ndarray (N, N)
     Force constant matrix of the initial supercell
    supercell: ASE Atoms object
      Initial supercell

    Returns:
    keep: ndarray of bool (N, N)
     pairs that pass the boundary check
    block: ndarray of int8 (N, N)
     sub-cell index s = sx + 2*sy + 4*sz of the matching image of j seen from i
    """
    print('Supercell chemical formula : ' + supercell.get_chemical_formula(mode='metal'))

    forceconstants = np.asarray(forceconstants, dtype=float)
    ions = supercell.get_positions()
    basisvec = np.array(supercell.cell, dtype=float)
    inv_basis = np.linalg.inv(basisvec)

    # keep and block are symmetric (D_ji = -D_ij), so only pairs j > i are needed
    keep = np.ones((natom, natom), dtype=bool)
    block = np.zeros((natom, natom), dtype=np.int8)
    for ii, jj, distances in _mic_upper_pairs(ions, basisvec):
        ok = check_dist(distances, basisvec)
        translation = ions[ii] + distances - ions[jj]
        ncell = np.rint(translation @ inv_basis).astype(np.int64)
        sub = ((ncell[:, 0] & 1) + 2 * (ncell[:, 1] & 1) + 4 * (ncell[:, 2] & 1)).astype(np.int8)
        keep[ii, jj] = ok
        keep[jj, ii] = ok
        block[ii, jj] = sub
        block[jj, ii] = sub

    # Warning: largest force constant between atoms close to the boundary
    MaxFC = forceconstants.max()
    rejected = forceconstants[~keep]
    Fmax = rejected.max() if rejected.size else 0.0
    if Fmax > 0 and (Fmax / MaxFC) > 0.2:
        print("""
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                    WARNING!                          !
 !       Hopefully, you know what are you doing.        !
 !       Supercell size is not big enough,              !
 !          results may not be reliable!!               !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 """)
    return keep, block


def doubled_fc(natom, forceconstants, supercell):
    """Generating doubled supercells force constants

    Only needed if the doubled matrix itself is wanted; the dimensionality
    scan (score.scanner) works without it.

    Parameters:
    natom: integer
     Number of atoms in the initial supercell
    forceconstants: ndarray
     Force constant matrix of the initial supercell
    supercell: ASE Atoms object
      Initial supercell

    Returns:
    forceconstants_doubled : scipy.sparse.csr_matrix (8N, 8N)
     Force constant matrix of the doubled supercell (zero entries not stored).
     Use .toarray() to get the dense matrix.
    """
    forceconstants = np.asarray(forceconstants, dtype=float)
    keep, block = doubled_cell_data(natom, forceconstants, supercell)
    ii, jj = np.nonzero(keep & (forceconstants != 0))
    values = forceconstants[ii, jj]
    sblock = block[ii, jj].astype(np.int64)
    # Replicate for the 8 sub-cells: row block t, column block (s XOR t)
    rows = np.concatenate([ii + t * natom for t in range(8)])
    cols = np.concatenate([jj + (sblock ^ t) * natom for t in range(8)])
    data = np.tile(values, 8)
    return coo_matrix((data, (rows, cols)), shape=(8 * natom, 8 * natom)).tocsr()


# ----------------------------------------------------------------------------
# Graphs at many thresholds
# ----------------------------------------------------------------------------
def _undirected_edges(forceconstants):
    """Edges (i < j) with weight max(F_ij, F_ji); zero and NaN entries skipped."""
    if issparse(forceconstants):
        m = forceconstants.tocoo()
        r, c, v = m.row, m.col, m.data
        good = ~np.isnan(v) & (r != c) & (v != 0)
        r, c, v = r[good], c[good], v[good]
        lo, hi = np.minimum(r, c), np.maximum(r, c)
        order = np.lexsort((-v, hi, lo))
        lo, hi, v = lo[order], hi[order], v[order]
        first = np.ones(len(v), dtype=bool)
        first[1:] = (lo[1:] != lo[:-1]) | (hi[1:] != hi[:-1])
        return lo[first], hi[first], v[first]
    A = np.asarray(forceconstants, dtype=float)
    if A.ndim == 3:  # old (N, N, 3) format of the doubled force constants
        A = A[:, :, 0]
    W = np.fmax(A, A.T)  # NaN means "no edge"
    lo, hi = np.nonzero(np.triu(W != 0, k=1) & ~np.isnan(W))
    return lo, hi, W[lo, hi]


class ThresholdGraph:
    """Connected components of {(i, j): F_ij >= t} for many thresholds t.

    An edge (i, j) exists when F_ij >= t or F_ji >= t (undirected graph, as
    in networkx.from_numpy_array). Components at any threshold t > 0 equal the
    components of the maximum spanning forest restricted to edges >= t, so the
    forest is computed once and each threshold costs only O(N).
    For t <= 0 all entries (including zeros) are edges: one component.
    """

    def __init__(self, forceconstants):
        self.n = forceconstants.shape[0]
        lo, hi, v = _undirected_edges(forceconstants)
        if len(v):
            # Rank transform (handles inf and ties): heaviest edge -> weight 1
            uniq, inv = np.unique(v, return_inverse=True)
            rank_weight = (len(uniq) - inv).astype(float)
            ranked = coo_matrix((rank_weight, (lo, hi)), shape=(self.n, self.n)).tocsr()
            del lo, hi, v, inv, rank_weight
            forest = minimum_spanning_tree(ranked).tocoo()
            f_vals = uniq[len(uniq) - np.rint(forest.data).astype(np.int64)]
            order = np.argsort(-f_vals, kind="stable")  # decreasing force constant
            self.f_rows = forest.row[order].astype(np.int64)
            self.f_cols = forest.col[order].astype(np.int64)
            self.f_vals = f_vals[order]
        else:
            self.f_rows = np.zeros(0, dtype=np.int64)
            self.f_cols = np.zeros(0, dtype=np.int64)
            self.f_vals = np.zeros(0)
        self._neg_vals = -self.f_vals  # ascending, for searchsorted
        self._cache = {}

    def labels(self, threshold):
        """Component label of every node for a given threshold."""
        threshold = float(threshold)
        if threshold <= 0:
            return np.zeros(self.n, dtype=np.int32)
        # forest edges with value >= threshold form a prefix
        nedge = int(np.searchsorted(self._neg_vals, -threshold, side="right"))
        if nedge not in self._cache:
            graph = coo_matrix((np.ones(nedge), (self.f_rows[:nedge], self.f_cols[:nedge])),
                               shape=(self.n, self.n))
            self._cache[nedge] = connected_components(graph, directed=False)[1]
        return self._cache[nedge]


def _span_sizes():
    """|span(S)| in (Z2)^3 for every subset S of {1..7} given as a 7-bit mask."""
    sizes = np.zeros(128, dtype=np.int64)
    for mask in range(128):
        group = {0}
        for v in range(1, 8):
            if mask >> (v - 1) & 1:
                group |= {g ^ v for g in group}
        sizes[mask] = len(group)
    return sizes


_SPAN_SIZE = _span_sizes()


class DoubledGraph:
    """Cluster sizes in the initial and in the doubled supercell for any threshold.

    The doubled-supercell graph is a (Z2)^3 "lift" of the graph G' made of the
    pairs that pass check_dist: edge (i, j) of G' becomes the 8 edges
    (i, t) - (j, t XOR s_ij), s_ij being the sub-cell index of the image.
    A cluster C' of G' then extends into a doubled-cell cluster of size
    |C'| * |H|, where H is the subgroup of (Z2)^3 generated by the sub-cell
    shifts accumulated around the loops of C'. Loops are measured against a
    fixed maximum spanning forest, so each threshold costs O(N) and the
    8N x 8N doubled matrix is never built.
    """

    def __init__(self, forceconstants, keep, block):
        A = np.asarray(forceconstants, dtype=float)
        n = A.shape[0]
        self.n = n
        W = np.fmax(A, A.T)
        self.graph = ThresholdGraph(W)                       # initial supercell
        Wp = np.where(keep & keep.T, W, np.nan)              # pairs kept in doubled cell
        self.graph_p = ThresholdGraph(Wp)

        # Sub-cell "potential" of every atom along the spanning forest of G'
        pot = np.zeros(n, dtype=np.int64)
        fr, fc = self.graph_p.f_rows, self.graph_p.f_cols
        ncomp, comp = connected_components(
            coo_matrix((np.ones(len(fr)), (fr, fc)), shape=(n, n)), directed=False)
        roots = np.unique(comp, return_index=True)[1]
        # virtual node n linked to one root per tree -> a single BFS
        rows = np.concatenate([fr, fc, np.full(ncomp, n), roots])
        cols = np.concatenate([fc, fr, roots, np.full(ncomp, n)])
        tree = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n + 1, n + 1)).tocsr()
        order, pred = breadth_first_order(tree, n, directed=True, return_predecessors=True)
        block64 = block.astype(np.int64)
        for node in order[1:]:
            parent = pred[node]
            if parent != n:
                pot[node] = pot[parent] ^ block64[parent, node]

        # Loop shift of every edge of G' and, for every atom, the largest force
        # constant of an incident edge carrying each nonzero shift v = 1..7
        lo, hi, w = _undirected_edges(Wp)
        shift = pot[lo] ^ pot[hi] ^ block64[lo, hi]
        nz = shift != 0
        lo, hi, w, shift = lo[nz], hi[nz], w[nz], shift[nz] - 1
        self.max_shift_fc = np.full((7, n), -np.inf)
        np.maximum.at(self.max_shift_fc, (shift, lo), w)
        np.maximum.at(self.max_shift_fc, (shift, hi), w)

    def cluster_sizes(self, threshold):
        """For each cluster of the initial supercell: (smallest atom index,
        size in initial supercell, size of the doubled-cell cluster of that atom)."""
        threshold = float(threshold)
        lab = self.graph.labels(threshold)
        _, first, size = np.unique(lab, return_index=True, return_counts=True)
        if threshold <= 0:
            # every pair is an edge: one cluster of N atoms, 8N in the doubled cell
            return first, size, np.full(len(first), 8 * self.n)

        labp = self.graph_p.labels(threshold)
        size_p = np.bincount(labp)
        # which loop shifts are present in each cluster of G'
        order = np.argsort(labp, kind="stable")
        starts = np.concatenate([[0], np.nonzero(np.diff(labp[order]))[0] + 1])
        comp_max = np.maximum.reduceat(self.max_shift_fc[:, order], starts, axis=1)
        present = ((comp_max >= threshold).astype(np.int64) << np.arange(7)[:, None]).sum(axis=0)
        group_size = _SPAN_SIZE[present]
        cp = labp[first]
        return first, size, size_p[cp] * group_size[cp]


def labels_to_sets(labels):
    """Convert component labels to a list of sets of node indices."""
    order = np.argsort(labels, kind="stable")
    splits = np.nonzero(np.diff(labels[order]))[0] + 1
    return [set(group.tolist()) for group in np.split(order, splits)]


def connectivity(forceconstants, forceconstants_doubled, thresholds):
    """Checking connectivity of atoms

    Parameters:
    forceconstants: ndarray
     Force constant matrix of the initial supercell
    forceconstants_doubled : ndarray or sparse matrix
     Force constant matrix of the doubled supercell
    thresholds : float
     selected threshold

    Returns:
    supercell_indices: list
     generate list of initial Supercell clusters indices
    supercell_doubled_indices: list
     generate list of doubled Supercell clusters indices
    """
    threshold = float(np.ravel(thresholds)[0])
    lab = ThresholdGraph(forceconstants).labels(threshold)
    lab_doubled = ThresholdGraph(forceconstants_doubled).labels(threshold)
    return labels_to_sets(lab), labels_to_sets(lab_doubled)
