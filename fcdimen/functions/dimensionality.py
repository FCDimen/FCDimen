import numpy as np


def _dimension_string(dimensionalities):
    """Build the dimensionality string, e.g. [0, 2, 2] -> '02D'."""
    digits = ''.join(str(i) for i in sorted(set(dimensionalities)))
    if digits == '':
        return "Unknown"
    return digits + "D"


def dimension_from_sizes(size, size_doubled):
    """Calculating dimensionality from cluster sizes

    Parameters:
    size: array of int
     size of each cluster in the initial supercell
    size_doubled: array of int
     size of the doubled-supercell cluster containing the same atoms

    Returns:
    dimensionality: String
     Dimensionality of structure (0D, 1D, 2D, 3D, mix or unknown)
    """
    size = np.asarray(size)
    size_doubled = np.asarray(size_doubled)
    dims = np.where(size_doubled == 2 * size, 1,
                    np.where(size_doubled == 4 * size, 2,
                             np.where(size_doubled == 8 * size, 3, 0)))
    return _dimension_string(dims.tolist())


def calc_dimension(supercell_indices, supercell_doubled_indices):
    """Calculating dimensionality

    Parameters:
    supercell_indices: list
     generate list of initial Supercell clusters indices
    supercell_doubled_indices: list
     generate list of doubled Supercell clusters indices

    Returns:
    dimensionality: String
     Dimensionality of structure (0D, 1D, 2D, 3D, mix or unknown)
    """
    supercell_clusters = {min(c): len(c) for c in supercell_indices}
    supercell_doubled_clusters = {min(c): len(c) for c in supercell_doubled_indices}

    keys = list(supercell_clusters.keys())
    size = [supercell_clusters[i] for i in keys]
    size_doubled = [supercell_doubled_clusters[i] for i in keys]
    return dimension_from_sizes(size, size_doubled)
