import numpy as np
from fcdimen.functions.cluster import doubled_cell_data, DoubledGraph
from fcdimen.functions.dimensionality import dimension_from_sizes


def scanner(natom, forceconstants_zero_diagonal, supercell, thresholds):
    """ Calculating dimensionality based on selected thresholds

    Parameters:
    natom:  Integer
     number of atoms
    forceconstants_zero_diagonal: ndarray
     force constants matrix where diagonal elements are zero
    supercell:  ASE Atoms object
      supercell
    thresholds: List
     list of thresholds for scanning

    Returns:
    dimensionalities_thresholds: Dictionary
     Predicted dimensionalities for each thresholds
    dimensionality_class1: string
     Predicted dimensionality with classification I
    dimensionality_class2: string
     Predicted dimensionality with classification II
    """
    forceconstants_zero_diagonal = np.asarray(forceconstants_zero_diagonal, dtype=float)
    keep, block = doubled_cell_data(natom, forceconstants_zero_diagonal, supercell)

    # Graphs are built once; every threshold then costs O(N)
    graph = DoubledGraph(forceconstants_zero_diagonal, keep, block)

    def dimension_at(threshold):
        _, size, size_doubled = graph.cluster_sizes(threshold)
        return dimension_from_sizes(size, size_doubled)

    if len(thresholds) != 1:
        dimensionality_class1 = dimension_at(0.5)
        MinFC = min(forceconstants_zero_diagonal.max(axis=0))
        dimensionality_class2 = dimension_at(MinFC - 0.01)

        dimensionalities_thresholds = {}
        for p in thresholds:
            dimensionalities_thresholds[str(p)] = dimension_at(p)
    else:
        dimensionality_class1 = dimension_at(thresholds[0])
        dimensionalities_thresholds = None
        dimensionality_class2 = None

    return dimensionalities_thresholds, dimensionality_class1, dimensionality_class2


def calc_score(dimensionalities_thresholds, maxforceconstant):
    """Calculating dimensionality scores

    Parameters:
    dimensionalities_thresholds: Dictionary
     calculated dimensionality of each threshold
    maxforceconstant: Float
     Maximum force constant in structure
    Returns:

    dimen_score: Dictionary
     dimensionalities and corresponding scores
    """
    # Same dimensionality labels and output order as before
    labels = ['0D', '1D', '2D', '3D', '01D', '02D', '03D', '012D', '013D',
              '023D', '12D', '13D', '123D', '23D']
    thresholds_of = {label: [] for label in labels}
    for threshold, dimensionality in dimensionalities_thresholds.items():
        if dimensionality in thresholds_of:
            thresholds_of[dimensionality].append(float(threshold))

    # difference of thresholds (maximum and minimum occurrence),
    # dimensionalities found at only one threshold are removed
    dT = {}
    for label in labels:
        values = thresholds_of[label]
        if len(values) > 1:
            # 0.0 and 0.1 are always 3D so we remove them scoring
            if label == '3D' and max(values) == 0.1:
                continue
            dT[label] = max(values) - min(values)

    # Normalize scores with Maximum force constatnt in structure
    dimen_score = {key: (value / maxforceconstant) for key, value in dT.items()}

    return dimen_score
