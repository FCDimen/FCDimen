import numpy as np


def force_matrix(forceconstants):
    """Read Force constants matrix and prepare requirements matrices

    Parameters:
    forceconstants: ndarray
      array of force constants matrix, shape (N, N, 3, 3)

    Returns:
    natom: Integer
     number of atoms
    forceconstants_zero_diagonal: ndarray
     force constants matrix where diagonal elements are zero
     (Frobenius norm of each 3x3 block)
    forceconstants_reshaped: ndarray
     Transpose of reshaped force constants matrix,
     forceconstants_reshaped[b, i, a, j] = forceconstants[i, j, a, b]
    """
    forceconstants = np.asarray(forceconstants, dtype=float)
    natom = forceconstants.shape[0]
    # only the square (N, N) block is used, as in the original implementation
    fc = forceconstants[:, :natom]

    forceconstants_reshaped = np.transpose(fc, (3, 0, 2, 1))
    forceconstants_zero_diagonal = np.sqrt(np.einsum('ijab,ijab->ij', fc, fc))
    np.fill_diagonal(forceconstants_zero_diagonal, 0)

    return natom, forceconstants_zero_diagonal, forceconstants_reshaped
