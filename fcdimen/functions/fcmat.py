import numpy as np


def force_matrix(forceconstatnts):
    """Read Force constants matrix and prepare requirements matrices

    Parameters:
    forceconstatnts: ndarray
      array of force constants matrix

    Returns:
    natom: Integer
     number of atoms
    forceconstatnts_zero: ndarray
     Zero matrix with size of number of atoms
    forceconstatnts_reshaped: ndarray
     Transpose of reshaped force constants matrix
    """
    force_mat = []
    for i in range(len(forceconstatnts)):
        for j in range(len(forceconstatnts)):
            for k in range(3):
                force_mat.append(forceconstatnts[i][j][k])

    natom = int(np.sqrt(len(force_mat) / 3))
    w = np.reshape(np.transpose(force_mat), [3, natom, natom, 3])
    forceconstatnts_reshaped = np.transpose(w, [0, 1, 3, 2])
    forceconstatnts_zero = np.zeros((natom, natom))

    for i in range(natom):
        for j in range(natom):
            forceconstatnts_zero[i, j] = np.linalg.norm(np.squeeze(forceconstatnts_reshaped[:, i, :, j]))

    for i in range(natom):
        forceconstatnts_zero[i, i] = 0

    return natom, forceconstatnts_zero, forceconstatnts_reshaped
