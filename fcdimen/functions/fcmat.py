import numpy as np

def force_matrix(forceconstatnts):
    """Read Force constants matrix

    Parameters:

    forceconstatnts: ndarray
      array of force constants matrix

    Returns:

    nat: Integer
     number of atoms
    Fp: ndarray
     Zero matrix with size of number of atoms
    Fp2: ndarray
     Transpose of reshaped force constants matrix
    """
    force_mat = []
    for i in range(len(forceconstatnts)):
        for j in range(len(forceconstatnts)):
            for k in range(3):
                force_mat.append(forceconstatnts[i][j][k])

    nat = int(np.sqrt(len(force_mat) / 3))
    w = np.reshape(np.transpose(force_mat), [3, nat, nat, 3])
    Fp2 = np.transpose(w, [0, 1, 3, 2])
    Fp = np.zeros((nat, nat))

    for i in range(nat):
        for j in range(nat):
            Fp[i, j] = np.linalg.norm(np.squeeze(Fp2[:, i, :, j]))

    for i in range(nat):
        Fp[i, i] = 0

    return nat, Fp, Fp2
