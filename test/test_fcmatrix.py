import pytest
import numpy as np
from fcdimen.functions.fcmat import force_matrix


fc = np.ndarray(shape=(2,2,3,3), dtype=float, order='F')

def test_fc_matrix():
    nat, Fp, Fp2 = force_matrix(fc)
    assert nat == 2
    assert len(Fp) == 2
    assert len(Fp2) == 3