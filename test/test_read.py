import pytest
import phonopy

from fcdimen.functions.reader import read_data



def test_read_phonon():
        ph, uca, fc, sc = read_data(filename="phonon.yaml")
        assert str(uca.symbols) == "S6Mo3"
        assert str(sc.symbols) == "S96Mo48"
        expected_result = 11.727323999864705

        assert abs(fc[0][0][0][0] - expected_result) < 1.0e-6

