import os
import re
import phonopy
from ase import Atoms


def _has_primitive_matrix(filename):
    """True if the phonopy yaml file defines a primitive_matrix."""
    pattern = re.compile(r"^primitive_matrix\s*:")
    try:
        with open(filename) as f:
            return any(pattern.match(line) for line in f)
    except OSError:
        return False


def _load(filename, **kwargs):
    """phonopy.load with settings that only matter for force constants.

    - Non-analytical term correction (Born charges) is not needed for the
      force constants, so it is switched off.
    - If the yaml file has no primitive_matrix, the primitive cell is the unit
      cell (default of phonopy < 4; phonopy >= 4 would guess one with 'auto').
    """
    if not _has_primitive_matrix(filename):
        kwargs.setdefault("primitive_matrix", "P")
    return phonopy.load(filename, is_nac=False, **kwargs)


def _to_ase(cell):
    return Atoms(cell.symbols, positions=cell.positions, cell=cell.cell, pbc=True)


def read_data(filename="phonon.yaml"):
    """Read yaml file from phonopy phonon calculations.

    Parameters:
    filename: string
       name of phonopy generated yaml file
       *Note: FORCE_SETS file should be present in same folder as yaml file

    Returns:
    phonon: Phonopy object
     force calculations results
    unitcell_structure: ASE Atoms object
     unitcell of loaded Phonopy object
    forceconstants: ndarray
      array of force constants matrix
    supercell_structure: ASE Atoms object
      supercell of loaded Phonopy object
    """
    force_sets = os.path.join(os.path.dirname(os.path.abspath(filename)), "FORCE_SETS")
    phonon = None
    if os.path.isfile(force_sets):
        try:
            # Load forces from phonopy object
            phonon = _load(filename, force_sets_filename=force_sets)
            phonon.produce_force_constants()
        except Exception as err:
            print("Reading with FORCE_SETS failed ({}), trying compact version".format(err))
            phonon = None

    if phonon is None:
        print("Compact version yaml file selected")
        # if FORCE_SETS not exist, yaml file should include forces or force constants
        phonon = _load(filename, is_compact_fc=True, log_level=1)
        phonon.produce_force_constants()

    unitcell_structure = _to_ase(phonon.unitcell)
    forceconstants = phonon.force_constants
    supercell_structure = _to_ase(phonon.supercell)

    return phonon, unitcell_structure, forceconstants, supercell_structure
