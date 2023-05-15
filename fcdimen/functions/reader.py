import os
import phonopy
from ase import Atoms

def read_data(filename="phonon.yaml"):
    """Read yaml file from phonopy phonon  calculations.

    Parameters:

    filename: string
       name of phonopy generated yaml file
       *Note: FORCE_SETS file should be present in same folder as yaml file

    Returns:
    ph: Phonopy object
    uca: ASE Atoms object
     unitcell of loaded Phonopy object
    fc: ndarray
      array of force constants matrix
    sc: ASE Atoms object    
      supercell of loaded Phonopy object    
    """
    try:
        # load forces from phonopy object
        ph = phonopy.load(filename, force_sets_filename="FORCE_SETS")
        uc = ph.unitcell

        #uca = Atoms(uc.get_chemical_symbols(), positions=uc.get_positions(), cell=uc.get_cell(), pbc=True) #deprecated
        uca = Atoms(uc.symbols, positions=uc.positions, cell=uc.cell, pbc=True)
        # Generate Force constants matrix
        ph.produce_force_constants()
        fc = ph.force_constants
        # Generate Supercell
        phsc = ph.supercell
        #sc = Atoms(phsc.get_chemical_symbols(), positions=phsc.get_positions(), cell=phsc.get_cell(), pbc=True) #d3precated
        sc = Atoms(phsc.symbols, positions=phsc.positions, cell=phsc.cell, pbc=True)
    except:
        print("Compact version yaml file selected")
        # if FORCE_SETS not exist, yaml file should include force constants compact constants
        #ph = phonopy.load(filename, is_compact_fc=False, log_level=1, produce_fc=True)
        ph = phonopy.load(filename, is_compact_fc=True, log_level=1)
        #Check freq negative
        #[freq, ev] = ph.get_frequencies_with_eigenvectors(q=[0,0,0])
        uc = ph.unitcell
        uca = Atoms(uc.symbols, positions=uc.positions, cell=uc.cell, pbc=True)
        ph.produce_force_constants()
        fc = ph.force_constants
        phsc = ph.supercell
        sc = Atoms(phsc.symbols, positions=phsc.positions, cell=phsc.cell, pbc=True)


    return ph, uca, fc, sc
