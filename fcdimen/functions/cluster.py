import numpy as np
import networkx as nx
from ase import Atoms
from fcdimen.functions.pbar import progressbar

def doubled_Fc (natom, FC, SC):
    """Generating doubled supercells force constants

    Parameters:
    natom: integer
     Number of atoms in the initial supercell
    FC: ndarray
     Force constant matrix of the initial supercell
    SC: ASE Atoms object
      Initial supercell

    Returns:
    FCNew: ndarray
     Force constant matrices array of the doubled supercell
    """
    print('Supercell chemical formula : ' + SC.get_chemical_formula(mode='metal'))
    
    # Get initial Supercell properties
    ions = SC.get_positions()
    basisvec = SC.cell
    distances = SC.get_all_distances(mic=True, vector=True)
    MaxFC = max(FC.max(axis=0))

    # Create the doubled Supercell
    natomNew = natom * 8
    FCNew = np.zeros([natomNew, natomNew, 3])
    basisvecNew = np.block([basisvec]) * 2

    ionsNew = np.block([[ions], [ions + basisvec[0, :]]])
    ionsNew = np.block([[ionsNew], [ionsNew + basisvec[1, :]]])
    ionsNew = np.block([[ionsNew], [ionsNew + basisvec[2, :]]])

    SCNew = Atoms(positions=ionsNew, cell=basisvecNew, pbc=True)
    distancesNew = SCNew.get_all_distances(mic=True, vector=True)
    # Inital value for Maximum force in new supercell
    Fmax = 0

    for k in progressbar(range(natomNew), "Progress: ", 40):
        nk = np.mod(k, natom)
        for m in range(natomNew):
            distances_diff = np.sum(((distancesNew[k, m, :] - distances[nk, :, :]) ** 2), axis=1)
            distances_diff = np.squeeze(distances_diff)
            mindiff = np.nonzero(distances_diff == min(distances_diff))[0][0]
            if (min(distances_diff) < 1e-8) and (check_dist(distances[nk, mindiff, :], basisvec)==True):
                FCNew[k, m, 0] = FC[nk, mindiff]
                FCNew[k, m, 1] = FC[nk, mindiff]
                FCNew[k, m, 2] = FC[nk, mindiff]

            elif (check_dist(distances[nk, mindiff, :], basisvec) == False) and (Fmax < FC[nk, mindiff]):
                Fmax = FC[nk, mindiff]
    
    if (Fmax / MaxFC) > 0.2:
         print("""
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                    WARNING!                          !
 !       Hopefully, you know what are you doing.        !
 !       Supercell size is not big enough,              !
 !          results may not be reliable!!               !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 """)
        
    return FCNew

def check_dist(distances, mindiff):
    """Check distances in doubled Supercell smaller
    than half of the value in the initial supercell and more than
    20% of MaxFC

    Parameters:
    distances: ndarray
     atomic distances in the initial supercell
    mindiff: integer
     Difference between atomic distances and minimum distance

    Returns:
    Boolean
     True if atomic distances are big enough
    """
    for i in range(3):
        if abs(np.sum(distances * mindiff[i]) / np.linalg.norm(mindiff[i])) > 0.48 * np.linalg.norm(mindiff[i]):
            return False
        else:
            return True


def connectivity(FC, FCNew, thresholds):
    """Checking connectivity of atoms

    Parameters:
    FC: ndarray
     Force constant matrix of the initial supercell
    FCNew: ndarray
     Force constant matrices array of the doubled supercell
    thresholds : float
     selected threshold(s)

    Returns:
    indicesSC: list
     list of initial Supercell clusters indices
    indicesSCNew: list
     list of doubeled Supercell clusters indices
    """
    # Check NetworkX version
    if int(nx.__version__.split(".")[0]) < 3:
        GraphSC = nx.from_numpy_matrix(FC >= thresholds)
        GraphSCNew = nx.from_numpy_matrix(FCNew[:,:,0] >= thresholds)
    else:
        GraphSC = nx.from_numpy_array(FC >= thresholds)
        GraphSCNew = nx.from_numpy_array(FCNew[:,:,0] >= thresholds)

    # Get the indices of the simple supercell
    indicesSC = [c for c in nx.connected_components(GraphSC)]
    indicesSCNew = [c for c in nx.connected_components(GraphSCNew)]

    return indicesSC, indicesSCNew
