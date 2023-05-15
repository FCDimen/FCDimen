#from tqdm import tqdm
from fcdimen.functions.cluster import doubled_fc, connectivity
from fcdimen.functions.dimensionality import calc_dimension
from fcdimen.functions.pbar import progressbar


def scanner(nat,Fp, Fp2, sc, ths):
    """
    Parameters:
        nat:  Integer
         number of atoms
        Fp: ndarray
         Zero matrix with size of number of atoms
        Fp2: ndarray
         Transpose of reshaped force constants matrix
        sc:  ASE Atoms object
         supercell
        ths: List
         list of thresholds for scanning
    Returns:
        fcdimlist: Dictionary
         Dictionary of predicted dimensionality with scores
    """
    
    fcdimlist = {}
    Fpb = doubled_fc (nat, Fp, sc)
    
    if len(ths) != 1:
        ind1, indices1 = connectivity(Fp, Fpb, 0.5)
        dimensionality1 = calc_dimension(ind1, indices1)
        minmaximumforce = min(Fp.max(axis=0))
        ind2, indices2 = connectivity(Fp, Fpb, (minmaximumforce-0.01))
        dimensionality2 = calc_dimension(ind2, indices2)
    

        #for p in progressbar(ths, "Progress: ", 40):
        for p in ths:
            ind, indices = connectivity(Fp, Fpb, p)
            dimensionality = calc_dimension(ind, indices)
            fcdimlist[str(p)] = dimensionality
    else:
        ind, indices = connectivity(Fp, Fpb, ths)
        dimensionality1 = calc_dimension(ind, indices)
        fcdimlist = None
        dimensionality2 = None


    return fcdimlist, dimensionality1, dimensionality2

def calc_score(dim_fca, maxforce):
    """Calculating dimensionality score

    Parameters:
    dim_fca: Dictionary
     calculated dimensionality of each threshold
    maxforce: Float
     Maximum force constant in structure
    Returns:

    dimen_score: Dictionary
     dimensionalities and corresponding scores
    """

    th0d = []
    th1d = []
    th2d = []
    th3d = []
    th01d = []
    th02d = []
    th03d = []
    th012d = []
    th013d = []
    th023d = []
    th12d = []
    th13d = []
    th123d = []
    th23d = []


    dT = {} # difference of thresholds (maximum and minimum occurrence)
    for i in dim_fca.keys():
            if dim_fca[i] == "0D":
                th0d.append(float(i))

            if dim_fca[i] == "1D":
                th1d.append(float(i))

            if dim_fca[i] == "2D":
                th2d.append(float(i))

            if dim_fca[i] == "3D":
                th3d.append(float(i))

            if dim_fca[i] == "01D":
                th01d.append(float(i))

            if dim_fca[i] == "02D":
                th02d.append(float(i))

            if dim_fca[i] == "03D":
                th03d.append(float(i))

            if dim_fca[i] == "012D":
                th012d.append(float(i))

            if dim_fca[i] == "013D":
                th013d.append(float(i))

            if dim_fca[i] == "023D":
                th023d.append(float(i))

            if dim_fca[i] == "12D":
                th12d.append(float(i))

            if dim_fca[i] == "13D":
                th13d.append(float(i))

            if dim_fca[i] == "123D":
                th123d.append(float(i))

            if dim_fca[i] == "23D":
                th23d.append(float(i))

    # removing empty threshold list that has only one threshold

    if len(th0d) > 1:
        # if (((max(th0d) - min(th0d))/maxforce))/2 > 0.5:
        # if (((max(th0d) - min(th0d))/maxforce)) > 0.9:
        dT['0D'] = max(th0d) - min(th0d)

    if len(th1d) > 1:
        dT['1D'] = max(th1d) - min(th1d)

    if len(th2d) > 1:
        dT['2D'] = max(th2d) - min(th2d)

    if len(th3d) > 1:
        # 0.0 and 0.1 are always 3D so we remove them scoring
        if max(th3d) != 0.1:
           dT['3D'] = max(th3d) - min(th3d)

    if len(th01d) > 1:
        dT['01D'] = max(th01d) - min(th01d)

    if len(th02d) > 1:
        dT['02D'] = max(th02d) - min(th02d)

    if len(th03d) > 1:
        dT['03D'] = max(th03d) - min(th03d)

    if len(th012d) > 1:
        dT['012D'] = max(th012d) - min(th012d)

    if len(th013d) > 1:
        dT['013D'] = max(th013d) - min(th013d)

    if len(th023d) > 1:
        dT['023D'] = max(th023d) - min(th023d)

    if len(th12d) > 1:
        dT['12D'] = max(th12d) - min(th12d)

    if len(th13d) > 1:
        dT['13D'] = max(th13d) - min(th13d)

    if len(th123d) > 1:
        dT['123D'] = max(th123d) - min(th123d)

    if len(th23d) > 1:
        dT['23D'] = max(th23d) - min(th23d)

    #Normalize scores with Maximum force constatnt in structure
    #dimen_score = {key: format((value/maxforce)*100, ".2f") for key, value in dT.items()}
    dimen_score = {key: (value/maxforce) for key, value in dT.items()}

    return dimen_score

