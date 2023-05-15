def calc_dimension(ind,indices):
    """Calculating dimensionality

    Parameters:

    ind: List
     list of index of cluster
    indices: List
     list of doubled clusters indexes 

    Returns:

    fcadim: String
     Dimensionality of structure (0D, 1D, 2D, 3D or mixture of them or unknown)
    """

    c0id = {}
    for i in range(len(ind)):
        x = list(ind[i])
        c0id[min(x)] = len(x)

    cxid = {}
    #for i in range(len(indices[0])):
    for i in range(len(indices)):
        #x = list(indices[0][i])
        x = list(indices[i])  # it is not 3 any more
        cxid[min(x)] = len(x)

    """cyid = {}
    for i in range(len(indices[1])):
        x = list(indices[1][i])
        cyid[min(x)] = len(x)

    czid = {}
    for i in range(len(indices[2])):
        x = list(indices[2][i])
        czid[min(x)] = len(x)"""

    
    dimenlist = []
    for i in c0id.keys():
        #dimenlist.append(sum([c0id[i] == cxid[i] / 2, c0id[i] == cyid[i] / 2, c0id[i] == czid[i] / 2]))
        #dimenlist.append(sum([c0id[i] == cxid[i] / 2, c0id[i] == cxid[i] / 4, c0id[i] == cxid[i] / 8]))
        if c0id[i] == cxid[i] / 2:
           dimenlist.append(1)     
        elif c0id[i] == cxid[i] / 4:
           dimenlist.append(2)
        elif c0id[i] == cxid[i] / 8:
           dimenlist.append(3)
        else:
           dimenlist.append(0)
           
           
    #print(dimenlist)       
    l = sorted(list(set(dimenlist)))

    ans = ' '
    for i in l:
        ans = ans + str(i)

    if ans == ' ':
        fcadim = "Unknown"
    else:
        fcadim = ans + "D"
        fcadim = fcadim.strip()
        
    return fcadim
