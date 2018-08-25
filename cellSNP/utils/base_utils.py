
def id_mapping(IDs1, IDs2, uniq_ref_only=True, IDs2_sorted=False):
    """
    Mapping IDs2 to IDs1. IDs1 (ref id) can have repeat values, but IDs2 need 
    to only contain unique ids.
    Therefore, IDs2[rv_idx] will be the same as IDs1.
    
    Parameters
    ----------
    IDs1 : array_like or list
        ids for reference.
    IDs2 : array_like or list
        ids waiting to map.
        
    Returns
    -------
    RV_idx : array_like, the same length of IDs1
        The index for IDs2 mapped to IDs1. If an id in IDs1 does not exist 
        in IDs2, then return a None for that id.
    """
    idx1 = sorted(range(len(IDs1)), key=IDs1.__getitem__)
    if IDs2_sorted:
        idx2 = range(len(IDs2))
    else:
        idx2 = sorted(range(len(IDs2)), key=IDs2.__getitem__)
    RV_idx1, RV_idx2 = [], []
    
    i, j = 0, 0
    while i < len(idx1):
        if j == len(idx2) or IDs1[idx1[i]] < IDs2[idx2[j]]:
            RV_idx1.append(idx1[i])
            RV_idx2.append(None)
            i += 1
        elif IDs1[idx1[i]] == IDs2[idx2[j]]:
            RV_idx1.append(idx1[i])
            RV_idx2.append(idx2[j])
            i += 1
            if uniq_ref_only: j += 1
        elif IDs1[idx1[i]] > IDs2[idx2[j]]:
            j += 1
            
    origin_idx = sorted(range(len(RV_idx1)), key=RV_idx1.__getitem__)
    RV_idx = [RV_idx2[i] for i in origin_idx]
    return RV_idx


def unique_list(X):
    """unique a list with index and count
    Example
    -------
    >>> unique_list([1,2,4,5,3,2,4])
    >>> ([1, 2, 3, 4, 5], [0, 1, 4, 2, 3], [1, 2, 1, 2, 1])
    """
    idx = sorted(range(len(X)), key=X.__getitem__)
    X_uniq = []
    X_count = []
    idx_uniq = []
    for i in idx:
        if len(X_uniq) == 0 or X[i] != X_uniq[-1]:
            X_uniq.append(X[i])
            X_count.append(1)
            idx_uniq.append(i)
        else:
            X_count[-1] += 1
    return X_uniq, idx_uniq, X_count
