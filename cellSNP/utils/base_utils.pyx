
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

'''
from libcpp.vector cimport vector as cpp_vector
from .cellsnp_utils cimport c_idxstr_t, c_idxstr_qsort, idxstr_arr_init, idxstr_arr_destroy, \
                            c_idxint_t, c_idxint_qsort, idxint_arr_init, idxint_arr_destroy, safe_strcmp

cdef str_id_mapping(IDs1, IDs2, bint uniq_ref_only=True, bint IDs2_sorted=False):
    """
    Mapping IDs2 to IDs1. IDs1 (ref id) can have repeat values, but IDs2 need 
    to only contain unique ids. 
    Therefore, IDs2[rv_idx] will be the same as IDs1.

    Parameters
    ----------
    IDs1 : array_like or list of strings.
        ids for reference.
    IDs2 : array_like or list of strings.
        ids waiting to map.
        
    Returns
    -------
    RV_idx : array_like, the same length of IDs1
        The index for IDs2 mapped to IDs1. If an id in IDs1 does not exist 
        in IDs2, then return a None for that id.

    @Note  it's faster than id_mapping() when IDs1 and IDs2 are small while slower when IDs1 or IDs2 are large,
        due to the overhead of translating the strings in IDs1 and IDs2 from <str> type to <bytes> type in Python3.
    """
    cdef int n1 = len(IDs1)
    cdef int n2 = len(IDs2)
    if n1 == 0: return []
    if n2 == 0: return [None] * n1
    IDs1 = [None if x is None else bytes(x, "utf8") for x in IDs1]
    IDs2 = [None if x is None else bytes(x, "utf8") for x in IDs2]
    cdef c_idxstr_t *a1 = idxstr_arr_init(IDs1, n1)
    c_idxstr_qsort(a1, n1)
    cdef c_idxstr_t *a2 = idxstr_arr_init(IDs2, n2)
    if not IDs2_sorted:
        c_idxstr_qsort(a2, n2)

    cdef int i, j, k
    i, j = 0, 0
    RV_idx1, RV_idx2 = [], []
    while i < n1:
        if j == n2:
            RV_idx1.append(a1[i].i)
            RV_idx2.append(None)
            i += 1
        else:
            k = safe_strcmp(a1[i].s, a2[j].s)
            if k < 0:
                RV_idx1.append(a1[i].i)
                RV_idx2.append(None)
                i += 1
            elif k == 0:
                RV_idx1.append(a1[i].i)
                RV_idx2.append(a2[j].i)
                i += 1
                if uniq_ref_only: j += 1
            else:
                j += 1
    cdef int m1 = len(RV_idx1)
    cdef c_idxint_t *origin_idx = idxint_arr_init(RV_idx1, m1)
    c_idxint_qsort(origin_idx, m1)
    RV_idx = [RV_idx2[origin_idx[i].i] for i in range(m1)]
    idxint_arr_destroy(origin_idx)
    idxstr_arr_destroy(a1)
    idxstr_arr_destroy(a2)
    return RV_idx

cdef str_unique_list(X):
    """unique a str list with index and count.
    Example
    -------
    >>> unique_str_list(['1','2','4','5','3','2','4'])
    >>> (['1', '2', '3', '4', '5'], [0, 1, 4, 2, 3], [1, 2, 1, 2, 1])

    @Note  it's faster than unique_list() when X is small while slower when X is large,
        due to the overhead of translating the strings in X from <str> type to <bytes> type in Python3.
    """
    cdef int n = len(X)
    if n == 0:
        return [], [], []
    X = [None if x is None else bytes(x, "utf8") for x in X]
    cdef c_idxstr_t *a = idxstr_arr_init(X, n)
    c_idxstr_qsort(a, n)
    cdef cpp_vector[char*] X_uniq
    cdef cpp_vector[int] X_count
    cdef cpp_vector[int] idx_uniq
    cdef int i, j, m
    cdef char *s;
    for j in range(n):
        i = a[j].i
        s = a[j].s
        m = X_uniq.size()
        if m == 0 or safe_strcmp(s, X_uniq[m - 1]):
            X_uniq.push_back(s)
            X_count.push_back(1)
            idx_uniq.push_back(i)
        else:
            X_count[X_count.size() - 1] += 1
    idxstr_arr_destroy(a)
    return X_uniq, idx_uniq, X_count
'''
