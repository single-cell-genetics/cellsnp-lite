
from libc.stdint cimport uint32_t
from pysam.libchtslib cimport BAM_CDIFF, BAM_CEQUAL, BAM_CINS, BAM_CMATCH, BAM_CSOFT_CLIP
from pysam.libcalignedsegment cimport AlignedSegment

cdef double c_max(double x, double y):
    return x if x > y else y

cdef double c_min(double x, double y):
    return x if x < y else y

cdef get_query_bases(AlignedSegment read, bint full_length=False):
    """
    @abstract            Return a list of bases in qurey sequence that are within the alignment.
    @param read          An AlignedSegment object. [AlignedSegment]
    @param full_length   If full_length is set, None values will be included for any soft-clipped or 
                         unaligned positions within the read. The returned list will thus be of the 
                         same length as the read. [bint]
    @return              A list of bases. [list]
    """
    cdef uint32_t k, i, l, pos
    cdef int op
    cdef bint _full = full_length

    cigar_tuples = read.cigartuples
    if not cigar_tuples:
        return []

    result = []
    pos = 0
    s = read.query_sequence

    for op, l in cigar_tuples:
        if op == BAM_CSOFT_CLIP or op == BAM_CINS:
            if _full:
                for i from 0 <= i < l:
                    result.append(None)
            pos += l
        elif op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
            for i from pos <= i < pos + l:
                result.append(s[i])
            pos += l
        # else: do nothing.
    return result

cdef get_query_qualities(AlignedSegment read, bint full_length=False):
    """
    @abstract            Return a list of qualities in qurey quality sequence that are within the alignment.
    @param read          An AlignedSegment object. [AlignedSegment]
    @param full_length   If full_length is set, None values will be included for any soft-clipped or 
                         unaligned positions within the read. The returned list will thus be of the 
                         same length as the read. [bint]
    @return              A list of qualities. [list]
    """
    cdef uint32_t k, i, l, pos
    cdef int op
    cdef bint _full = full_length

    cigar_tuples = read.cigartuples
    if not cigar_tuples:
        return []

    result = []
    pos = 0
    s = read.query_qualities

    for op, l in cigar_tuples:
        if op == BAM_CSOFT_CLIP or op == BAM_CINS:
            if _full:
                for i from 0 <= i < l:
                    result.append(None)
            pos += l
        elif op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
            for i from pos <= i < pos + l:
                result.append(chr(s[i] + 33))
            pos += l
        # else: do nothing.
    return result

'''
cimport libc.stdlib as c_stdlib
cimport libc.string as c_string

cdef c_idxstr_t* idxstr_arr_init(X, const int n):
    """
    @abstract     Initialize an array of c_idxstr_t structures based on the input list of strings.
    @param X      Python list of strings (could include None in the list).
    @param n      Num of strings in the list X.
    @return       Pointer to the array if success, NULL if input list is empty.
    """
    if n == 0: return NULL
    cdef c_idxstr_t *a = <c_idxstr_t*> c_stdlib.calloc(n, sizeof(c_idxstr_t))
    if a is NULL:
        raise MemoryError
    cdef int i
    for i in range(n):
        a[i].s = NULL if X[i] is None else <char*> X[i]
        a[i].i = i
    return a

cdef void idxstr_arr_destroy(c_idxstr_t *a):
    c_stdlib.free(a)   # check if the char* in the c_idxstr_t should be freed.

cdef int safe_strcmp(const char *x, const char *y) nogil:
    if x is not NULL:
        if y is not NULL: return c_string.strcmp(x, y)
        else: return 1
    elif y is not NULL: return -1
    else: return 0

cdef int c_idxstr_cmp(const void *x, const void *y) nogil:
    return safe_strcmp((<c_idxstr_t*>x).s, (<c_idxstr_t*>y).s)

cdef void c_idxstr_qsort(c_idxstr_t *a, const int n):
    c_stdlib.qsort(a, n, sizeof(c_idxstr_t), c_idxstr_cmp)

cdef c_idxint_t* idxint_arr_init(X, const int n):
    """
    @abstract     Initialize an array of c_idxint_t structures based on the input list of integers.
    @param X      Python list of integers (could include None in the list).
    @param n      Num of integers in the list X.
    @return       Pointer to the list if success, NULL if input list is empty.
    """
    if n == 0: return NULL
    cdef c_idxint_t *a = <c_idxint_t*> c_stdlib.calloc(n, sizeof(c_idxint_t))
    if a is NULL:
        raise MemoryError
    cdef int i
    for i in range(n):
        a[i].i = i
        if X[i] is None:
            a[i].d = -1
            a[i].not_none = False
        else:
            a[i].d = X[i]
            a[i].not_none = True
    return a

cdef void idxint_arr_destroy(c_idxint_t *a):
    c_stdlib.free(a)

cdef int safe_intcmp(const int x, bint x_not_none, const int y, bint y_not_none) nogil:
    if x_not_none:
        if y_not_none: return x - y
        else: return 1
    elif y_not_none: return -1
    else: return 0

cdef int c_idxint_cmp(const void *x, const void *y) nogil:
    cdef c_idxint_t *tx = <c_idxint_t*>x
    cdef c_idxint_t *ty = <c_idxint_t*>y
    return safe_intcmp(tx.d, tx.not_none, ty.d, ty.not_none)

cdef void c_idxint_qsort(c_idxint_t *a, const int n):
    c_stdlib.qsort(a, n, sizeof(c_idxint_t), c_idxint_cmp)
'''
