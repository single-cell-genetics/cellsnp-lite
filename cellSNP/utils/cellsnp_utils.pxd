
from pysam.libcalignedsegment cimport AlignedSegment

cdef double c_max(double x, double y)
cdef double c_min(double x, double y)

cdef get_query_bases(AlignedSegment read, bint full_length=*)
cdef get_query_qualities(AlignedSegment read, bint full_length=*)

"""
ctypedef struct c_idxstr_t:
    char *s
    int i
cdef c_idxstr_t* idxstr_arr_init(X, const int n)
cdef void idxstr_arr_destroy(c_idxstr_t *a)
cdef int safe_strcmp(const char *x, const char *y) nogil
cdef int c_idxstr_cmp(const void *x, const void *y) nogil
cdef void c_idxstr_qsort(c_idxstr_t *a, const int n)

ctypedef struct c_idxint_t:
    int d
    int i
    bint not_none
cdef c_idxint_t* idxint_arr_init(X, const int n)
cdef void idxint_arr_destroy(c_idxint_t *a)
cdef int safe_intcmp(const int x, bint x_not_none, const int y, bint y_not_none) nogil
cdef int c_idxint_cmp(const void *x, const void *y) nogil
cdef void c_idxint_qsort(c_idxint_t *a, const int n)
"""