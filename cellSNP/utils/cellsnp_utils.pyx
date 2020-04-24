
from libc.stdint cimport uint32_t
from pysam.libchtslib cimport BAM_CDIFF, BAM_CEQUAL, BAM_CINS, BAM_CMATCH, BAM_CSOFT_CLIP
from pysam.libcalignedsegment cimport AlignedSegment

def get_query_bases(AlignedSegment read, bint full_length=False):
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

def get_query_qualities(AlignedSegment read, bint full_length=False):
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
