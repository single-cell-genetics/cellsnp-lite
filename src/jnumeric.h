/* Numeric operatoins API/routine
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef SZ_JNUMERIC_H
#define SZ_JNUMERIC_H

#include "htslib/kstring.h"        // do not use "kstring.h" as it's different from "htslib/kstring.h"

/*
* Number Functions
 */
#define min2(a, b) ((a) < (b) ? (a) : (b))
#define max2(a, b) ((a) > (b) ? (a) : (b))

/* JNUMERIC functions
@abstract     Provide common functions for numeric operations.
@example

#include <stdio.h>
#include "htslib/kstring.h"
#include "general_util.h"

JNUMERIC_INIT(test, int)

int main(void) {
    int a[5] = {3, 6, 2, 1, 10};
    int sum, idx;
    kstring_t s = KS_INITIALIZE;
    sum = get_sum_of_arr(test, a, 5);
    idx = get_idx_of_max(test, a, 5);
    if (join_arr_to_str(test, a, 5, ':', "%d", &s) < 5) {
        fprintf(stderr, "Error: cannot join elements to string.\n");
    } else {
        printf("sum = %d, idx_of_max = %d, join_str = %s\n", sum, idx, s.s);
    }
    ks_free(&s);
    return 0;
}
 */
#define JNUMERIC_INIT2(SCOPE, name, type)							        \
    SCOPE type get_sum_of_arr_##name(type *a, const int n) {						\
        int i;												\
        type sum;												\
        for (i = 0, sum = 0; i < n; i++) sum += a[i];							\
        return sum;											\
    }														\
    SCOPE int get_idx_of_max_##name(type *a, const int n) {						\
        int i, j;												\
        type max;												\
        for (i = 1, j = 0, max = a[0]; i < n; i++) {							\
            if (a[i] > max) { max = a[i]; j = i; }							\
        }												\
        return j;										\
    }													\
    SCOPE int join_arr_to_str_##name(type *a, const int n, char c, char *fmt, kstring_t *s) {		\
        int i;	/* TODO: to test if the ret of ksxxx functions is right. */							\
        for (i = 0; i < n - 1; i++) { ksprintf(s, fmt, a[i]); kputc_(c, s); }					\
        if (n >= 1) { ksprintf(s, fmt, a[i]); }										\
        return n;											\
    }

/*@abstract    Macro to declare the JNUMERIC functions.
@param name    Name of JNUMERIC.
@param type    The numeric type: int, double, size_t etc.

@example       JNUMERIC_INIT(test, int).
 */
#define JNUMERIC_INIT(name, type) JNUMERIC_INIT2(static inline, name, type)

/*@abstract   Calculate sum of a numeric array.
@param name   Name of the JNUMERIC.
@param a      Pointer of the numeric array [type*].
@param n      Size of the array [int].
@return       Sum of the array [type].
 */
#define get_sum_of_arr(name, a, n) get_sum_of_arr_##name(a, n)

/*@abstract   Get the index of the maximum element in a numeric array.
@param name   Name of the JNUMERIC.
@param a      Pointer of the numeric array [type*].
@param n      Size of the array [int].
@return       Index of the maximum element [int].
 */
#define get_idx_of_max(name, a, n) get_idx_of_max_##name(a, n)

/*@abstract   Join elements of a numeric array by a delimiter into string.
@param name   Name of the JNUMERIC.
@param a      Pointer of the numeric array [type*].
@param n      Size of the array [int].
@param c      The delimiter char [int].
@param fmt    The format string for ksprintf/printf [char*].
@param s      Pointer of kstring_t structure.
@return       Num of elements that are successfully joined.
 */
#define join_arr_to_str(name, a, n, c, fmt, s) join_arr_to_str_##name(a, n, c, fmt, s)

#endif
