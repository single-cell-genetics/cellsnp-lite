
#ifndef SZ_GENERAL_UTIL_H
#define SZ_GENERAL_UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "htslib/kstring.h"        // do not use "kstring.h" as it's different from "htslib/kstring.h"

/*
* Number Functions
 */
#define min2(a, b) ((a) < (b) ? (a) : (b))
#define max2(a, b) ((a) > (b) ? (a) : (b))

/* SZ_NUMERIC_OP functions
@abstract     Provide common functions for numeric operations.
@example

#include <stdio.h>
#include "htslib/kstring.h"
#include "general_util.h"

SZ_NUMERIC_OP_INIT(test, int)

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
#define SZ_NUMERIC_OP_INIT2(SCOPE, name, type)																			\
    SCOPE type get_sum_of_arr_##name(type *a, const int n) {															\
        int i;																											\
        type sum;																										\
        for (i = 0, sum = 0; i < n; i++) sum += a[i];																	\
        return sum;																										\
    }																													\
    SCOPE int get_idx_of_max_##name(type *a, const int n) {																\
        int i, j;																											\
        type max;																									\
        for (i = 1, j = 0, max = a[0]; i < n; i++) {																	\
            if (a[i] > max) j = i;																					\
        }																												\
        return j;																									\
    }																												\
    SCOPE int join_arr_to_str_##name(type *a, const int n, char c, char *fmt, kstring_t *s) {							\
        int i;	/* TODO: to test if the ret of ksxxx functions is right. */												\
        for (i = 0; i < n - 1; i++) { ksprintf(s, fmt, a[i]); kputc_(c, s); }											\
        if (n >= 1) { ksprintf(s, fmt, a[i]); }																		\
        return n;																									\
    }

/*@abstract    Macro to declare the SZ_NUMERIC_OP functions.
@param name    Name of SZ_NUMERIC_OP.
@param type    The numeric type: int, double, size_t etc.

@example       SZ_NUMERIC_OP_INIT(test, int).
 */
#define SZ_NUMERIC_OP_INIT(name, type) SZ_NUMERIC_OP_INIT2(static inline, name, type)

/*@abstract   Calculate sum of a numeric array.
@param name   Name of the SZ_NUMERIC_OP.
@param a      Pointer of the numeric array [type*].
@param n      Size of the array [int].
@return       Sum of the array [type].
 */
#define get_sum_of_arr(name, a, n) get_sum_of_arr_##name(a, n)

/*@abstract   Get the index of the maximum element in a numeric array.
@param name   Name of the SZ_NUMERIC_OP.
@param a      Pointer of the numeric array [type*].
@param n      Size of the array [int].
@return       Index of the maximum element [int].
 */
#define get_idx_of_max(name, a, n) get_idx_of_max_##name(a, n)

/*@abstract   Join elements of a numeric array by a delimiter into string.
@param name   Name of the SZ_NUMERIC_OP.
@param a      Pointer of the numeric array [type*].
@param n      Size of the array [int].
@param c      The delimiter char [int].
@param fmt    The format string for ksprintf/printf [char*].
@param s      Pointer of kstring_t structure.
@return       Num of elements that are successfully joined.
 */
#define join_arr_to_str(name, a, n, c, fmt, s) join_arr_to_str_##name(a, n, c, fmt, s)

/*
* String Functions
 */

/*@abstract   A safe version of strdup() as the parameter of strdup() cannot be NULL.
@param s      Pointer of string to be duplicated [char*].
@return       Pointer of newly duplicated string if success, NULL if the original pointer is NULL [char*].
 */
#define safe_strdup(s) ((s) ? strdup(s) : NULL)

/*@abstract    Free space of a char** array.
@param a       Pointer to the char* array.
@param n       Num of char* elements in the array.
@return        Void.
 */
static inline void str_arr_destroy(char **a, const int n) {
    int i;      
    for (i = 0; i < n; i++) free(a[i]);
    free(a);
}

/*abstract    Join an array of strings (char*) by a delim char.
@param a      Pointer to the char* array.
@param n      Num of char* elements in the array.
@param c      The delimiter char.
@param s      Pointer of the kstring_t structure to save the result.
@return       Num of elements that are successfully joined.
 */
static inline int str_arr_join(char **a, const int n, int c, kstring_t *s) {
    int i;      /* TODO: to test if the ret of ksxxx functions is right. */	
    for (i = 0; i < n - 1; i++) { kputs(a[i], s); kputc_(c, s); }
    if (n >= 1) { kputs(a[i], s); }
    return n;
}

/* 
* File Functions
 */

/*@abstract   Join together two pathes.
@param p1     Pointer to the first path.
@param p2     Pointer to the second path.
@return       Pointer to the joined path if success, NULL otherwise.

@note         Only works for Unix system as the path seperator used in this function is '/'.
 */
static inline char* join_path(const char *p1, const char *p2) {
    kstring_t ks = KS_INITIALIZE;
    kstring_t *s = &ks;
    char *p = NULL;
    int n1;
    if (NULL == p1 || (n1 = strlen(p1)) <= 0 || NULL == p2) { ks_free(s); return NULL; }
    kputs(p1, s);
    if (p1[n1 - 1] != '/') { kputc('/', s); }
    kputs(p2, s); 
    p = strdup(ks_str(s));
    ks_free(s);
    return p;
}

#define TMP_BUFSIZE 1048576
/*@abstract      Merge several files into one.
@param in_fn     Names of input files to be merged.
@param n         Num of input files.
@param out_fn    Name of output file.
@return          Num of files have been merged.
 */
static int merge_files(char **in_fn, const int n, const char *out_fn) {
    char buf[TMP_BUFSIZE];
    FILE *in = NULL, *out = NULL;
    int i = 0;
    size_t m;
    if (NULL == in_fn || NULL == out_fn) { return 0; }
    if (NULL == (out = fopen(out_fn, "w"))) { fprintf(stderr, "[E::%s] could not open '%s'\n", __func__, out_fn); return 0; }
    for (i = 0; i < n; i++) {
        if (NULL == (in = fopen(in_fn[i], "r"))) { fprintf(stderr, "[E::%s] could not open '%s'\n", __func__, in_fn[i]); goto fail; }
        while ((m = fread(buf, 1, TMP_BUFSIZE, in)) > 0) { fwrite(buf, 1, m, out); }
        if (ferror(in)) { goto fail; }
        else { fclose(in); in = NULL; }
    }
    fclose(out);
    return i;
  fail:
    if (in) fclose(in);
    if (out) fclose(out);
    return i;
}

/*@abstract    Output contents of several files to stream.
@param in_fn   Names of input files to be merged.
@param n       Num of input files.
@param out     Pointer of output stream.
@return        Num of files have been merged.

@note          Do not close the out stream inside this function! It should be closed by the caller.
 */
static int merge_files_to_fp(char **in_fn, const int n, FILE *out) {
    char buf[TMP_BUFSIZE];
    FILE *in = NULL;
    int i = 0;
    size_t m;
    if (NULL == in_fn || NULL == out) { return 0; }
    for (i = 0; i < n; i++) {
        if (NULL == (in = fopen(in_fn[i], "r"))) { fprintf(stderr, "[E::%s] could not open '%s'\n", __func__, in_fn[i]); goto fail; }
        while ((m = fread(buf, 1, TMP_BUFSIZE, in)) > 0) { fwrite(buf, 1, m, out); }
        if (ferror(in)) { goto fail; }
        else { fclose(in); in = NULL; }
    }
    return i;
  fail:
    if (in) fclose(in);
    return i;
}

/*@abstract    Output contents of several files to zipped stream.
@param in_fn   Names of input files to be merged.
@param n       Num of input files.
@param out     Pointer of zipped output stream.
@return        Num of files have been merged.

@note          Do not close the out stream inside this function! It should be closed by the caller.
 */
static int merge_files_to_zip_fp(char **in_fn, const int n, gzFile out) {
    char buf[TMP_BUFSIZE];
    FILE *in = NULL;
    int i = 0;
    size_t m;
    if (NULL == in_fn || NULL == out) { return 0; }
    for (i = 0; i < n; i++) {
        if (NULL == (in = fopen(in_fn[i], "r"))) { fprintf(stderr, "[E::%s] could not open '%s'\n", __func__, in_fn[i]); goto fail; }
        while ((m = fread(buf, 1, TMP_BUFSIZE, in)) > 0) { gzwrite(out, buf, m); }
        if (ferror(in)) { goto fail; }
        else { fclose(in); in = NULL; }
    }
    return i;
  fail:
    if (in) fclose(in);
    return i;
}

/*@abstract    Output contents of several zipped files to zipped stream.
@param in_fn   Names of input zipped files to be merged.
@param n       Num of input files.
@param out     Pointer of zipped output stream.
@return        Num of files have been merged.

@note          Do not close the out stream inside this function! It should be closed by the caller.
 */
static int merge_zip_files_to_zip_fp(char **in_fn, const int n, gzFile out) {
    char buf[TMP_BUFSIZE];
    gzFile in = NULL;
    int i = 0;
    size_t m;
    if (NULL == in_fn || NULL == out) { return 0; }
    for (i = 0; i < n; i++) {
        if (NULL == (in = gzopen(in_fn[i], "rb"))) { fprintf(stderr, "[E::%s] could not open '%s'\n", __func__, in_fn[i]); goto fail; }
        while ((m = gzread(in, buf, TMP_BUFSIZE)) > 0) { gzwrite(out, buf, m); }
        gzclose(in); in = NULL;
    }
    return i;
  fail:
    if (in) gzclose(in);
    return i;
}
#undef TMP_BUFSIZE

/*@abstract      Remove several files.
@param fn        Names of files to be removed.
@param n         Num of files.
@return          Num of files have been removed.
 */
static inline int remove_files(char **fn, const int n) {
    int i;
    for (i = 0; i < n; i++) {
        if (0 != access(fn[i], F_OK) || remove(fn[i]) < 0) break;
    }
    return i;
}

/* 
*POOL
 */

/* The SZ_POOL structure is something like kv_t in klib. But the SZ_POOL is much lighter as it provides less APIs.
It's often used when need dynamically allocate and free/reset memory. The main feature (diff from kv_t) of SZ_POOL are:
1. It can automately allocate memory for elements in the pool.
2. It can automately reset the eleements in the pool.

An example:
#include <stdio.h>
#include "general_util.h"

typedef struct _pair { int a, b; } pair;
#define free_pair(p) (free(p))
static inline void reset_pair(pair *p) { p->a = p->b = -1; }

SZ_POOL_INIT(tp, pair, free_pair, reset_pair)

int main(void) {
    sz_pool_t(tp) *pool = sz_pool_init(tp);
    pair *p = sz_pool_get(tp, pool);
    p->a = 1;
    p->b = 2;
    int sum = p->a + p->b;
    printf("sum is %d\n", sum);
    sz_pool_reset(tp, pool);
    p = sz_pool_get(tp, pool);
    sum = p->a + p->b;
    printf("reset: sum is %d\n", sum);
    p->a = 10;
    p->b = 20;
    // do something.
    sz_pool_destroy(tp, pool);
    return 0;
}
 */

/* 
The SZ_POOL_INIT2() macro
@abstract            Define APIs of the SZ_POOL
@param SCOPE         Decoration of the API functions. e.g. static inline.
@param name          Name of the pool.
@param base_type     Basic type of the elements in the pool. The real type of elements would be base_type*.
@param base_free_f   The function used to free base_type*.
@param base_reset_f  The function used to reset base_type*.

The sz_pool_##name##_t
@abstract         Pool structure that can only get elements from while cannot push elements into.
                  The elements in the pool would be pointers of base_type, i.e. base_type*.
@param l          Pos of next element that can be used.
@param n          Size of elements that exist in the pool.
@param m          Total size of the pool.
@param a          Pointer to the array of base_type*.
@param is_reset   If the pool has been reset. 1/0.
*/
#define SZ_POOL_INIT2(SCOPE, name, base_type, base_free_f, base_reset_f)                               			\
    typedef struct {                                                                                                        \
        size_t l, n, m;                                                                                             \
        base_type **a;																								\
        int is_reset;  		                                                                                            \
    } sz_pool_##name##_t;                                                                                                       \
    SCOPE sz_pool_##name##_t* sz_pool_init_##name(void) {                                           							\
        return (sz_pool_##name##_t*) calloc(1, sizeof(sz_pool_##name##_t));                                               \
    }                                                                                                                                       \
    SCOPE void sz_pool_destroy_##name(sz_pool_##name##_t *p) {                                     											\
        size_t k;                                                                                                               \
        for (k = 0; k < p->n; k++) { base_free_f(p->a[k]); }							                         			\
        free(p->a); free(p);                                                                               						\
    }                                                                                                                                       \
    SCOPE base_type* sz_pool_get_##name(sz_pool_##name##_t *p) {														\
        if (p->l < p->n) { 																									\
            base_type *t = p->a[p->l++]; 																					\
            if (p->is_reset) base_reset_f(t); 																					\
            return t; 																											\
        } else if (p->n >= p->m) { 																								\
            if (p->m) { p->m++; kroundup32(p->m); }																				\
            else { p->m = 16; } 																										\
            p->a = (base_type**) realloc(p->a, sizeof(base_type*) * p->m); 														\
        }   																											\
        p->a[p->n++] = (base_type*) calloc(1, sizeof(base_type));    /* after kroundup32: whatif p->m <= p->n */							\
        return p->a[p->l++];    /* assert p->l = p->n */																\
    }																													\
    SCOPE void sz_pool_reset_##name(sz_pool_##name##_t *p) { p->l = 0; p->is_reset = 1; }										\
    SCOPE size_t sz_pool_size_##name(sz_pool_##name##_t *p) { return p->n; }												\
    SCOPE size_t sz_pool_max_##name(sz_pool_##name##_t *p) { return p->m; }												\
    SCOPE size_t sz_pool_used_##name(sz_pool_##name##_t *p) { return p->l; }												\
    SCOPE base_type* sz_pool_A_##name(sz_pool_##name##_t *p, size_t i) { return p->a[i]; }

#define SZ_POOL_INIT(name, base_type, base_free_f, base_reset_f) 														\
        SZ_POOL_INIT2(static inline, name, base_type, base_free_f, base_reset_f)

/*@abstract    Declare a SZ_POOL. 
@param name    Name of the pool.
@return        Void.
*/
#define sz_pool_t(name) sz_pool_##name##_t

/*@abstract    Initialize a named pool
@param name    Name of the pool.
@return        Pointer to the pool if success, NULL otherwise.
*/
#define sz_pool_init(name) sz_pool_init_##name()

/*@abstract    Destroy a named pool. 
@param name    Name of the pool.
@param p       Pointer to the pool.
@return        Void.
*/
#define sz_pool_destroy(name, p) sz_pool_destroy_##name(p)

/*@abstract    Get an available element from the pool.
@param name    Name of the pool.
@param p       Pointer to the pool.
@return        An available element in the pool.
 */
#define sz_pool_get(name, p) sz_pool_get_##name(p)

/*@abstract    Reset the pool without reallocating memory.
@param name    Name of the pool.
@param p       Pointer to the pool.
@return        Void.
 */
#define sz_pool_reset(name, p) sz_pool_reset_##name(p)

/*@abstract  Return the number of elements have been used in the pool.
@param name  Name of the pool.
@param p     Pointer to the pool.
@return      The number of elements have been used in the pool.
 */
#define sz_pool_used(name, p) sz_pool_used_##name(p)

/*@abstract  Return the number of elements exist in the pool.
@param name  Name of the pool.
@param p     Pointer to the pool.
@return      The number of elements exist in the pool.
 */
#define sz_pool_size(name, p) sz_pool_size_##name(p)

/*@abstract  Return total size of the pool.
@param name  Name of the pool.
@param p     Pointer to the pool.
@return      Total size of the pool.
 */
#define sz_pool_max(name, p) sz_pool_max_##name(p)

/*@abstract  Return the element of certain index in the pool.
@param name  Name of the pool.
@param p     Pointer to the pool.
@param i     Index of the element in the pool [size_t].
@return      The element [base_type*].

@note        1. Be careful that the index i must be less than the total size of pool.
             2. Besides, the returned element could be NULL or element that have not been reset, so there are some steps 
                to be taken before the element can be used. If you want to get an element that can be used immediately, 
                please use sz_pool_get() instead.
 */
#define sz_pool_A(name, p, i) sz_pool_A_##name(p, i)
#endif