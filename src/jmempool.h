/* jmempool.h - Memory management API/routine
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef SZ_JMEMPOOL_H
#define SZ_JMEMPOOL_H

#include <stdlib.h>
#include "htslib/kstring.h"        // do not use "kstring.h" as it's different from "htslib/kstring.h"

/* 
* MEMPOOL
 */

/* The JMEMPOOL structure is something like kv_t in klib. But the JMEMPOOL is much lighter as it provides less APIs.
It's often used when need dynamically allocate and free/reset memory. The main feature (diff from kv_t) of JMEMPOOL are:
1. It can automately allocate memory for elements in the pool.
2. It can automately reset the eleements in the pool.
@TODO: add init function as parameter.

An example:
#include <stdio.h>
#include "general_util.h"

typedef struct _pair { int a, b; } pair;
#define free_pair(p) (free(p))
static inline void reset_pair(pair *p) { p->a = p->b = -1; }

JMEMPOOL_INIT(tp, pair, free_pair, reset_pair)

int main(void) {
    jmempool_t(tp) *pool = jmempool_init(tp);
    pair *p = jmempool_get(tp, pool);
    p->a = 1;
    p->b = 2;
    int sum = p->a + p->b;
    printf("sum is %d\n", sum);
    jmempool_reset(tp, pool);
    p = jmempool_get(tp, pool);
    sum = p->a + p->b;
    printf("reset: sum is %d\n", sum);
    p->a = 10;
    p->b = 20;
    // do something.
    jmempool_destroy(tp, pool);
    return 0;
}
 */

/* 
The JMEMPOOL_INIT2() macro
@abstract            Define APIs of the JMEMPOOL
@param SCOPE         Decoration of the API functions. e.g. static inline.
@param name          Name of the pool.
@param base_type     Basic type of the elements in the pool. The real type of elements would be base_type*.
@param base_free_f   The function used to free base_type*.
@param base_reset_f  The function used to reset base_type*.

The jmempool_##name##_t
@abstract         Pool structure that can only get elements from while cannot push elements into.
                  The elements in the pool would be pointers of base_type, i.e. base_type*.
@param l          Pos of next element that can be used.
@param n          Size of elements that exist in the pool.
@param m          Total size of the pool.
@param a          Pointer to the array of base_type*.
@param is_reset   If the pool has been reset. 1/0.
*/
#define JMEMPOOL_INIT2(SCOPE, name, base_type, base_free_f, base_reset_f)    \
    typedef struct {                                   \
        size_t l, n, m;                                \
        base_type **a;					  \
        int is_reset;  		                       \
    } jmempool_##name##_t;                             \
    SCOPE jmempool_##name##_t* jmempool_init_##name(void) {                     \
        return (jmempool_##name##_t*) calloc(1, sizeof(jmempool_##name##_t));    \
    }                                                                           \
    SCOPE void jmempool_destroy_##name(jmempool_##name##_t *p) {                 \
        size_t k;                                                                \
        for (k = 0; k < p->n; k++) { base_free_f(p->a[k]); }			 \
        free(p->a); free(p);                                                     \
    }                                                                            \
    SCOPE base_type* jmempool_get_##name(jmempool_##name##_t *p) {		\
        if (p->l < p->n) { 							\
            base_type *t = p->a[p->l++]; 					\
            if (p->is_reset) base_reset_f(t); 					\
            return t; 								\
        } else if (p->n >= p->m) { 						\
            if (p->m) { p->m++; kroundup32(p->m); }				\
            else { p->m = 16; } 						\
            p->a = (base_type**) realloc(p->a, sizeof(base_type*) * p->m); 	\
        }   												\
        p->a[p->n++] = (base_type*) calloc(1, sizeof(base_type));    /* whatif still p->m <= p->n */   \
        return p->a[p->l++];    /* assert p->l = p->n */					\
    }												\
    SCOPE void jmempool_reset_##name(jmempool_##name##_t *p) { p->l = 0; p->is_reset = 1; }	\
    SCOPE size_t jmempool_size_##name(jmempool_##name##_t *p) { return p->n; }			\
    SCOPE size_t jmempool_max_##name(jmempool_##name##_t *p) { return p->m; }			\
    SCOPE size_t jmempool_used_##name(jmempool_##name##_t *p) { return p->l; }			\
    SCOPE base_type* jmempool_A_##name(jmempool_##name##_t *p, size_t i) { return p->a[i]; }

#define JMEMPOOL_INIT(name, base_type, base_free_f, base_reset_f) 			\
        JMEMPOOL_INIT2(static inline, name, base_type, base_free_f, base_reset_f)

#define jmempool_t(name) jmempool_##name##_t

//@return Pointer to the pool if success, NULL otherwise.
#define jmempool_init(name) jmempool_init_##name()

//@return Void.
#define jmempool_destroy(name, p) jmempool_destroy_##name(p)

//@return An available element in the pool.
#define jmempool_get(name, p) jmempool_get_##name(p)

//@return Void.
#define jmempool_reset(name, p) jmempool_reset_##name(p)

//@abstract  Return the number of elements have been used in the pool.
#define jmempool_used(name, p) jmempool_used_##name(p)

//@abstract  Return the number of elements exist in the pool.
#define jmempool_size(name, p) jmempool_size_##name(p)

//@abstract  Return total size of the pool.
#define jmempool_max(name, p) jmempool_max_##name(p)

/*@abstract  Return the element of certain index in the pool.
@param name  Name of the pool.
@param p     Pointer to the pool.
@param i     Index of the element in the pool [size_t].
@return      The element [base_type*].
@note        1. Be careful that the index i must be less than the total size of pool.
             2. Besides, the returned element could be NULL or element that have not been reset, so there are some steps 
                to be taken before the element can be used. If you want to get an element that can be used immediately, 
                please use jmempool_get() instead.
 */
#define jmempool_A(name, p, i) jmempool_A_##name(p, i)

#endif

