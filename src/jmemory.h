/* Memory management API/routine
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef SZ_JMEMORY_H
#define SZ_JMEMORY_H

#include <stdlib.h>
#include "htslib/kstring.h"        // do not use "kstring.h" as it's different from "htslib/kstring.h"

/* 
*POOL
 */

/* The SZ_POOL structure is something like kv_t in klib. But the SZ_POOL is much lighter as it provides less APIs.
It's often used when need dynamically allocate and free/reset memory. The main feature (diff from kv_t) of SZ_POOL are:
1. It can automately allocate memory for elements in the pool.
2. It can automately reset the eleements in the pool.
@TODO: add init function as parameter.

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
#define SZ_POOL_INIT2(SCOPE, name, base_type, base_free_f, base_reset_f)                      \
    typedef struct {                                                                            \
        size_t l, n, m;                                                                         \
        base_type **a;										  \
        int is_reset;  		                                                                      \
    } sz_pool_##name##_t;                                                                       \
    SCOPE sz_pool_##name##_t* sz_pool_init_##name(void) {                                        \
        return (sz_pool_##name##_t*) calloc(1, sizeof(sz_pool_##name##_t));               \
    }                                                                                         \
    SCOPE void sz_pool_destroy_##name(sz_pool_##name##_t *p) {                                  \
        size_t k;                                                                             \
        for (k = 0; k < p->n; k++) { base_free_f(p->a[k]); }					     \
        free(p->a); free(p);                                                                   \
    }                                                                                       \
    SCOPE base_type* sz_pool_get_##name(sz_pool_##name##_t *p) {				\
        if (p->l < p->n) { 									\
            base_type *t = p->a[p->l++]; 							\
            if (p->is_reset) base_reset_f(t); 							\
            return t; 										\
        } else if (p->n >= p->m) { 								\
            if (p->m) { p->m++; kroundup32(p->m); }						\
            else { p->m = 16; } 									\
            p->a = (base_type**) realloc(p->a, sizeof(base_type*) * p->m); 			\
        }   												\
        p->a[p->n++] = (base_type*) calloc(1, sizeof(base_type));    /* whatif still p->m <= p->n */ \
        return p->a[p->l++];    /* assert p->l = p->n */						\
    }													\
    SCOPE void sz_pool_reset_##name(sz_pool_##name##_t *p) { p->l = 0; p->is_reset = 1; }		\
    SCOPE size_t sz_pool_size_##name(sz_pool_##name##_t *p) { return p->n; }				\
    SCOPE size_t sz_pool_max_##name(sz_pool_##name##_t *p) { return p->m; }				\
    SCOPE size_t sz_pool_used_##name(sz_pool_##name##_t *p) { return p->l; }				\
    SCOPE base_type* sz_pool_A_##name(sz_pool_##name##_t *p, size_t i) { return p->a[i]; }

#define SZ_POOL_INIT(name, base_type, base_free_f, base_reset_f) 					\
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
