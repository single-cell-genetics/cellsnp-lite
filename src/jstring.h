/* String operatoins API/routine
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef SZ_JSTRING_H
#define SZ_JSTRING_H

#include <stdlib.h>
#include <string.h>
#include "htslib/kstring.h"        // do not use "kstring.h" as it's different from "htslib/kstring.h"

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

#endif
