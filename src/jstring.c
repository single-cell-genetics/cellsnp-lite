/* String operatoins API/routine
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#include <stdlib.h>
#include <string.h>
#include "htslib/kstring.h"        // do not use "kstring.h" as it's different from "htslib/kstring.h"
#include "jstring.h"

/*
* String Functions
 */

void str_arr_destroy(char **a, const int n) {
    int i;      
    for (i = 0; i < n; i++) free(a[i]);
    free(a);
}

int str_arr_join(char **a, const int n, int c, kstring_t *s) {
    int i;      /* TODO: to test if the ret of ksxxx functions is right. */	
    for (i = 0; i < n - 1; i++) { kputs(a[i], s); kputc_(c, s); }
    if (n >= 1) { kputs(a[i], s); }
    return n;
}

