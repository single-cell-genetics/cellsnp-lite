/* Thread operatoins API/routine
 * Author: Xianjie Huang <hxj5@hku.hk>
 */

#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "jfile.h"
#include "thread.h"

/* 
* Thread API
*/

/*@abstract  Create the thread_data structure.
@return      Pointer to the structure if success, NULL otherwise.
@note        The pointer returned successfully by thdata_init() should be freed
             by thdata_destroy() when no longer used.
 */
inline thread_data* thdata_init(void) { return (thread_data*) calloc(1, sizeof(thread_data)); }

inline void thdata_destroy(thread_data *p) { free(p); }

inline void thdata_print(FILE *fp, thread_data *p) {
    fprintf(fp, "\tm = %ld, n = %ld\n", p->m, p->n);
    fprintf(fp, "\ti = %d, ret = %d\n", p->i, p->ret);
}

