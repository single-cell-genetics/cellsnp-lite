/* Thread operatoins API/routine
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef CSP_THREAD_H
#define CSP_THREAD_H

#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "jfile.h"

/* 
* Thread API
*/
/*@abstract    The data structure used as thread-func parameter.
@param gs      Pointer to the global_settings structure.
@param n       Pos of next element in the snp-list/chrom-list to be used by certain thread.
@param m       Total size of elements to be used by certain thread, must not be changed.
@param i       Id of the thread data.
@param ret     Running state of the thread.
@param ns      Num of SNPs that passed all filters.
@param nr_*    Num of records for each output matrix file. 
@param out_*   Pointers of output files.
 */
typedef struct {
    global_settings *gs;
    size_t m, n;   // for snp-list or chrom-list.
    int i;
    int ret;
    size_t ns, nr_ad, nr_dp, nr_oth;
    jfile_t *out_mtx_ad, *out_mtx_dp, *out_mtx_oth, *out_vcf_base, *out_vcf_cells;
} thread_data;

/*@abstract  Create the thread_data structure.
@return      Pointer to the structure if success, NULL otherwise.
@note        The pointer returned successfully by thdata_init() should be freed
             by thdata_destroy() when no longer used.
 */
static inline thread_data* thdata_init(void);
static inline void thdata_destroy(thread_data *p);
static inline void thdata_print(FILE *fp, thread_data *p);

#endif
