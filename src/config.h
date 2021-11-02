/* config.h - Global configure
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef CSP_CONFIG_H
#define CSP_CONFIG_H

#define DEBUG 0
#define VERBOSE 1
/* DEVELOP defined to 1 means some codes for future version of cellsnp will be included. */
#define DEVELOP 0

#define CSP_NAME "cellsnp-lite"
#define CSP_VERSION "1.2.1"
#define CSP_AUTHOR "hxj5<hxj5@hku.hk>"

#define JF_ZIP_TYPE         2  // 1, gzip; 2, bgzip.
#define CSP_FIT_MULTI_SMP   1  // nthread auto fit multi samples? 0, no; 1, yes.

#endif

