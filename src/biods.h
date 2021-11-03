//biods.h

#ifndef __BIODS_H__
#define __BIODS_H__

#include "jmempool.h"

#define pool_str_free(s) free(*(s))
#define pool_str_reset(s) free(*(s))
JMEMPOOL_INIT(str, char*, pool_str_free, pool_str_reset)
typedef jmempool_t(str) pool_str_t;

#endif

