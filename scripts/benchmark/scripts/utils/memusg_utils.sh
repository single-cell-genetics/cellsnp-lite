#!/bin/bash 
#this script contains several functions related to memusg.
#hxj5<hxj5@hku.hk>

#@abstract  Extract time usage from performance file.
#@param $1  Name of performance file [STR]
#@return    RetContent: Extracted time usage in seconds if success, "" otherwise [STR]
#           RetCode: 0 if success, 1 otherwise.
#@example   get_memusage_time perf.txt
function get_memusg_time() {
    if [ $# -lt 1 ]; then
        echo ""
        return 1
    else
        tail $1 | grep -E '^elapsed time: [0-9.]+$' | awk '{print $NF}'
        return 0
    fi
}

#@abstract  Extract memory usage from performance file.
#@param $1  Name of performance file [STR]
#@return    RetContent: Extracted memory usage in KB if success, "" otherwise [STR]
#           RetCode: 0 if success, 1 otherwise.
#@example   get_memusage_mem perf.txt
function get_memusg_mem() {
    if [ $# -lt 1 ]; then
        echo ""
        return 1
    else
        tail $1 | grep -E '^peak rss: [0-9.]+$' | awk '{print $NF}'
        return 0
    fi
}

