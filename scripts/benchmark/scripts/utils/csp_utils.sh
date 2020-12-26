#!/bin/bash
#this script contains several functions for cellSNP data processing.
#hxj5<hxj5@hku.hk>

#abstract   Compare output files, in two dirs, of different runs of cellSNP.
#@param $1  The first dir [STR]
#@param $2  The second dir [STR]
#@return    RetCode: 0 if no difference, 1 if has difference, -1 if error.
#           No RetContent
#example    csp_diff_dir "c-cellSNP-res" "py-cellSNP-res"
function csp_diff_dir() {
    if [ $# -lt 2 -o ! -d "$1" -o ! -d "$2" ]; then
        return -1
    fi
    local dir1=`get_abspath_dir "$1"`
    local dir2=`get_abspath_dir "$2"`
    local flist1=`ls "$dir1" | tr ' ' '\n' | grep '^cellSNP' | grep -v 'cells.vcf'`
    local flist2=`ls "$dir2" | tr ' ' '\n' | grep '^cellSNP' | grep -v 'cells.vcf'`
    local nf1=`echo -e "$flist1" | wc -l`
    local nf2=`echo -e "$flist2" | wc -l`
    if [ $nf1 -ne $nf2 ]; then
        return 1
    fi
    for fn in $flist1; do
        if [ ! -f "$dir2/$fn" ]; then
            return 1
        else
            is_gzip "$fn"
            if [ $? -eq 1 ]; then
                local cmd="diff <(zcat $dir1/$fn) <(zcat $dir2/$fn) &> /dev/null"
            else
                local cmd="diff $dir1/$fn $dir2/$fn &> /dev/null"
            fi
            echo $cmd && eval $cmd
            if [ $? -ne 0 ]; then
                return 1
            fi
        fi
    done
    return 0
}

