#!/bin/bash
#this script is aimed to compare results of mode 2 and mode 3
#hxj5<hxj5@hku.hk>

# print usage message of this script. e.g. print_usage test.sh
function print_usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options:"
    echo "  --dir2 DIR         Dir of cellsnp-lite mode 2 results"
    echo "  --dir3 DIR         Dir of cellsnp-lite mode 3 results"
    echo "  --variants FILE    Original input variants file"
    echo "  --ncores INT       Number of cores"
    echo "  -O, --out-dir DIR  Directory of outputing files."
    echo "  --rootdir DIR      Path to root dir of this project."
    echo "  -h, --help         This message."
    echo
}

# parse command line args
script_name=$0
if [ $# -lt 1 ]; then
    print_usage $script_name
    exit 1
fi

cmdline=`echo $0 $*`
ARGS=`getopt -o O:h --long dir2:,dir3:,variants:,ncores:,out-dir:,rootdir:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        --dir2) dir2=$2; shift 2;;
        --dir3) dir3=$2; shift 2;;
        --variants) var_file=$2; shift 2;;
        --ncores) ncores=$2; shift 2;;
        -O|--out-dir) out_dir=$2; shift 2;;
        --rootdir) root_dir=$2; shift 2;;
        -h|--help) print_usage $script_name; shift; exit 0;;
        --) shift; break;;
        *) echo "Internal error!" >&2; exit 1;;
    esac
done

# check args.
if [ -z "$root_dir" ] || [ ! -d "$root_dir" ]; then
    echo "Error: root_dir invalid!" >&2
    exit 1
fi
script_dir=$root_dir/scripts
source $script_dir/utils/base_utils.sh
check_path_exist $dir2 "dir2"
check_path_exist $dir3 "dir3"
check_path_exist $var_file "original variant file"
check_arg_null $ncores "number of cores"
check_arg_null $out_dir "out dir"
safe_mkdir $out_dir

dir2=`get_abspath_dir $dir2`
dir3=`get_abspath_dir $dir3`
out_dir=`get_abspath_dir $out_dir`
util_dir=$script_dir/utils

bin_rscript=$script_dir/bin/Rscript
bin_get_ref=$script_dir/bm/get_ref_mtx.r
bin_accu=$script_dir/bm/analysis_accu.sh
bin_csp2mtx=$root_dir/run/bm2_bm3/snp/csp2mtx.sh
bin_commit_ver=$script_dir/utils/get_git_last_commit.sh

### preprocess 
log_dir=$out_dir/log
safe_mkdir $log_dir
out_log=$log_dir/`basename $script_name`.out.log
err_log=$log_dir/`basename $script_name`.err.log

# print the command line.
echo "=> START @`get_now_str`"
echo "=> ABSTRACT this script is aimed to compare results of mode 2 and mode 3"
echo "=> COMMAND $cmdline"
echo "=> VERSION Rscript `$bin_rscript --version`"
echo "=> VERSION data dir"
$bin_commit_ver -d $root_dir 2> /dev/null
echo
echo "=> OUTLOG $out_log"
echo "=> ERRLOG $err_log"
echo "=> OUTPUT"
echo

# core part
cat /dev/null > $out_log
cat /dev/null > $err_log

part_aim="for mode 2: convert mode 2 vcf to cellSNP-style matrices"
res_dir2=$out_dir/m2
safe_mkdir $res_dir2
log_file=$log_dir/csp2mtx.log
cmd="$bin_csp2mtx --vcf0 $var_file --vcf1 $dir2/cellSNP.cells.vcf.gz \\
       --out-dir $res_dir2 --rootdir $root_dir --ncores $ncores &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="for mode 3: get ref mtx for cellsnp-lite"
res_dir3=$out_dir/m3
safe_mkdir $res_dir3
ref_mtx3=$res_dir3/cellSNP.ref.m3.mtx
log_file=$log_dir/get_csp_ref3.log
cmd="$bin_rscript $bin_get_ref --ad $dir3/cellSNP.tag.AD.mtx --dp $dir3/cellSNP.tag.DP.mtx \\
       -o $ref_mtx3 --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="diff on ref & alt & alt matrices of mode 2 and mode 3"
ref_log_file=$log_dir/diff.ref.log
alt_log_file=$log_dir/diff.alt.log
oth_log_file=$log_dir/diff.oth.log
cmd="diff $res_dir2/cellsnp.ref.mtx $ref_mtx3 &> $ref_log_file && \\
     diff $res_dir2/cellsnp.alt.mtx $dir3/cellSNP.tag.AD.mtx &> $alt_log_file && \\
     diff $res_dir2/cellsnp.oth.mtx $dir3/cellSNP.tag.OTH.mtx &> $oth_log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="plot accuracy"
res_dir=$out_dir/accu
safe_mkdir $res_dir
log_file=$log_dir/plot.accu.log
cmd="$bin_accu --ref1 $res_dir2/cellsnp.ref.mtx --alt1 $res_dir2/cellsnp.alt.mtx --name1 mode2 \\
       --variant1 $res_dir2/cellsnp.variants.tsv --ref2 $ref_mtx3 --alt2 $dir3/cellSNP.tag.AD.mtx --name2 mode3 \\
       --variant2 $dir3/cellSNP.base.vcf.gz -O $res_dir --rootdir $root_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

echo "=> END @`get_now_str`"
