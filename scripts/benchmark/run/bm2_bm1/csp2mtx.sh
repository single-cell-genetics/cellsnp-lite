#!/bin/bash
#this script is aimed to fix ref & alt of cellsnp mode 2 vcf and output to cellSNP-style matrices
#hxj5<hxj5@hku.hk>

# print usage message of this script. e.g. print_usage test.sh
function print_usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options:"
    echo "  --vcf0 FILE      Vcf containing target SNPs"
    echo "  --vcf1 FILE      Output vcf of app"
    echo "  --out-dir DIR    Directory of outputing files."
    echo "  --rootdir DIR    Path to root dir of this project."
    echo "  --ncores INT     Number of cores"
    echo "  -h, --help       This message."
    echo
}

# parse command line args
script_name=$0
if [ $# -lt 1 ]; then
    print_usage $script_name
    exit 1
fi

cmdline=`echo $0 $*`
ARGS=`getopt -o h --long vcf0:,vcf1:,out-dir:,rootdir:,ncores:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        --vcf0) vcf0=$2; shift 2;;
        --vcf1) vcf1=$2; shift 2;;
        --out-dir) out_dir=$2; shift 2;;
        --rootdir) root_dir=$2; shift 2;;
        --ncores) ncores=$2; shift 2;;
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
check_path_exist $vcf0 "vcf of target SNPs"
check_path_exist $vcf1 "output vcf of app"
check_arg_null $ncores "number of cores"
check_arg_null $out_dir "out dir"
safe_mkdir $out_dir

vcf0_sort=$out_dir/`basename ${vcf0%.vcf.gz}`.sort.vcf
vcf1_sort=$out_dir/`basename ${vcf1%.vcf.gz}`.sort.vcf

bin_csp2mtx=$root_dir/run/bm2_bm1/csp2mtx2.py
bin_python=$script_dir/bin/python
bin_bcftools=$script_dir/bin/bcftools
bin_commit_ver=$script_dir/utils/get_git_last_commit.sh

### preprocess
log_dir=$out_dir/log
safe_mkdir $log_dir
out_log=$log_dir/`basename $script_name`.out.log
err_log=$log_dir/`basename $script_name`.err.log

# print the command line.
echo "=> START @`get_now_str`"
echo "=> ABSTRACT this script is aimed to fix ref & alt of cellsnp mode 2 vcf and output to cellSNP-style matrices"
echo "=> COMMAND $cmdline"
echo "=> VERSION python `$bin_python --version`"
echo "=> VERSION bcftools `$bin_bcftools --version`"
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

# core part
part_aim="sort vcf0"
cmd="$bin_bcftools sort -m 80G $vcf0 > $vcf0_sort"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="csp2mtx (include filtering)"
cmd="$bin_python $bin_csp2mtx --vcf0 $vcf0_sort --vcf1 $vcf1 --outdir $out_dir"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

echo "=> END @`get_now_str`"
