#!/bin/bash
#this script is aimed to analysis results of different apps for mode 1
#hxj5<hxj5@hku.hk>

# print usage message of this script. e.g. print_usage test.sh
function print_usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options:"
    echo "  --cellsnp-dir DIR  Dir of cellsnp results"
    echo "  --vartrix-dir DIR  Dir of vartrix results"
    echo "  --variants FILE    Original input variants file"
    echo "  --perf DIR         Performance Dir"
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
ARGS=`getopt -o O:h --long cellsnp-dir:,vartrix-dir:,variants:,perf:,out-dir:,rootdir:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        --cellsnp-dir) csp_dir=$2; shift 2;;
        --vartrix-dir) vtx_dir=$2; shift 2;;
        --variants) var_file=$2; shift 2;;
        --perf) perf_dir=$2; shift 2;;
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
check_path_exist $csp_dir "cellsnp dir"
check_path_exist $vtx_dir "vartrix dir"
check_path_exist $var_file "original variant file"
check_path_exist $perf_dir "performance dir"
check_arg_null $out_dir "out dir"
safe_mkdir $out_dir

csp_dir=`get_abspath_dir $csp_dir`
vtx_dir=`get_abspath_dir $vtx_dir`
out_dir=`get_abspath_dir $out_dir`
util_dir=$script_dir/utils

bin_rscript=$script_dir/bin/Rscript
bin_get_ref=$script_dir/bm/get_ref_mtx.r
bin_perf=$script_dir/bm/plot_performance.r
bin_accu=$script_dir/bm/analysis_accu.sh
bin_commit_ver=$script_dir/utils/get_git_last_commit.sh

### preprocess 
log_dir=$out_dir/log
safe_mkdir $log_dir
out_log=$log_dir/`basename $script_name`.out.log
err_log=$log_dir/`basename $script_name`.err.log

# print the command line.
echo "=> START @`get_now_str`"
echo "=> ABSTRACT this script is aimed to analysis results of different apps for mode 1"
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

csp_alt_mtx=$csp_dir/cellSNP.tag.AD.mtx
csp_dp_mtx=$csp_dir/cellSNP.tag.DP.mtx
csp_var_file=$csp_dir/cellSNP.base.vcf.gz
vtx_ref_mtx=$vtx_dir/ref.mtx
vtx_alt_mtx=$vtx_dir/alt.mtx
vtx_var_file=$vtx_dir/variants.tsv

part_aim="merge perf files"
res_dir=$out_dir/perf
safe_mkdir $res_dir
perf_file=$res_dir/perf.txt
cmd="cat $perf_dir/perf_*.txt | awk 'NR == 1 || ! (\$0 ~ /^app/)' > $perf_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="plot performance for all tools"
res_dir=$out_dir/perf
safe_mkdir $res_dir
log_file=$log_dir/plot.perf.log
cmd="$bin_rscript $bin_perf -i $perf_file -o $res_dir/perf.summary.tsv \\
       -f $res_dir/perf.tiff --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="plot performance for cellsnp-lite and vartrix"
res_dir=$out_dir/perf2
safe_mkdir $res_dir
log_file=$log_dir/plot.perf2.log
cmd="grep -v '^cellSNP' $perf_file > $res_dir/perf2.txt && \\
     $bin_rscript $bin_perf -i $res_dir/perf2.txt -o $res_dir/perf.summary.tsv \\
       -f $res_dir/perf.tiff --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="prepare for plot accuracy: get ref mtx for cellsnp-lite"
res_dir=$out_dir/accu
safe_mkdir $res_dir
csp_ref_mtx=$res_dir/cellSNP.tag.ref.mtx
log_file=$log_dir/get_csp_ref.log
cmd="$bin_rscript $bin_get_ref --ad $csp_alt_mtx --dp $csp_dp_mtx \\
       -o $csp_ref_mtx --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="compare vartrix variants list with original input variants"
log_file=$log_dir/vartrix.variants.diff.log
cmd="diff <(cat $vtx_var_file | awk -F'_' '{printf(\"%s\t%s\n\", \$1, \$2 + 1)}') \\
          <(zcat $var_file | grep -v '^#' | awk '{printf(\"%s\t%s\n\", \$1, \$2)}') &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

vtx_var_file=$var_file

part_aim="plot accuracy"
res_dir=$out_dir/accu
safe_mkdir $res_dir
log_file=$log_dir/plot.accu.log
cmd="$bin_accu --ref1 $csp_ref_mtx --alt1 $csp_alt_mtx --name1 cellsnp-lite \\
       --variant1 $csp_var_file --ref2 $vtx_ref_mtx --alt2 $vtx_alt_mtx --name2 vartrix \\
       --variant2 $vtx_var_file -O $res_dir --rootdir $root_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

echo "=> END @`get_now_str`"
