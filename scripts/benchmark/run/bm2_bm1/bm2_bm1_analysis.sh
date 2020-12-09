#bin/bash 
#this script is aimed to compare results for mode 1 and mode 2
#hxj5<hxj5@hku.hk>

# print usage message of this script. e.g. print_usage test.sh
function print_usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options:"
    echo "  --dir1 DIR         Dir of cellsnp mode 1 results"
    echo "  --dir2 DIR         Dir of cellsnp mode 2 results"
    echo "  --variants FILE    Original target-region variants file"
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
ARGS=`getopt -o O:h --long dir1:,dir2:,variants:,ncores:,out-dir:,rootdir:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        --dir1) dir1=$2; shift 2;;
        --dir2) dir2=$2; shift 2;;
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
check_path_exist $dir1 "cellsnp mode 1 dir"
check_path_exist $dir2 "cellsnp mode 2 dir"
check_path_exist $var_file "original variant file"
check_arg_null $ncores "number of cores"
check_arg_null $out_dir "out dir"
safe_mkdir $out_dir

out_dir=`get_abspath_dir $out_dir`
dir1=`get_abspath_dir $dir1`
dir2=`get_abspath_dir $dir2`
script_dir=$root_dir/scripts
util_dir=$script_dir/utils

bin_commit_ver=$script_dir/utils/get_git_last_commit.sh
bin_python=$script_dir/bin/python
bin_rscript=$script_dir/bin/Rscript
bin_get_ref=$script_dir/bm/get_ref_mtx.r
bin_accu=$script_dir/bm/analysis_accu.sh
bin_csp2mtx=$root_dir/run/bm2_bm1/csp2mtx.sh

log_dir=$out_dir/log
safe_mkdir $log_dir
out_log=$log_dir/`basename $script_name`.ncores${ncores//,/-}.out.log
err_log=$log_dir/`basename $script_name`.ncores${ncores//,/-}.err.log

# print the command line.
echo "=> START @`get_now_str`"
echo "=> ABSTRACT compare results for mode 1 and mode 2"
echo "=> COMMAND $cmdline"
echo "=> VERSION python `$bin_python --version`"
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

part_aim="prepare for plot accuracy: get ref mtx for cellsnp-lite mode 1"
res_dir=$out_dir/accu
safe_mkdir $res_dir
ref_mtx1=$res_dir/cellsnp.m1.tag.ref.mtx
log_file=$log_dir/get_csp_ref1.log
cmd="$bin_rscript $bin_get_ref --ad $dir1/cellSNP.tag.AD.mtx --dp $dir1/cellSNP.tag.DP.mtx \\
       -o $ref_mtx1 --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="prepare for plot accuracy: [csp2mtx] get fixed ref & alt mtx for cellsnp-lite mode 2"
res_dir2=$out_dir/accu/cellsnp-m2
safe_mkdir $res_dir2
log_file=$log_dir/csp2mtx.log
cmd="$bin_csp2mtx --vcf0 $var_file --vcf1 $dir2/cellSNP.cells.vcf.gz \\
       --out-dir $res_dir2 --rootdir $root_dir --ncores $ncores &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="diff ref, alt & oth matrices"
ref_diff_log=$log_dir/ref.diff.log
alt_diff_log=$log_dir/alt.diff.log
oth_diff_log=$log_dir/oth.diff.log
cmd="diff $ref_mtx1 $res_dir2/cellsnp.ref.mtx &> $ref_diff_log && \\
     diff $dir1/cellSNP.tag.AD.mtx $res_dir2/cellsnp.alt.mtx &> $alt_diff_log && \\
     diff $dir1/cellSNP.tag.OTH.mtx $res_dir2/cellsnp.oth.mtx &> $oth_diff_log"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="plot accuracy"
res_dir=$out_dir/accu
safe_mkdir $res_dir
log_file=$log_dir/plot.accu.log
cmd="$bin_accu --ref1 $ref_mtx1 --alt1 $dir1/cellSNP.tag.AD.mtx --name1 cellsnp-lite-m1 \\
       --variant1 $dir1/cellSNP.base.vcf.gz --ref2 $res_dir2/cellsnp.ref.mtx --alt2 $res_dir2/cellsnp.alt.mtx \\
       --name2 cellsnp-lite-m2 --variant2 $res_dir2/cellsnp.variants.tsv -O $res_dir --rootdir $root_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

echo "=> END @`get_now_str`"
