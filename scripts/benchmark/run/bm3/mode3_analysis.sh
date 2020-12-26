#!/bin/bash
#this script is aimed to analysis results of different apps for mode 3
#hxj5<hxj5@hku.hk>

# print usage message of this script. e.g. print_usage test.sh
function print_usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options:"
    echo "  --cellsnp-dir DIR  Dir of cellsnp results"
    echo "  --bcftools-dir DIR Dir of bcftools results"
    echo "  --variants FILE    Original input variants file"
    echo "  --perf DIR         Performance Dir"
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
ARGS=`getopt -o O:h --long cellsnp-dir:,bcftools-dir:,variants:,perf:,ncores:,out-dir:,rootdir:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        --cellsnp-dir) csp_dir=$2; shift 2;;
        --bcftools-dir) bcft_dir=$2; shift 2;;
        --variants) var_file=$2; shift 2;;
        --perf) perf_dir=$2; shift 2;;
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
check_path_exist $csp_dir "cellsnp dir"
check_path_exist $bcft_dir "bcftools dir"
check_path_exist $var_file "original variant file"
check_path_exist $perf_dir "performance dir"
check_arg_null $ncores "number of cores"
check_arg_null $out_dir "out dir"
safe_mkdir $out_dir

csp_dir=`get_abspath_dir $csp_dir`
bcft_dir=`get_abspath_dir $bcft_dir`
out_dir=`get_abspath_dir $out_dir`
util_dir=$script_dir/utils

bin_rscript=$script_dir/bin/Rscript
bin_get_ref=$script_dir/bm/get_ref_mtx.r
bin_perf=$script_dir/bm/plot_performance.r
bin_accu=$root_dir/run/bm3/analysis_accu.sh
bin_bcf2mtx=$root_dir/run/bm3/bcf2mtx.sh
bin_python=$script_dir/bin/python
bin_csp2mtx=$root_dir/run/bm3/csp2mtx_noN.py
bin_commit_ver=$script_dir/utils/get_git_last_commit.sh

### preprocess 
log_dir=$out_dir/log
safe_mkdir $log_dir
out_log=$log_dir/`basename $script_name`.out.log
err_log=$log_dir/`basename $script_name`.err.log

# print the command line.
echo "=> START @`get_now_str`"
echo "=> ABSTRACT this script is aimed to analysis results of different apps for mode 3"
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

csp_alt_mtx=$csp_dir/cellSNP.tag.AD.mtx
csp_dp_mtx=$csp_dir/cellSNP.tag.DP.mtx
csp_oth_mtx=$csp_dir/cellSNP.tag.OTH.mtx
csp_var_file=$csp_dir/cellSNP.base.vcf.gz
csp_vcf=$csp_dir/cellSNP.cells.vcf.gz

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

part_aim="plot performance for bcftools and cellsnp-lite"
res_dir=$out_dir/perf2
safe_mkdir $res_dir
log_file=$log_dir/plot.perf2.log
cmd="grep -v '^cellSNP' $perf_file > $res_dir/perf2.txt && \\
     $bin_rscript $bin_perf -i $res_dir/perf2.txt -o $res_dir/perf.summary.tsv \\
       -f $res_dir/perf.tiff --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="prepare for plot accuracy: get ref mtx for cellsnp-lite"
res_dir=$out_dir/cellsnp
safe_mkdir $res_dir
csp_ref_mtx=$res_dir/cellSNP.tag.ref.mtx
log_file=$log_dir/get_csp_ref.log
cmd="$bin_rscript $bin_get_ref --ad $csp_alt_mtx --dp $csp_dp_mtx \\
       -o $csp_ref_mtx --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="convert bcftools output to cellSNP-style matrices"
res_dir=$out_dir/bcftools
safe_mkdir $res_dir
bcft_ref_mtx=$res_dir/bcftools.ref.mtx
bcft_alt_mtx=$res_dir/bcftools.alt.mtx
bcft_oth_mtx=$res_dir/bcftools.oth.mtx
bcft_var_file=$res_dir/bcftools.variants.tsv
log_file=$log_dir/bcf2mtx.log
cmd="$bin_bcf2mtx --vcf0 $var_file --vcf1 $bcft_dir/bcftools.vcf.gz \\
       --out-dir $res_dir --rootdir $root_dir --ncores $ncores &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="plot accuracy"
res_dir=$out_dir/accu
safe_mkdir $res_dir
log_file=$log_dir/plot.accu.log
cmd="$bin_accu --ref1 $csp_ref_mtx --alt1 $csp_alt_mtx --oth1 $csp_oth_mtx --name1 cellsnp-lite \\
       --variant1 $csp_var_file --ref2 $bcft_ref_mtx --alt2 $bcft_alt_mtx --oth2 $bcft_oth_mtx --name2 bcftools \\
       --variant2 $bcft_var_file -O $res_dir --rootdir $root_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="not count N for cellsnp-lite"
res_dir=$out_dir/cellsnp-noN
safe_mkdir $res_dir
noN_dir=$res_dir
log_file=$log_dir/csp-noN.log
cmd="$bin_python $bin_csp2mtx --vcf0 $csp_var_file --vcf1 $csp_vcf --outdir $res_dir"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="plot accuracy for noN"
res_dir=$out_dir/accu-noN
safe_mkdir $res_dir
log_file=$log_dir/plot.accu-noN.log
cmd="$bin_accu --ref1 $noN_dir/cellsnp.ref.mtx --alt1 $noN_dir/cellsnp.alt.mtx --oth1 $noN_dir/cellsnp.oth.mtx \\
       --name1 cellsnp-lite-noN --variant1 $noN_dir/cellsnp.variants.tsv --ref2 $bcft_ref_mtx --alt2 $bcft_alt_mtx \\
       --oth2 $bcft_oth_mtx --name2 bcftools --variant2 $bcft_var_file -O $res_dir --rootdir $root_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

#part_aim="diff ref, alt & oth matrices"

echo "=> END @`get_now_str`"
