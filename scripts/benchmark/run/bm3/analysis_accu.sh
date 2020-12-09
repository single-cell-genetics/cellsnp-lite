#!/bin/bash 
#this script is aimed to compare the accuracy of different apps
#hxj5<hxj5@hku.hk>

# print usage message of this script. e.g. print_usage test.sh
function print_usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options:"
    echo "  --name1 STR        Name of app1"
    echo "  --variant1 FILE    Variant file of app1: vcf or region file"
    echo "  --ref1 FILE        Ref mtx of app1"
    echo "  --alt1 FILE        Alt mtx of app1"
    echo "  --oth1 FILE        OTH mtx of app1"
    echo "  --name2 STR        Name of app2"
    echo "  --variant2 FILE    Variant file of app2: vcf or region file"
    echo "  --ref2 FILE        Ref mtx of app2"
    echo "  --alt2 FILE        Alt mtx of app2"
    echo "  --oth2 FILE        OTH mtx of app2"
    echo "  -O, --out-dir DIR  Directory of outputing files."
    echo "  --rootdir DIR      Path to root dir of this project."
    echo "  -h, --help         This message."
    echo
    echo "Notes:"
    echo "  Region file's first two columns should be <chrom> <pos> and"
    echo "  pos is 1-based."
    echo
}

# parse command line args
script_name=$0
if [ $# -lt 1 ]; then
    print_usage $script_name
    exit 1
fi

cmdline=`echo $0 $*`
ARGS=`getopt -o O:h --long name1:,variant1:,ref1:,alt1:,oth1:,name2:,variant2:,ref2:,alt2:,oth2:,out-dir:,rootdir:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        --name1) name1=$2; shift 2;;
        --variant1) var_file1=$2; shift 2;;
        --ref1) ref_file1=$2; shift 2;;
        --alt1) alt_file1=$2; shift 2;;
        --oth1) oth_file1=$2; shift 2;;
        --name2) name2=$2; shift 2;;
        --variant2) var_file2=$2; shift 2;;
        --ref2) ref_file2=$2; shift 2;;
        --alt2) alt_file2=$2; shift 2;;
        --oth2) oth_file2=$2; shift 2;;
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
check_arg_null "$name1" "name of app1"
check_path_exist $var_file1 "variant file of app1"
check_path_exist $ref_file1 "ref mtx of app1"
check_path_exist $alt_file1 "alt mtx of app1"
check_path_exist $oth_file1 "oth mtx of app1"
check_arg_null "$name2" "name of app2"
check_path_exist $var_file2 "variant file of app2"
check_path_exist $ref_file2 "ref mtx of app2"
check_path_exist $alt_file2 "alt mtx of app2"
check_path_exist $oth_file2 "oth mtx of app2"
check_arg_null $out_dir "out dir"
safe_mkdir $out_dir

root_dir=`get_abspath_dir $root_dir`
script_dir=`get_abspath_dir $script_dir`
out_dir=`get_abspath_dir $out_dir`
util_dir=$script_dir/utils

# global settings
bin_commit_ver=$script_dir/utils/get_git_last_commit.sh
bin_rscript=$script_dir/bin/Rscript
bin_get_ref=$script_dir/bm/get_ref_mtx.r
bin_match_idx=$script_dir/bm/match_snp_index.r
bin_merge_snp=$script_dir/bm/merge_snp.r
bin_plot_accu=$script_dir/bm/plot_accuracy.r
bin_update_idx=$script_dir/bm/update_mtx_snp_index.r

log_dir=$out_dir/log
safe_mkdir $log_dir
out_log=$log_dir/`basename $script_name`.out.log
err_log=$log_dir/`basename $script_name`.err.log

# print the command line.
echo "=> START @`get_now_str`"
echo "=> ABSTRACT this script is aimed to compare the accuracy of different apps"
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

# merge variants
part_aim="Merge variants"
merged_var_file=$out_dir/merged.variants.tsv
log_file=$log_dir/merge_snp.log
cmd="$bin_rscript $bin_merge_snp -1 $var_file1 -2 $var_file2 -o $merged_var_file --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

nvar=`cat $merged_var_file | wc -l`

# match variant index
part_aim="match variant index for app1"
idx_file1=$out_dir/${name1}.match.idx
log_file=$log_dir/${name1}.match.idx.log
cmd="$bin_rscript $bin_match_idx -1 $var_file1 -2 $merged_var_file -o $idx_file1 --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="match variant index for app2"
idx_file2=$out_dir/${name2}.match.idx
log_file=$log_dir/${name2}.match.idx.log
cmd="$bin_rscript $bin_match_idx -1 $var_file2 -2 $merged_var_file -o $idx_file2 --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

# update mtx snp index
part_aim="Update snp index of ref mtx for app1"
new_ref_file1=$out_dir/${name1}.update.ref.mtx
log_file=$log_dir/${name1}.update.ref.mtx.log
cmd="$bin_rscript $bin_update_idx --mtx $ref_file1 --index $idx_file1 -N $nvar -o $new_ref_file1 --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="Update snp index of alt mtx for app1"
new_alt_file1=$out_dir/${name1}.update.alt.mtx
log_file=$log_dir/${name1}.update.alt.mtx.log
cmd="$bin_rscript $bin_update_idx --mtx $alt_file1 --index $idx_file1 -N $nvar -o $new_alt_file1 --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="Update snp index of oth mtx for app1"
new_oth_file1=$out_dir/${name1}.update.oth.mtx
log_file=$log_dir/${name1}.update.oth.mtx.log
cmd="$bin_rscript $bin_update_idx --mtx $oth_file1 --index $idx_file1 -N $nvar -o $new_oth_file1 --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="Update snp index of ref mtx for app2"
new_ref_file2=$out_dir/${name2}.update.ref.mtx
log_file=$log_dir/${name2}.update.ref.mtx.log
cmd="$bin_rscript $bin_update_idx --mtx $ref_file2 --index $idx_file2 -N $nvar -o $new_ref_file2 --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="Update snp index of alt mtx for app2"
new_alt_file2=$out_dir/${name2}.update.alt.mtx
log_file=$log_dir/${name2}.update.alt.mtx.log
cmd="$bin_rscript $bin_update_idx --mtx $alt_file2 --index $idx_file2 -N $nvar -o $new_alt_file2 --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="Update snp index of oth mtx for app2"
new_oth_file2=$out_dir/${name2}.update.oth.mtx
log_file=$log_dir/${name2}.update.oth.mtx.log
cmd="$bin_rscript $bin_update_idx --mtx $oth_file2 --index $idx_file2 -N $nvar -o $new_oth_file2 --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

# plot accuracy
part_aim="plot accuracy"
accu_dir=$out_dir/accu
safe_mkdir $accu_dir
log_file=$log_dir/accu.log
cmd="$bin_rscript $bin_plot_accu --ref1 $new_ref_file1 --alt1 $new_alt_file1 --name1 $name1 --ref2 $new_ref_file2 --alt2 $new_alt_file2 --name2 $name2 -O $accu_dir --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

echo "=> END @`get_now_str`"
