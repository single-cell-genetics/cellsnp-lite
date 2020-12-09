#!/bin/bash
#this script is aimed to compare genotyping results of cellsnp-mode 2 and freebayes on souporcell dataset with ground truth.
#hxj5<hxj5@hku.hk>

function count_snp() {
    local infile=$1
    local bin_cat=zcat
    if [ "${infile%.gz}" == "$infile" ]; then
        bin_cat=cat
    fi
    $bin_cat $infile | grep -v '^#' | wc -l
}

# print usage message of this script. e.g. print_usage test.sh
function print_usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options:"
    echo "  --cellsnp FILE        Vcf of cellsnp-lite genotyping results."
    echo "  --freebayes FILE      Vcf of freebayes genotyping results."
    echo "  --array FILE          Vcf of array genotyping."
    echo "  --perf DIR            Performance Dir"
    echo "  --fix-het-GP          If use, GP of SNPs having different REF/ALT would be set to 0."
    echo "  -O, --out-dir DIR     Directory of outputing files."
    echo "  --rootdir DIR         Project Dir"
    echo "  -h, --help            This message."
    echo
}

# parse command line args
script_name=$0
if [ $# -lt 1 ]; then
    print_usage $script_name
    exit 1
fi

fix_het_gp=0

cmdline=`echo $0 $*`
ARGS=`getopt -o O:h --long cellsnp:,freebayes:,array:,perf:,fix-het-GP,out-dir:,rootdir:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        --cellsnp) csp_vcf=$2; shift 2;;
        --freebayes) fby_vcf=$2; shift 2;;
        --array) arr_vcf=$2; shift 2;;
        --perf) perf_dir=$2; shift 2;;
        --fix-het-GP) fix_het_gp=1; shift;;
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
check_path_exist $csp_vcf "cellsnp vcf"
check_path_exist $fby_vcf "freebayes vcf"
check_path_exist $arr_vcf "array vcf"
check_path_exist $perf_dir "performance dir"
check_arg_null $out_dir "out dir"
safe_mkdir $out_dir

root_dir=`get_abspath_dir $root_dir`
script_dir=`get_abspath_dir $script_dir`
util_dir=$script_dir/utils
out_dir=`get_abspath_dir $out_dir`

# global settings
bin_commit_ver=$script_dir/utils/get_git_last_commit.sh
bin_rscript=$script_dir/bin/Rscript
bin_perf=$script_dir/bm/plot_performance.r
bin_pre=$root_dir/run/bm2/pr/mode2a_pre.sh
bin_python=$script_dir/bin/python
bin_pr=$root_dir/run/bm2/pr/pr.py
bin_pr2=$root_dir/run/bm2/pr/pr2.py
bin_bcftools=$script_dir/bin/bcftools

log_dir=$out_dir/log
safe_mkdir $log_dir
out_log=$log_dir/`basename $script_name`.ncores${ncores//,/-}.out.log
err_log=$log_dir/`basename $script_name`.ncores${ncores//,/-}.err.log

# print the command line.
echo "=> START @`get_now_str`"
echo "=> ABSTRACT this script is aimed to compare genotyping results of cellsnp-mode 2 and freebayes on souporcell dataset with ground truth."
echo "=> COMMAND $cmdline"
echo "=> VERSION python `$bin_python --version`"
echo "=> VERSION Rscript `$bin_rscript --version`"
echo "=> VERSION bcftools `$bin_bcftools --version`"
echo "=> VERSION data dir"
$bin_commit_ver -d $root_dir 2> /dev/null
echo
echo "=> OUTLOG $out_log"
echo "=> ERRLOG $err_log"
echo "=> OUTPUT"
echo

# run each software. 
cat /dev/null > $out_log
cat /dev/null > $err_log

# core part
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

part_aim="plot performance for cellsnp-lite and freebayes"
res_dir=$out_dir/perf2
safe_mkdir $res_dir
log_file=$log_dir/plot.perf2.log
cmd="grep -v '^cellSNP' $perf_file > $res_dir/perf2.txt && \\
     $bin_rscript $bin_perf -i $res_dir/perf2.txt -o $res_dir/perf.summary.tsv \\
       -f $res_dir/perf.tiff --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

echo "Num of SNPs of cellsnp-lite GT vcf: `count_snp $csp_vcf` for $csp_vcf"
echo "Num of SNPs of freebayes GT vcf: `count_snp $fby_vcf` for $fby_vcf"
echo "Num of SNPs of array GT vcf: `count_snp $arr_vcf` for $arr_vcf"
echo "Num of Het SNPs of array GT vcf: `$bin_bcftools query -f '[%GT]\n' $arr_vcf | awk '\$0 == \"0/1\" || \$0 == \"1/0\"' | wc -l` for $arr_vcf"
echo

#chrom  pos  ref_alt  DP|AD  GT  TGT  $gp_type  arr_ref_alt  arr_GT  **arr_TGT**  ra_same  bi_same  het_app  het_arr  GQ  all_GP  het_GP

### Precision-Recall Curve of cellsnp-lite & array (chrom+pos)
part_aim="preprocess genotyping vcf of cellsnp-lite"
csp_dir=$out_dir/cellsnp-lite
safe_mkdir $csp_dir
csp_tsv0=$csp_dir/cellsnp-lite.tsv
if [ $fix_het_gp -eq 0 ]; then
    cmd="$bin_pre --app-name cellsnp-lite --app-vcf $csp_vcf --array $arr_vcf -O $csp_dir --rootdir $root_dir"
else
    cmd="$bin_pre --app-name cellsnp-lite --app-vcf $csp_vcf --array $arr_vcf --fix-het-GP -O $csp_dir --rootdir $root_dir"
fi
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

echo "Num of SNPs shared by cellsnp-lite & array: `count_snp $csp_tsv0` for $csp_tsv0"
echo "Num of Het SNPs having different REF/ALT with cellsnp-lite: \
      `cat $csp_tsv0 | grep -v '^#' | awk '$14 && ! $11' | wc -l` for $csp_tsv0"     # het_arr && ! ra_same
echo

part_aim="filter vcf of cellsnp-lite"
csp_tsv=${csp_tsv0/.tsv/.filter.tsv}
cmd="cat $csp_tsv0 | grep -v '^#' > $csp_tsv"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="plot Precision-Recall curve of all SNPs for cellsnp-lite"
fig=$csp_dir/cellsnp-lite.array.pr.tiff
cmd="cat $csp_tsv | awk '{printf(\"%s\t%s\n\", \$12, \$15)}' > ${csp_tsv%.tsv}.pr.tsv && \\
     $bin_python $bin_pr --name cellsnp-lite --infile ${csp_tsv%.tsv}.pr.tsv --outfig $fig --color red"  # awk: bi_same GQ
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

### Precision-Recall Curve of freebayes & array (chrom+pos)
part_aim="preprocess genotyping vcf of freebayes"
fby_dir=$out_dir/freebayes
safe_mkdir $fby_dir
fby_tsv0=$fby_dir/freebayes.tsv
if [ $fix_het_gp -eq 0 ]; then
    cmd="$bin_pre --app-name freebayes --app-vcf $fby_vcf --array $arr_vcf -O $fby_dir --rootdir $root_dir"
else
    cmd="$bin_pre --app-name freebayes --app-vcf $fby_vcf --array $arr_vcf --fix-het-GP -O $fby_dir --rootdir $root_dir"
fi
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

echo "Num of SNPs shared by freebayes & array: `count_snp $fby_tsv0` for $fby_tsv0"
echo "Num of Het SNPs having different REF/ALT with freebayes: \
      `cat $fby_tsv0 | grep -v '^#' | awk '$14 && ! $11' | wc -l` for $fby_tsv0"    # het_arr && ! ra_same
echo

part_aim="filter freebayes vcf"
fby_tsv=${fby_tsv0/.tsv/.filter.tsv}
cmd=`cat $fby_tsv0 | grep -v '^#' > $fby_tsv`
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="plot Precision-Recall curve of all SNPs for freebayes"
fig=$fby_dir/freebayes.array.pr.tiff
cmd="cat $fby_tsv | awk '{printf(\"%s\t%s\n\", \$12, \$15)}' > ${fby_tsv%.tsv}.pr.tsv && \\
     $bin_python $bin_pr --name freebayes --infile ${fby_tsv%.tsv}.pr.tsv --outfig $fig --color blue"    # awk: bi_same GQ
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

### Precision-Recall Curve of cellsnp-lite & freebayes & array (chrom+pos)
part_aim="get intersect part (chrom+pos) of celsnp, freebayes and array"
share_dir3=$out_dir/share3
safe_mkdir $share_dir3
csp_tsv3=$share_dir3/cellsnp-lite.share.tsv
fby_tsv3=$share_dir3/freebayes.share.tsv
cmd="awk 'ARGIND == 1 {a[\$1\">\"\$2]=1; next} ARGIND == 2 && a[\$1\">\"\$2]' $fby_tsv $csp_tsv > $csp_tsv3 && \\
     awk 'ARGIND == 1 {a[\$1\">\"\$2]=1; next} ARGIND == 2 && a[\$1\">\"\$2]' $csp_tsv $fby_tsv > $fby_tsv3"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

echo "Num of SNPs shared (chrom+pos) by cellsnp-lite, freebayes & array (after filtering): `count_snp $csp_tsv3` for $csp_tsv3"
echo "Num of SNPs shared (chrom+pos) by cellsnp-lite, freebayes & array (after filtering): `count_snp $fby_tsv3` for $fby_tsv3"
echo "Num of Het SNPs (of array share3 (chrom+pos)) having different REF/ALT with cellsnp-lite: \
      `cat $csp_tsv3 | grep -v '^#' | awk '$14 && ! $11' | wc -l` for $csp_tsv3"     # het_arr && ! ra_same
echo "Num of Het SNPs (of array share3 (chrom+pos)) having different REF/ALT with freebayes: \
      `cat $fby_tsv3 | grep -v '^#' | awk '$14 && ! $11' | wc -l` for $fby_tsv3"
echo

part_aim="plot Precision-Recall curve of all SNPs for share part (chrom+pos) of two apps"
fig=$share_dir3/share3.pr.tiff
cmd="cat $csp_tsv3 | awk '{printf(\"%s\t%s\n\", \$12, \$15)}' > ${csp_tsv3%.tsv}.pr.tsv && \\
     cat $fby_tsv3 | awk '{printf(\"%s\t%s\n\", \$12, \$15)}' > ${fby_tsv3%.tsv}.pr.tsv && \\
     $bin_python $bin_pr2 --name1 cellsnp-lite --infile1 ${csp_tsv3%.tsv}.pr.tsv --name2 freebayes \
       --infile2 ${fby_tsv3%.tsv}.pr.tsv --outfig $fig"     # awk: bi_same GQ
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

echo "=> END @`get_now_str`"
