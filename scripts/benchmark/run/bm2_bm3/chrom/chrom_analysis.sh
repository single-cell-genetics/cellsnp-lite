#!/bin/bash
#this script is aimed to analysis results of cellsnp-lite mode 2 and bcftools (without input region file).
#hxj5<hxj5@hku.hk>

# print usage message of this script. e.g. print_usage test.sh
function print_usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options:"
    echo "  --cellsnp-dir DIR  Dir of cellsnp results"
    echo "  --bcftools-dir DIR Dir of bcftools results"
    echo "  --perf DIR         Performance dir"
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
ARGS=`getopt -o O:h --long cellsnp-dir:,bcftools-dir:,perf:,ncores:,out-dir:,rootdir:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        --cellsnp-dir) csp_dir=$2; shift 2;;
        --bcftools-dir) bcft_dir=$2; shift 2;;
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
check_path_exist $perf_dir "performance dir"
check_arg_null $ncores "number of cores"
check_arg_null $out_dir "out dir"
safe_mkdir $out_dir

csp_dir=`get_abspath_dir $csp_dir`
bcft_dir=`get_abspath_dir $bcft_dir`
out_dir=`get_abspath_dir $out_dir`
util_dir=$script_dir/utils

bin_bcftools=$script_dir/bin/bcftools
bin_python=$script_dir/bin/python
bin_bcf2depth=$script_dir/bm/bcf2depth.py
bin_csp2depth=$script_dir/bm/csp2depth.py
bin_commit_ver=$script_dir/utils/get_git_last_commit.sh
bin_rscript=$script_dir/bin/Rscript
bin_perf=$script_dir/bm/plot_performance.r

### preprocess 
log_dir=$out_dir/log
safe_mkdir $log_dir
out_log=$log_dir/`basename $script_name`.out.log
err_log=$log_dir/`basename $script_name`.err.log

# print the command line.
echo "=> START @`get_now_str`"
echo "=> ABSTRACT this script is aimed to analysis results of cellsnp-lite mode 2 and bcftools (without input region file)"
echo "=> COMMAND $cmdline"
echo "=> VERSION bcftools `$bin_bcftools --version`"
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
res_dir=$out_dir/result
safe_mkdir $res_dir

chrom_file=$res_dir/chroms.grep.pattern.lst
for ch in `seq 1 22` X Y; do
    echo "^$ch" >> $chrom_file
done

part_aim="merge perf files"
pf_dir=$out_dir/perf
safe_mkdir $pf_dir
perf_file=$pf_dir/perf.txt
cmd="cat $perf_dir/perf_*.txt | awk 'NR == 1 || ! (\$0 ~ /^app/)' > $perf_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="plot performance for all tools"
log_file=$log_dir/plot.perf.log
cmd="$bin_rscript $bin_perf -i $perf_file -o $pf_dir/perf.summary.tsv \\
       -f $pf_dir/perf.tiff --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="plot performance for bcftools and cellsnp-lite"
pf_dir=$out_dir/perf2
safe_mkdir $pf_dir
log_file=$log_dir/plot.perf2.log
cmd="grep -v '^cellSNP' $perf_file > $pf_dir/perf2.txt && \\
     $bin_rscript $bin_perf -i $pf_dir/perf2.txt -o $pf_dir/perf2.summary.tsv \\
       -f $pf_dir/perf2.tiff --utilDir $util_dir &> $log_file"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

# cellsnp-lite target chroms have already been 1-22, X, Y
part_aim="convert cellsnp-lite mode 2 output to allele depth of A,C,G,T for each sample"
csp_in=$res_dir/cellsnp-lite.vcf
csp_out=$res_dir/cellsnp-lite.depth.ACGT.sort.tsv
cmd="zcat $csp_dir/cellSNP.cells.vcf.gz > $csp_in && \\
     $bin_python $bin_csp2depth --vcf $csp_in --outfile ${csp_out}.tmp && \\
     grep -wf $chrom_file ${csp_out}.tmp | sort -k1,1V -k2,2n > $csp_out && \\
     rm ${csp_out}.tmp" 
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="convert cellsnp-lite mode 2 output to allele depth of A,C,G,T,N for each sample"
csp_out2=$res_dir/cellsnp-lite.depth.ACGTN.sort.tsv
cmd="$bin_python $bin_csp2depth --vcf $csp_in --countN --outfile ${csp_out2}.tmp && \\
     grep -wf $chrom_file ${csp_out2}.tmp | sort -k1,1V -k2,2n > $csp_out2 && \\
     rm ${csp_out2}.tmp" 
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="filter bcftools output"
bcft_flt=$res_dir/bcftools.flt.vcf
cmd="$bin_bcftools view -i 'INFO/DP > 0' -V indels --threads $ncores $bcft_dir/bcftools.vcf.gz > $bcft_flt"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="convert bcftools output to allele depth of A,C,G,T for each sample"
bcft_out=$res_dir/bcftools.depth.ACGT.sort.tsv
cmd="$bin_python $bin_bcf2depth --vcf $bcft_flt --outfile ${bcft_out}.tmp && \\
     grep -wf $chrom_file ${bcft_out}.tmp | sort -k1,1V -k2,2n > $bcft_out && \\
     rm ${bcft_out}.tmp"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="convert bcftools output to allele depth of A,C,G,T,N for each sample"
bcft_out2=$res_dir/bcftools.depth.ACGTN.sort.tsv
cmd="$bin_python $bin_bcf2depth --vcf $bcft_flt --countN --outfile ${bcft_out2}.tmp && \\
     grep -wf $chrom_file ${bcft_out2}.tmp | sort -k1,1V -k2,2n > $bcft_out2 && \\
     rm ${bcft_out2}.tmp"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

part_aim="compare depth files of cellsnp-lite mode 2 and bcftools"
diff_log=$res_dir/diff.log
cmd="diff $csp_out $bcft_out > $diff_log"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

echo "=> END @`get_now_str`"
