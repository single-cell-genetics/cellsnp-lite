#!/bin/bash
#this script is aimed to preprocess genotyping results of cellsnp-lite and freebayes.
#hxj5<hxj5@hku.hk>

#@abstract  Preprocess cellsnp-lite vcf or freebayes vcf
#@param $1  App name [str]
#@param $2  Vcf of cellsnp-lite or freebayes [str]
#@param $3  Vcf of array [str]
#@param $4  App info tsv [str]
#@param $5  Out Dir [str]
#@param $6  If to fix het GP [int]
#@return    No RetCode or RetContent
#@example   preprocess csp.vcf cellsnp-lite arr.vcf csp.tsv ./test
function preprocess() {
    # todo: no checking inputs 
    local app=$1
    local app_vcf=$2
    local arr_vcf=$3
    local app_tsv=$4
    local out_dir=$5
    local fix_het_gp=$6

    local bin_cat=zcat
    if [ "${app_vcf%.gz}" == "$app_vcf" ]; then
        bin_cat=cat
    fi

    echo "+ basic stat for $app vcf $app_vcf"
    echo "+ total number of SNPs is `$bin_cat $app_vcf | grep -v '^#' | wc -l`"
    echo "+ number of SNPs for each chrom is"
    $bin_cat $app_vcf | grep -v '^#' | awk '{print $1}' | uniq -c
    echo
    
    local part_aim="intersect $app genotyping vcf with array vcf"
    local app_vcf1=$out_dir/${app}.intersect.vcf.gz
    local cmd=
    if [ -f "${app_vcf}.csi" ] || [ -f "${app_vcf}.tbi" ]; then
        cmd="$bin_bcftools view -R $arr_vcf -Oz $app_vcf > $app_vcf1"
    else
        cmd="$bin_bcftools index $app_vcf && \\
             $bin_bcftools view -R $arr_vcf -Oz $app_vcf > $app_vcf1"
    fi
    eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"
    
    part_aim="intersect array vcf with $app vcf"
    local arr_vcf1=$out_dir/array.${app}.vcf.gz
    if [ -f "${arr_vcf}.csi" ] || [ -f "${arr_vcf}.tbi" ]; then
        cmd="$bin_bcftools view -R $app_vcf1 -Oz $arr_vcf > $arr_vcf1"
    else
        cmd="$bin_bcftools index $arr_vcf && \\
             $bin_bcftools view -R $app_vcf1 -Oz $arr_vcf > $arr_vcf1"
    fi
    eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"
    
    part_aim="extract info from $app vcf"
    local app_tsv0=${app_tsv}.tmp
    if [ "$app" == "cellsnp-lite" ]; then
        cmd="$bin_bcftools query -f '%CHROM\t%POS\t%REF|%ALT[\t%DP|%AD\t%GT\t%TGT\t%PL]\n' $app_vcf1 \
               > $app_tsv0"
        gp_type=PL
    elif [ "$app" == "freebayes" ]; then
        cmd="$bin_bcftools query -f '%CHROM\t%POS\t%REF|%ALT[\t%DP|%AD\t%GT\t%TGT\t%GL]\n' $app_vcf1 \
               > $app_tsv0"
        gp_type=GL
    else
        error_exit "Error: wrong app" 1
    fi
    eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"
    
    part_aim="extract info from array vcf"
    local arr_tsv=$out_dir/array.${app}.tsv
    cmd="$bin_bcftools query -f '%CHROM\t%POS\t%REF|%ALT[\t%GT\t%TGT\t%GC]\n' $arr_vcf1 > $arr_tsv"
    eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"
    
    part_aim="check $app info and array info"
    cmd="diff <(awk '{print \$1, \$2}' $app_tsv0) <(awk '{print \$1, \$2}' $arr_tsv)"
    eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"
    
    part_aim="merge $app info and array info and compare GT"
    bin_gq=$bin_gl2gq
    if [ "$app" == "cellsnp-lite" ]; then
        bin_gq=$bin_pl2gq
    fi
    cmd="paste $app_tsv0 <(awk '{printf(\"%s\t%s\t%s\n\", \$3, \$4, \$5)}' $arr_tsv) | \
           $bin_cmp_gt | $bin_gq -v fix_het_gp=$fix_het_gp > $app_tsv"
    eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"
    #chrom  pos  ref_alt  DP|AD  GT  TGT  $gp_type  arr_ref_alt  arr_GT  arr_TGT  ra_same  bi_same  het_app  het_arr  GQ  all_GP  het_GP
}

# print usage message of this script. e.g. print_usage test.sh
function print_usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options:"
    echo "  --app-name STR        App name (cellsnp-lite or freebayes)."
    echo "  --app-vcf FILE        Vcf of app (cellsnp-lite or freebayes) genotyping."
    echo "  --array FILE          Vcf of array genotyping."
    echo "  --fix-het-GP          If use, the GP of SNPs having different REF/ALT would be set to 0."
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
ARGS=`getopt -o O:h --long app-name:,app-vcf:,array:,fix-het-GP,out-dir:,rootdir:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        --app-name) app_name=$2; shift 2;;
        --app-vcf) app_vcf=$2; shift 2;;
        --array) arr_vcf=$2; shift 2;;
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
check_arg_null $app_name "app name"
check_path_exist $app_vcf "app vcf"
check_path_exist $arr_vcf "array vcf"
check_arg_null $out_dir "out dir"
safe_mkdir $out_dir

root_dir=`get_abspath_dir $root_dir`
script_dir=`get_abspath_dir $script_dir`
out_dir=`get_abspath_dir $out_dir`

# global settings
bin_commit_ver=$script_dir/utils/get_git_last_commit.sh
bin_bcftools=$script_dir/bin/bcftools
bin_bgzip=$script_dir/bin/bgzip
bin_gl2gq=$root_dir/run/bm2/pr/gl2gq.awk
bin_pl2gq=$root_dir/run/bm2/pr/pl2gq.awk
bin_cmp_gt=$root_dir/run/bm2/pr/cmp_gt.awk

log_dir=$out_dir/log
safe_mkdir $log_dir
out_log=$log_dir/`basename $script_name`.ncores${ncores//,/-}.out.log
err_log=$log_dir/`basename $script_name`.ncores${ncores//,/-}.err.log

# print the command line.
echo "=> START @`get_now_str`"
echo "=> ABSTRACT this script is aimed to preprocess genotyping results of cellsnp-lite and freebayes."
echo "=> COMMAND $cmdline"
echo "=> VERSION bcftools `$bin_bcftools --version`"
echo "=> VERSION bgzip `$bin_bgzip --version`"
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

# preprocess genotyping vcf
app_vcf2=$app_vcf
if [ "$app_name" == "freebayes" ]; then
    part_aim="filter non-snp for freebayes"
    app_vcf2=$out_dir/freebayes.filter.vcf.gz
    cmd="cat $app_vcf | awk '\$0 ~ /^#/ || (length(\$4) == 1 && length(\$5) == 1)' | $bin_bgzip -c > $app_vcf2"
    eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"
fi

app_tsv=$out_dir/${app_name}.tsv
preprocess $app_name $app_vcf2 $arr_vcf $app_tsv $out_dir $fix_het_gp

echo "=> END @`get_now_str`"
