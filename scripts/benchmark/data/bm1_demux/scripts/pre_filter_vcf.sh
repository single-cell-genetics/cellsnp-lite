#!/bin/bash 
#this script is aimed to filter the input SNP vcf before it is loaded by cellSNP or vartrix.
#the indel, biallelic and duplicates will be filtered.
#author:hxj5

script_name=$0
work_dir=`cd $(dirname $0); pwd`
path_base_utils=$work_dir/base_utils.sh

source $path_base_utils

# print usage message of this script. e.g. print_usage test.sh
function print_usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options"
    echo "  -i, --input FILE      The input SNP vcf to be filtered."
    echo "  -o, --output FILE     The output SNP vcf that the filtered records will be written to."
    echo "  --indel FILE          The vcf that the filtered indel SNPs will be written to."
    echo "  --biallele FILE       The vcf that the filtered biallele SNPs will be written to."
    echo "  --duplicate FILE      The vcf that the filtered duplicate SNPs will be written to."
    echo "  -h, --help            This message."
    echo
}

function get_cat_cmd() {
    local cmd=cat
    is_gzip "$1"
    if [ $? -eq 1 ]; then
        cmd=zcat
    fi
    echo $cmd
}

function format_gzip_file() {
    is_gzip "$1"
    if [ $? -eq 1 ]; then
        local tmp_file=$1.$$.tmp
        mv $1 $tmp_file
        gzip $tmp_file
        mv $tmp_file.gz $1
    fi
}

# parse command line args
if [ $# -lt 1 ]; then
    print_usage $script_name
    exit 1
fi

cmdline="$0 $*"
ARGS=`getopt -o i:o:h --long input:,output:,indel:,biallele:,duplicate:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        -i|--input) infile=$2; shift 2;;
        -o|--output) outfile=$2; shift 2;;
        --indel) indel_file=$2; shift 2;;
        --biallele) biallele_file=$2; shift 2;;
        --duplicate) duplicate_file=$2; shift 2;;
        -h|--help) print_usage $script_name; shift; exit 0;;
        --) shift; break;;
        *) echo "Internal error!" >&2; exit 1;;
    esac
done

# check args.
check_path_exist $infile "input file"

# print the command line.
echo COMMAND $cmdline
echo
echo "================BEGIN @`date '+%Y-%m-%d %H:%M:%S'`================"

path_cat_in=`get_cat_cmd $infile`
path_cat_out=`get_cat_cmd $outfile`

echo "There are `$path_cat_in $infile | grep -v '^#' | wc -l` total SNPs in input vcf"

# output the vcf header
echo "Output vcf header ..."
echo
header=`$path_cat_in $infile | sed '/^#/!q' | grep '^#'`
echo -e "$header" > $outfile
if [ -n "$indel_file" ]; then echo -e "$header" > $indel_file; fi
if [ -n "$biallele_file" ]; then echo -e "$header" > $biallele_file; fi
if [ -n "$duplicate_file" ]; then echo -e "$header" > $duplicate_file; fi

# filter indel, biallelic and duplicates.
if [ -n "$indel_file" ]; then
    echo "Filtering indel SNPs ..."
    cmd="$path_cat_in $infile | grep -v '^#' | awk '(length(\$4) > 1 && \$4 !~ /,/) || (length(\$5) > 1 && \$5 !~ /,/)' >> $indel_file"
    echo COMMAND "$cmd" && eval $cmd
    format_gzip_file $indel_file
    path_cat_indel=`get_cat_cmd $indel_file`
    echo "There are `$path_cat_indel $indel_file | grep -v '^#' | wc -l` indel SNPs filtered."
    echo
fi

if [ -n "$biallele_file" ]; then
    echo "Filtering biallele SNPs ..."
    cmd="$path_cat_in $infile | grep -v '^#' | awk 'length(\$5) > 1 && \$5 ~ /,/' >> $biallele_file"
    echo COMMAND "$cmd" && eval $cmd
    format_gzip_file $biallele_file
    path_cat_biallele=`get_cat_cmd $biallele_file`
    echo "There are `$path_cat_biallele $biallele_file | grep -v '^#' | wc -l` biallele SNPs filtered."
    echo
fi

if [ -n "$duplicate_file" ]; then
    echo "Filtering duplicate SNPs ..."
    cmd="$path_cat_in $infile | grep -v '^#' | awk '{k=\$1\"_\"\$2\"_\"\$4\"_\"\$5; if(k in a){print \$0} else{a[k]=1}}' >> $duplicate_file"
    echo COMMAND "$cmd" && eval $cmd
    format_gzip_file $duplicate_file
    path_cat_duplicate=`get_cat_cmd $duplicate_file`
    echo "There are `$path_cat_duplicate $duplicate_file | grep -v '^#' | wc -l` duplicate SNPs filtered."
    echo
fi

echo "Filtering indel, biallelic and duplicates ..."
cmd="$path_cat_in $infile | grep -v '^#' | awk 'length(\$4) == 1 && length(\$5) == 1' | \
            awk '{k=\$1\"_\"\$2\"_\"\$4\"_\"\$5; if(! (k in a)){a[k]=\$0}} END {for(k in a){print a[k]}}' | \
            sort -k1,1V -k2,2n -k4,4 -k5,5 >> $outfile"
echo COMMAND "$cmd" && eval $cmd
format_gzip_file $outfile
echo "There are `$path_cat_out $outfile | grep -v '^#' | wc -l` SNPs left."
echo

echo "================END @`date '+%Y-%m-%d %H:%M:%S'`================"
echo "All Done!"
