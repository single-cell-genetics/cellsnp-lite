#!/bin/bash 
#this script is aimed to run different softwares for mode 2.
#hxj5<hxj5@hku.hk>

# default settings
## default global settings
nrep=1                # repeat times.
ncores=1,2            # number of cores to be used, seperated by comma
intval=1              # interval time between repeats in seconds.
all_tools="cellSNP cellSNP-dedup cellSNP-umi cellsnp-lite cellsnp-lite-dedup cellsnp-lite-umi freebayes freebayes-dedup"
all_tool_str=`echo $all_tools | tr ' ' '|'`
def_tools=`echo $all_tools | tr ' ' ','`

# default settings
min_count=1
min_maf=0
min_len=0
min_mapq=20

conda_env=CSP   # for freebayes scripts like freebayes-parallele etc.

#@abstract  Format the performance results to string and write it to performance file.
#@param $1  Name of log file [STR]
#@param $2  Name of performance file [STR]
#@param $3  Name of app [STR]
#@param $4  ncores [INT]
#@param $5  Index of repeat [INT]
#@return    RetCode: 0 if success, 1 otherwise. No RetContent.
#@example   write_perf log.txt perf.txt ccellSNP 10 2
function write_perf() {
    if [ $# -lt 5 ]; then
        return 1
    fi
    local time_used=`get_memusg_time $1`
    if [ -z "$time_used" ]; then
        return 1
    fi
    local mem_used=`get_memusg_mem $1`
    if [ -z "$mem_used" ]; then
        return 1
    fi
    echo -e "$3\t$4\t$5\t$time_used\t$mem_used" >> $2
    return 0
}

#@abstract  Run app
#@param $1  Full name of app [STR]
#@param $2  Command to run [STR]
#@param $3  Number of cores [INT]
#@param $4  Id of repeat [INT]
#@param $5  Result Dir [STR]
#@param $6  Path to perf file [STR]
#@return    No RetCode or RetContent
function run_app() {
    if [ $# -lt 6 ]; then
        error_exit "Error: too few parameters for run_app!" 1
    fi
    local app=$1
    local cmd=$2
    local nc=$3
    local i=$4
    local res_dir=$5
    local perf_file=$6
    part_aim="run $app REP $i for NCORE $nc"
    res_out_log=$res_dir/${app}_ncores${nc}_rep${i}.out.log
    res_err_log=$res_dir/${app}_ncores${nc}_rep${i}.err.log
    cmd="$cmd > $res_out_log 2> $res_err_log"
    eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"
    write_perf $res_err_log $perf_file $app $nc $i
    if [ $? -ne 0 ]; then
        error_exit "Error: failed to write perf results for ${app}."
    fi
}

# print usage message of this script. e.g. print_usage test.sh
function print_usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Input Data:"
    echo "  --bam FILE         Input bam file."
    echo "  --fasta FILE       Input fasta file."
    echo
    echo "Tool Settings:"
    echo "  -t, --tools STR    Choose from $all_tool_str, separated by comma"
    echo "                     [$def_tools]"
    echo "  --min-mapq FLOAT   Min MAPQ [$min_mapq]"
    echo
    echo "Output Options:"
    echo "  -O, --out-dir DIR  Directory of outputing files."
    echo "  -f, --perf FILE    File for performance results."
    echo
    echo "Global Settings:"
    echo "  -n, --repeat INT   Repeat times [$nrep]"
    echo "  -p, --ncores INT   Number of CPUs to be used, seperated by comma [$ncores]"
    echo "  -s, --interval INT Interval for each repeat in seconds [$intval]"
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
ARGS=`getopt -o t:O:f:n:s:p:h --long bam:,fasta:,tools:,min-mapq:,out-dir:,perf:,repeat:,interval:,ncores:,rootdir:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        --bam) bam_file=$2; shift 2;;
        --fasta) fasta_file=$2; shift 2;;
        -t|--tools) tools=$2; shift 2;;
        --min-mapq) min_mapq=$2; shift 2;;
        -O|--out-dir) out_dir=$2; shift 2;;
        -f|--perf) perf_file=$2; shift 2;;
        -n|--repeat) nrep=$2; shift 2;;
        -p|--ncores) ncores=$2; shift 2;;
        -s|--interval) intval=$2; shift 2;;
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
source $script_dir/utils/memusg_utils.sh
check_path_exist $bam_file "bam file"
check_path_exist $fasta_file "fasta file"
check_arg_null $out_dir "out dir"
safe_mkdir $out_dir
check_arg_null "$perf_file" "perf_file"

if [ -z "$tools" ]; then
    tools=$def_tools
fi
tools=`echo $tools | tr ',' ' '`
for tl in $tools; do
    is_in $tl "${all_tools}"
    if [ $? -eq 0 ]; then
        error_exit "Error: invalid tool name '$tl'" 1
    fi
done

root_dir=`get_abspath_dir $root_dir`
script_dir=`get_abspath_dir $script_dir`
out_dir=`get_abspath_dir $out_dir`

# global settings
work_dir=`cd $(dirname $0); pwd`
bin_commit_ver=$script_dir/utils/get_git_last_commit.sh
bin_cellsnp=$script_dir/bin/cellSNP
bin_cellsnp_lite=$script_dir/bin/cellsnp-lite
bin_freebayes=$script_dir/bin/freebayes
bin_memusg=$script_dir/utils/memusg

log_dir=$out_dir/log
safe_mkdir $log_dir
out_log=$log_dir/`basename $script_name`.ncores${ncores//,/-}.out.log
err_log=$log_dir/`basename $script_name`.ncores${ncores//,/-}.err.log

# print the command line.
echo "=> START @`get_now_str`"
echo "=> ABSTRACT this script is aimed to run different softwares for mode 2."
echo "=> COMMAND $cmdline"
echo "=> VERSION cellSNP `$bin_cellsnp`"
echo "=> VERSION cellsnp-lite `$bin_cellsnp_lite -V`"
echo "=> VERSION freebayes: `$bin_freebayes --version`"
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
echo -e "app\tncore\trep\ttime\tmem" > $perf_file
multi_cores=`echo "$ncores" | tr ',' ' '`
chroms="`seq 1 22` X Y"
chroms=`echo $chroms | tr ' ' ',' | sed 's/,$//'`
for nc in $multi_cores; do
    for i in `seq 1 $nrep`; do
        is_in cellSNP-dedup "$tools"
        if [ $? -eq 1 ]; then 
            app=cellSNP-dedup
            res_dir=$out_dir/${app}_ncores${nc}_rep${i}
            safe_mkdir $res_dir
            cmd="$bin_memusg -t -H $bin_cellsnp -s $bam_file -O $res_dir -I Sample0 --chrom $chroms -p $nc \
                     --cellTAG None --UMItag None --minCOUNT $min_count --minMAF $min_maf \
                     --minLEN $min_len --minMAPQ $min_mapq --maxFLAG 255"
            run_app $app "$cmd" $nc $i $res_dir $perf_file
            sleep $intval
        fi

        is_in cellSNP "$tools"
        if [ $? -eq 1 ]; then 
            app=cellSNP
            res_dir=$out_dir/${app}_ncores${nc}_rep${i}
            safe_mkdir $res_dir
            cmd="$bin_memusg -t -H $bin_cellsnp -s $bam_file -O $res_dir -I Sample0 --chrom $chroms -p $nc \
                     --cellTAG None --UMItag None --minCOUNT $min_count --minMAF $min_maf \
                     --minLEN $min_len --minMAPQ $min_mapq --maxFLAG 4096"
            run_app $app "$cmd" $nc $i $res_dir $perf_file
            sleep $intval
        fi

        is_in cellSNP-umi "$tools"
        if [ $? -eq 1 ]; then 
            app=cellSNP-umi
            res_dir=$out_dir/${app}_ncores${nc}_rep${i}
            safe_mkdir $res_dir
            cmd="$bin_memusg -t -H $bin_cellsnp -s $bam_file -O $res_dir -I Sample0 --chrom $chroms -p $nc \
                     --cellTAG None --UMItag UR --minCOUNT $min_count --minMAF $min_maf \
                     --minLEN $min_len --minMAPQ $min_mapq --maxFLAG 4096"
            run_app $app "$cmd" $nc $i $res_dir $perf_file
            sleep $intval
        fi

        is_in cellsnp-lite-dedup "$tools"
        if [ $? -eq 1 ]; then 
            app=cellsnp-lite-dedup
            res_dir=$out_dir/${app}_ncores${nc}_rep${i}
            safe_mkdir $res_dir
            cmd="$bin_memusg -t -H $bin_cellsnp_lite -s $bam_file -O $res_dir -I Sample0 --chrom $chroms -p $nc \
                     --cellTAG None --UMItag None --minCOUNT $min_count --minMAF $min_maf \
                     --minLEN $min_len --minMAPQ $min_mapq --exclFLAG 1796 --inclFLAG 0 --gzip --genotype"
            run_app $app "$cmd" $nc $i $res_dir $perf_file
            sleep $intval
        fi

        is_in cellsnp-lite "$tools"
        if [ $? -eq 1 ]; then 
            app=cellsnp-lite
            res_dir=$out_dir/${app}_ncores${nc}_rep${i}
            safe_mkdir $res_dir
            cmd="$bin_memusg -t -H $bin_cellsnp_lite -s $bam_file -O $res_dir -I Sample0 --chrom $chroms -p $nc \
                     --cellTAG None --UMItag None --minCOUNT $min_count --minMAF $min_maf \
                     --minLEN $min_len --minMAPQ $min_mapq --exclFLAG 772 --inclFLAG 0 --gzip --genotype"
            run_app $app "$cmd" $nc $i $res_dir $perf_file
            sleep $intval
        fi

        is_in cellsnp-lite-umi "$tools"
        if [ $? -eq 1 ]; then 
            app=cellsnp-lite-umi
            res_dir=$out_dir/${app}_ncores${nc}_rep${i}
            safe_mkdir $res_dir
            cmd="$bin_memusg -t -H $bin_cellsnp_lite -s $bam_file -O $res_dir -I Sample0 --chrom $chroms -p $nc \
                     --cellTAG None --UMItag UR --minCOUNT $min_count --minMAF $min_maf \
                     --minLEN $min_len --minMAPQ $min_mapq --exclFLAG 772 --inclFLAG 0 --gzip --genotype"
            run_app $app "$cmd" $nc $i $res_dir $perf_file
            sleep $intval
        fi

        is_in freebayes-dedup "$tools"
        if [ $? -eq 1 ]; then
            app=freebayes-dedup
            res_dir=$out_dir/${app}_ncores${nc}_rep${i}
            safe_mkdir $res_dir
            out_vcf=$res_dir/freebayes.vcf
            script_run=$res_dir/${app}_ncores${nc}_rep${i}.run.sh
            echo "#!/bin/bash" > $script_run
            echo "source ~/.bashrc" >> $script_run
            echo "eval \"\$(conda shell.bash hook)\"" >> $script_run    
            echo "conda activate $conda_env" >> $script_run
            cmd="freebayes-parallel <(fasta_generate_regions.py ${fasta_file}.fai 100000) $nc \
                   -f $fasta_file --genotype-qualities $bam_file > $out_vcf"
            echo "$cmd" >> $script_run
            echo "conda deactivate" >> $script_run
            chmod u+x $script_run
            cmd="$bin_memusg -t -H $script_run"
            run_app $app "$cmd" $nc $i $res_dir $perf_file
            sleep $intval
        fi

        is_in freebayes "$tools"
        if [ $? -eq 1 ]; then
            app=freebayes
            res_dir=$out_dir/${app}_ncores${nc}_rep${i}
            safe_mkdir $res_dir
            out_vcf=$res_dir/freebayes.vcf
            script_run=$res_dir/${app}_ncores${nc}_rep${i}.run.sh
            echo "#!/bin/bash" > $script_run
            echo "source ~/.bashrc" >> $script_run
            echo "eval \"\$(conda shell.bash hook)\"" >> $script_run    
            echo "conda activate $conda_env" >> $script_run
            cmd="freebayes-parallel <(fasta_generate_regions.py ${fasta_file}.fai 100000) $nc \
                   -f $fasta_file --genotype-qualities --use-duplicate-reads $bam_file > $out_vcf"
            echo "$cmd" >> $script_run
            echo "conda deactivate" >> $script_run
            chmod u+x $script_run
            cmd="$bin_memusg -t -H $script_run"
            run_app $app "$cmd" $nc $i $res_dir $perf_file
            sleep $intval
        fi

    done
done

echo "=> END @`get_now_str`"

