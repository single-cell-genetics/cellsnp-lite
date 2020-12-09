#!/bin/bash 
#this script is aimed to compare different softwares for mode 1 with different cores.
#hxj5<hxj5@hku.hk>

all_tools="cellSNP cellsnp-lite vartrix"
all_tool_str=`echo $all_tools | tr ' ' '|'`

#@abstract  Print usage message of this program
#@param $1  Path to this program [str]
#@return    No RetCode
#@example   print_usage ./test.sh
function print_usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options:"
    echo "  --bam FILE         Input bam file."
    echo "  --snp FILE         Input snp file."
    echo "  --barcode FILE     Input barcode file."
    echo "  --fasta FILE       Input fasta file."
    echo "  -t, --tools STR    Tools, choose from ${all_tool_str}."
    echo "                     Separated by comma."
    echo "  -p, --ncores STR   Number of cores, separatd by comma"
    echo "  -O, --outdir DIR   Path to output dir."
    echo "  --rootdir DIR      Path to root dir of this project."
    echo "  -h, --help         This message."
    echo 
}

### parse command line args
script_name=$0
if [ $# -lt 1 ]; then
    print_usage $script_name
    exit 1
fi

# parse args
cmdline=`echo $0 $*`
ARGS=`getopt -o t:p:O:h --long bam:,snp:,barcode:,fasta:,tools:,ncores:,rootdir:,outdir:,remove-tmp,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        --bam) bam_file=$2; shift 2;;
        --snp) snp_file=$2; shift 2;;
        --barcode) barcode_file=$2; shift 2;;
        --fasta) fasta_file=$2; shift 2;;
        -t|--tools) tools=$2; shift 2;;
        -p|--ncores) ncores=$2; shift 2;;
        -O|--outdir) out_dir=$2; shift 2;;
        --rootdir) root_dir=$2; shift 2;;
        -h|--help) print_usage $script_name; shift; exit 0;;
        --) shift; break;;
        *) echo "Internal error!" >&2; exit 1;;
    esac
done

# check cmdline args
if [ -z $root_dir ] || [ ! -d $root_dir ]; then 
    echo "Error: root_dir invalid!" >&2
    exit 1
fi
script_dir=$root_dir/scripts
source $script_dir/utils/base_utils.sh
check_path_exist $bam_file "bam file"
check_path_exist $snp_file "snp file"
check_path_exist $barcode_file "barcode file"
check_path_exist $fasta_file "fasta file"
check_arg_null $tools "tools"
check_arg_null $ncores "number of cores"
check_arg_null $out_dir "output dir"
safe_mkdir $out_dir

out_dir=`get_abspath_dir $out_dir`
root_dir=`get_abspath_dir $root_dir`
script_dir=`get_abspath_dir $script_dir`

res_dir=$out_dir/result
mkdir $res_dir &> /dev/null
log_dir=$out_dir/log
mkdir $log_dir &> /dev/null
perf_dir=$out_dir/perf
mkdir $perf_dir &> /dev/null
#data_dir=$root_dir/data/bm1

run_name=bm1
bin_mode1=$root_dir/run/bm1/mode1.sh
script_run=$log_dir/${run_name}_`date '+%Y%m%d%H%M%S'`_$$_${RANDOM}.sh
cmd_log=${script_run%.sh}.cmdline.txt

# settings of run_mode1.sh
nrep=3                # repeat times.
intval=3              # interval time between repeats in seconds.
min_mapq=20
no_dup=0
primary_aln=1

# The codes below should not be changed!
# qsub settings
qsub_qname=cgsd          # queue name
qsub_nodes=1
qsub_mem=200          # gb
qsub_walltime=100     # hours

# output command line
echo "command dir = `pwd`" > $cmd_log
echo "command line = $cmdline" >> $cmd_log

# run APPs with different number of cores.
multi_cores=`echo $ncores | tr ',' ' '`
for nc in $multi_cores; do
    perf_file=$perf_dir/perf_ncores${nc}.txt
    out_log=$log_dir/`basename $bin_mode1`.ncores${nc}.out.log
    err_log=$log_dir/`basename $bin_mode1`.ncores${nc}.err.log
    no_dup_opt=""
    if [ $no_dup -eq 1 ]; then no_dup_opt="--no-dup"; fi
    primary_aln_opt=""
    if [ $primary_aln -eq 1 ]; then primary_aln_opt="--primary"; fi
    cmd="qsub -q $qsub_qname -N bm1_ncores${nc} \
        -l nodes=$qsub_nodes:ppn=$[nc + 1],mem=${qsub_mem}gb,walltime=${qsub_walltime}:00:00 \
        -o $out_log -e $err_log -- \\
        $bin_mode1 -O $res_dir -f $perf_file -n $nrep -p $nc -s $intval --bam $bam_file --snp $snp_file \
            --barcode $barcode_file --fasta $fasta_file -t $tools --min-mapq $min_mapq --rootdir $root_dir \
            $no_dup_opt $primary_aln_opt"
    echo >> $script_run
    echo "$cmd" >> $script_run
done

### submit each sample to queue
chmod a+x $script_run
$script_run
