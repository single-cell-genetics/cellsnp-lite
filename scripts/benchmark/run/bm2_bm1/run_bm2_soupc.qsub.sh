#bin/bash 
#this script is aimed to run cellSNP mode 2 on souporcell dataset.
#hxj5<hxj5@hku.hk>
# declare a name for this job
#PBS -N bm2_bm1_pre
# request the queue for this job
#PBS -q cgsd
# request a total of x processors for this job (y nodes and z processors per node)
#PBS -l nodes=1:ppn=25
# request memory
#PBS -l mem=80gb
# specify walltime
#PBS -l walltime=200:00:00
# out log file
#PBS -o run.out.log
# err log file
#PBS -e run.err.log

cd $PBS_O_WORKDIR

root_dir=./portal
if [ ! -d $root_dir ]; then
    echo "Error: $root_dir does not exist." >&2
    exit 1
fi

out_dir=~/projects/csp-bm/result/submit2/bm2_bm1/run
mkdir -p $out_dir &> /dev/null
ncores=25

source $root_dir/scripts/utils/base_utils.sh
root_dir=`get_abspath_dir $root_dir`
safe_mkdir $out_dir
script_name=$0

data_dir=$root_dir/data/bm2_soupc
bam_file=$data_dir/soupc.bam
barcode_file=$data_dir/soupc.barcodes.tsv

out_dir=`get_abspath_dir $out_dir`
script_dir=$root_dir/scripts
bin_cellsnp_lite=$script_dir/bin/cellsnp-lite
bin_commit_ver=$script_dir/utils/get_git_last_commit.sh

log_dir=$out_dir/log
safe_mkdir $log_dir
out_log=$log_dir/`basename $script_name`.ncores${ncores//,/-}.out.log
err_log=$log_dir/`basename $script_name`.ncores${ncores//,/-}.err.log

# print the command line.
cmdline="$0 $*"
echo "=> START @`get_now_str`"
echo "=> ABSTRACT run cellSNP mode 2 on souporcell dataset"
echo "=> COMMAND $cmdline"
echo "=> VERSION cellsnp-lite `$bin_cellsnp_lite -V`"
echo "=> VERSION data dir"
$bin_commit_ver -d $root_dir 2> /dev/null
echo
echo "=> OUTLOG $out_log"
echo "=> ERRLOG $err_log"
echo "=> OUTPUT"
echo

# core part
part_aim="run cellSNP mode 2 on souporcell dataset"
target_chroms="`seq 1 22` X Y"
target_chroms=`echo $target_chroms | tr ' ' ',' | sed 's/,$//'`
m2_out_log=$log_dir/cellsnp-lite.mode2.soupc.out.log
m2_err_log=$log_dir/cellsnp-lite.mode2.soupc.err.log
cmd="$bin_cellsnp_lite -s $bam_file -O $out_dir -b $barcode_file --chrom $target_chroms \\
       --cellTAG CB --UMItag UB --minCOUNT 1 --minMAF 0 --minLEN 0 --minMAPQ 20 \\
       --exclFLAG 772 --inclFLAG 0 -p $ncores --gzip --genotype \\
       > $m2_out_log 2> $m2_err_log"
eval_cmd "$cmd" "$part_aim" "$out_log" "$err_log"

echo "=> END @`get_now_str`"
