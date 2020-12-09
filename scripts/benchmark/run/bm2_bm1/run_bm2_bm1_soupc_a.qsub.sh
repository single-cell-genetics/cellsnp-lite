#bin/bash 
#this script is aimed to compare results of souporcell dataset for mode 1 and mode 2
#hxj5<hxj5@hku.hk>
# declare a name for this job
#PBS -N bm2_bm1a
# request the queue for this job
#PBS -q cgsd
# request a total of x processors for this job (y nodes and z processors per node)
#PBS -l nodes=1:ppn=10
# request memory
#PBS -l mem=200gb
# specify walltime
#PBS -l walltime=300:00:00
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

ncores=10
dir1=~/projects/csp-bm/result/submit2/bm1_soupc/run/result/cellsnp-lite_ncores8_rep1
dir2=~/projects/csp-bm/result/submit2/bm2_bm1/run
vcf0=~/data/cellsnp-bm/bm1_demux/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.noindel.nobiallele.nodup.vcf.gz
out_dir=~/projects/csp-bm/result/submit2/bm2_bm1/analysis

source $root_dir/scripts/utils/base_utils.sh
root_dir=`get_abspath_dir $root_dir`
safe_mkdir $out_dir
out_dir=`get_abspath_dir $out_dir`

cmd="$root_dir/run/bm2_bm1/bm2_bm1_analysis.sh \
  --dir1 $dir1      \
  --dir2 $dir2       \
  --variants $vcf0   \
  --ncores $ncores      \
  -O $out_dir          \
  --rootdir $root_dir"
echo "$cmd" > $out_dir/cmdline.txt && eval $cmd

