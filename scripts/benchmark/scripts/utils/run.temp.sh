#!/bin/bash
#this script is aimed to <run ...>.
#hxj5<hxj5@hku.hk>

# if conda environment is needed
#source /home/xianjie/.bashrc       # this cmd will activate conda base env.
#eval "$(conda shell.bash hook)"    # this will allow run conda inside shell scripts. from https://github.com/conda/conda/issues/7980
#bin_conda=conda
#conda_env=
#$bin_conda activate $conda_env

work_dir=`cd $(dirname $0); pwd`
root_dir=$work_dir/portal             # path to 'portal' dir
script_dir=$root_dir/scripts
source $script_dir/utils/base_utils.sh

# pathes to data
data_dir=$root_dir/data
bam_file=$data_dir/raw/cnv_GX109-T1c-CNV/outs/possorted_bam.bam
out_file=

# pathes to scripts
bin_commit_ver=$script_dir/utils/get_git_last_commit.sh 

# default settings
region=chr1:1-10000000

# should not change codes below.
cmdline="$0 $*"
echo "=> START @`get_now_str`"
echo "=> ABSTRACT this script is aimed to <run ...>."
echo "=> COMMAND $cmdline"
echo "=> VERSION cellsnp-lite `cellsnp-lite -V`"
# if use conda env
#conda_dir=$work_dir/conda_env
#safe_mkdir $conda_dir
#conda_env_yml=$conda_dir/${conda_env}.environment.yml
#conda_spec_file=$conda_dir/${conda_env}.spec-file.txt
#conda env export > $conda_env_yml
#conda list --explicit > $conda_spec_file
#echo "=> VERSION conda env $conda_env_yml $conda_spec_file"
echo "=> VERSION data dir"
$bin_commit_ver -d $root_dir 2> /dev/null
echo "=> OUTPUT"
echo

# <the first part>
echo "=> PART <the first part> ..."
cmd="<...>"    # e.g. cmd="samtools > $out_log 2> $err_log"
echo "$cmd" && eval $cmd
if [ $? -ne 0 ]; then
    echo "Error: failed to run PART <the first part> S1 ..." >&2
    exit 1
fi
echo
cmd="<...>"
echo "$cmd" && eval $cmd
if [ $? -ne 0 ]; then
    echo "Error: failed to run PART <the first part> S2 ..." >&2
    exit 1
fi
echo

# <the second part>
echo "=> PART <the second part> ..."
cmd="<...>"
echo "$cmd" && eval $cmd
if [ $? -ne 0 ]; then
    echo "Error: failed to run PART <the second part> ..." >&2
    exit 1
fi
echo

#$bin_conda deactivate
echo "=> END @`get_now_str`"
