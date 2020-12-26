#!/bin/bash 
#this script is aimed to compare output results among different versions of cellSNP.
#exit_code: 0 if no difference, 1 if error or has difference.
#hxj5<hxj5@hku.hk>

script_name=$0
work_dir=`cd $(dirname $0); pwd`
path_base_utils=$work_dir/../utils/base_utils.sh
path_csp_utils=$work_dir/../utils/csp_utils.sh
source $path_base_utils
source $path_csp_utils     # import csp_diff_dir

# print usage message of this script. e.g. print_usage test.sh
function print_usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options"
    echo "  -d, --dir DIR         Dir that contains the subdirs of output results."
    echo "  -h, --help            This message."
    echo
}

# parse command line args
if [ $# -lt 1 ]; then
    print_usage $script_name
    exit 1
fi

cmdline="$0 $*"
ARGS=`getopt -o d:h --long dir:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        -d|--dir) dir0=$2; shift 2;;
        -h|--help) print_usage $script_name; shift; exit 0;;
        --) shift; break;;
        *) echo "Internal error!" >&2; exit 1;;
    esac
done

# check args.
check_path_exist $dir0 "dir"

# print the command line.
echo "================BEGIN @`date '+%Y-%m-%d %H:%M:%S'`================"
echo "=> COMMAND $cmdline"
echo "=> OUTPUT"
echo

cd $dir0
exit_code=0
idx=1
dir1=
for dir2 in `ls | tr ' ' '\n' | grep 'cellSNP'`; do
    if [ ! -d "$dir2" ]; then
        continue
    fi
    if [ $idx -gt 1 ]; then
        echo "=> Comparing $dir1 and $dir2 ..."
        csp_diff_dir "$dir1" "$dir2"
        if [ $? -ne 0 ]; then
            echo "=> WARNING $dir1 and $dir2 have different output!" >&2
            exit_code=1
        else
            echo "=> INFO $dir1 and $dir2 have the same output!"
        fi
        echo
    fi
    dir1=$dir2
    idx=`expr $idx + 1`
done

echo "================END @`date '+%Y-%m-%d %H:%M:%S'`================"
echo "=> All Done!"

exit $exit_code
