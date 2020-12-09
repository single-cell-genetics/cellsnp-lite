#!/bin/bash
#this script is aimed to output the last git commit of certain branch dir to file.
#hxj5<hxj5@hku.hk>

exec 3>&1
exec >&2

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
    echo "  -d, --dir DIR         Dir of the branch."
    echo "  -o, --outfile FILE    The output file. If not set, use stdout."
    echo "  -h, --help            This message."
    echo
}

# parse command line args
if [ $# -lt 1 ]; then
    print_usage $script_name
    exit 1
fi

cmdline="$0 $*"
ARGS=`getopt -o d:o:h --long dir:outfile:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        -d|--dir) br_dir=$2; shift 2;;
        -o|--outfile) out_file=$2; shift 2;;
        -h|--help) print_usage $script_name; shift; exit 0;;
        --) shift; break;;
        *) echo "Internal error!" >&2; exit 1;;
    esac
done

# check args.
check_path_exist $br_dir "branch dir"
out_opt=">&3"
if [ -n "$out_file" ]; then
    out_opt="> $out_file"
fi

# print the command line.
echo "================BEGIN @`date '+%Y-%m-%d %H:%M:%S'`================"
echo "=> COMMAND $cmdline"
echo "=> OUTPUT"
echo

cmd="echo -e \"\`cd $br_dir; git log | awk '{if (\$0 ~ /^commit/) n++; if (n > 1) exit; print \$0;}'\`\" $out_opt"
echo "$cmd" && eval $cmd
if [ $? -ne 0 ]; then
    echo "Error: failed to print git log." >&2
    exit 1
fi
echo

echo "=> All Done!"
echo "================END @`date '+%Y-%m-%d %H:%M:%S'`================"
