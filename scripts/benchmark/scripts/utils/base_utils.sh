#!/bin/bash 
#this script contains some common-used functions.
#hxj5<hxj5@hku.hk>

#@abstract  Check if an arg is not empty.
#@param $1  The arg [STR]
#@param $2  Name of the arg [STR]
#@return    RetContent: error message if file does not exist; No RetCode.
#@note      If arg is empty, then exit the program.
#@example   check_arg_null "123" "test"
function check_arg_null() {
    if [ $# -lt 2 ]; then
        error_exit "Error: $1 should be provided!" 1
    fi
}

#@abstract  Check if a path exists.
#@param $1  Path to the file/dir [STR]
#@param $2  Name of the file/dir [STR]
#@return    RetContent: error message if file does not exist; No RetCode.
#@note      If the path does not exist, then exit the program.
#@example   check_path_exist "~/demo/test.txt" "test"
function check_path_exist() {
    if [ $# -lt 2 ]; then
        error_exit "Error: $1 needed!" 1
    elif [ -z "$1" -o ! -e "$1" ]; then
        error_exit "Error: valid $2 needed!" 1
    fi
}

#@abstract  Echo error message and exit
#@param $1  Error message [str]
#@param $2  Exit code [int]
#@return    Void
#@example   error_exit "Error: invalid input!" 1
function error_exit() {
    local error_code=1
    if [ -n $2 ]; then error_code=$2; fi
    echo "[`get_now_str`] $1" >&2
    exit $error_code
}

#@abstract  Execute command and check running status
#@param $1  Command [str]
#@param $2  Aim of this command [str]
#@param $3  Path to out log for the command [str]
#@param $4  Path to err log for the command [str]
#@return    No RetCode
#@note      Will exit the program if error.
#@example   eval_cmd "echo hello world" "test this function" "./out.log" "./err.log"
function eval_cmd() {
    if [ $# -lt 4 ]; then
        error_exit "Error: too few arguments for function eval_cmd" 1
    fi
    local cmd="$1"
    local aim="$2"
    local out_log=$3
    local err_log=$4

    # print info to script out log
    echo "=> $aim"
    echo "==================== START `get_now_str` ===================="
    echo "    $cmd"

    # print info to command out & err log
    echo "=> $aim" >> $out_log
    echo "==================== START `get_now_str` ====================" >> $out_log
    echo "    $cmd" >> $out_log
    echo "=> $aim" >> $err_log

    exec 3>&1
    exec 4>&2
    exec 1>>$out_log
    exec 2>>$err_log

    eval "$cmd"
    local ret=$?

    exec 1>&3
    exec 2>&4

    if [ $ret -ne 0 ]; then
        error_exit "Error: failed to run PART $aim." 1
    fi

    # print info to script out log
    echo "===================== END `get_now_str` ====================="
    echo

    # print info to command out & err log
    echo "===================== END `get_now_str` =====================" >> $out_log
    echo >> $out_log
    echo >> $err_log
}

#@abstract  Get abs path of an existing dir.
#@param $1  Relative path of dir [STR]
#@return    RetContent: abspath of the dir if success, "" otherwise.
#           RetCode: 0 if success, 1 otherwise.
#@example   get_abspath_dir ./test_dir
function get_abspath_dir() {
    if [ $# -lt 1 -o ! -d $1 ]; then
        echo ""
        return 1
    else
        cd "$1"; pwd;
        return 0
    fi
}

#@abstract  Get absolute path to file
#@param $1  Path to file [str]
#@return    RetContent: Absolute path to file if success, otherwise return the original path [str]
#           RetCode: 0 if success, 1 otherwise [int]
#@example   get_abspath_file ./test.sh
function get_abspath_file() {
    if [ $# -lt 1 -o ! -f $1 ]; then 
        echo "" 
        return 1; 
    else 
        echo `realpath $1`
    fi
}

#@abstract  Return the date string of today.
#@param     Void.
#@return    RetContent: today str [STR].
#           RetCode: 0
#@example   get_today_str
function get_today_str() {
    echo `date '+%Y-%m-%d'`
}

#@abstract  Return the time string of now.
#@param     Void.
#@return    RetContent: now str [STR].
#           RetCode: 0
#@example   get_now_str
function get_now_str() {
    echo `date '+%Y-%m-%d %H:%M:%S'`
}

#@abstract  Is an element in a characterized vector.
#@param $1  The element [STR]
#@param $2  The characterized vector, separated by space [STR]
#@return    RetCode: 1/0: 1 for yes, 0 for no [INT]
#           No RetContent.
#@example   is_in hello "hello world"
function is_in() {
    for e in $2; do
        if [ "$1" = "$e" ]; then
            return 1
        fi
    done
    return 0
}

#@abstract  Is the file in gzip format.
#@param $1  Path of the file [STR]
#@return    RetCode: 1 if gzip, 0 otherwise [INT]; No RetContent.
#@example   is_gzip "test.gz"
function is_gzip() {
    if [ -z "$1" -o ${1##*.} != "gz" ]; then
        return 0
    else
        return 1
    fi
}

function safe_mkdir() {
    if [ ! -e "$1" ]; then
        mkdir -p "$1" &> /dev/null
        return 0
    elif [ ! -d "$1" ]; then
        echo "Error: $1 has already existed" >&2
        return 1
    else
        return 0
    fi
}

function to_lower() {
    echo ${1,,}      # support bash 4.+
}

function to_upper() {
    echo ${1^^}      # support bash 4.+
}

