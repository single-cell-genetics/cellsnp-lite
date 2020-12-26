#this script contains several common used APIs.
#hxj5<hxj5@hku.hk>

#@abstract    Check if the arg is not NULL.
#@param path  The arg.
#@param name  Name of the arg [STR]
#@return      Void.
check_arg_null <- function(arg, name) {
  if (is.null(arg)) {
    msg <- paste0("Error: ", name, " is needed!\n")
    error_exit(msg, 1)
  }
}

#@abstract    Check if a file/dir exists.
#@param path  The file path [STR]
#@param name  Name of the file [STR]
#@return      Void.
#@note        If file/dir does not exist, then quit the program.
check_path_exists <- function(path, name) {
  if (is.null(path) || 0 == length(path) || ! file.exists(path)) {
    msg <- paste0("Error: ", name, " is invalid!\n")
    error_exit(msg, 1)
  }
}

#@abstract    Echo error message and exit
#@param msg   Error message [str]
#@param code  Exit code [int]
#@return      Void
error_exit <- function(msg, code = 1) {
  write(paste0("[", get_now_str(), "] ", msg), file = stderr())
  quit("no", code)
}

#@abstract  Get abs path of an existing dir.
#@param d   Relative path of dir [STR]
#@return    Abspath of the dir [STR]
get_abspath_dir <- function(d) {
  return(normalizePath(d))
}

#@abstract  Get dirname of current script from command line.
#@param     Void.
#@return    The dirname if succuess, NULL otherwise [STR]
get_dir_of_cmd_script <- function() {
  args <- commandArgs()
  idx <- grep("--file", args)
  if (0 == length(idx)) { return(NULL) }
  res <- strsplit(args[idx], "=", fixed = T)
  if (0 == length(res) || length(res[[1]]) < 2) { return(NULL) }
  fpath <- res[[1]][2]
  return(dirname(fpath))
}

#@abstract  Return the date string of today.
#@param fmt The date format [STR]
#@return    The date string [STR]
get_today_str <- function(fmt = "%Y-%m-%d") {
    return(strftime(Sys.Date(), fmt))
}

#@abstract  Return the time string of now.
#@param fmt The time format [STR]
#@return    The time string [STR]
get_now_str <- function(fmt = "%Y-%m-%d %H:%M:%S") {
    return(strftime(Sys.time(), fmt))
}

#@abstract  Return the command line str.
#@param     Void.
#@return    Command line str [STR]
get_cmdline <- function() {
    return(paste(commandArgs(), collapse = " "))
}

#@abstract  Join path of dir and file.
#@param dir Path of dir [STR]
#@param fn  Name of file [STR]
#@return    Joined path if success, NULL otherwise [STR]
join_path <- function(dir, fn) {
  if (0 == (n <- nchar(dir)) || 0 == nchar(fn)) { return(NULL) }
  sep <- "/"
  if ('/' == substr(dir, n, n)) { sep <- "" }
  return(paste0(dir, sep, fn))
}

safe_mkdir <- function(d) {
  if (! dir.exists(d)) dir.create(d)
}

