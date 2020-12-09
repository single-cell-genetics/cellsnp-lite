#!/usr/bin/env Rscript
#this script is aimed to output the REF mtx based on AD and DP mtx.
#hxj5<hxj5@hku.hk>

# parse command line args.
args <- commandArgs(trailingOnly = TRUE)
if (0 == length(args)) {
  print("Welcome!")
  print("use -h or --help for help on argument.")
  quit("no", 1)
}

options(warn = -1)
library(argparse)
options(warn = 0)

parser <- ArgumentParser(
  description = "", 
  formatter_class = "argparse.RawTextHelpFormatter"
)
parser$add_argument("--ad", type = "character", help = "The AD mtx file.")
parser$add_argument("--dp", type = "character", help = "The DP mtx file.")
parser$add_argument("-o", "--outfile", type = "character", help = "The REF mtx file.")
parser$add_argument("--utilDir", type = "character", default = ".", 
                    help = "The util dir [.]")

args <- parser$parse_args()

# check args.
if (is.null(args$utilDir) || ! dir.exists(args$utilDir)) {
  write("Error: the valid util Dir needed!", file = stderr())
  quit("no", 1) 
}

old_wd <- getwd()
setwd(args$utilDir)
source("base_utils.r")
source("mtx_utils.r")
setwd(old_wd)

check_path_exists(args$ad, "AD mtx")
check_path_exists(args$dp, "DP mtx")
check_arg_null(args$outfile, "the output REF mtx file")

# core part
# load ad mtx and dp mtx.
ad <- parse_mtx_file(args$ad)
if (is.null(ad)) {
  error_exit("Error: the input AD file is invalid!")
}

dp <- parse_mtx_file(args$dp)
if (is.null(dp)) {
  error_exit("Error: the input DP file is invalid!")
}

if (ad$nrow != dp$nrow || ad$ncol != dp$ncol || ad$nval > dp$nval) {
  error_exit("Error: the headers of AD and DP mtx are not compatible!")
}

# substract ad from dp to get values of ref.
mdata <- merge(ad$data, dp$data, by = c("row", "col"), all = T, 
               suffixes = c("_ad", "_dp"), sort = T)

nna <- sum(is.na(mdata$value_dp))   # Each records in ad mtx should have a corresponding record in dp mtx. 
if (nna > 0) {
  msg <- paste0("Error: AD mtx has ", nna, " records that are not in DP mtx!")
  error_exit(msg)
}

mdata$value_ad[is.na(mdata$value_ad)] <- 0
mdata$value_ref <- mdata$value_dp - mdata$value_ad
ref_data <- mdata[mdata$value_ref > 0, c("row", "col", "value_ref")]

# output to file.
mtx_ref <- list(
  nrow = ad$nrow,
  ncol = ad$ncol,
  nval = nrow(ref_data),
  data = ref_data
)
write_mtx_file(mtx_ref, args$outfile)

print("All Done!")
