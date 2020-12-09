#!/usr/bin/env Rscript
#this script is aimed to merge SNPs of two files.
#both two files should be in the format: the first two columns are CHROM and POS.
#Using VCFs as inputs also works well.
#hxj5<hxj5@hku.hk>

#@abstract  Get the first two columns
#@param fn  Filename [STR]
#@return    A dataframe containing the target columns if success, NULL otherwise.
get_cols <- function(fn) {
  t1 <- read.table(fn, header = F, sep = "\t", comment.char = "#", nrows = 3)
  if (nrow(t1) < 1) { return(NULL) }
  nc <- ncol(t1)
  if (nc < 2) { return(NULL) }
  if (nc <= 5) {
    d1 <- read.table(fn, header = F, sep = "\t", comment.char = "#")
  } else {
    col_class <- c("character", "integer", "NULL", "character", "character", rep("NULL", nc - 5))
    d1 <- read.table(fn, header = F, sep = "\t", comment.char = "#", colClasses = col_class)
  }
  colnames(d1) <- c("chrom", "pos", "ref", "alt")
  return(d1)
}

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
parser$add_argument("-1", "--file1", type = "character", 
                    help = "The first file.")
parser$add_argument("-2", "--file2", type = "character", 
                    help = "The second file.")
parser$add_argument("-o", "--outfile", type = "character", 
                    help = "The merged file.")
parser$add_argument("--utilDir", type = "character", default = ".", 
                    help = "The dir containing base utils scripts [.]")

args <- parser$parse_args()

# check args.
if (is.null(args$utilDir) || ! dir.exists(args$utilDir)) {
  write("Error: the valid utils Dir needed!", file = stderr())
  quit("no", 1) 
}

old_wd <- getwd()
setwd(args$utilDir)
source("base_utils.r")
setwd(old_wd)

check_path_exists(args$file1, "file1")
check_path_exists(args$file2, "file2")
check_arg_null(args$outfile, "out file")

# core part
# load two files.
d1 <- get_cols(args$file1)
d2 <- get_cols(args$file2)

# merge positions of file1 and file2.
dm <- merge(d1, d2, by = c("chrom", "pos", "ref", "alt"), 
            all = T, sort = T)  

# output merged positions
write.table(dm, file = args$outfile, sep = "\t", quote = F, row.names = F, col.names = F, append = F)
print("All Done!")
