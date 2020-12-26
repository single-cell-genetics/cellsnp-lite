#!/usr/bin/env Rscript
#this script is aimed to update mtx with new SNP indexes.
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
library(dplyr)
options(warn = 0)

parser <- ArgumentParser(
  description = "", 
  formatter_class = "argparse.RawTextHelpFormatter"
)
parser$add_argument("--mtx", type = "character", 
                    help = "The mtx file to be updated.")
parser$add_argument("--index", type = "character", 
                    help = "The new SNP indexes.")
parser$add_argument("-N", "--nsnp", type = "integer",
                    help = "Total number of new SNPs, which would be written to mtx header.")
parser$add_argument("-o", "--outfile", type = "character", 
                    help = "Output updated mtx file.")
parser$add_argument("--utilDir", type = "character", default = ".", 
                    help = "The dir containing utils scripts [.]")

args <- parser$parse_args()

# check args.
if (is.null(args$utilDir) || ! dir.exists(args$utilDir)) {
  write("Error: the valid utils Dir needed!", file = stderr())
  quit("no", 1) 
}

old_wd <- getwd()
setwd(args$utilDir)
source("base_utils.r")
source("mtx_utils.r")
setwd(old_wd)

check_path_exists(args$mtx, "mtx file")
check_path_exists(args$index, "index file")
check_arg_null(args$nsnp, "total new SNPs")
check_arg_null(args$outfile, "out file")

# core part
# load data
lmtx <- parse_mtx_file(args$mtx)
dindex <- read.table(args$index, header = F)
colnames(dindex) <- c("new_row")
total_snp <- args$nsnp

if (lmtx$nrow != nrow(dindex)) {
  msg <- "Error: number of SNPs in two files are not the same!"
  msg <- paste0(msg, "\n", paste0("Error: nsnp = ", lmtx$nrow, "; nindex = ", nrow(dindex), ";"))
  error_exit(msg)
}

dindex$row <- 1:nrow(dindex)
dm <- merge(lmtx$data, dindex, by = c("row"), all.x = T, all.y = F, sort = T)

nna <- sum(is.na(dm$new_row))   # Each snp_idx should be within the matching indexes.
if (nna > 0) {
  msg <- paste0("Error: mtx has ", nna, " records that are not in original vcf!")
  error_exit(msg)
}

# output the mtx whose snp_idxes have been updated.
dm <- dm[, c("new_row", "col", "value")]
colnames(dm)[1] <- "row"
lmtx$data <- dm %>%
  arrange(row, col)
lmtx$nrow <- total_snp
write_mtx_file(lmtx, args$outfile)

print("All Done!")
