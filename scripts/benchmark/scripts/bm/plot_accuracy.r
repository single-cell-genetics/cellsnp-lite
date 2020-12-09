#!/usr/bin/env Rscript
#this script is aimed to visualize the results of comparing the output coverage mtx of two Apps.
#hxj5<hxj5@hku.hk>

#ref1:tbl
#  row col value
#
#ref_bs:list
#  $union:tbl
#    row col value.x value.y
#  $uniq_x_idx:tbl
#    row
#  $uniq_y_idx:tbl
#    row
#  $overlap_idx:tbl
#    row
#
#mcmp:tbl
#  row st_type st_value mtx_type
#
#tags:tbl
#  row tag_type tag_value app
#
#aly_overlap:tbl
#  row st_type st_value mtx_type tag_type tag_value app


COR_OF_SD0 <- -2    # special correlation value when one of the vector has sd = 0

print2 <- function(msg) {
  print(paste0("[", get_now_str(), "] ", msg))
}

preprocess_mtx <- function(mtx, grp_col = 1) {
  mtx$data <- mtx$data %>% 
    as_tibble() %>%
    filter(value > 0)   # necessary!
  mtx$nsnp <- mtx$nrow
  mtx$nsmp <- mtx$ncol
  mtx$nrec <- mtx$nval
  if (grp_col != 1) {
    mtx$nrow <- mtx$nsmp
    mtx$ncol <- mtx$nsnp
    mtx$data <- mtx$data %>%
      rename(row = col, col = row)
  }
  return(mtx)
}

#@abstract    Basic stat for two mtx.
#@param mtx1  The first mtx to be compared, returned by parse_mtx_file() [list]
#@param mtx2  The second mtx to be compared, returned by parse_mtx_file() [list]
#@return      Void
basic_stat <- function(mtx1, mtx2) {
  if (mtx1$nsnp != mtx2$nsnp || mtx1$nsmp != mtx2$nsmp) {
    error_exit("Error: invalid headers of mtx files!")
  }
  print(paste0("nsnp = ", mtx1$nsnp, "; nsmp = ", mtx1$nsmp, "; nrec1 = ", 
               mtx1$nrec, "; nrec2 = ", mtx2$nrec))
  
  union <- mtx1$data %>%
    full_join(mtx2$data, by = c("row", "col")) %>%
    mutate(value.x = ifelse(is.na(value.x), 0, value.x),
           value.y = ifelse(is.na(value.y), 0, value.y))
  print(paste0("uniq rows for union of mtx1 and mtx2 = ", 
               nrow(union %>% select(row) %>% distinct())))
  
  uniq_x_idx <- mtx1$data %>%
    select(row) %>%
    distinct() %>%
    anti_join(mtx2$data %>% select(row), by = "row")
  print(paste0("uniq rows only for mtx1 = ", nrow(uniq_x_idx)))
  
  uniq_y_idx <- mtx2$data %>%
    select(row) %>%
    distinct() %>%
    anti_join(mtx1$data %>% select(row), by = "row")  
  print(paste0("uniq rows only for mtx2 = ", nrow(uniq_y_idx)))
  
  overlap_idx <- mtx1$data %>%
    select(row) %>%
    distinct() %>%
    semi_join(mtx2$data %>% select(row), by = "row")
  print(paste0("uniq rows for overlap of mtx1 and mtx2 = ", nrow(overlap_idx)))
  
  return(list(
    union = union,
    uniq_x_idx = uniq_x_idx,
    uniq_y_idx = uniq_y_idx,
    overlap_idx = overlap_idx
  ))
}

#@abstract  A safe version of cor()
#@param v1  Vector 1 [vector]
#@param v2  Vector 2 with the same length of vector 1 [vector]
#@return    The correlation value [DOUBLE]
safe_cor <- function(v1, v2) {
  sd1 <- sd(v1)
  sd2 <- sd(v2)
  if (sd1 == 0 || sd2 == 0) {
    i <- sample(1:length(v1), 1)
    if (sd1 == 0) {
      v1[i] <- v1[i] + 1e-6
    }
    if (sd2 == 0) {
      v2[i] <- v2[i] + 1e-6
    }
    return(cor(v1, v2))
  } else {
    return(cov(v1, v2) / (sd1 * sd2))
  }
}

#@abstract  A safe version of cor()
#@param v1  Vector 1 [vector]
#@param v2  Vector 2 with the same length of vector 1 [vector]
#@return    The correlation value [DOUBLE]
#@note      When either vector's sd is 0, then return -2.
safe_cor2 <- function(v1, v2) {
  sd1 <- sd(v1)
  sd2 <- sd(v2)
  if (sd1 == 0 || sd2 == 0) {
    return(COR_OF_SD0)
  } else {
    return(cov(v1, v2) / (sd1 * sd2))
  }
}

safe_mae <- function(v1, v2) {
  ave <- (v1 + v2) / 2
  ave[ave == 0] <- 1     # here if ave = 0, then v1 = v2 = 0.
  res <- mean(abs(v1 - v2) / ave)
  return(res)
}

#@abstract    Compare each SNPs within two mtx.
#@param mmtx  Merged two ref mtx or merged two alt mtx [tbl]
#@param lc    Length of col elements of each SNP [INT]
#@param type  Name of mtx type [STR]
#@return      The stat results [tbl]
cmp_mtx <- function(mmtx, lc, type) {
  tb <- tibble(col = 1:lc)
  res <- mmtx %>%
    group_by(row) %>%
    group_modify(~ {
      .x %>%
        right_join(tb, by = "col") %>%
        mutate(value.x = ifelse(is.na(value.x), 0, value.x),
               value.y = ifelse(is.na(value.y), 0, value.y)) %>%
        summarise(cor = safe_cor2(value.x, value.y), mae = safe_mae(value.x, value.y))
    }) %>%
    ungroup() %>%
    gather(st_type, st_value, -row) %>%
    mutate(mtx_type = type)
  return(res)
}

#@abstract   Get values of ref, ad, dp values of each SNP.
#@param ref  The ref mtx [tbl]
#@param alt  The alt mtx [tbl]
#@return     The mtx with a SNP per row [tbl]
get_tag_values <- function(ref, alt) {
  dp <- ref %>%
    full_join(alt, by = c("row", "col")) %>%
    mutate(value.x = ifelse(is.na(value.x), 0, value.x),
           value.y = ifelse(is.na(value.y), 0, value.y)) %>%
    group_by(row) %>%
    summarise(r = sum(value.x), a = sum(value.y)) %>%
    ungroup() %>%
    mutate(d = r + a) %>%
    gather(tag_type, tag_value, -row)
  return(dp)
}

#@abstract    Left merge mtx with tag values.
#@param mtx   Mtx [tbl]
#@param tags  Values of tags, one SNP per row [tbl]
#@return      Merged results [tbl]
left_merge_tags <- function(mtx, tags) {
  res <- mtx %>%
    left_join(tags, by = "row")
  return(res)
}

options(warn = -1)
library(argparse)
options(warn = 0)

# default settings
stat_types <- c("cor", "mae")
all_fig_types <- c("boxplot", "density")
def_fig_types <- paste(all_fig_types, collapse = ",")
all_range <- c("overlap", "union")
def_range <- "overlap"
def_grp_col <- 1
def_width <- 8
def_height <- 6
def_dpi <- 300

min_cor <- 0.9
max_mae <- 0.1
cor_unit_range <- 0.1
mae_unit_range <- 0.1

# parse command line args.
args <- commandArgs(trailingOnly = TRUE)
if (0 == length(args)) {
  print("Welcome!")
  print("use -h or --help for help on argument.")
  quit("no", 1)
}

parser <- ArgumentParser(
  description = "", 
  formatter_class = "argparse.RawTextHelpFormatter"
)
parser$add_argument("--ref1", type = "character", help = "Ref matrix file of app 1.")
parser$add_argument("--alt1", type = "character", help = "Alt matrix file of app 1.")
parser$add_argument("--name1", type = "character", help = "Name of app 1.")
parser$add_argument("--ref2", type = "character", help = "Ref matrix file of app 2.")
parser$add_argument("--alt2", type = "character", help = "Alt matrix file of app 2.")
parser$add_argument("--name2", type = "character", help = "Name of app 2.")
parser$add_argument("-O", "--outdir", type = "character",
                    help = "Outdir for result summary files.")
parser$add_argument("-f", "--outfig", type = "character", default = def_fig_types, 
                    help = paste0("Result figure file types: boxplot|density, separated by comma [", def_fig_types, "]"))
parser$add_argument("--groupcol", type = "integer", default = def_grp_col,
                    help = paste0("Col for grouping: 1|2 [", def_grp_col, "]"))
parser$add_argument("--range", type = "character", default = def_range,
                    help = paste0("Range of merged mtx: overlap|union [", def_range, "]"))
parser$add_argument("--title", help = "If set, will add title to plots.")
parser$add_argument("--width", type = "double", default = def_width, 
                    help = paste0("Result file width [", def_width, "]"))
parser$add_argument("--height", type = "double", default = def_height, 
                    help = paste0("Result file height [", def_height, "]"))
parser$add_argument("--dpi", type = "integer", default = def_dpi, 
                    help = paste0("DPI [", def_dpi, "]"))
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

check_path_exists(args$ref1, "ref1 file")
check_path_exists(args$alt1, "alt1 file")
check_path_exists(args$ref2, "ref2 file")
check_path_exists(args$alt2, "alt2 file")
check_arg_null(args$outdir, "out dir")
if (! dir.exists(args$outdir)) {
  dir.create(args$outdir)
}
check_arg_null(args$name1, "name1")
check_arg_null(args$name2, "name2")

fig_types <- args$outfig
if (is.null(args$outfig) || 0 == length(args$outfig) || "" == args$outfig) {
  fig_types <- def_fig_types
}
fig_types <- strsplit(fig_types, ",")[[1]]
if (0 == length(fig_types) || any(! fig_types %in% all_fig_types)) {
  error_exit("Error: invalid fig types!")
}

grp_col <- args$groupcol
if (is.null(grp_col) || ! grp_col %in% c(1, 2)) {
  error_exit("Error: wrong group col!")
}

if (is.null(args$range) || ! args$range %in% all_range) {
  error_exit("Error: invalid range!")
}

# load data
options(warn = -1)
library(stringr)
library(dplyr, warn.conflicts = F)
library(tidyr)
library(ggplot2)
options(warn = 0)

print2("loading data ...")
ref1 <- parse_mtx_file(args$ref1)
check_arg_null(ref1, "ref1")

ref2 <- parse_mtx_file(args$ref2)
check_arg_null(ref2, "ref2")

alt1 <- parse_mtx_file(args$alt1)
check_arg_null(alt1, "alt1")

alt2 <- parse_mtx_file(args$alt2)
check_arg_null(alt2, "alt2")

ref1 <- preprocess_mtx(ref1, grp_col)
ref2 <- preprocess_mtx(ref2, grp_col)
alt1 <- preprocess_mtx(alt1, grp_col)
alt2 <- preprocess_mtx(alt2, grp_col)
write("", file = stdout())

# basic statistics
print2("basic stat for ref mtx files ...")
ref_bs <- basic_stat(ref1, ref2)
write("", file = stdout())

print2("basic stat for alt mtx files ...")
alt_bs <- basic_stat(alt1, alt2)
write("", file = stdout())

print2(paste0("Query range is '", args$range, "'"))
write("", file = stdout())

# compare two mtx.
print2("Comparing ref & alt mtx files ...")
if (args$range == "overlap") {
  mcmp <- rbind(
    cmp_mtx(ref_bs$union %>% 
              semi_join(ref_bs$overlap_idx, by = "row"), ref1$ncol, "r"),
    cmp_mtx(alt_bs$union %>%
              semi_join(alt_bs$overlap_idx, by = "row"), alt1$ncol, "a")
  )
} else if (args$range == "union") {
  mcmp <- rbind(
    cmp_mtx(ref_bs$union, ref1$ncol, "r"),
    cmp_mtx(alt_bs$union, alt1$ncol, "a")
  )
} else {
  error_exit("Error: invalid range!")
}
mcmp
mcmp_file <- join_path(args$outdir, "mcmp.tsv")
write.table(mcmp, mcmp_file, quote = F, sep = "\t", row.names = F)
print2(paste0("The mcmp file is saved to ", mcmp_file))
write("", file = stdout())

print2("Total uniq rows: ")
total_uniq <- mcmp %>% 
  group_by(mtx_type) %>%
  distinct(row) %>%
  summarise(total_rows = n()) %>%
  ungroup()
total_uniq
write("", file = stdout())

print2("The ratio of uniq rows in each range for cor:")
cor_flt_std <- tibble(
  st_flt = rep(c(0.99, 0.95, 0.9, 0.8), 2),
  mtx_type = rep(c("r", "a"), each = 4)
)
cor_range_ratio <- cor_flt_std %>%
  group_by(st_flt) %>%
  group_modify(~ {
    mcmp %>%
      group_by(mtx_type) %>%
      filter(st_type == "cor" & st_value >= .y$st_flt[1]) %>%
      distinct(row) %>%
      summarise(nrows = n()) %>%
      ungroup()
  }) %>%
  ungroup() %>%
  right_join(cor_flt_std, by = c("mtx_type", "st_flt")) %>%
  mutate(nrows = ifelse(is.na(nrows), 0, nrows)) %>%
  full_join(total_uniq, by = "mtx_type") %>%
  mutate(row_ratio = nrows / total_rows)
cor_range_ratio
cor_range_ratio_file <- join_path(args$outdir, "cor_range_ratio.tsv")
write.table(cor_range_ratio, cor_range_ratio_file, quote = F, sep = "\t", row.names = F)
print2(paste0("The cor_range_ratio is saved to ", cor_range_ratio_file))
write("", file = stdout())

print2("The ratio of uniq rows in each range for mae:")
mae_flt_std <- tibble(
  st_flt = rep(c(0.01, 0.05, 0.1, 0.2), 2),
  mtx_type = rep(c("r", "a"), each = 4)
)
mae_range_ratio <- mae_flt_std %>%
  group_by(st_flt) %>%
  group_modify(~ {
    mcmp %>%
      group_by(mtx_type) %>%
      filter(st_type == "mae" & st_value < .y$st_flt[1]) %>%
      distinct(row) %>%
      summarise(nrows = n()) %>%
      ungroup()
  }) %>%
  ungroup() %>%
  right_join(mae_flt_std, by = c("mtx_type", "st_flt")) %>%
  mutate(nrows = ifelse(is.na(nrows), 0, nrows)) %>%
  full_join(total_uniq, by = "mtx_type") %>%
  mutate(row_ratio = nrows / total_rows)
mae_range_ratio
mae_range_ratio_file <- join_path(args$outdir, "mae_range_ratio.tsv")
write.table(mae_range_ratio, mae_range_ratio_file, quote = F, sep = "\t", row.names = F)
print2(paste0("The mae_range_ratio is saved to ", mae_range_ratio_file))
write("", file = stdout())

# visualize the result of comparing
for (st in stat_types) {
  pdata <- mcmp %>% filter(st_type == st)
  st_name <- ifelse(st == "cor", "Correlation", "Mean Absolute Error")
  for (ft in fig_types) {
    fig_path <- join_path(args$outdir, paste0(ft, "_grpcol_", grp_col, "_", st, 
                    "_", args$range, ".tiff"))
    if ("boxplot" == ft) {
      p <- pdata %>%
        ggplot(aes(x = mtx_type, y = st_value)) +
        geom_boxplot() +
        scale_x_discrete(labels = c("r" = "Ref", "a" = "Alt"),
                         limits = c("a", "r")) +
        labs(x = "Mtx Type", y = st_name)
    } else if ("density" == ft) {
      p <- pdata %>%
        ggplot(aes(x = st_value, color = mtx_type)) +
        geom_density(fill = "transparent") +
        scale_colour_discrete(labels = c("r" = "Ref", "a" = "Alt"),
                              limits = c("a", "r")) +
        labs(x = st_name, y = "Density", color = "Mtx Type")
    }
    if (! is.null(args$title)) {
        p <- p + 
          labs(title = paste0("Comparison between the output results of ", args$name1, " and ", args$name2))
    }
    if (grepl("tiff?$", fig_path, perl = TRUE, ignore.case = TRUE)) {
      ggsave(fig_path, p, width = args$width, height = args$height, dpi = args$dpi, compress="lzw")
    } else {
      ggsave(fig_path, p, width = args$width, height = args$height, dpi = args$dpi)
    }
    msg <- paste0("the result of stat:", st, "; fig:", ft, "; is saved to ", fig_path)
    print2(msg)
  }
}
write("", file = stdout())

# the ratio of SD = 0
print2("The ratio of SD = 0")
sd0_ratio <- mcmp %>%
  filter(st_type == "cor" & st_value == COR_OF_SD0) %>%
  group_by(mtx_type) %>%
  distinct(row) %>%
  summarise(nrows = n()) %>%
  ungroup() %>%
  full_join(total_uniq, by = "mtx_type") %>%
  mutate(nrows = ifelse(is.na(nrows), 0, nrows)) %>%
  mutate(row_ratio = nrows / total_rows)
sd0_ratio
sd0_ratio_file <- join_path(args$outdir, "sd0_ratio.tsv")
write.table(sd0_ratio, sd0_ratio_file, quote = F, sep = "\t", row.names = F)
print2(paste0("The sd0_ratio is saved to ", sd0_ratio_file))
write("", file = stdout())

#if (args$range == "union") {
#  print("All Done!")
#  quit("no", 0)
#}
#
## analyse the bad points.
#print("Analyse bad points ...")
#tags <- rbind(
#  get_tag_values(ref1$data, alt1$data) %>% mutate(app = "1"),
#  get_tag_values(ref2$data, alt2$data) %>% mutate(app = "2")
#)
#aly_overlap <- left_merge_tags(mcmp, tags)
#
## detailed analysis
## x-y:mae_range-dp
#st_range <- tibble(
#  st_type = c("cor", "mae"),
#  st_unit = c(cor_unit_range, mae_unit_range)
#)
#aly_range_overlap <- aly_overlap %>%
#  left_join(st_range, by = "st_type") %>%
#  mutate(st_upper = (floor(st_value / st_unit) + 1) * st_unit) %>%  # will generate [) ranges
#  mutate(st_upper = ifelse(st_upper > 1, 1, st_upper)) %>%
#  mutate(st_upper = as.character(st_upper)) %>%
#  select(-st_value, -st_unit)
#
#for (st in stat_types) {
#  fig_path <- join_path(args$outdir, paste0(st, "_grpcol_", grp_col, "_aly_",
#                                            args$range, ".tiff"))
#  p <- aly_range_overlap %>%
#    filter(st_type == st) %>%
#    ggplot(aes(x = st_upper, y = tag_value, color = app)) +
#    geom_boxplot() +
#    facet_grid(tag_type ~ mtx_type)
#  if (grepl("tiff?$", fig_path, perl = TRUE, ignore.case = TRUE)) {
#    ggsave(fig_path, p, width = args$width, height = args$height, dpi = args$dpi, compress="lzw")
#  } else {
#    ggsave(fig_path, p, width = args$width, height = args$height, dpi = args$dpi)
#  }
#  msg <- paste0("the detailed result of stat:", st, "; is saved as ", fig_path)
#  print(msg)
#}
#write("", file = stdout())

print2("All Done!")
