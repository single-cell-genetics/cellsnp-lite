#!/usr/bin/env Rscript
#this script is aimed to visualize the performance of benchmarking for cellsnp-lite.
#modified from https://github.com/shenwei356/seqkit/blob/master/benchmark/plot.R
#hxj5<hxj5@hku.hk>

# tools
#"#b79d1a", "#ff7600", "#d20015", "#00824c", "#00518b"
#gray yellow, orange red, red, green, blue
all_colors <- c("#00518b", "#ff7600", "#00824c", "#d20015", "#b79d1a")

# default settings
def_width <- 8
def_height <- 6
def_dpi <- 300

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
parser$add_argument("-i", "--infile", type = "character",
                    help = "Performance log file.")
parser$add_argument("-o", "--outfile", type = "character",
                    default = "", help = "[Optional] Result summary file.")
parser$add_argument("-f", "--outfig", type = "character",
                    default = "", help = "Result figure file.")
parser$add_argument("--title", type = "character",
                    help = "Title of output figure.")
parser$add_argument("--width", type = "double",
                    default = def_width, help = paste0("Result file width [", def_width, "]"))
parser$add_argument("--height", type = "double",
                    default = def_height, help = paste0("Result file height [", def_height, "]"))
parser$add_argument("--dpi", type = "integer",
                    default = def_dpi, help = paste0("DPI [", def_dpi, "]"))
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

check_path_exists(args$infile, "input file")
#check_arg_null(args$outfile, "output file")
check_arg_null(args$outfig, "output figure")
#check_arg_null(args$title, "title of output figure")
check_arg_null(args$width, "figure width")
check_arg_null(args$height, "figure height")
check_arg_null(args$dpi, "figure dpi")

outfig <- args$outfig
title <- args$title

# load data
options(warn = -1)
library(dplyr, warn.conflicts = F)
library(ggplot2)
library(ggrepel)
options(warn = 0)

raw_dat <- read.table(args$infile, header = T, sep = "\t")

# grouping and calc mean and sd values of each group
dat <- raw_dat %>% 
  as_tibble() %>%
  mutate(ncore_f = factor(ncore, levels = sort(unique(ncore)))) %>%
  group_by(app, ncore_f) %>%
  summarise(time_mean = mean(time), time_sd = sd(time),
            mem_mean = mean(mem), mem_sd = sd(mem)) %>%
  ungroup()

if (length(args$outfile) > 0 && args$outfile != "") {
  dat_out <- dat %>%
    mutate(time_mean = round(time_mean, digits = 0)) %>%
    mutate(time_sd = round(time_sd, digits = 0)) %>%
    mutate(mem_mean = round(mem_mean, digits = 0)) %>%
    mutate(mem_sd = round(mem_sd, digits = 0)) %>%
    arrange(ncore_f, app)
  write.table(dat_out, args$outfile, quote = F, sep = "\t", row.names = F)
  print(paste0("The output summary file is saved as '", args$outfile, "'"))
}

# humanize time unit
max_time <- max(dat$time_mean)
time_unit <- "s"
if (max_time > 3600) {
  dat <- dat %>% mutate(time_mean2 = time_mean / 3600)
  time_unit <- "h"
} else if (max_time > 60) {
  dat <- dat %>% mutate(time_mean2 = time_mean / 60)
  time_unit <- "m"
} else {
  dat <- dat %>% mutate(time_mean2 = time_mean / 1)
  time_unit <- "s"
}

# humanize mem unit
max_mem <- max(dat$mem_mean)
mem_unit <- "KB"
if (max_mem > 1024 * 1024) {
  dat <- dat %>% mutate(mem_mean2 = mem_mean / 1024 / 1024)
  mem_unit <- "GB"
} else if (max_mem > 1024) {
  dat <- dat %>% mutate(mem_mean2 = mem_mean / 1024)
  mem_unit <- "MB"
} else {
  dat <- dat %>% mutate(mem_mean2 = mem_mean / 1)
  mem_unit <- "KB"
}

# visualize the mean values
all_tools <- sort(unique(dat$app))
tool_colors <- all_colors[1:length(all_tools)]
p <- dat %>%
  ggplot(aes(x = mem_mean2, y = time_mean2, color = app, 
             shape = ncore_f, label = app)) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = time_mean2, color = app), size = 0.1, alpha = 0.4) +
  geom_vline(aes(xintercept = mem_mean2, color = app), size = 0.1, alpha = 0.4) +
  geom_text_repel(size = 4, max.iter = 400000) +
  scale_color_manual(values = tool_colors) +
  xlim(0, max(dat$mem_mean2)) +
  ylim(0, max(dat$time_mean2)) +
  xlab(paste0("Peak Memory (", mem_unit, ")")) +
  ylab(paste0("Time (", time_unit, ")")) +
  labs(color = "Tools", shape = "Ncores")

p <- p +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", size = 1.2),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_line(size = 0.8),
    axis.ticks.x = element_line(size = 0.8),
    
    strip.background = element_rect(
      colour = "white", fill = "white",
      size = 0.2
    ),
    
    legend.text = element_text(size = 12),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.6, "cm"),
    legend.key = element_blank(),
    legend.text.align = 0,
    legend.box.just = "left",
    strip.text.x = element_text(angle = 0, hjust = 0),
    
    text = element_text(
      #size = 14, family = "arial", face = "bold"
      size = 12, face = "bold"
    )
  )

if (! is.null(title) && title != "") {
  p <- p + labs(title = title) +
    theme(plot.title = element_text(size = 15))
}

if (grepl("tiff?$", outfig, perl = TRUE, ignore.case = TRUE)) {
  ggsave(outfig, p, width = args$width, height = args$height, dpi = args$dpi, compress="lzw")
} else {
  ggsave(outfig, p, width = args$width, height = args$height, dpi = args$dpi)
}
print(paste0("The output figure file is saved as '", outfig, "'"))

print("All Done!")
