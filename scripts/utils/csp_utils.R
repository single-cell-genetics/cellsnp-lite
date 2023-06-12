# csp_utils.R - functions for post-processing cellsnp output.


# library(Matrix)


#' Update three sparse matrices
#'
#' Update three sparse matrices `AD`, `DP`, and `OTH` given a VCF file storing
#' a list of SNPs that passed filtering.
#'
#' @param in_dir A string. The input cellSNP dir.
#' @param out_dir A string. The output dir, should not be the same with `in_dir`.
#' @param filtered_vcf_file A string. The VCF file storing SNPs that passed filtering.
#' @param is_vcf_gzipped A bool. Whether the `base.vcf` in `in_dir` is gzipped.
#' @param verbose A bool. Whether to output detailed log information.
#' @return Void.
#' @examples
#' update_cellsnp_matrices(
#'   in_dir = "./csp_mode1",
#'   out_dir = "./updated_matrix",
#'   filtered_vcf_file = "./csp_mode1.base.chr2.vcf.gz"
#' )
update_cellsnp_matrices <- function(
  in_dir,
  out_dir,
  filtered_vcf_file,
  is_vcf_gzipped = TRUE,
  verbose = TRUE)
{
  if (in_dir == out_dir)
    stop("Error: output dir should not be the input dir!")

  if (! dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

  assert_e(filtered_vcf_file)
  flt_snp <- read.delim(filtered_vcf_file, header = FALSE, comment.char = "#",
                        stringsAsFactors = FALSE)
  flt_snp_id <- paste0(flt_snp$V1, "_", flt_snp$V2)     # chrom_pos

  if (verbose)
    print(sprintf("The %d post-filtering SNPs loaded.", nrow(flt_snp)))

  suffix <- ifelse(is_vcf_gzipped, ".gz", "")
  raw_vcf_file <- paste0(in_dir, "/cellSNP.base.vcf", suffix)
  assert_e(raw_vcf_file)

  raw_snp <- read.delim(raw_vcf_file, header = FALSE, comment.char = "#",
                        stringsAsFactors = FALSE)
  raw_snp_id <- paste0(raw_snp$V1, "_", raw_snp$V2)     # chrom_pos

  if (verbose)
    print(sprintf("The %d raw SNPs loaded.", nrow(raw_snp)))

  flt_snp_index <- raw_snp_id %in% flt_snp_id
  
  for (mtx_name in c("AD", "DP", "OTH")) {
    mtx_file <- sprintf("%s/cellSNP.tag.%s.mtx", in_dir, mtx_name)
    assert_e(mtx_file)

    if (verbose)
      print(sprintf("Load %s matrix ...", mtx_name))

    mtx <- Matrix::readMM(file = mtx_file)

    if (verbose)
      print(sprintf("n_snp=%d; n_cell=%d.", nrow(mtx), ncol(mtx)))

    if (verbose)
      print(sprintf("Update %s matrix ...", mtx_name))

    flt_mtx <- mtx[flt_snp_index, ]

    if (verbose)
      print(sprintf("n_snp=%d; n_cell=%d.", nrow(flt_mtx), ncol(flt_mtx)))

    new_matrix_file <- sprintf("%s/cellSNP.tag.%s.mtx", out_dir, mtx_name)
    output_cellsnp_matrix(flt_mtx, new_matrix_file)
  }

  system(sprintf("cp %s %s/cellSNP.base.vcf%s", filtered_vcf_file, out_dir, suffix))
}


assert_e <- function(path) {
  if (is.null(path) || (! file.exists(path)))
    stop(sprintf("Error: '%s' does not exist!", path))
}


output_cellsnp_matrix <- function(mtx, out_file)
{
  # Do not use Matrix::writeMM, as it will omit the "1"s when
  # all the values of the third column are 1s.

  n_snp <- nrow(mtx)
  n_cell <- ncol(mtx)
  n_count <- length(mtx@i)

  dat <- data.frame(
    snp_idx = mtx@i + 1,
    cell_idx = mtx@j + 1,
    count = mtx@x
  )
  dat <- dat[order(dat$snp_idx, dat$cell_idx), ]

  write(c("%%MatrixMarket matrix coordinate integer general", 
          "%",
          sprintf("%d\t%d\t%d", n_snp, n_cell, n_count)), 
        file = out_file,
        sep = "\n")
  write.table(dat, out_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE, append = TRUE)
}


