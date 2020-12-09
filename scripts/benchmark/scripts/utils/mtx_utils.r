#this script contains several mtx utils.
#hxj5<hxj5@hku.hk>

library(tibble)

#@abstract  Parse mtx file to get mtx info.
#@param fn  Name of mtx file [STR]
#@return    A list with four elements if success, NULL otherwise [list]
#           The four elements are:
#             $nrow     Number of rows [INT]
#             $ncol     Number of columns [INT]
#             $nval     Number of values [INT]
#             $data     The matrix data without header [dataframe]
parse_mtx_file <- function(fn) {
  if (! file.exists(fn)) { return(NULL) }
  df <- read.table(fn, header = F, comment.char = "%")
  if (nrow(df) < 2) { return(NULL) }
  colnames(df)[1:3] <- c("row", "col", "value")
  return(list(
    nrow = df[1, 1], 
    ncol = df[1, 2],
    nval = df[1, 3],
    data = df[-1, ] %>% as_tibble()
  ))
}

#@abstract   Write the mtx data to the mtx file.
#@param mtx  The mtx data [list]
#@param fn   Name of mtx file [STR]
#@return     Void.
write_mtx_file <- function(mtx, fn) {
  write("%%MatrixMarket matrix coordinate integer general", file = fn)
  write("%", file = fn, append = T)
  write(paste(c(mtx$nrow, mtx$ncol, mtx$nval), collapse = "\t"), 
        file = fn, append = T)
  write.table(mtx$data, fn, append = T, sep = "\t", row.names = F, 
              col.names = F)
}

