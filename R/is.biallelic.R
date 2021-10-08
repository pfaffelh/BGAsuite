#' is.biallelic
#' 
#' A function to check if the data is biallelic
#'
#' @usage is.biallelic(sample)
#' @param sample An N x M matrix with entrie AA,...,TT.
#'

#' @return A binary variable.
#' @export
#'
#' @examples
#' sample = forenseq.snipper.reference$sample
#' is.biallelic(sample)

# sample must be a matrix
# output is TRUE or FALSE, and a list of FALSE/TRUE and non-biallelic SNPs
is.biallelic<-function(sample) {
  M = ncol(sample)
  res = TRUE
  nonbi.SNPs = NULL
  for(m in 1:M) {
    loc = as.vector(sample[,m])
    loc = loc[!is.na(loc)]
    if(length(loc)) {
      loc = unlist(strsplit(loc, ","))
      loc = unlist(strsplit(loc, ""))
    }
    if(length(unique(loc))>2) {
      res = FALSE
      nonbi.SNPs = c(nonbi.SNPs, m)
    }
  }
#  if(SNP.return) res = list(res = res, nonbi.SNPs = nonbi.SNPs)
  res
}

