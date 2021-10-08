#' get.y
#' 
#' Obtain allele counts for a sample
#'
#' get.y(sample)
#' @param sample A N x M matrix with entries AA,..., TT. If N = 1, sample can also be a vector.
#' @return A N x 4 x M array, where N is the sample size, and M is the number of markers.
#'
#' @examples
#' get.y(forenseq.snipper.reference$sample)
#' @export

get.y<-function(sample) {
  if(is.vector(sample)) sample = matrix(sample, nrow = 1, dimnames = list("1", names(sample)))
  M = ncol(sample)
  N = nrow(sample)
  
  res = array(0, dim = c(N, 4, M), dimnames=list(rownames(sample), c("A", "C", "G", "T"), colnames(sample)))
  for(i in 1:N) {
    for(m in 1:M) {
      loc = as.character(sample[i,m])
      if(!is.na(loc)) {
        loc = unlist(strsplit(loc, ","))
        loc = unlist(strsplit(loc, ""))
        res[i,loc[1],m] = res[i,loc[1],m] +1
        res[i,loc[2],m] = res[i,loc[2],m] +1
      } else {
        res[i,,m]=NA
      }
    }
  }
  res
}

