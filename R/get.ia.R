#' get.ia
#' 
#' Obtaining individual admixture from a sample, using a reference/training database
#' @usage get.ia(test.sample, training.sample, training.samLoc, tol=1e-6, verbose = TRUE)
#' @param test.sample An N x M matrix with entries AA,..., TT.
#' @param training.sample An N x M matrix with entries AA,..., TT.
#' If N = 1, test.sample can also be a vector.
#' From training.sample, allele frequencies are computed, and from test.data the individual admixture
#' is computed by iteration of F.ia.iterate until convergence. 
#' @param training.samLoc A vector of length N giving the locations of all N samples. unique(training.samLoc) gives all different geographical locations of the data. We set K = length(unique(training.samLoc)).
#' @param tol The iteration stops if two consecutive calls of F.ia.iterate lead to vectors with L1-distance below tol.
#' @param verbose If TRUE, the number of iterations is reported, FALSE is silent mode
#' 
#' @return An N x K matrix, where entry (i,k) is the estimated contribution of ancestry from population k for sample i.
#' @export
#'
#' @examples
#' training.sample = VISAGE.BASIC.TOOL.snipper.reference$sample
#' training.samLoc = VISAGE.BASIC.TOOL.snipper.reference$samLoc
#' test.sample = VISAGE.BASIC.TOOL.EGDP$sample[1:10,]
#' get.ia(test.sample, training.sample, training.samLoc, verbose = TRUE)

get.ia<-function(test.sample, training.sample, training.samLoc, tol=1e-6, verbose = TRUE) {

  test.y = get.y(test.sample)
  freqs = get.freqs(training.sample, training.samLoc)
  N = length(test.y[,1,1])
  M = length(test.y[1,1,])
  K = length(freqs[,1,1])
  
  res = matrix(0, nrow = N, ncol = K, dimnames = list(names(test.y[,1,1]), names(freqs[,1,1])))

  for(k in 1:N) {
    loc.y = test.y[k,,]
    loc.freqs = freqs[,,!is.na(loc.y[1,])]
    loc.y = loc.y[,!is.na(loc.y)[1,]]

    # initialize admixture proportions
    res[k,] = rep(1/K, K)

    err = 1
    j=1 
    while(err>tol) {
      j = j+1
      loc = F.ia.iterate(res[k,], loc.freqs, loc.y)
      err = sum(abs(res[k,] - loc))
      res[k,] = loc
    }
    if(verbose) cat("\t Iterations: ", j, "\n")
  }
  res    
}

