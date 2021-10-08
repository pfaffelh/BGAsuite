#' F.ia.iterate
#' 
#' A function iterated in order to estimate individual admixture.
#' @param q A vector of length K giving the last estimate for individual admixture of the test data.
#' @param p An K x 4 x M array containing the allele frequencies of the training data.
#' @param loc.y A 4 x K matrix containing the allele counts of the test data.

#' @return A vector of length K, the new estimate for individual admixture of the test data.
#' @export
#'

F.ia.iterate<-function(q, p, loc.y) {
  K = length(q)
  M = length(p[1,1,])

  E = matrix(0, nrow=K, ncol = M)
  loc = NULL
  for(i in 1:4) loc = rbind(loc, q %*% p[,i,]) # loc is a 4xM matrix
  loc[loc==0] = 1e-16
  loc[loc==1] = 1-1e-16
  for(k in 1:K) {
    E[k,] = loc.y[1,] * p[k,1,] / loc[1,] + loc.y[2,] * p[k,2,] / loc[2,] + loc.y[3,] * p[k,3,] / loc[3,] + loc.y[4,] * p[k,4,] / loc[4,]
  }
  res = rowSums(E)/(2*M) * q
  res/(sum(res))
}

