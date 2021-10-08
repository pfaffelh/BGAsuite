#' loglik.q
#' 
#' A function to obtain the log-likelihood of the individual admixture of a single sample
#'
#' @usage loglik.q(sample, p, q, sum=TRUE)
#' @param sample A vector of length M giving the sample.
#' @param p A K x 4 x M array of allele frequencies.
#' @param q A vector of length K of the individual admixture estimates for the sample.
#' @param sum A binary variable indicating if the sum over all markers should be reported or all contributions.

#' @return Either a single value (if sum=FALSE) or a vector of length M (if sum=TRUE), giving the contributions to the log-likelihood.
#' @export
#'
#' @examples
#' training.sample = forenseq.snipper.reference$sample
#' training.samLoc = forenseq.snipper.reference$samLoc
#' test.sample = forenseq.EGDP$sample[1,]
#' p = get.freqs(training.sample, training.samLoc)
#' q = get.ia(test.sample, training.sample, training.samLoc)
#' loglik.q(test.sample, p, q)

loglik.q<-function(sample, p, q, sum=TRUE) {
  y = get.y(sample)[1,,]
  loc.p = p[,,!is.na(y[1,])]
  loc.y = y[,!is.na(y)[1,]]

  res = 0*loc.y[1,]
  
  loc = NULL # loc[i,m] will give the probability to obtain allele i at marker m
  for(i in 1:4) loc = rbind(loc, q %*% loc.p[,i,])
  for(i in 1:4) res = res + log(choose(2, loc.y[i,])) + log(loc[i,]^loc.y[i,])
  
  if(sum) res = sum(res)
  res
}

