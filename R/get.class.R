#' get.class
#' 
#' A function to obtain probabilities that a test sample belongs to a population/class. At the momenta, only a naive Bayesian classifier is implemented.
#' @usage get.class(test.sample, training.sample, training.samLoc, method = "naiveBayes")
#' @param test.sample An N_te x M matrix with entries AA,..., TT.
#' @param training.sample An N_Tr x M matrix with entries AA,..., TT.
#' From training.sample, allele frequencies are computed, and from test.data the individual admixture
#' is computed by iteration of F.ia.iterate until convergence. 
#' @param training.samLoc A vector of length N giving the locations of all N samples. unique(training.samLoc) give
#' @param method The method used for classification. Until now, only a naive Bayes approach is implemented.
#' 
#' @return An N x K matrix, where entry (i,k) is the probability that test-sample belongs to population k.
#' @export
#'
#' @examples
#' training.sample = forenseq.snipper.reference$sample
#' training.samLoc = forenseq.snipper.reference$samLoc
#' test.sample = forenseq.EGDP$sample
#' get.class(test.sample, training.sample, training.samLoc)

get.class<-function(test.sample, training.sample, training.samLoc, method = "naiveBayes") {
  pops = unique(training.samLoc)
  K = length(pops)
  
  test.y = get.y(test.sample)
  training.y = get.y(training.sample)
  N = length(test.y[,1,1])
  M = length(test.y[1,1,])
  if(method == "naiveBayes") {
    # determine number of occurrences 
    training.allelecounts = array(1, dim = c(length(pops), 4, length(training.y[1,1,])), dimnames=list(pops, c("A", "C", "G", "T"), names(training.y[1,1,])))
    training.allelecounts[,"A",] = t(simplify2array(by(training.y[,"A",],training.samLoc,colSums,na.rm=TRUE)))
    training.allelecounts[,"C",] = t(simplify2array(by(training.y[,"C",],training.samLoc,colSums,na.rm=TRUE)))
    training.allelecounts[,"G",] = t(simplify2array(by(training.y[,"G",],training.samLoc,colSums,na.rm=TRUE)))
    training.allelecounts[,"T",] = t(simplify2array(by(training.y[,"T",],training.samLoc,colSums,na.rm=TRUE)))

    res = matrix(0, nrow = length(test.y[,1,1]), ncol = K, dimnames = list(rownames(test.y),pops))

    for(i in 1:N) {
      loc.y = test.y[i,,]
      loc.allelecounts = training.allelecounts[,,!is.na(loc.y[1,])]
      loc.y = loc.y[,!is.na(loc.y)[1,]]

      for(k in 1:K) {
        a = sum(log((1 + loc.allelecounts[k,"A",])) * loc.y["A",])
        a = a + sum(log((1 + loc.allelecounts[k,"C",])) * loc.y["C",])
        a = a + sum(log((1 + loc.allelecounts[k,"G",])) * loc.y["G",])
        a = a + sum(log((1 + loc.allelecounts[k,"T",])) * loc.y["T",])
        res[i,k] = a
      }
    } 

    for(i in 1:nrow(res)) {
      res[i,] = exp(res[i,] - max(res[i,]))
    }
    res / rowSums(res)
  } else {
    cat("Sorry, method", method, "is not implemented.\n")
  }
}

