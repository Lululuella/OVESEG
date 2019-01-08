#' p-value by weighted permutation scheme
#'
#' This function estimates p-values by aggregating weighted permutations.
#' @param tt a vector of test statistics.
#' @param ttperm a matrix of test statistics from permutaitons.
#'     Rows correspond to probes and columns to one permutation.
#' @param W a matrix containing weights for each spot in \code{ttperm}.
#'     Provided by \code{\link{postProbNull}}.
#' @details P-values are estimated by weightedly accumulating test statistics
#' from permutations that are larger than observations
#' @return  p-values.
#' @export
#' @examples
#' #generate some example data#'
#' t.obs <- rnorm(100)
#' t.perm <- matrix(rnorm(100*1000),nrow=100)
#' w <- matrix(runif(100*1000),nrow=100)
#'
#' pv <- pvalueWeightedEst(t.obs, t.perm, w)
pvalueWeightedEst <- function(tt, ttperm, w)
{
    return(pvalue_weighted_est(c(tt), c(ttperm), c(w)))
    # tt <- c(tt)
    # ttperm <- c(ttperm)
    # w0 <- rep(1, length(tt))
    #
    # r0 <- rank(-c(tt, ttperm), ties.method = 'last')
    # w1 <- c(w0, w)
    # w1[r0] <- w1
    # wrank <- cumsum(w1)
    #
    # r <- r0[seq_along(tt)]
    # r2 <- rank(-tt)
    # r3 <- wrank[r] - r2
    # pv <- r3/sum(w)
    # return(pv)
}
