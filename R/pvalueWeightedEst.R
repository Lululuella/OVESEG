#' pvalueWeightedEst
#'
#' This function reduces data dimension by loading matrix and then project
#' dimension-reduced data to the hyperplane orthogonal to c(1,0,...,0),
#' i.e., the first axis in the new coordinate system..
#' @param tt A data set that will be internally coerced into a matrix.
#'     Each row is a gene and each column is a sample. Missing values
#'     are not supported.
#'     All-zero rows will be removed internally.
#' @param ttstar The matrix whose rows are loading vectors
#' @param W The matrix whose rows are loading vectors
#' @details This function can project gene expression vectors to simplex
#' filtering. This function helps observe all genes in simplex plot.
#' @return  The data after perspective projection.
#' @export
#' @examples
#' #obtain data
#' data(RocheBT)
#' data <- RocheBT$y
#'
#' #preprocess data
pvalueWeightedEst <-function (tt, ttstar, w)
{
    tt <- c(tt)
    ttstar <- c(ttstar)
    w0 <- rep(1, length(tt))

    r0 <- rank(-c(tt, ttstar), ties.method = 'last')
    w1 <- c(w0, w)
    w1[r0] <- w1
    wrank <- cumsum(w1)

    r <- r0[seq_along(tt)]
    r2 <- rank(-tt)
    r3 <- wrank[r] - r2
    pv <- r3/sum(w)
    return(pv)
}
