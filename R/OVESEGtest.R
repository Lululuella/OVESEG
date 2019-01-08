#' OVESEG-test
#'
#' This function performs OVESEG-test to assess significance of molecular
#' markers.
#' @param y a numeric matrix containing log-expression or logCPM
#'     (log2-counts per million) values.
#'     Data frame, or SummarizedExperiment object will be
#'     internally coerced into a matrix.
#'     Rows correspond to probes and columns to samples.
#'     Missing values are not permitted.
#' @param group categorical vector or factor giving group membership of
#'     columns of y. At least two categories need to be presented.
#' @param weights optional numeric matrix containing prior weights
#'     for each spot.
#' @param alpha parameter specifying within-group variance estimator to be used.
#'     'moderated': empirical Bayes moderated variance estimator as used in
#'     \code{\link[limma]{eBayes}}.
#'     numberic value: a constant value added to pooled variance estimator
#'     (\eqn{\alpha + \sigma}).
#'     NULL: no estimator; all variances are set to be 1.
#' @param NumPerm an integer specifying the number of permutation resamplings
#'     (default 999).
#' @param seed an integer seed for the random number generator.
#' @param detail.return a logical indicating whether more details about
#'     posterior probability estimation will be returned.
#' @param BPPARAM a BiocParallelParam object indicating whether parallelization
#'     should be used for permutation resamplings. The default is bpparam().
#' @details OVESEG-test is a statistically-principled method that can detect
#' tissue/cell-specific marker genes among many subtypes.
#' OVESEG-test statistics are designed to mathematically match the definition
#' of molecular markers, and a novel permutation scheme are employed to
#' estimate the corresponding distribution under null hypotheses where the
#' expression patterns of non-markers can be highly complex.
#' @return a list containing the following components:
#' \item{tstat}{a vector of OVESEG-test statistics for probes.}
#' \item{groupOrder}{If \code{order.return} is TRUE, a matrix with each row
#' being group indexes ordered based on decreasing expression levels. If
#' \code{order.return} is FALSE, a vector with each element being a probe's
#' highest expressed group index.}
#' \item{fit}{a \code{MArrayLM} fitted model object produced by \code{lmFit}.
#' This is returned only when \code{lmfit.return} is TRUE.}
#' @export
#' @examples
#' data(RocheBT)
#' rtest <- OVESEGtest(RocheBT$y, RocheBT$group, NumPerm=99,
#'                     BPPARAM=BiocParallel::SerialParam())
#' \dontrun{
#' #parallel computing
#' rtest <- OVESEGtest(RocheBT$y, RocheBT$group, NumPerm=99,
#'                     BPPARAM=BiocParallel::SnowParam())
#'}
OVESEGtest <- function(y, group, weights = NULL, alpha = 'moderated',
                       NumPerm = 999, seed = 111, detail.return=TRUE,
                       BPPARAM=bpparam()){
    K <- length(unique(group))
    group <- factor(as.character(group))

    message('Calculating posterior probabilities of null hypotheses')
    ppnull <- postProbNull(y, group, weights = weights, alpha = alpha,
                            detail.return = detail.return)

    ## permutation among top M groups
    tstat.perm <- c()
    topidx.perm <- c()
    for (M in seq.int(2,K)) {
        message('Permuting top ', M, ' groups')
        ovet.perm <- OVEtstatPermTopM(y, group, ppnull$groupOrder, M=M,
                            weights = weights, alpha = alpha,
                            NumPerm = NumPerm, seed = seed, BPPARAM = BPPARAM)
        tstat.perm <- cbind(tstat.perm, ovet.perm$tstat.perm)
        topidx.perm <- cbind(topidx.perm, ovet.perm$topidx.perm)
    }

    message('Calculating p-values')
    wPerm <- ppnull$W[,rep(1:(K-1),each=NumPerm+1)]
    testing <- tstat.perm[,1]
    pv.overall <- pvalueWeightedEst(testing, tstat.perm, wPerm)

    pv.oneside <- pv.overall
    pv.oneside.max <- rep(0, K)
    for (M in seq_len(K)) {
        idx <- which(ppnull$groupOrder[,1] == M)
        tstat.perm.oneside <- tstat.perm
        tstat.perm.oneside[topidx.perm != M] <- -1
        pvtmp <- pvalueWeightedEst(c(0,testing[idx]), tstat.perm.oneside, wPerm)
        pv.oneside[idx] <- pvtmp[-1]
        pv.oneside.max[M] <- pvtmp[1]
    }

    minpv <- min(rowSums(ppnull$W)) / (NumPerm + 1) / sum(ppnull$W)
    pv.multiside <- pv.oneside * K
    pv.multiside[pv.multiside < 0] <- 0
    pv.multiside[pv.multiside > 1] <- 1

    return(c(list(pv.oneside=pv.oneside, pv.oneside.max=pv.oneside.max,
                pv.multiside=pv.multiside,
                pv.overall=pv.overall), ppnull))
}
