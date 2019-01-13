#' OVESEG-test
#'
#' This function performs OVESEG-test to assess significance of molecular
#' markers.
#' @param y a numeric matrix containing log-expression or logCPM
#'     (log2-counts per million) values.
#'     Data frame or SummarizedExperiment object will be
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
#'     Numeric value: a constant value added to pooled variance estimator
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
#' \item{pv.overall}{a vector of p-values calculated by all permutations
#' regardless of upregulated subtypes.}
#' \item{pv.oneside}{a vector of subtype-specific p-values calculated
#' specifically for the upregulated subtype of each probe.}
#' \item{pv.oneside.max}{subtype-specific p-values when observed test
#' statistic equal to zero.}
#' \item{pv.multiside}{pv.oneside*K (K-time comparison correction) and
#' truncated at 1.}
#' \item{W}{a matrix of posterior probabilities for each component null
#' hypothesis given an observed probe. Rows correspond to probes and columns
#' to one hypothesis.}
#' \item{label}{a vector of group labels.}
#' \item{groupOrder}{a matrix with each row being group indexes ordered based
#' on decreasing expression levels.
#' Group indexes are positions in \code{label}.}
#' \item{F.p.value}{a matrix with each column giving p-values corresponding
#' to F-statistics on certain groups.}
#' \item{lfdr}{a matrix with each column being local false discovery rates
#' estimated based on one column of weighted F.p.value matrix.}
#' \item{fit}{a \code{MArrayLM} fitted model object produced by \code{lmFit}.}
#' \code{F.p.value}, \code{lfdr} and \code{fit} are returned only when
#' \code{detail.return} is TRUE.
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

    message('Calculating p-values for each group specifically')
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

    return(c(list(pv.overall=pv.overall,
                pv.oneside=pv.oneside,
                pv.oneside.max=pv.oneside.max,
                pv.multiside=pv.multiside), ppnull))
}

