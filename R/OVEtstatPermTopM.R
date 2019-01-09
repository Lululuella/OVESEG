#' OVESEG-test statistics after permuting top M groups
#'
#' This function permutes group labels among highest expressed M groups and
#' then computes new OVESEG-test statistics.
#' @param y a numeric matrix containing log-expression or logCPM
#'     (log2-counts per million) values.
#'     Data frame, or SummarizedExperiment object will be
#'     internally coerced into a matrix.
#'     Rows correspond to probes and columns to samples.
#'     Missing values are not permitted.
#' @param group categorical vector or factor giving group membership of
#'     columns of y. At least two categories need to be presented.
#' @param groupOrder an integer matrix with each row giving group indexes
#'     ordered based on decreasing expression levels .
#' @param M an integer indicating the number of groups being permutated.
#'     The range is \eqn{[2,K]}, where K is the total number of groups.
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
#' @param BPPARAM a BiocParallelParam object indicating whether parallelization
#'     should be used for permutation resamplings. The default is bpparam().
#' @details Top M expressed groups will be involved in permutation. There are
#' /eqn{C_K^M} probe patterns which are highly expressed in certain M
#' groups among the total K groups. Probes of the same pattern share the same
#' shuffled labels.
#'
#' To improve the time efficiency, some functions within permutation loop are
#' implemented using c++.
#' @return a list containing the following components:
#' \item{tstat.perm}{a numeric matrix with each column giving OVESEG-test
#' statistics over the expressions after one permutation.}
#' \item{topidx.perm}{a integer matrix with each colulmn giving the highest
#' expressed group index over the expressions after one permutation.}
#' @export
#' @examples
#' data(RocheBT)
#' ppnull <- postProbNull(RocheBT$y, RocheBT$group, detail.return = TRUE)
#' rperm <- OVEtstatPermTopM(RocheBT$y, RocheBT$group, ppnull$groupOrder, M=2,
#'                          NumPerm=99, BPPARAM=BiocParallel::SerialParam())
#' \dontrun{
#' #parallel computing
#' rperm <- OVEtstatPermTopM(RocheBT$y, RocheBT$group, ppnull$groupOrder, M=2,
#'                          NumPerm=99, BPPARAM=BiocParallel::SnowParam())
#'}
OVEtstatPermTopM <- function(y, group, groupOrder, M, weights=NULL,
                        alpha='moderated', NumPerm=999, seed=111,
                        BPPARAM=bpparam()) {

    K <- length(unique(group))
    if (K < M) {
        stop("M should be larger than then number of gorups!")
    }
    group <- factor(as.character(group))

    seeds <- unlist(lapply(seq_len(NumPerm), function(x) seed + x - 1))

    combM <- t(combn(K, M))
    ncombM <- nrow(combM)
    geneSubset <- integer(nrow(y))
    for(j in 1:ncombM) {
        geneidx <- apply(groupOrder[,1:M], 1,
                        function(x) setequal(x, combM[j,]))
        geneSubset[geneidx] <- j
    }

    NPerm <- NumPerm + 1

    if (is.null(BPPARAM)) {
        oveseg.perm <- lapply(seq_len(NPerm), permfunc,
                        y, group, weights, alpha, combM, geneSubset, seeds)
    } else {
        oveseg.perm <- bplapply(seq_len(NPerm), permfunc,
                        y, group, weights, alpha, combM, geneSubset, seeds,
                        BPPARAM = BPPARAM)
    }

    tstat.perm <- vapply(oveseg.perm, function(x) x$tstat,
                                FUN.VALUE= numeric(nrow(y)))

    topidx.perm <- vapply(oveseg.perm, function(x) x$groupOrder,
                                FUN.VALUE= integer(nrow(y)))


    return(list(tstat.perm=tstat.perm, topidx.perm=topidx.perm))
}




#' Internal function for one permutation task
#'
#' @param p an integer indicating permutation index
#' @param y an expressions matrix
#' @param group a integer vector indicating group labels
#' @param weights optional numeric matrix containing prior weights
#' @param alpha parameter specifying within-group variance estimator to be used
#' @param combM a integer matrix with each row giving one choice of M groups
#' @param geneSubset a integer vector indicating the probe pattern of combM
#' @param seed an integer seed for the random number generator
#' @return test statistics and upregulated group indexes after one permutaion
permfunc <- function(p, y, group, weights, alpha, combM, geneSubset, seeds) {
    if(p > 1){
        groupid <- unlist(lapply(group,
                                 function(x) which(levels(group)==x)))
        data.shuffled <- shuffle_topm(y, groupid, weights,
                                      combM, geneSubset, seeds[p-1])
        y.p <- data.shuffled$y
        weights.p <- data.shuffled$weights
    } else {
        y.p <- y
        weights.p <- weights
    }

    ovesegt <- OVESEGtstat(y.p, group, weights = weights.p, alpha = alpha,
                           order.return = FALSE, lmfit.return = FALSE)

    return(ovesegt)
}

