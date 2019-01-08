#' Posterior probabilities of each component null hypothesis
#'
#' This function computes posterior probabilities of each component null
#' hypothesis given observed probes. Such probe-wise probabilities
#' will be used as weights for aggregating permutations.
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
#' @param detail.return a logical indicating whether more details (e.g. lfdr)
#' will be returned.
#' @details Posterior probabilities of each component null hypothesis given
#' observed probes are estimated from ANOVA test on certain groups and local
#' fdr. There are totally (\eqn{K-1}) null hypotheses, where \eqn{K} is the
#' number of groups.
#' @return a list containing the following components:
#' \item{W}{a matrix of posterior probabilities for each component null
#' hypothesis given an observed probe. Rows correspond to probes and columns
#' to one hypothesis.}
#' \item{groupOrder}{a matrix with each row being group indexes ordered based
#' on decreasing expression levels.}
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
#' ppnull <- postProbNull(RocheBT$y, RocheBT$group, alpha='moderated')
#' ppnull <- postProbNull(RocheBT$y, RocheBT$group, alpha=0.1)
postProbNull <- function(y, group, weights=NULL, alpha='moderated',
                        detail.return=TRUE) {
    K <- length(unique(group))
    group <- factor(as.character(group))

    ###OVESEG-test statistics
    ovesegt <- OVESEGtstat(y, group, weights = weights, alpha = alpha,
                order.return = TRUE, lmfit.return = TRUE)

    groupOrder <- ovesegt$groupOrder

    groupMean <- ovesegt$fit$coefficients
    if (!is.null(alpha)) {
        if (!is.numeric(alpha)) {
            varWithinGroup <- ovesegt$fit$s2.post
            df <- ovesegt$fit$df.total
        } else {
            varWithinGroup <- ovesegt$fit$sigma^2
            df <- ncol(y) - K
        }
    } else {
        varWithinGroup <- 1
        df <- Inf
    }

    if (is.null(weights)) {
        weights <- matrix(1, nrow(y), ncol(y))
    }

    ### compute between-group variance for top M groups
    varBtwnGroup <- matrix(0, nrow(y), K-1)
    for(i in seq_len(nrow(y))) {
        numSample <- unlist(lapply(levels(group),
                            function(x) sum(weights[i,group==x])))
        for(M in seq.int(2,K)) {
            gpidx <- groupOrder[i,seq_len(M)]
            gsampleidx <- group %in% levels(group)[gpidx]
            wmean <- c(y[i,gsampleidx] %*% weights[i,gsampleidx]) /
                sum(weights[i,gsampleidx])
            varBtwnGroup[i,M-1] <- 1 / (M-1) *
                sum((groupMean[i,gpidx] - wmean)^2 * numSample[gpidx])
        }
    }

    ### f-statistics for ANOVA
    fstat <- varBtwnGroup / varWithinGroup

    ### compute genewise posterior weights by local fdr
    F.p.value <- matrix(0, nrow(y), K-1)
    lfdr <- matrix(0, nrow(y), K-1)
    W <- matrix(0, nrow(y), K-1)
    Wsum <- matrix(0, nrow(y), K)
    for(M in seq.int(K,2)) {
        df1 <- M - 1
        df2 <- df
        if (df2[1] > 1e+06) {
            F.p.value[,M-1] <-
                pchisq(df1 * fstat[,M-1], df1, lower.tail = FALSE)
        }
        else {
            F.p.value[,M-1] <- pf(fstat[,M-1], df1, df2, lower.tail = FALSE)
        }
        Wpvalue <- rep(seq_len(nrow(F.p.value)), ceiling(1000*(1-Wsum[,M])))
        lfdrweighted <- fdrtool(F.p.value[Wpvalue,M-1],
                                statistic="pvalue",
                                plot=FALSE, verbose=FALSE)$lfdr
        lfdr[Wsum[,M]<1,M-1] <- lfdrweighted[!duplicated(Wpvalue)]

        W[,M-1] <- (1 - Wsum[,M]) * lfdr[,M-1]
        Wsum[,M-1] <- W[,M-1] + Wsum[,M]
    }

    if (detail.return) {
        return(list(W=W, groupOrder=groupOrder, F.p.value=F.p.value, lfdr=lfdr,
                    fit=ovesegt$fit))
    } else {
        return(list(W=W, groupOrder=groupOrder))
    }
}
