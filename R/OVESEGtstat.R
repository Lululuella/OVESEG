#' OVESEG-test statistics
#'
#' This function computes OVESEG-test statistics.
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
#' @param order.return a logical indicating whether the order of groups
#'     will be returned. If FALSE, only the highest expressed group index is
#'     return for each probe.
#' @param lmfit.return a logical indicating whether a \code{MArrayLM} fitted
#'     model object produced by \code{\link[limma]{lmFit}}
#'     should be returned.
#' @details OVESEG-test statistics are designed to mathematically match the
#' definition of molecular markers:
#' \deqn{\max_{k=1,...,K}\left\{min_{l \neq k}\left\{ \frac{\mu_k(j)-\mu_l(j)}
#' {\sigma(j)\sqrt{\frac{1}{N_k}+\frac{1}{N_l}}} \right\}\right\}}
#' where \eqn{j} is probe index, \eqn{K} is the number of groups, and
#' \eqn{\mu_k} is the mean expression of group \eqn{k}.
#' @return a list containing the following components:
#' \item{tstat}{a vector of OVESEG-test statistics for probes.}
#' \item{label}{a vector of group labels.}
#' \item{groupOrder}{If \code{order.return} is TRUE, a matrix with each row
#' being group indexes ordered based on decreasing expression levels. If
#' \code{order.return} is FALSE, a vector with each element being a probe's
#' highest expressed group index. Group indexes are positions in \code{label}.}
#' \item{fit}{a \code{MArrayLM} fitted model object produced by \code{lmFit}.
#' Returned only when \code{lmfit.return} is TRUE.}
#' @export
#' @examples
#' data(RocheBT)
#' rtstat <- OVESEGtstat(RocheBT$y, RocheBT$group, alpha='moderated')
#' rtstat <- OVESEGtstat(RocheBT$y, RocheBT$group, alpha=0.1)
OVESEGtstat <- function(y, group, weights=NULL, alpha='moderated',
                        order.return=FALSE, lmfit.return=FALSE) {
    K <- length(unique(group))
    if (K < 2) {
        stop("At least two groups are needed.")
    }

    if (is(y, "data.frame")) {
        data <- as.matrix(y)
    } else if (is(y, "SummarizedExperiment")) {
        data <- SummarizedExperiment::assay(y)
    } else if (is(y, "matrix") == FALSE) {
        stop("Only matrix, data frame and SummarizedExperiment
            object are supported for expression data!")
    }

    if (length(group) != ncol(y)) {
        stop('length(group) should be the same as ncol(y)!')
    }

    if (sum(table(group) < 2) > 0) {
        warning('At least one group has less than 2 samples!')
    }

    group <- factor(as.character(group))
    label <- levels(group)

    fit <- lmFit(y, model.matrix(~0+group), weights=weights)
    coeff.stded <- pairwise_tstat_unscaled(fit$coefficients,fit$stdev.unscaled)
    tstat <- row_min(coeff.stded)

    if (!is.null(alpha)) {
        if (!is.numeric(alpha)) {
            fit <- eBayes(fit, trend = FALSE)
            tstat <- tstat / sqrt(fit$s2.post)
        } else {
            tstat <- tstat / sqrt(alpha + (fit$sigma)^2)
        }
    }

    if (order.return == FALSE) {
        groupOrder <- row_which_max(fit$coefficients)
    } else {
        group.o <- apply(coeff.stded, 1, order) + 1
        group.o <- rbind(1, group.o)
        group.o <- t(group.o)
        groupOrder <- t(apply(cbind(fit$coefficients, group.o), 1, function(x)
            order(x[seq_len(K)], decreasing=TRUE)[x[K + seq_len(K)]]))
    }

    if (lmfit.return) {
        return(list(tstat=tstat, label=label, groupOrder=groupOrder, fit=fit))
    } else {
        return(list(tstat=tstat, label=label, groupOrder=groupOrder))
    }
}
