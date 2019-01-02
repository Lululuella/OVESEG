#' posteriorProbNull
#'
#' This function reduces data dimension by loading matrix and then project
#' dimension-reduced data to the hyperplane orthogonal to c(1,0,...,0),
#' i.e., the first axis in the new coordinate system..
#' @param y A data set that will be internally coerced into a matrix.
#'     Each row is a gene and each column is a sample. Missing values
#'     are not supported.
#'     All-zero rows will be removed internally.
#' @param group The matrix whose rows are loading vectors;
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
posteriorProbNull <- function(y, group, weights = NULL, alpha = 'moderated',
                              detail.return = TRUE) {
    K <- length(unique(group))


    ovesegt <- OVESEGtstat(y, group, weights = weights, alpha = alpha,
                order.return = 'all', lmfit.return = TRUE)

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

    if(is.null(weights)) {
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
                                plot = FALSE)$lfdr
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
