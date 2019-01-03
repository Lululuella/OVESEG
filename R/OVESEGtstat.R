#' OVESEG test statistics
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
OVESEGtstat <- function(y, group, weights = NULL, alpha = 'moderated',
                        order.return = 'top1', lmfit.return = FALSE, cplusplus = TRUE) {
    K <- length(unique(group))
    if (K < 2) {
        stop("At least two groups are needed.")
    }
    #if (table(group))


    group <- factor(group)
    fit <- lmFit(y, model.matrix(~0+group), weights=weights)

    if (cplusplus) {
        coeff.stded = pairwise_tstat_unscaled(fit$coefficients,fit$stdev.unscaled);
        tstat = row_min(coeff.stded)
    } else {
        coeff.stded <- t(apply(cbind(fit$coefficients,fit$stdev.unscaled), 1,
                               function(x) {o <- order(x[1:K], decreasing = T);
                               t <- rep(NA, K-1);
                               for(i in 2:K) {
                                   t[i-1] <- (x[o[1]]- x[o[i]])/sqrt(x[o[1]+K]^2 + x[o[i]+K]^2);
                               }
                               t}))

        tstat <- apply(coeff.stded, 1, min)
    }



    if (!is.null(alpha)) {
        if (!is.numeric(alpha)) {
            fit <- eBayes(fit, trend = FALSE)
            tstat <- tstat / sqrt(fit$s2.post)
        } else {
            tstat <- tstat / sqrt(alpha + (fit$sigma)^2)
        }
    }


    if (order.return == 'top1') {
        if (cplusplus) {
            groupOrder <- row_which_max(fit$coefficients)
        } else {
            groupOrder <- apply(fit$coefficients, 1, which.max)
        }
    } else {
        group.o <- t(apply(coeff.stded, 1, order)) + 1
        group.o <- cbind(1, group.o)
        groupOrder <- t(apply(cbind(fit$coefficients, group.o), 1,
            function(x) order(x[1:K], decreasing = T)[x[K + 1:K]]))
    }

    if (lmfit.return) {
        return(list(tstat=tstat, groupOrder=groupOrder, fit=fit))
    } else {
        return(list(tstat=tstat, groupOrder=groupOrder))
    }
}
