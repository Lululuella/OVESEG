#' Probability of one group being upregulated under null
#'
#' This function estimates probabilities of any one group being upregulated
#' than other groups under null hypotheses.
#' @param ppnull a list returned by \code{\link{postProbNull}} or
#'     \code{\link{OVESEGtest}}.
#' @details The probability of one group being upregulated under null
#' hypotheses is calculated by accumulating and normalizing genewise posterior
#' probability of null hypotheses. The group with higher probability tends to
#' get more False Positive MGs.
#' @return a numeric vector indicating probabilities of each group being
#' upregulated than others under null hypotheses.
#' @export
#' @examples
#' data(RocheBT)
#' ppnull <- postProbNull(RocheBT$y, RocheBT$group, alpha='moderated')
#' pk <- nullDistri(ppnull)
nullDistri <- function(ppnull)
{
    W <- ppnull$W
    groupOrder <- ppnull$groupOrder
    K <- ncol(groupOrder)
    ngene <- nrow(groupOrder)

    prob <- matrix(0, ngene, K)
    for(M in seq.int(2, K)) {
        for(i in seq_len(ngene)) {
            prob[i,groupOrder[i,seq_len(M)]] <-
                prob[i,groupOrder[i,seq_len(M)]] + W[i,M-1] / M
        }
    }
    propK <- colSums(prob)/sum(prob)
    names(propK) <- ppnull$label
    return(propK)
}
