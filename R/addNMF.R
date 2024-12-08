#' Non-negative Matrix Factorization
#'
#' These functions perform Non-negative Matrix Factorization on data stored in a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
#' 
#' @param x a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
#'   
#' @param k \code{numeric vector}. A number of latent vectors/topics. 
#' (Default: \code{2})
#' 
#' @param name \code{Character scalar}. The name to be used to store the result 
#' in the reducedDims of the output. (Default: \code{"NMF"})
#'  
#' @param assay.type \code{Character scalar}. Specifies which assay to use for 
#' NMF ordination. (Default: \code{"counts"})
#'  
#' @param eval.metric \code{Character scalar}. Specifies the evaluation metric
#' that will be used to select the model with the best fit. Must be one of the
#' following options: \code{"evar"} (explained variance; maximized),
#' \code{"sparseness.basis"} (degree of sparsity in the basis matrix;
#' maximized), \code{"sparseness.coef"} (degree of sparsity in the coefficient
#' matrix; maximized), \code{"rss"} (residual sum of squares; minimized),
#' \code{"silhouette.coef"} (quality of clustering based on the coefficient
#' matrix; maximized), \code{"silhouette.basis"} (quality of clustering based
#' on the basis matrix; maximized), \code{"cophenetic"} (correlation between
#' cophenetic distances and original distances; maximized), \code{"dispersion"}
#' (spread of data points within clusters; minimized). (Default: \code{"evar"})
#' 
#' @param ... optional arguments passed to \code{nmf::NMF}.
#' 
#' @return 
#' For \code{getNMF}, the ordination matrix with feature loadings matrix
#' as attribute \code{"loadings"}.
#'  
#' For \code{addNMF}, a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object is returned containing the ordination matrix in
#' \code{reducedDims(x, name)} with the following attributes:
#' \itemize{
#'   \item "loadings" which is a matrix containing the feature loadings
#'   \item "NMF_output" which is the output of function \code{nmf::NMF} 
#'   \item "best_fit" which is the result of the best fit if k is a vector of
#'   integers
#' }
#' 
#' @details 
#' The functions \code{getNMF} and \code{addNMF} internally use \code{nmf::NMF} 
#' compute the ordination matrix and 
#' feature loadings.
#'  
#' If k is a vector of integers, NMF output is calculated for all the rank 
#' values contained in k, and the best fit is selected based on 
#' \code{eval.metric} value.
#'  
#' @name addNMF
#' 
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' 
#' # Reduce the number of features
#' tse <- agglomerateByPrevalence(tse, rank = "Phylum")
#' 
#' # Run NMF and add the result to reducedDim(tse, "NMF").
#' tse <- addNMF(tse, k = 2, name = "NMF")
#' 
#' # Extract feature loadings
#' loadings_NMF <- attr(reducedDim(tse, "NMF"), "loadings")
#' head(loadings_NMF)
#' 
#' # Estimate models with number of topics from 2 to 4. Perform 2 runs.
#' tse <- addNMF(tse, k = c(2, 3, 4), name = "NMF_4", nrun = 2)
#' 
#' # Extract feature loadings
#' loadings_NMF_4 <- attr(reducedDim(tse, "NMF_4"), "loadings")
#' head(loadings_NMF_4)
#' 
NULL

#' @export
#' @rdname addNMF
setMethod("getNMF", "SummarizedExperiment",
    function(
        x, k = 2, assay.type = "counts", eval.metric = "evar", ...){
    .require_package("NMF")
    # Both NMF and DelayedArray have method seed(). When running
    # NMF::nmf() an error occurs due to wrong method. That is why NMF
    # is first loaded into the session.
    if( "NMF" %in% (.packages()) ){
        detach("package:NMF", unload = TRUE)
    }
    library("NMF")
    .check_assay_present(assay.type, x)
    # Calculate NMF ordination
    mat <- t(assay(x, assay.type))
    res <- NMF::nmf(mat, k, ...)
    # Check oif the output includes multiple ordination with different k values.
    # If it includes multiple, get the best fit based on certain evaluation
    # metric.
    if( is(res, "NMF.rank") ){
        best_fit <- .get_best_nmf_fit(res, eval.metric)
    } else{
        best_fit <- res
    }
    # Get scores and loadings, add loadings and NMF output to attributes of
    # scores
    scores <- best_fit@fit@W
    loadings <- best_fit@fit@H
    attr(scores, "loadings") <- t(loadings)
    attr(scores, "NMF_output") <- res
    # Add best fit if multiple k values
    if( is(res, "NMF.rank") ){
        attr(scores, "best_fit") <- best_fit
    }
    # The NMF package is unloaded
    detach("package:NMF", unload = TRUE)
    # Return scores with loadings, metrics and model as attribute
    return(scores)
}
)

#' @export
#' @rdname addNMF
setMethod("addNMF", "SummarizedExperiment",
    function(
        x, k = 2, assay.type = "counts", eval.metric = "evar", name = "NMF",
        ...){
    # Input checks
    if( !.is_a_string(name) ){
        stop("'name' must be a non-empty single character value.",
            call. = FALSE)
    }
    # Fit the model
    nmf <- getNMF(x, k = k, assay.type = assay.type, ...)
    # Add scores matrix with loadings as attribute to reducedDims
    x <- .add_values_to_reducedDims(x, name = name, values = nmf)
    return(x)
    }
)

################################ HELP FUNCTIONS ################################
# This function is for evaluating a fit of NMF models
.get_best_nmf_fit <- function(res, eval.metric){
    # Get whether the metric is maximized or minimized
    maximize <- c(
        "sparseness.basis" = TRUE,
        "sparseness.coef" = TRUE,
        "rss" = FALSE,
        "evar" = TRUE,
        "silhouette.coef" = TRUE,
        "silhouette.basis" = TRUE,
        "cophenetic" = TRUE,
        "dispersion" = FALSE
        )
    maximize <- maximize[ eval.metric ]
    if( maximize ){
        FUN <- which.max
    } else{
        FUN <- which.min
    }
    # Get the index of best fit
    measures <- res[["measures"]]
    values <- measures[[eval.metric]]
    ind <- FUN(values)
    # Get the model of best fit
    model <- res[["fit"]][[ind]]
    return(model)
}
