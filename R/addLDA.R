#' Latent Dirichlet Allocation
#'
#' These functions perform Latent Dirichlet Allocation on data stored in a 
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
#' 
#' @param x a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
#'   
#' @param k \code{Integer vector}. A number of latent vectors/topics. 
#'  (Default: \code{2})
#' 
#' @param name \code{Character scalar}. The name to be used to store the result 
#'  in the reducedDims of the output. (Default: \code{"LDA"})
#'  
#' @param assay.type \code{Character scalar}. Specifies which assay to use for 
#'  LDA ordination. (Default: \code{"counts"})
#' 
#' @param eval.metric \code{Character scalar}. Specifies evaluation metric that
#' will be used to select the model with the best fit. Must be either
#' \code{"perplexity"} (\code{topicmodels::perplexity}) or \code{"coherence"}
#' (\code{topicdoc::topic_coherence}, the best model is selected based on mean
#' coherence). (Default: \code{"perplexity"})
#' 
#' @param ... optional arguments passed to \code{\link[topicmodels:LDA]{LDA}}
#' 
#' @return 
#' For \code{getLDA}, the ordination matrix with feature loadings matrix
#' as attribute \code{"loadings"}.
#'  
#' For \code{addLDA}, a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object is returned containing the ordination matrix in
#' \code{reducedDim(..., name)} with feature loadings matrix as attribute
#' \code{"loadings"}.
#'  
#' @details 
#' The functions \code{getLDA} and \code{addLDA} internally use 
#' \code{\link[topicmodels:LDA]{LDA}} to compute the ordination matrix and 
#' feature loadings.
#'  
#' @name addLDA
#' 
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' 
#' # Reduce the number of features 
#' tse <- agglomerateByPrevalence(tse, rank="Phylum")
#' 
#' # Run LDA and add the result to reducedDim(tse, "LDA")
#' tse <- addLDA(tse)
#' 
#' # Extract feature loadings
#' loadings <- attr(reducedDim(tse, "LDA"), "loadings")
#' head(loadings)
#' 
#' # Estimate models with number of topics from 2 to 10
#' tse <- addLDA(tse, k = c(2, 3, 4, 5, 6, 7, 8, 9, 10), name = "LDA_10")
#' # Get the evaluation metrics
#' tab <- attr(reducedDim(tse, "LDA_10"),"eval_metrics")
#' # Plot
#' plot(tab[["k"]], tab[["perplexity"]], xlab = "k", ylab = "perplexity")
NULL

#' @export
#' @rdname addLDA
setMethod("getLDA", "SummarizedExperiment",
    function(x, k = 2, assay.type = "counts", eval.metric = "perplexity", ...){
        .require_package("topicmodels")
        # Input checks
        if( !.is_integer(k) ){
            stop("'k' must be an integer.", call. = FALSE)
        }
        .check_assay_present(assay.type, x)
        if( !(.is_a_string(eval.metric) &&
                eval.metric %in% c("perplexity", "coherence")) ){
            stop(
                "'eval.metric' must be 'perplexity' or 'coherence'.",
                call. = FALSE)
        }
        # If eval.metric is coherence, tpicdoc must be installed
        if( eval.metric == "coherence" ){
            .require_package("topicdoc")
        }
        #
        df <- as.data.frame(t(assay(x, assay.type)))
        # Fit LDA models with different k values
        models <- lapply(k, function(i){
            topicmodels::LDA(df, k = i, ...)
        })
        names(models) <- k
        # Evaluate the goodness of fit of models
        metrics <- .calculate_lda_metrics(models, df)
        # Get the model with the best fit
        which_FUN <- switch (eval.metric,
            "perplexity" = which.min,
            "coherence" = which.max
        )
        model <- models[[ which_FUN(metrics[[eval.metric]]) ]]
        # Calculate scores and loadings
        posteriors <- topicmodels::posterior(model, df)
        scores <- t(as.data.frame(posteriors$topics))
        loadings <- t(as.data.frame(posteriors$terms)) 
        # Add loadings as attribute of the scores matrix
        attr(scores, "loadings") <- loadings
        # Add LDA model as attribute of the scores matrix
        attr(scores, "model") <- model
        # Add evaluation metrics table
        attr(scores, "eval_metrics") <- metrics
        # Return scores with loadings, metrics and model as attribute
        return(scores)
    }
)

#' @export
#' @rdname addLDA
setMethod("addLDA", "SummarizedExperiment",
    function(x, k = 2, assay.type = "counts", name = "LDA", ...){
        # Input checks
        if( !.is_a_string(name) ){
            stop("'name' must be a non-empty single character value.",
                call. = FALSE)
        }
        # Fit the model
        scores <- t(getLDA(x, k, assay.type, ...))
        # Add scores matrix with loadings as attribute to reducedDims
        x <- .add_values_to_reducedDims(x, name, values = scores)
        return(x)
    }
)

################################ HELP FUNCTIONS ################################
# This function is for evaluating a fit of LDA models
#' @importFrom dplyr bind_rows
.calculate_lda_metrics <- function(models, df){
    # Loop over each model
    metrics <- lapply(models, function(model){
        # Calculate perplexity
        res <- topicmodels::perplexity(model)
        names(res) <- "perplexity"
        # If topicdoc package is available, calculate also coherence. Each topic
        # has own coherence, so calculate mean.
        if( require("topicdoc", quietly = TRUE) ){
            coherence <- topicdoc::topic_coherence(model, df)
            names(coherence) <- paste0("coherence_", seq_len(length(coherence)))
            coherence <- c(coherence = mean(coherence), coherence)
            res <- c(res, coherence)
        }
        return(res)
    })
    # Create a data.frame from metrics
    metrics <- as.data.frame( bind_rows(metrics) )
    metrics[["k"]] <- as.numeric( names(models) )
    return(metrics)
}
