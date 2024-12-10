#' @name
#' getAbundant
#' 
#' @title
#' Determine abundant and rare taxa. Rare taxa can be further classified
#' to conditionally rare and permanently rare.
#' 
#' @description
#' These functions determine abundant and rare taxa based on the abundances
#' of taxa. Compared to \code{\link[=getPrevalence]{getPrevalent}} and
#' \code{\link[=getPrevalence]{getRare}}, these functions determine abundant
#' and rare taxa based on abundance while the first mentioned are based on
#' prevalence.
#' 
#' @details
#' Abundant and rare taxa 
#' 
#' @return
#' Vector of taxa.
#' 
#' @inheritParams runCCA
#' 
#' @param threshold
#'
#' @param ... additional arguments.
#' \itemize{
#'   \item \code{by}: \code{Character scalar}. Specifies how significance is
#'   calculated. (Default: \code{"margin"})
#' }
#'
#' @examples
#' 
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' 
#' # Agglomerate to family level
#' tse <- agglomerateByRank(tse, rank = "Family")
#' # Tranform to relative abundances
#' tse <- transformAssay(tse, method = "relabundance", pseudocount = TRUE)
#' 
#' # Get abundant taxa
#' abundant <- getAbundant(tse, assay.type = "relabundance")
#' abundant |> head()
#' 
#' # Get all rare taxa that have average relartive abundance below 10%
#' rare <- getLowAbundant(tse, assay.type = "relabundance", abundance.th = 10/100)
#' rare |> head()
#' 
#' # Get rare taxa that are not permanently or conditionally rare
#' rare <- getLowAbundant(tse, assay.type = "relabundance", permanent.th = 5, condition.th = 100)
#' rare |> head()
#' 
#' # Get permanently rare taxa
#' prt <- getPermanentlyRare(tse, assay.type = "relabundance", permanent.th = 5)
#' prt |> head()
#' 
#' # Get conditionally rare taxa
#' prt <- getConditionallyRare(tse, assay.type = "relabundance", condition.th = 100)
#' prt |> head()
#' 
#' @seealso
#' \code{\link[=getPrevalence]{getPrevalent}} and
#' \code{\link[=getPrevalence]{getRare}}
#'
NULL

#' @rdname getAbundant
#' @export
setGeneric("getAbundant", signature = "x", function(x, ...)
    standardGeneric("getAbundant"))

#' @rdname getAbundant
#' @export
setGeneric("getLowAbundant", signature = "x", function(x, ...)
    standardGeneric("getLowAbundant"))

#' @rdname getAbundant
#' @export
setGeneric("getConditionallyRare", signature = "x", function(x, ...)
    standardGeneric("getConditionallyRare"))

#' @rdname getAbundant
#' @export
setGeneric("getPermanentlyRare", signature = "x", function(x, ...)
    standardGeneric("getPermanentlyRare"))

################################# getAbundant ##################################

#' @rdname getAbundant
#' @export
setMethod("getAbundant", signature = c(x = "TreeSummarizedExperiment"),
    function(x, ...){
        x <- .check_and_get_altExp(x, ...)
        res <- callNextMethod(x, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getAbundant", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts", ...){
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        res <- getAbundant(mat, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getAbundant", signature = c(x = "ANY"),
    function(x, abundance.th = 1/100, ...){
        if( is.null(rownames(x)) ){
            rownames(x) <- seq_len(nrow(x))
        }
        res <- .classify_by_abundance(x, abundance.th, ...)
        res <- rownames(x)[ which( res == "abundant" ) ]
        return(res)
    }
)

################################ getLowAbundant ################################

#' @rdname getAbundant
#' @export
setMethod("getLowAbundant", signature = c(x = "TreeSummarizedExperiment"),
    function(x, ...){
        x <- .check_and_get_altExp(x, ...)
        res <- callNextMethod(x, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getLowAbundant", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts", abundance.th = 1/100, ...){
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        res <- getLowAbundant(mat, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getLowAbundant", signature = c(x = "ANY"),
    function(x, abundance.th = 1/100, ...){
        if( is.null(rownames(x)) ){
            rownames(x) <- seq_len(nrow(x))
        }
        res <- .classify_by_abundance(x, abundance.th, ...)
        res <- rownames(x)[ which( res == "rare" ) ]
        return(res)
    }
)

############################# getConditionallyRare #############################

#' @rdname getAbundant
#' @export
setMethod("getConditionallyRare", signature = c(x = "TreeSummarizedExperiment"),
    function(x, ...){
        x <- .check_and_get_altExp(x, ...)
        res <- callNextMethod(x, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getConditionallyRare", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts", threshold = 1/100, ...){
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        res <- getConditionallyRare(mat, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getConditionallyRare", signature = c(x = "ANY"),
    function(
        x, abundance.th = 1/100, condition.th = 100, ...){
        if( !.is_a_numeric(condition.th) && condition.th > 0 ) stop(".")
        #
        if( is.null(rownames(x)) ){
            rownames(x) <- seq_len(nrow(x))
        }
        res <- .classify_by_abundance(x, abundance.th, condition.th = condition.th, ...)
        res <- rownames(x)[ which( res == "crt" ) ]
        return(res)
    }
)

############################## getPermanentlyRare ##############################

#' @rdname getAbundant
#' @export
setMethod("getPermanentlyRare", signature = c(x = "TreeSummarizedExperiment"),
    function(x, ...){
        x <- .check_and_get_altExp(x, ...)
        res <- callNextMethod(x, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getPermanentlyRare", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts", threshold = 1/100, ...){
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        res <- getPermanentlyRare(mat, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getPermanentlyRare", signature = c(x = "ANY"),
    function(
        x, abundance.th = 1/100, permanent.th = 5, ...){
        if( !.is_a_numeric(permanent.th) && permanent.th > 0) stop(".")
        #
        if( is.null(rownames(x)) ){
            rownames(x) <- seq_len(nrow(x))
        }
        res <- .classify_by_abundance(x, abundance.th, permanent.th = permanent.th, ...)
        res <- rownames(x)[ which( res == "prt" ) ]
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

.classify_by_abundance <- function(mat, abundance.th, condition.th = NULL, permanent.th = NULL, ...){
    if( (!is.null(condition.th) || !is.null(permanent.th)) && any(mat==0) ) stop("CRT or PRT cannot be calculated when zeroes.")
    #
    means <- rowMeans2(mat, ...)
    #
    res <- rep("abundant", length(means))
    res[ means <= abundance.th ] <- "rare"
    #
    ratio <- NULL
    if( !is.null(condition.th) || !is.null(permanent.th) ){
        ratio <- rowMaxs(mat, ...)/rowMins(mat, ...)
    }
    if( !is.null(condition.th) ){
        res[ res == "rare" & ratio > condition.th ] <- "crt"
    }
    if( !is.null(permanent.th) ){
        res[ res == "rare" & ratio <= permanent.th ] <- "prt"
    }
    return(res)
}
#############################################


