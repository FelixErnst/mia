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
#' These functions identify abundant and rare taxa in a dataset.
#' Abundant taxa are characterized by high average abundance across the dataset,
#' while rare taxa are characterized by consistently low abundance.
#'
#' Conditionally rare taxa exhibit variable abundance, being abundant in some
#' samples and rare in others. In contrast, permanently rare taxa consistently
#' maintain low abundance across all samples.
#'
#' \itemize{
#'   \item Abundant taxa: Taxa with an average abundance exceeding
#'   \code{abundant.th}.
#'
#'   \item Low abundant / rare taxa: Taxa with an average abundance not
#'   exceeding \code{abundant.th}. Optionally, if specified, they must also
#'   satisfy the condition \eqn{crt.th >=
#'   \frac{abundance_{max}}{abundance_{min}} > prt.th}.
#'
#'   \item Conditionally rare or low abundant taxa (CRT): Taxa with an average
#'   abundance not exceeding \code{abundant.th} and with a maximum-to-minimum
#'   abundance ratio (\eqn{\frac{abundance_{max}}{abundance_{min}}}) greater
#'   than \code{crt.th}.
#'
#'   \item Permanently rare or low abundant taxa (PRT): Taxa with an average
#'   abundance not exceeding \code{abundant.th} and with a maximum-to-minimum
#'   abundance ratio (\eqn{\frac{abundance_{max}}{abundance_{min}}}) less than
#'   or equal to \code{prt.th}.
#' }
#'
#' @return
#' For \code{getAbundant}, \code{getLowAbundant},
#' \code{getConditionallyLowAbundant}, and \code{getPermanentlyLowAbundant} a
#' \code{vector} of taxa. For \code{getAbudanceClass} a vector of abundance
#' classes for each feature. For \code{addAbudanceClass}, a
#' \code{SummarizedExperiment} object.
#'
#' @inheritParams addAlpha
#'
#' @param assay.type \code{Character scalar}. Specifies the name of assay
#' used in calculation. (Default: \code{"relabundance"})
#'
#' @param abundant.th \code{Numeric scalar}. Specifies threshold that is used
#' to separate abundant features from rare. (Default: \code{1/100})
#'
#' @param crt.th \code{Numeric scalar}. Specifies threshold that is used
#' to separate conditionally rare features from other rare features.
#' (Default: \code{100})
#'
#' @param prt.th \code{Numeric scalar}. Specifies threshold that is used
#' to separate permanently rare features from other rare features.
#' (Default: \code{5})
#'
#' @param name \code{Character scalar}. Specifies name of column in
#' \code{rowData} where the results will be stored.
#' (Default: \code{"abundance_class"})
#'
#' @param ... additional arguments.
#'
#' @examples
#'
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # Agglomerate to family level
#' tse <- agglomerateByRank(tse, rank = "Family")
#' # Transform to relative abundances. Note that we add pseudocount. This is
#' # because otherwise we cannot calculate CRT and PRT due to zeroes and
#' # zero division in calculating abundance ratio.
#' tse <- transformAssay(tse, method = "relabundance", pseudocount = TRUE)
#'
#' # Get abundant taxa
#' abundant <- getAbundant(tse, assay.type = "relabundance")
#' abundant |> head()
#'
#' # Get all rare taxa that have average relative abundance below 10%
#' rare <- getLowAbundant(
#'     tse, assay.type = "relabundance", abundant.th = 10/100)
#' rare |> head()
#'
#' # Get rare taxa that are not permanently or conditionally rare
#' rare <- getLowAbundant(
#'     tse, assay.type = "relabundance", prt.th = 5, crt.th = 100)
#' rare |> head()
#'
#' # Get permanently rare taxa
#' prt <- getPermanentlyLowAbundant(
#'     tse, assay.type = "relabundance", prt.th = 5)
#' prt |> head()
#'
#' # Get conditionally rare taxa
#' prt <- getConditionallyLowAbundant(
#'     tse, assay.type = "relabundance", crt.th = 100)
#' prt |> head()
#'
#' # To classify all features, one can use *AbundantClass function
#' tse <- addAbundanceClass(tse)
#' # When one uses add* function, the results are stored to rowData
#' rowData(tse)
#'
#' @seealso
#' \code{\link[=getPrevalence]{getPrevalent}} and
#' \code{\link[=getPrevalence]{getRare}}
#'
#' @references
#'
#' Sizhong Y. et al. (2017) Community structure of rare methanogenic archaea:
#' insight from a single functional group- _FEMS Microbiol. Ecol._ 93(11).
#' \url{https://doi.org/10.1093/femsec/fix126}
#'
NULL

################################# getAbundant ##################################

#' @rdname getAbundant
#' @export
setMethod("getAbundant", signature = c(x = "SingleCellExperiment"),
    function(x, ...){
        x <- .check_and_get_altExp(x, ...)
        res <- callNextMethod(x, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getAbundant", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "relabundance", ...){
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        res <- getAbundant(mat, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getAbundant", signature = c(x = "ANY"),
    function(x, abundant.th = 1/100, ...){
        if( is.null(rownames(x)) ){
            rownames(x) <- seq_len(nrow(x))
        }
        #
        res <- .classify_by_abundance(x, abundant.th, ...)
        res <- rownames(x)[ which( res == "abundant" ) ]
        return(res)
    }
)

################################ getLowAbundant ################################

#' @rdname getAbundant
#' @export
setMethod("getLowAbundant", signature = c(x = "SingleCellExperiment"),
    function(x, ...){
        x <- .check_and_get_altExp(x, ...)
        res <- callNextMethod(x, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getLowAbundant", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "relabundance", abundant.th = 1/100, ...){
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        res <- getLowAbundant(mat, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getLowAbundant", signature = c(x = "ANY"),
    function(x, abundant.th = 1/100, ...){
        if( is.null(rownames(x)) ){
            rownames(x) <- seq_len(nrow(x))
        }
        #
        res <- .classify_by_abundance(x, abundant.th, ...)
        res <- rownames(x)[ which( res == "rare" ) ]
        return(res)
    }
)

######################### getConditionallyLowAbundant ##########################

#' @rdname getAbundant
#' @export
setMethod("getConditionallyLowAbundant",
    signature = c(x = "SingleCellExperiment"),
    function(x, ...){
        x <- .check_and_get_altExp(x, ...)
        res <- callNextMethod(x, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getConditionallyLowAbundant",
    signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "relabundance", ...){
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        res <- getConditionallyLowAbundant(mat, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getConditionallyLowAbundant", signature = c(x = "ANY"),
    function(x, abundant.th = 1/100, crt.th = 100, ...){
        if( !.is_a_numeric(crt.th) && crt.th > 0 ){
            stop("'crt.th' must be positive numeric value.", call. = FALSE)
        }
        if( is.null(rownames(x)) ){
            rownames(x) <- seq_len(nrow(x))
        }
        #
        res <- .classify_by_abundance(x, abundant.th, crt.th = crt.th, ...)
        res <- rownames(x)[ which( res == "crt" ) ]
        return(res)
    }
)

########################## getPermanentlyLowAbundant ###########################

#' @rdname getAbundant
#' @export
setMethod("getPermanentlyLowAbundant",
    signature = c(x = "SingleCellExperiment"),
    function(x, ...){
        x <- .check_and_get_altExp(x, ...)
        res <- callNextMethod(x, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getPermanentlyLowAbundant",
    signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "relabundance", ...){
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        res <- getPermanentlyLowAbundant(mat, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getPermanentlyLowAbundant", signature = c(x = "ANY"),
    function(x, abundant.th = 1/100, prt.th = 5, ...){
        if( !.is_a_numeric(prt.th) && prt.th > 0 ){
            stop("'prt.th' must be positive numeric value.", call. = FALSE)
        }
        if( is.null(rownames(x)) ){
            rownames(x) <- seq_len(nrow(x))
        }
        #
        res <- .classify_by_abundance(x, abundant.th, prt.th = prt.th, ...)
        res <- rownames(x)[ which( res == "prt" ) ]
        return(res)
    }
)

############################### getAbundanceClass ##############################

#' @rdname getAbundant
#' @export
setMethod("getAbundanceClass", signature = c(x = "SingleCellExperiment"),
    function(x, ...){
        x <- .check_and_get_altExp(x, ...)
        res <- callNextMethod(x, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getAbundanceClass", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "relabundance", ...){
        .check_assay_present(assay.type, x)
        mat <- assay(x, assay.type)
        res <- getAbundanceClass(mat, ...)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("getAbundanceClass", signature = c(x = "ANY"),
    function(x, abundant.th = 1/100, crt.th = 100, prt.th = 5, ...){
        if( is.null(rownames(x)) ){
            rownames(x) <- seq_len(nrow(x))
        }
        #
        res <- .classify_by_abundance(x, abundant.th, crt.th, prt.th, ...)
        return(res)
        }
)

############################### addAbundanceClass ##############################

#' @rdname getAbundant
#' @export
setMethod("addAbundanceClass", signature = c(x = "SingleCellExperiment"),
    function(x, ...){
        x <- .check_and_get_altExp(x, ...)
        args <- c(list(x = x), list(...))
        args <- args[ !names(args) %in% "altexp" ]
        res <- do.call(callNextMethod, args)
        return(res)
    }
)

#' @rdname getAbundant
#' @export
setMethod("addAbundanceClass", signature = c(x = "SummarizedExperiment"),
    function(x, name = "abundance_class", ...){
        if( !.is_a_string(name) ){
            stop("'name' must be a single character value.", call. = FALSE)
        }
        res <- getAbundanceClass(x, ...) |> list()
        x <- .add_values_to_colData(x, res, name, MARGIN = 1L)
        return(x)
    }
)

################################ HELP FUNCTIONS ################################

# This function classifies the features based on their abundance. The result
# is a vector, an abundance class for each feature.
#' @importFrom DelayedMatrixStats rowMeans2 rowMaxs rowMins
.classify_by_abundance <- function(
        mat, abundant.th, crt.th = NULL, prt.th = NULL, ...){
    if( !(.is_a_numeric(abundant.th) && abundant.th > 0) ){
        stop("'abundant.th' must be positive numeric value.", call. = FALSE)
    }
    if( !( is.null(crt.th) || (.is_a_numeric(crt.th) && crt.th > 0)) ){
        stop("'crt.th' must be positive numeric value or NULL.", call. = FALSE)
    }
    if( !( is.null(prt.th) || (.is_a_numeric(prt.th) && prt.th > 0) ) ){
        stop("'prt.th' must be positive numeric value or NULL.", call. = FALSE)
    }
    if( (!is.null(crt.th) || !is.null(prt.th)) && any(mat==0) ){
        stop("CRT or PRT cannot be calculated when abundance contains zeroes. ",
            "Please consider adding pseudocount.", call. = FALSE)
    }
    # Classify to abundant and rare based on average abundance
    means <- rowMeans2(mat, ...)
    res <- rep("abundant", length(means))
    res[ means <= abundant.th ] <- "rare"
    # Rare features are further assigned based on abundance MAX/MIN ratio
    ratio <- rowMaxs(mat, ...)/rowMins(mat, ...)
    if( !is.null(crt.th) ){
        res[ res == "rare" & ratio > crt.th ] <- "crt"
    }
    if( !is.null(prt.th) ){
        res[ res == "rare" & ratio <= prt.th ] <- "prt"
    }
    return(res)
}
