#' Estimate divergence
#'
#' Estimate divergence against a given reference sample.
#'
#' @inheritParams addDissimilarity
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#'
#' @param assay_name Deprecated. Use \code{assay.type} instead.
#'
#' @param reference \code{Character scalar}. A column name from
#' \code{colData(x)} or either \code{"mean"} or \code{"median"}.
#' (Default: \code{"median"})
#'
#' @param ... optional arguments passed to
#' \code{\link[=addDissimilarity]{addDissimilarity}}. Additionally:
#' \itemize{
#'   \item \code{dimred}: \code{Character scalar}. Specifies the name of
#'   dimension reduction result from \code{reducedDim(x)}. If used, these
#'   values are used to calculate divergence instead of the assay. Can be
#'   disabled with \code{NULL}. (Default: \code{NULL})
#' }
#'
#' @return \code{x} with additional
#' \code{\link[SummarizedExperiment:colData]{colData}} named \code{name}
#'
#' @details
#'
#' Microbiota divergence (heterogeneity / spread) within a given sample
#' set can be quantified by the average sample dissimilarity or beta
#' diversity with respect to a given reference sample.
#'
#' The calculation makes use of the function \code{getDissimilarity()}. The
#' divergence
#' measure is sensitive to sample size. Subsampling or bootstrapping can be
#' applied to equalize sample sizes between comparisons.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[=addAlpha]{addAlpha}}
#'   \item \code{\link[=addDissimilarity]{addDissimilarity}}
#'   \item \code{\link[scater:plotColData]{plotColData}}
#' }
#'
#' @name addDivergence
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # By default, reference is median of all samples. The name of column where
#' # results is "divergence" by default, but it can be specified.
#' tse <- addDivergence(tse)
#'
#' # The method that are used to calculate distance in divergence and
#' # reference can be specified. Here, euclidean distance is used. Reference is
#' # the first sample. It is recommended # to add reference to colData.
#' tse[["reference"]] <- rep(colnames(tse)[[1]], ncol(tse))
#' tse <- addDivergence(
#'     tse, name = "divergence_first_sample",
#'     reference = "reference",
#'     method = "euclidean")
#'
#' # Here we compare samples to global mean
#' tse <- addDivergence(tse, name = "divergence_average", reference = "mean")
#'
#' # All three divergence results are stored in colData.
#' colData(tse)
#'
NULL

#' @rdname addDivergence
#' @export
setMethod("addDivergence", signature = c(x="SummarizedExperiment"),
    function(x, name = "divergence", ...){
        ################### Input check ###############
        # Check name
        if( !.is_a_string(name) ){
            stop("'name' must be a non-empty character value.", call. = FALSE)
        }
        ################# Input check end #############
        # Calculate values
        res <- getDivergence(x, ...)
        # Add them to colData
        x <- .add_values_to_colData(x, list(res), name)
        return(x)
    }
)

#' @rdname addDivergence
#' @export
setMethod("getDivergence", signature = c(x="SummarizedExperiment"),
    function(
        x, assay.type = assay_name, assay_name = "counts",
        reference = "median", method = "bray", ...){
        ################### Input check ###############
        # Get altExp if user has specified
        x <- .check_and_get_altExp(x, ...)
        # Check assay.type
        .check_assay_present(assay.type, x)
        # Check reference
        ref_type <- .get_reference_type(reference, x)
        if( is.null(ref_type) ){
            stop(
                "'reference' must be a column from colData or either 'mean' ",
                "or 'median'.", call. = FALSE)
        }
        # If there are no colnames, add them. They are not added to returned
        # values; they are used just in calculation.
        if( is.null(colnames(x)) ){
            colnames(x) <- paste0("sample_", seq_len(ncol(x)))
        }
        ################# Input check end #############
        # Get assay and references
        args <- .get_matrix_and_reference(
            x, assay.type, reference, ref_type, ...)
        # Calculate sample-wise divergence
        args <- c(args, list(method = method), list(...))
        res <- do.call(.calc_divergence, args)
        # Get only values and ensure that their order is correct
        res <- res[match(colnames(x), res[["sample"]]), "value"]
        return(res)
        }
)
############################## HELP FUNCTIONS ##############################

# This function returns reference type.
# reference must be a column from colData, or either "median" or "mean".
# We also support providing a numeric vector or single sample name, but
# those are not recommended for user to not make the function too complex to
# use (too many options).
.get_reference_type <- function(reference, x){
    is_col <- .is_a_string(reference) && reference %in% colnames(colData(x)) &&
        all(!is.na(x[[reference]]) %in% colnames(x))
    is_mean_or_median <- .is_a_string(reference) && reference %in% c(
        "mean", "median")
    is_num_vector <- is.numeric(reference) && length(reference) == nrow(x)
    is_char_vector <- is.character(reference) && length(reference) == ncol(x) &&
        all(reference %in% colnames(x))
    is_sample <- .is_a_string(reference) && reference %in% colnames(x)
    #
    res <- NULL
    if( is_col ){
        res <- "colData_column"
    } else if(is_mean_or_median){
        res <- reference
    } else if( is_num_vector ){
        res <- "num_vector"
    } else if( is_char_vector ){
        res <- "char_vector"
    } else if( is_sample ){
        res <- "sample"
    }
    return(res)
}

# This function gets the abundance table along with reference information
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
.get_matrix_and_reference <- function(
        x, assay.type, reference, ref_type, dimred = NULL,
        ref.name = "temporal_reference_for_divergence", ...){
    #
    if( !.is_a_string(ref.name) ){
        stop("'ref.name' must be a single character value.", call. = FALSE)
    }
    #
    if( !is.null(dimred) && !is(x, "SingleCellExperiment") ){
        stop("If 'dimred' is specified, 'x' must be SingleCellExperiment.",
            call. = FALSE)
    }
    #
    if( !(is.null(dimred) || (
        (.is_a_string(dimred) && dimred %in% reducedDimNames(x)) ||
        .is_integer(dimred) && dimred > 0 && dimred <= length(reducedDims(x))
        ))){
        stop("'dimred' must be NULL or specify a name from reducedDimNames(x).",
            call.= FALSE)
    }
    # Get assay or reducedDim
    if( is.null(dimred) ){
        mat <- assay(x, assay.type)
    } else{
        mat <- t(reducedDim(x, dimred))
    }

    # If reference type is median or mean, calculate it
    if( ref_type %in% c("median", "mean") ){
        reference <- apply(mat, 1, ref_type)
    }
    # In case of numeric values, add them to matrix
    if( ref_type %in% c("num_vector", "median", "mean") ){
        reference <- matrix(reference)
        colnames(reference) <- ref.name
        mat <- cbind(mat, reference)
        reference <- ref.name
    }
    # In case of colData variable, get name reference samples from there
    if( ref_type %in% c("colData_column") ){
        reference <- x[[reference]]
    }
    # If the reference is only one sample, replicate it to cover all samples
    if( .is_a_string(reference) ){
        reference <- rep(reference, ncol(mat))
    }
    # Check that all reference samples are included in the data
    if( !all(reference %in% colnames(mat) | is.na(reference)) ){
        stop("All reference samples must be included in the data.",
            call. = FALSE)
    }
    # Return a list with matrix and reference samples for each sample
    res <- list(mat, reference)
    return(res)
}

# For each sample-pair, this function calculates dissimilarity.
#' @importFrom dplyr mutate
.calc_divergence <- function(mat, reference, method, ...){
    # Create sample-pair data.frame
    reference <- data.frame(sample = colnames(mat), reference = reference)
    # Exclude NA values
    reference <- reference[!is.na(reference$reference), ]
    # For dissimilarity calculation, the samples must be in rows
    mat <- t(mat)
    # Loop through sample-pairs
    temp <- t(reference) |> as.data.frame()
    temp <- lapply(temp, function(sample_pair){
        # Calculate dissimilarity between a sample pair
        temp <- mat[ sample_pair, ]
        temp <- getDissimilarity(temp, method, ...)
        # Get only the single value
        temp <- temp[[1]]
        return(temp)
    })
    # Add values to data.frame that holds sample pairs
    temp <- unlist(temp)
    reference[["value"]] <- temp
    return(reference)
}
