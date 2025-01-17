#' Transform assay
#'
#' Variety of transformations for abundance data, stored in \code{assay}.
#' See details for options.
#'
#' @inheritParams getDominant
#' @inheritParams getDissimilarity
#'   
#' @param method \code{Character scalar}. Specifies the transformation
#'   method.
#' 
#' @param MARGIN \code{Character scalar}. Determines whether the
#'   transformation is applied sample (column) or feature (row) wise.
#'   (Default: \code{"samples"})
#' 
#' @param pseudocount \code{Logical scalar} or \code{numeric scalar}. 
#'   When \code{TRUE}, automatically adds half of the minimum positive 
#'   value of \code{assay.type} (missing values ignored by default: 
#'   \code{na.rm = TRUE}).
#'   When FALSE, does not add any pseudocount (pseudocount = 0).
#'   Alternatively, a user-specified numeric value can be added as pseudocount.
#'   (Default: \code{FALSE}).
#' 
#' @param name \code{Character scalar}. A name for the column of the 
#'   \code{colData} where results will be stored. (Default: \code{"method"})
#'   
#' @param altexp \code{Character vector} or \code{NULL}. Specifies the names 
#' of alternative experiments to which the transformation should also be 
#' applied. If \code{NULL}, the transformation is only applied to the main 
#' experiment. (Default: \code{altExpNames(x)}).
#'
#' @param ... additional arguments passed e.g. on to \code{vegan:decostand}.
#' \itemize{
#'   \item \code{reference}: \code{Character scalar}. Used to
#'   to fill reference sample's column in returned assay when calculating alr. 
#'   (Default: \code{NA})
#'   \item \code{ref_vals} Deprecated. Use \code{reference} instead.
#'   \item \code{percentile}: \code{Numeric scalar} or \code{NULL} (css). Used
#'   to set the  percentile value that calculates the scaling factors in the css
#'   normalization. If \code{NULL}, percentile is estimated from the data by
#'   calculating the portion of samples that exceed the \code{threshold}.
#'   (Default: \code{NULL})
#'   \item \code{scaling}: \code{Numeric scalar}. Adjusts the normalization
#'   scale  by dividing the calculated scaling factors, effectively changing 
#'   the magnitude of the normalized counts. (Default: \code{1000}).
#'   \item \code{threshold}: \code{Numeric scalar}. Specifies relative
#'   difference threshold and determines the first point where the relative
#'   change in  differences between consecutive quantiles exceeds this
#'   threshold. (Default: \code{0.1}).
#' }
#' @details
#'
#' \code{transformAssay} function provides a variety of options for
#' transforming abundance data. The transformed data is calculated and stored
#' in a new \code{assay}.
#'
#' The \code{transformAssay} provides sample-wise (column-wise) or feature-wise
#' (row-wise) transformation to the abundance table
#' (assay) based on specified \code{MARGIN}.
#' 
#' The available transformation methods include:
#'
#' \itemize{
#' 
#' \item 'alr', 'chi.square', 'clr', 'frequency', 'hellinger', 'log',
#' 'normalize', 'pa', 'rank', 'rclr' relabundance', 'rrank', 'standardize',
#' 'total': please refer to 
#' \code{\link[vegan:decostand]{decostand}} for details.
#' 
#' \item 'css': Cumulative Sum Scaling (CSS) can be used to normalize count data
#' by accounting for differences in library sizes. By default, the function
#' determines the normalization percentile for summing and scaling 
#' counts. If you want to specify the percentile value, good default value
#' might be \code{0.5}.The method is inspired by the CSS methods in
#' \code{\link[https://www.bioconductor.org/packages/metagenomeSeq/]{metagenomeSeq}}
#' package.
#' 
#' \item 'log10': log10 transformation can be used for reducing the skewness
#' of the data.
#' \deqn{log10 = \log_{10} x}{%
#' log10 = log10(x)}
#' where \eqn{x} is a single value of data.
#' 
#' \item 'log2': log2 transformation can be used for reducing the skewness of
#' the data.
#' \deqn{log2 = \log_{2} x}{%
#' log2 = log2(x)}
#' where \eqn{x} is a single value of data.
#'
#' }
#'
#' @return
#' \code{transformAssay} returns the input object \code{x}, with a new 
#' transformed abundance table named \code{name} added in the
#' \code{\link{assay}}.
#' 
#' @references
#'  
#' Paulson, J., Stine, O., Bravo, H. et al. (2013)
#' Differential abundance analysis for microbial marker-gene surveys
#' _Nature Methods_ 10, 1200–1202.
#' doi:10.1038/nmeth.2658
#'
#' @name transformAssay
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # By specifying 'method', it is possible to apply different transformations, 
#' # e.g. compositional transformation.
#' tse <- transformAssay(tse, method = "relabundance")
#' 
#' # The target of transformation can be specified with "assay.type"
#' # Pseudocount can be added by specifying 'pseudocount'.
#' 
#' # Perform CLR with smallest positive value as pseudocount
#' tse <- transformAssay(
#'     tse, assay.type = "relabundance", method = "clr", 
#'     pseudocount = TRUE
#'     )
#'                       
#' head(assay(tse, "clr"))
#' 
#' # Perform CSS normalization.
#' tse <- transformAssay(tse, method = "css")
#' head(assay(tse, "css"))
#' 
#' # With MARGIN, you can specify the if transformation is done for samples or
#' # for features. Here Z-transformation is done feature-wise.
#' tse <- transformAssay(tse, method = "standardize", MARGIN = "features")
#' head(assay(tse, "standardize"))
#' 
#' # Name of the stored table can be specified.
#' tse <- transformAssay(tse, method="hellinger", name="test")
#' head(assay(tse, "test"))
#'
#' # pa returns presence absence table.
#' tse <- transformAssay(tse, method = "pa")
#' head(assay(tse, "pa"))
#' 
#' # rank returns ranks of taxa.
#' tse <- transformAssay(tse, method = "rank")
#' head(assay(tse, "rank"))
#'
#' # In order to use other ranking variants, modify the chosen assay directly:
#' assay(tse, "rank_average", withDimnames = FALSE) <- colRanks(
#'     assay(tse, "counts"), ties.method = "average", preserveShape = TRUE)  
#' 
#' # Using altexp parameter. First agglomerate the data and then apply
#' # transformation.
#' tse <- GlobalPatterns
#' tse <- agglomerateByRanks(tse)
#' tse <- transformAssay(tse, method = "relabundance")
#' # The transformation is applied to all alternative experiments
#' altExp(tse, "Species")
#' 
NULL

#' @rdname transformAssay
#' @export
setMethod("transformAssay", signature = c(x = "SummarizedExperiment"),
    function(x,
        assay.type = "counts", assay_name = NULL,
        method = c("alr", "chi.square", "clr", "css", "frequency", 
            "hellinger", "log", "log10", "log2", "max", "normalize", 
            "pa", "range", "rank", "rclr", "relabundance", "rrank",
            "standardize", "total", "z"),
        MARGIN = "samples",
        name = method,
        pseudocount = FALSE,
        ...){
        #
        x <- .transform_assay(
            x = x, method = method, name = name, assay.type = assay.type,
            MARGIN = MARGIN, pseudocount =  pseudocount, ...)
        return(x)
    }
)

#' @rdname transformAssay
#' @export
setMethod("transformAssay", signature = c(x = "SingleCellExperiment"),
    function(x, altexp = altExpNames(x), ...){
        # Check altexp
        if( !(is.null(altexp) || all(altexp %in% altExpNames(x))) ){
            stop("'altexp' should be NULL or specify names from ",
                "altExpNames(x).", call. = FALSE)
        }
        #
        # Transform the main object
        x <- .transform_assay(x, ...)
        # Transform alternative experiments
        altExps(x)[altexp] <- lapply(altExps(x)[altexp], function(y){
            .transform_assay(y, ...)
        })
        return(x)
    }
)

########################### HELP FUNCTIONS #####################################

############################### .transform_assay ###############################
# A generic functon that takes SE as input and transforms the specified assay
# with specified method. By calling this generic function, we can apply same
# methods for SE and altExps of TreeSE.
.transform_assay <- function(
        x, assay.type = "counts", assay_name = NULL,
        method = c(
            "alr", "chi.square", "clr", "css", "frequency", 
            "hellinger", "log", "log10", "log2", "max", "normalize", 
            "pa", "range", "rank", "rclr", "relabundance", "rrank",
            "standardize", "total", "z"),
        MARGIN = "samples",
        name = method,
        pseudocount = FALSE,
        ...){
    # Input check
    if(!is.null(assay_name)){
        .Deprecated(old="assay_name", new="assay.type", "Now assay_name is 
                deprecated. Use assay.type instead.")
        assay.type <- assay_name
    }
    # Check assay.type
    .check_assay_present(assay.type, x)
    # Check name
    if(!.is_non_empty_string(name) || name == assay.type){
        stop("'name' must be a non-empty single character value and be ",
            "different from `assay.type`.", call. = FALSE)
    }
    # If method is not single string, user has not specified transform method,
    # or has given e.g. a vector
    if(!.is_non_empty_string(method)){
        stop("'method' must be a non-empty single character value.",
            call. = FALSE)
    }
    method <- match.arg(method, several.ok = FALSE)
    # Check that MARGIN is 1 or 2
    MARGIN <- .check_MARGIN(MARGIN)
    # Check pseudocount
    if( !.is_a_bool(pseudocount) && !(is.numeric(pseudocount) && 
            length(pseudocount) == 1 && pseudocount >= 0) ){
        stop("'pseudocount' must be TRUE, FALSE or a number equal to or ",
            "greater than 0.", call. = FALSE)
    }
    # Input check end
    # Get the method and abundance table
    method <- match.arg(method)
    assay <- assay(x, assay.type)
    # Apply pseudocount, if it is not 0 or FALSE
    assay <- .apply_pseudocount(assay, pseudocount, ...)
    # Store pseudocount value and set attr equal to NULL. The function above,
    # add the used pseudocount to attributes.
    pseudocount <- attr(assay, "pseudocount")
    attr(assay, "pseudocount") <- NULL
    # Calls help function that does the transformation
    # Help function is different for mia and vegan transformations
    if( method %in% c("log10", "log2", "css") ){
        transformed_table <- .apply_transformation(
            assay, method, MARGIN, ...)
    } else {
        transformed_table <- .apply_transformation_from_vegan(
            assay, method, MARGIN, ...)
    }
    # Add pseudocount info to transformed table
    attr(transformed_table, "parameters")$pseudocount <- pseudocount
    # Assign transformed table to assays
    assay(x, name, withDimnames = FALSE) <- transformed_table
    return(x)
}

##############################.apply_transformation#############################
# Help function for transformAssay, takes abundance table
# as input and returns transformed table. This function utilizes mia's
# transformation functions.
.apply_transformation <- function(assay, method, MARGIN, ...){
    # Transpose if MARGIN is row
    if( MARGIN == 1L ){
        assay <- t(assay)
    }
    # Function is selected based on the "method" variable
    FUN <- switch(
        method,
        log10 = .calc_log,
        log2 = .calc_log,
        css = .calc_css
    )
    # Get transformed table
    transformed_table <- do.call(
        FUN, list(mat = assay, method = method, ...) )
    # Transpose back to normal if MARGIN is row
    if( MARGIN == 1L ){
        transformed_table <- t(transformed_table)
    }
    # Add method and margin to attributes
    attr(transformed_table, "mia") <- method
    attr(transformed_table, "parameters")$margin <- MARGIN
    return(transformed_table)
}

########################.apply_transformation_from_vegan########################
# Help function for transformAssay, takes abundance
# table as input and returns transformed table. This function utilizes vegan's
# transformation functions.
.apply_transformation_from_vegan <- function(
        mat, method, MARGIN, reference = ref_vals, ref_vals = NA, ...){
    # Input check
    # Check reference
    if( length(reference) != 1 ){
        stop("'reference' must be a single value specifying the ",
            "values of the reference sample.", call. = FALSE)
    }
    # Input check end
    # Ensure that the matrix has proper dimnames
    if (is.null(rownames(mat))) {
        rownames(mat) <- paste0("feature", seq_len(nrow(mat)))
    }
    if (is.null(colnames(mat))) {
        colnames(mat) <- paste0("sample", seq_len(ncol(mat)))
    }
    # Adjust method if mia-specific alias was used
    method <- ifelse(method == "relabundance", "total", method)
    if (method == "z") {
        .Deprecated(old="z", new="standardize")
    }
    method <- ifelse(method == "z", "standardize", method)
    
    # If method is ALR, vegan drops one column/sample, because it is used
    # as a reference. To work with TreeSE, reference sample must be added back.
    # Get the original order of samples/features
    orig_dimnames <- dimnames(mat)

    # Call vegan::decostand and apply transformation
    transformed_table <- vegan::decostand(
        mat, method = method, MARGIN = MARGIN, ...)
    
    # Add reference sample back if ALR
    if( method %in% c("alr") ){
        transformed_table <- .adjust_alr_table(
            mat = transformed_table, orig_dimnames = orig_dimnames,
            reference = reference)
    }
    # If table is transposed (like in chi.square), transpose back
    if(identical(rownames(transformed_table), colnames(mat)) &&
        identical(colnames(transformed_table), rownames(mat)) &&
        ncol(transformed_table) != ncol(mat) &&
        nrow(transformed_table != nrow(mat))){
        transformed_table <- t(transformed_table)
    }
    return(transformed_table)
}

####################################.calc_log###################################
# This function applies log transformation to abundance table.
.calc_log <- function(mat, method, ...){
    # If abundance table contains zeros or negative values, gives an error,
    # because it is not possible to calculate log from zeros. Otherwise,
    # calculates log.
    if ( any(mat < 0, na.rm = TRUE) ){
        stop("The assay contains negative values and ", method,
            " transformation is being applied without pseudocount.",
            "`pseudocount` must be specified manually.", call. = FALSE)
    } else if ( any(mat == 0, na.rm = TRUE) ){
        stop("The assay contains zeroes and ", method,
            " transformation is being applied without pseudocount.",
            "`pseudocount` must be set to TRUE.", call. = FALSE)
    }
    # Calculate log2 or log10 abundances
    if(method == "log2"){
        mat <- log2(mat)
    } else{
        mat <- log10(mat)
    }
    # Add parameter to attributes
    attr(mat, "parameters") <- list()
    return(mat)
}

################################### .calc_css ##################################
# This function applies cumulative sum scaling (CSS) to the abundance table.
.calc_css <- function(mat, percentile = NULL, scaling = 1000, ...) {
    # Input check
    if( !(is.null(percentile) || .is_a_numeric(percentile)) ){
        stop("'percentile' must be a numeric value or NULL.", call. = FALSE)
    }
    if( !.is_a_numeric(scaling) ){
        stop("'scaling' must be an integer value.", call. = FALSE)
    }
    #
    # If 'percentile' is not provided, calculate it dynamically based on the
    # abundance matrix. This allows flexibility in the function by determining
    # where the low-abundance features start to diverge,
    # which is crucial for CSS normalization.
    if( is.null(percentile) ){
        percentile <- .calc_css_percentile(mat, ...)
        message("'percentile' set to: ", percentile)
    }
    # Calculate the scaling factors for each sample based on the provided or
    # calculated percentile. These factors determine how much of the data
    # (below the percentile threshold) is used to normalize each sample.
    scaling_factors <- .calc_scaling_factors(mat, percentile)
    # Normalize the count data by dividing each sample's counts by its
    # corresponding scaling factor. This step adjusts for varying distributions
    # of low-abundance features across samples,
    # making the counts more comparable by removing biases from varying library
    # sizes or sample depth.
    scaled_factors <- scaling_factors / scaling
    normalized_data <- sweep(mat, 2, scaled_factors, "/")
    # Add information on used parameters to normalized matrix
    attr(normalized_data, "parameters") <- list(
        percentile = percentile,
        scaling_factors = scaling_factors,
        scaling = scaling
        )
    return(normalized_data)
}

############################# .calc_css_percentile #############################
# Calculates the cumulative sum scaling (CSS) percentile from the input data.
# The percentile determines the point at which a sample's abundance values start
# to deviate significantly from the average abundance profile of all samples.
#' @importFrom DelayedMatrixStats colSums2 colQuantiles rowMeans2 rowMedians
.calc_css_percentile <- function(mat, threshold = 0.1, ...) {
    # Input check
    if( !(.is_a_numeric(threshold) && threshold > 0 && threshold < 1)  ){
        stop("'threshold' must be a numeric value between 0 and 1.",
            call. = FALSE)
    }
    #
    # Replace zero values to NA, i.e. not detected. That ensures that they are
    # not included in the percentile calculation; only non-zero values are.
    mat[ mat == 0 ] <- NA
    # Calculate number of features detected
    found_features <- colSums2(!is.na(mat))
    # Stop if there are columns with only 1 or 0 found features. The percentile
    # estimation is not reliable with so few features.
    if( any(found_features <= 1) ){
        stop(
            "There are samples that contain 1 or less found features.",
            "'percentile' cannot be estimated from the data. Specify it ",
            "manually.", call. = FALSE)
    }
    # Calculate quantiles for the features present in each sample.
    # Use a 'probs' sequence from 0 to 1, with increments based on the maximum
    # number of features found. This allows quantile calculation for every
    # feature increment.
    # The quantiles gives us a detailed view of the abundance distribution in
    # each sample, helping to locate low-abundance regions.
    quantiles <- colQuantiles(
        mat, probs = seq(0, 1, length.out = max(found_features)), na.rm = TRUE)
    quantiles <- t(quantiles)
    # Sort the columns based on abundance. This means that each sample is sorted
    # so that smallest values come first. Ultimately, this means that values
    # should increase from first to row to last. Rows do not represent single
    # taxa anymore but abundance level.
    # This helps compare the accumulation of low-abundance features across
    # samples.
    mat <- as.matrix(mat)
    mat <- apply(mat, 2, function(col){
        col[ is.na(col) ] <- 0
        col <- sort(col)
        return(col)
    })
    # Compute a reference profile as the average abundance at each abundance
    # level across samples. This creates a baseline for comparison, representing
    # the "average" behavior of abundance in the dataset.
    ref <- rowMeans2(mat)
    ref <- ref[ ref > 0 ]
    # Calculate the differences between the reference profile and each sample's
    # quantiles. This allows us to detect how much a sample's abundance deviates
    # from the average abundance profile.
    difference <- ref - quantiles
    # Calculate the median absolute differences between the reference and sample
    # profiles. This helps to quantify the typical deviation in abundance
    # accumulation across samples.
    difference <- rowMedians(abs(difference))
    # Compute the relative change in absolute differences between the reference
    # profile and each sample's quantiles. Identify the first index where this
    # relative change exceeds the specified threshold. This index represents
    # the count value at which the deviation between sample abundance levels and
    # the reference profile becomes significant, highlighting a notable
    # difference in abundance between samples.
    res <- abs(diff(difference)) / difference[-1]
    res <- which(res > threshold)
    res <- res[[1]]
    # Convert the count value to a relative proportion by dividing by the
    # maximum number of features. This scales the count to a proportion of the
    # total feature count, normalizing it against the library size.
    res <- res/max(found_features)
    return(res)
}

############################### .calc_scaling_factors ##########################
# Calculates the cumulative sum scaling (CSS) factors for each sample based on
# the given percentile. The scaling factors determine the cumulative sum of
# counts below the percentile for each sample, which is used for normalizing
# the count data.
#' @importFrom DelayedMatrixStats colQuantiles
.calc_scaling_factors <- function(mat, percentile) {
    # Replace zero values in the matrix with NA to indicate missing or not
    # detected values.
    mat_tmp <- mat
    mat_tmp[mat_tmp == 0] <- NA
    # Calculate the abundance threshold for each sample at the specified
    # percentile. This threshold indicates the abundance value below which the
    # given proportion of counts falls in the sample.
    quantiles <- colQuantiles(mat_tmp, probs = percentile, na.rm = TRUE)
    # For each sample, sum the counts that are less than or equal to the
    # percentile threshold. This cumulative sum is used as the scaling factor,
    # adjusting the sample's abundance to account for low-abundance features.
    # Improve precision by subtracting the smallest positive floating-point
    # number from each value before comparison.
    scaling_factors <- vapply(seq_len(ncol(mat_tmp)), function(i){
        col_values <- mat[, i] - .Machine$double.eps
        sum(col_values[col_values <= quantiles[i]], na.rm = TRUE)
    }, numeric(1))
    # Give scaling_factors names corresponding to the original matrix columns
    names(scaling_factors) <- colnames(mat)
    return(scaling_factors)
}

#################################.calc_rel_abund################################
# This function is for other functions to use internally.
.calc_rel_abund <- function(mat, ...){
    mat <- .apply_transformation_from_vegan(
        mat, method = "relabundance", MARGIN = 2)
    return(mat)
}

###############################.adjust_alr_table################################
# vegan::decostand returns ALR transformed abundance table without reference
# sample. Because in TreeSE all assays must have same row and column names,
# the reference sample is assigned back to transformed abundance table.
.adjust_alr_table <- function(mat, orig_dimnames, reference){
    # Store attributes
    attributes <- attributes(mat)
    # Get original and current sample/feature names and dimensions of reference
    # based on what dimensions misses names
    MARGIN <- ifelse(length(orig_dimnames[[1]]) == nrow(mat), 2, 1)
    orig_names <- if(MARGIN == 1){orig_dimnames[[1]]} else {orig_dimnames[[2]]}
    current_names <- if(MARGIN == 1){rownames(mat)} else {colnames(mat)}
    nrow <- ifelse(MARGIN == 1, 1, nrow(mat))
    ncol <- ifelse(MARGIN == 1, ncol(mat), 1)
    # Get the name of reference sample/feature and names of other dimension
    reference_name <- setdiff(orig_names, current_names)
    var_names <- if(MARGIN == 1){colnames(mat)} else {rownames(mat)}
    if(MARGIN == 1){
        ref_dimnames <- list(reference_name, var_names)
    } else {
        ref_dimnames <- list(var_names, reference_name)
        }
    # Reference sample as NAs or with symbols that are specified by user
    reference_sample <- matrix(
        reference, nrow = nrow, ncol = ncol, dimnames = ref_dimnames)
    # Add reference sample/feature
    if(MARGIN == 1){
        mat <- rbind(mat, reference_sample)
        # Preserve the original order
        mat <- mat[orig_names, ]
    } else {
        mat <- cbind(mat, reference_sample)
        # Preserve the original order
        mat <- mat[, orig_names]
    }
    # Add those attributes that were related to calculation
    attributes(mat) <- c(
        attributes(mat),
        attributes[ !names(attributes) %in% c("dim", "dimnames") ])
    return(mat)
}


###############################.apply_pseudocount###############################
# This function applies pseudocount to abundance table.
.apply_pseudocount <- function(mat, pseudocount, na.rm = TRUE, ...){
    if( .is_a_bool(pseudocount) ){
        # If pseudocount TRUE and some NAs, a warning is issued
        if ( pseudocount && any(is.na(mat)) ){
            warning("The assay contains missing values (NAs). These will be ",
                "ignored in pseudocount calculation.", call. = FALSE)
        }    
        # If pseudocount TRUE but some negative values, numerical pseudocount
        # needed
        if ( pseudocount && any(mat < 0, na.rm = TRUE) ){
            stop("The assay contains negative values. ",
                "'pseudocount' must be specified manually.", call. = FALSE)
        }
        # If pseudocount TRUE, set it to half of non-zero minimum value
        # else set it to zero.
        # Get min value
        value <- min(mat[mat > 0], na.rm = na.rm)
        value <- value/2
        pseudocount <- ifelse(pseudocount, value, 0)
        # Report pseudocount if positive value
        if ( pseudocount > 0 ){
            message("A pseudocount of ", pseudocount, " was applied.")
        }
    }
    # Give warning if pseudocount should not be added
    # Case 1: only positive values
    if( pseudocount != 0 && all(mat > 0, na.rm = TRUE) ){
        warning(
            "The assay contains only positive values. ",
            "Applying a pseudocount may be unnecessary.", call. = FALSE)
    }
    # Case 2: some negative values
    if( pseudocount != 0 && any(mat < 0, na.rm = TRUE) ){
        warning(
            "The assay contains some negative values. ",
            "Applying a pseudocount may produce meaningless data.",
            call. = FALSE)
    }
    # Add pseudocount
    mat <- mat + pseudocount
    # Set attr equal to pseudocount
    attr(mat, "pseudocount") <- pseudocount
    return(mat)
}
