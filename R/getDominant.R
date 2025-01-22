#' Get dominant taxa
#'
#' These functions return information about the most dominant taxa in a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#'
#' @inheritParams agglomerateByVariable
#' @inheritParams getPrevalence
#'
#' @param name \code{Character scalar}. A name for the column of the
#' \code{colData} where results will be stored.
#' (Default: \code{"dominant_taxa"})
#'
#' @param other.name \code{Character scalar}. A name for features that are not
#' included in n the most frequent dominant features in the data.
#' (Default: \code{"Other"})
#'
#' @param group \code{Character scalar}. Defines a group. Must be one of the
#' columns from \code{rowData(x)}. (Default: \code{NULL})
#'
#' @param rank Deprecated. Use \code{group} instead.
#'
#' @param n \code{Numeric scalar}. The number of features that are the most
#' frequent
#' dominant features. Default is NULL, which defaults that each sample is
#' assigned
#' a dominant taxon. (Default: \code{NULL})
#'
#' @param complete \code{Logical scalar}. A value to manage multiple dominant
#' taxa for a sample.
#' Default for getDominant is TRUE to include all equally dominant taxa
#' for each sample. complete = FALSE samples one taxa for the samples that have
#' multiple.
#' Default for addDominant is FALSE to add a column with only one
#' dominant taxon assigned for each sample into colData. complete = TRUE adds a
#' list that includes all dominant taxa for each sample into colData.
#'
#' @param ... Additional arguments passed on to \code{agglomerateByRank()} when
#' \code{rank} is specified.
#'
#' @details
#' \code{addDominant} extracts the most abundant taxa in a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object, and stores the information in the \code{colData}. It is a wrapper for
#' \code{getDominant}.
#'
#' With \code{group} parameter, it is possible to agglomerate rows based on
#' groups. If the value is one of the columns in \code{taxonomyRanks()},
#' \code{agglomerateByRank()} is applied. Otherwise,
#' \code{agglomerateByVariable()} is utilized.
#' E.g. if 'Genus' rank is used, all abundances of same Genus
#' are added together, and agglomerated features are returned.
#' See corresponding functions for additional arguments to deal with
#' missing values or special characters.
#'
#' @return \code{getDominant} returns a named character vector \code{x}
#' while \code{addDominant} returns
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' with additional column in
#' \code{\link[SummarizedExperiment:colData]{colData}} named \code{*name*}.
#'
#' @name getDominant
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' x <- GlobalPatterns
#'
#' # Finds the dominant taxa.
#' sim.dom <- getDominant(x, group = "Genus")
#'
#' # Add information to colData
#' x <- addDominant(x, group = "Genus", name ="dominant_genera")
#' colData(x)
NULL

#' @rdname getDominant
#' @importFrom IRanges relist
#' @export
setMethod("getDominant", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = assay_name, assay_name = "counts", group = rank,
            rank = NULL, other.name = "Other", n = NULL, complete = TRUE, ...){
        # Input check
        # Check assay.type
        .check_assay_present(assay.type, x)
        # group check
        if(!is.null(group)){
            if(!.is_a_string(group)){
                stop("'group' must be an single character value.",
                    call. = FALSE)
            }
        }
        # If "group" is not NULL, species are aggregated according to the
        # taxonomic rank that is specified by user.
        if (!is.null(group) && group %in% taxonomyRanks(x)) {
            x <- agglomerateByRank(x, rank = group, ...)
        # or factor that is specified by user
        } else if (!is.null(group)) {
            x <- agglomerateByVariable(x, by = "rows", group = group, ...)
        }
        # Get assay
        mat <- assay(x, assay.type)
        # apply() function finds the indices of taxa's that has the highest
        # abundance.
        # rownames() returns the names of taxa that are the most abundant.
        idx <- as.list(apply(t(mat) == colMaxs(mat),1L,which))
        # Get rownames based on indices
        taxa <- rownames(mat)[unlist(idx)]

        # If multiple dominant taxa were found, names contain taxa in addition
        # to
        # sample name. Names are converted so that they include only sample
        # names.
        names(taxa) <- rep( names(idx), times = lengths(idx) )

        # If individual sample contains multiple dominant taxa (they have equal
        # counts) and if complete is FALSE, the an arbitrarily chosen dominant
        # taxa is returned
        if( length(taxa)>ncol(x) && !complete){
            # Store order
            order <- unique(names(taxa))
            # there are multiple dominant taxa in one sample (counts are equal),
            # length
            # of dominant is greater than rows in colData.
            taxa <- split(taxa, rep(names(taxa), lengths(taxa)) )
            # Order the data
            taxa <- taxa[order]
            names <- names(taxa)
            # If complete is set FALSE, and there are multiple dominant taxa,
            # one of them is arbitrarily chosen
            taxa <- lapply(taxa, function(item) {
                        return(sample(item, 1)) })
            taxa <- unname(unlist(lapply(taxa, function (x) {
                        unlist(x)})))
            names(taxa) <- names
            warning(
                "Multiple dominant taxa were found for some samples. ",
                "Use complete = TRUE for details.", call. = FALSE)
        }

        # Name "Other" the features that are not included in n the most abundant
        # in the data
        if(!is.null(n)){
            flat_taxa <- unlist(taxa, recursive = TRUE)
            top <- .top_n(flat_taxa, n=n)
            top <- names(top)
            # Group the rest into the "Other" category
            taxa <- lapply(taxa, function(x){
                ind <- vapply(x, function(y) y %in% top, FUN.VALUE = logical(1))
                if( any(ind) ){
                    res <- x[ ind ]
                } else{
                    res <- other.name
                }
                return(res)
            })
            if ( all(lengths(taxa) == 1 ) ){
                taxa <- unlist(taxa)
            }
        }
        return(taxa)
    }
)

#' @rdname getDominant
#' @export
setMethod("addDominant", signature = c(x = "SummarizedExperiment"),
    function(x, name = "dominant_taxa", other.name = "Other", n = NULL,
            complete = FALSE, ...) {
        # name check
        if(!.is_non_empty_string(name)){
            stop("'name' must be a non-empty single character value.",
                call. = FALSE)
        }
        # other.name check
        if(!.is_non_empty_string(other.name)){
            stop("'other.name' must be a non-empty single character value.",
                call. = FALSE)
        }
        dom.taxa <- getDominant(
            x, other.name = other.name, n = n, complete = complete, ...)
        # Add list into colData if there are multiple dominant taxa
        if(length(unique(names(dom.taxa))) < length(names(dom.taxa))) {
            # Store order
            order <- unique(names(dom.taxa))
            grouped <- split(dom.taxa, rep(names(dom.taxa)), lengths(dom.taxa))
            grouped <- grouped[order]
            dom.taxa <- grouped
            warning("A new column that was added in colData(x) is a list",
                    call. = FALSE)
        }
        colData(x)[[name]] <- dom.taxa
        return(x)
    }
)

########################## HELP FUNCTIONS summary ##############################

# top entries in a vector or given field in a data frame
.top_n <- function(x, n = NULL, na.rm = FALSE) {
    if (na.rm){
        inds <- which(x == "NA")
        if (length(inds) > 0){
            x[inds] <- NA
            warning("Interpreting NA string as missing value NA. Removing",
                length(inds), "entries", call. = FALSE)
        }
        x <- x[!is.na(x)]
    }
    # Create a frequency table of unique values of the dominant taxa for each
    # sample
    s <- rev(sort(table(x)))
    # Include only n the most frequent taxa
    if (!is.null(n)){
        s <- s[seq_len(min(n, length(s)))]
    }
    return(s)
}
