#' Calculation prevalence information for features across samples
#'
#' These functions calculate the population prevalence for taxonomic ranks in a
#' \code{\link{SummarizedExperiment-class}} object.
#'
#' @inheritParams getDissimilarity
#' 
#' @param x
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}.
#' 
#' @param assay_name Deprecated. Use \code{assay.type} instead.
#'
#' @param detection \code{Numeric scalar}. Detection threshold for
#' absence/presence. If \code{as_relative = FALSE},
#' it sets the counts threshold for a taxon to be considered present.
#' If \code{as_relative = TRUE}, it sets the relative abundance threshold
#' for a taxon to be considered present. (Default: \code{0})
#'
#' @param include.lowest \code{Logical scalar}. Should the lower boundary of the
#'   detection and prevalence cutoffs be included? (Default: \code{FALSE})
#' 
#' @param include_lowest Deprecated. Use \code{include.lowest} instead.
#'
#' @param sort \code{Logical scalar}. Should the result be sorted by prevalence?
#'   (Default: \code{FALSE})
#'
#' @param rank \code{Character scalar}. Defines a taxonomic rank. Must be a
#' value of \code{taxonomyRanks()} function.
#'
#' @param na.rm \code{Logical scalar}. Should NA values be omitted?
#' (Default: \code{TRUE})
#' 
#' @param update.tree \code{Logical scalar}. Should
#' \code{rowTree()} also be agglomerated? (Default: \code{FALSE})
#'
#' @param ... additional arguments
#' \itemize{
#'   \item If \code{!is.null(rank)} arguments are passed on to
#'   \code{\link[=agglomerate-methods]{agglomerateByRank}}. See
#'   \code{\link[=agglomerate-methods]{?agglomerateByRank}} for more details.
#'
#'   \item for \code{getPrevalent}, \code{getRare}, \code{subsetByPrevalent}
#'   and \code{subsetByRare} additional parameters passed to
#'   \code{getPrevalence}
#'
#'   \item for \code{getPrevalentAbundance} additional parameters passed to
#'   \code{getPrevalent}
#' }
#'
#' @details
#' \code{getPrevalence} calculates the frequency of samples that exceed
#' the detection threshold. For \code{SummarizedExperiment} objects, the
#' prevalence is calculated for the selected taxonomic rank, otherwise for the
#' rows. The absolute population prevalence can be obtained by multiplying the
#' prevalence by the number of samples (\code{ncol(x)}). 
#'
#' The core abundance index from \code{getPrevalentAbundance} gives the relative
#' proportion of the core species (in between 0 and 1). The core taxa are
#' defined as those that exceed the given population prevalence threshold at the
#' given detection level as set for \code{getPrevalent}.
#'
#' \code{subsetPrevalent} and \code{subsetRareFeatures} return a subset of
#' \code{x}.
#' The subset includes the most prevalent or rare taxa that are calculated with
#' \code{getPrevalent} or \code{getRare} respectively.
#'
#' @return
#' \code{subsetPrevalent} and \code{subsetRareFeatures} return subset of
#' \code{x}.
#'
#' All other functions return a named vectors:
#' \itemize{
#'   \item \code{getPrevalence} returns a \code{numeric} vector with the
#'     names being set to either the row names of \code{x} or the names after
#'     agglomeration.
#'
#'   \item \code{getPrevalentAbundance} returns a \code{numeric} vector with
#'     the names corresponding to the column name of \code{x} and include the
#'     joint abundance of prevalent taxa.
#'
#'   \item \code{getPrevalent} and \code{getRare} return a
#'     \code{character} vector with only the names exceeding the threshold set
#'     by \code{prevalence}, if the \code{rownames} of \code{x} is set.
#'     Otherwise an \code{integer} vector is returned matching the rows in
#'     \code{x}.
#' }
#'
#' @seealso
#' \code{\link[=agglomerate-methods]{agglomerateByRank}},
#' \code{\link[=getTop]{getTop}}
#'
#'
#' @name getPrevalence
#' @export
#'
#' @references
#' A Salonen et al. The adult intestinal core microbiota is determined by
#' analysis depth and health status. Clinical Microbiology and Infection
#' 18(S4):16 20, 2012.
#' To cite the R package, see citation('mia')
#'
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' # Get prevalence estimates for individual ASV/OTU
#' prevalence.frequency <- getPrevalence(tse,
#'                                       detection = 0,
#'                                       sort = TRUE)
#' head(prevalence.frequency)
#'
#' # Get prevalence estimates for phyla
#' # - the getPrevalence function itself always returns population frequencies
#' prevalence.frequency <- getPrevalence(tse,
#'                                       rank = "Phylum",
#'                                       detection = 0,
#'                                       sort = TRUE)
#' head(prevalence.frequency)
#'
#' # - to obtain population counts, multiply frequencies with the sample size,
#' # which answers the question "In how many samples is this phylum detectable"
#' prevalence.count <- prevalence.frequency * ncol(tse)
#' head(prevalence.count)
#'
#' # Detection threshold 1 (strictly greater by default);
#' # Note that the data (GlobalPatterns) is here in absolute counts
#' # (and not compositional, relative abundances)
#' # Prevalence threshold 50 percent (strictly greater by default)
#' prevalent <- getPrevalent(
#'     tse,
#'     rank = "Phylum",
#'     detection = 10,
#'     prevalence = 50/100)
#' head(prevalent)
#'
#' # Add relative aundance data
#' tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")
#'
#' # Gets a subset of object that includes prevalent taxa
#' altExp(tse, "prevalent") <- subsetByPrevalent(tse,
#'                                              rank = "Family",
#'                                              assay.type = "relabundance",
#'                                              detection = 0.001,
#'                                              prevalence = 0.55)
#' altExp(tse, "prevalent")
#'
#' # getRare returns the inverse
#' rare <- getRare(tse,
#'     rank = "Phylum",
#'     assay.type = "relabundance",
#'     detection = 1/100,
#'     prevalence = 50/100)
#' head(rare)
#'
#' # Gets a subset of object that includes rare taxa
#' altExp(tse, "rare") <- subsetByRare(
#'     tse,
#'     rank = "Class",
#'     assay.type = "relabundance",
#'     detection = 0.001,
#'     prevalence = 0.001)
#' altExp(tse, "rare")
#'
#' # Names of both experiments, prevalent and rare, can be found from slot
#' # altExpNames
#' tse
#'
#' data(esophagus)
#' getPrevalentAbundance(esophagus, assay.type = "counts")
#'
NULL

#' @rdname getPrevalence
#' @export
setMethod("getPrevalence", signature = c(x = "ANY"), function(
    x, detection = 0, include.lowest = include_lowest, include_lowest = FALSE, 
    sort = FALSE, na.rm = TRUE, ...){
        # input check
        if (!.is_numeric_string(detection)) {
            stop("'detection' must be a single numeric value or coercible to ",
                "one.", call. = FALSE)
        }
        #
        if(!.is_a_bool(na.rm)){
            stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        detection <- as.numeric(detection)
        if(!.is_a_bool(include.lowest)){
            stop("'include.lowest' must be TRUE or FALSE.", call. = FALSE)
        }
        if(!.is_a_bool(sort)){
            stop("'sort' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        # Give warning if there are taxa with NA values
        if( any( is.na(x) ) ){
            msg <- paste0(
                "The abundance table contains NA values and they are",
                ifelse(na.rm, " ", " not "), "excluded (see 'na.rm').")
            warning(msg, call. = FALSE)
        }
        #
        if (include.lowest) {
            prev <- x >= detection
        } else {
            prev <- x > detection
        }
        # Calculate prevalence for each taxa
        prev <- rowSums(prev, na.rm = na.rm)
        # Always return prevalence as a relative frequency.
        # This helps to avoid confusion with detection limit
        prev <- prev / ncol(x)
        if (sort) {
            prev <- rev(sort(prev))
        }
        prev
    }
)

#' @rdname getPrevalence
#' @export
setMethod("getPrevalence", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = assay_name, assay_name = "counts",
            rank = NULL, ...){
        # check assay
        .check_assay_present(assay.type, x)
        x <- .merge_features(x, rank = rank, ...)
        mat <- assay(x, assay.type)
        # Calculate abundance
        mat <- .to_rel_abund(mat, ...)
        getPrevalence(mat, ...)
    }
)
############################# getPrevalent #################################
#' @name getPrevalence
#'
#' @param prevalence Prevalence threshold (in 0 to 1). The
#'   required prevalence is strictly greater by default. To include the
#'   limit, set \code{include.lowest} to \code{TRUE}.
#'
#' @details
#' \code{getPrevalent} returns taxa that are more prevalent with the
#' given detection threshold for the selected taxonomic rank.
#'
#' @aliases getPrevalent
#'
#' @export
NULL

.norm_rownames <- function(x){
    if(is.null(rownames(x))){
        rownames(x) <- seq_len(nrow(x))
    } else if(anyDuplicated(rownames(x))) {
        rownames(x) <- make.unique(rownames(x))
    }
    x
}

.get_prevalent_indices <- function(x, prevalence = 50/100,
                                include.lowest = FALSE, ...){
    # input check
    if (!.is_numeric_string(prevalence)) {
        stop("'prevalence' must be a single numeric value or coercible to ",
            "one.", call. = FALSE)
    }

    prevalence <- as.numeric(prevalence)
    if(!.is_a_bool(include.lowest)){
        stop("'include.lowest' must be TRUE or FALSE.", call. = FALSE)
    }
    # rownames must bet set and unique, because if sort = TRUE, the order is
    # not preserved
    x <- .norm_rownames(x)
    pr <- getPrevalence(x, rank = NULL, ...)

    # get logical vector which row does exceed threshold
    if (include.lowest) {
        f <- pr >= prevalence
    } else {
        f <- pr > prevalence
    }
    # get it back into order of x
    m <- match(rownames(x),names(f))
    taxa <- f[m]
    # Gets indices of most prevalent taxa
    indices <- which(taxa)
    # revert the order based on f
    m <- match(names(f),names(indices))
    m <- m[!is.na(m)]
    indices <- indices[m]
    #
    indices
}

.get_prevalent_taxa <- function(x, rank = NULL, ...){
    if(is(x,"SummarizedExperiment")){
        x <- .merge_features(x, rank = rank, ...)
    }
    indices <- .get_prevalent_indices(x, ...)
    # If named input return named output
    if( !is.null(rownames(x)) ){
        # Gets the names
        taxa <- rownames(x)[indices]
    } else {
        # Otherwise indices are returned
        taxa <- unname(indices)
    }
    unique(taxa)
}

#' @rdname getPrevalence
#' @export
setMethod("getPrevalent", signature = c(x = "ANY"),
    function(x, prevalence = 50/100, include.lowest = include_lowest, 
        include_lowest = FALSE, ...){
            .get_prevalent_taxa(x, rank = NULL, prevalence = prevalence,
                include.lowest = include.lowest, ...)
    }
)

#' @rdname getPrevalence
#' @export
setMethod("getPrevalent", signature = c(x = "SummarizedExperiment"),
    function(x, rank = NULL, prevalence = 50/100,
        include.lowest = include_lowest, include_lowest = FALSE, ...){
            .get_prevalent_taxa(x, rank = rank, prevalence = prevalence,
                include.lowest = include.lowest, ...)
    }
)

############################# getRare ######################################

#' @name getPrevalence
#'
#' @details
#' \code{getRare} returns complement of \code{getPrevalent}.
#'
#' @export
NULL

.get_rare_indices <- function(x, ...){
    indices <- .get_prevalent_indices(x = x, ...)
    # reverse the indices
    indices_x <- seq_len(nrow(x))
    f <- !(indices_x %in% indices)
    indices_new <- indices_x[f]
    indices_new
}

.get_rare_taxa <- function(x, rank = NULL, ...){
    if(is(x,"SummarizedExperiment")){
        x <- .merge_features(x, rank = rank, ...)
    }
    indices <- .get_rare_indices(x, ...)
    #
    if( !is.null(rownames(x)) ){
        # Gets the names
        taxa <- rownames(x)[indices]
    } else {
        # Otherwise indices are returned
        taxa <- indices
    }
    unique(taxa)
}

#' @rdname getPrevalence
#' @export
setMethod("getRare", signature = c(x = "ANY"),
    function(x, prevalence = 50/100, include.lowest = include_lowest, 
        include_lowest = FALSE, ...){
            .get_rare_taxa(x, rank = NULL, prevalence = prevalence,
                include.lowest = include.lowest, ...)
    }
)

#' @rdname getPrevalence
#' @export
setMethod("getRare", signature = c(x = "SummarizedExperiment"),
    function(x, rank = NULL, prevalence = 50/100,
        include.lowest = include_lowest, include_lowest = FALSE, ...){
            .get_rare_taxa(x, rank = rank, prevalence = prevalence,
                include.lowest = include.lowest, ...)
    }
)

############################# subsetByPrevalent ############################

#' @rdname getPrevalence
#' @export
setMethod("subsetByPrevalent", signature = c(x = "SummarizedExperiment"),
    function(x, rank = NULL, ...){
        x <- .merge_features(x, rank = rank, ...)
        prevalent_indices <- .get_prevalent_indices(x, ...)
        x[prevalent_indices, ]
    }
)

#' @rdname getPrevalence
#' @export
setMethod("subsetByPrevalent", signature = c(x = "TreeSummarizedExperiment"),
    function(x, update.tree = FALSE, ...){
        # Check that update.tree is logical value
        if( !.is_a_bool(update.tree) ){
            stop("'update.tree' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        x <- callNextMethod(x, ...)
        # Agglomerate tree if specified
        if( update.tree ){
            x <- .agglomerate_trees(x, ...)
        }
    return(x)
    }
)

############################# subsetByRare #################################

#' @rdname getPrevalence
#' @export
setMethod("subsetByRare", signature = c(x = "SummarizedExperiment"),
    function(x, rank = NULL, ...){
        x <- .merge_features(x, rank = rank, ...)
        rare_indices <- .get_rare_indices(x, ...)
        x[rare_indices, ]
    }
)

#' @rdname getPrevalence
#' @export
setMethod("subsetByRare", signature = c(x = "TreeSummarizedExperiment"),
    function(x, update.tree = FALSE, ...){
        # Check that update.tree is logical value
        if( !.is_a_bool(update.tree) ){
            stop("'update.tree' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        x <- callNextMethod(x, ...)
        # Agglomerate tree if specified
        if( update.tree ){
            x <- .agglomerate_trees(x, ...)
        }
        return(x)
    }
)

############################# getPrevalentAbundance ############################

#' @rdname getPrevalence
#' @export
setMethod("getPrevalentAbundance", signature = c(x = "ANY"),
    function(x, ...){
        x <- .calc_rel_abund(x)
        cm <- getPrevalent(x, ...)
        if (length(cm) == 0) {
            stop("With the given abundance and prevalence thresholds, no taxa ",
                "were found. Try to change detection and prevalence ",
                "parameters.", call. = FALSE)
        }
        colSums(x[cm, ,drop=FALSE])
    }
)

#' @rdname getPrevalence
#' @export
setMethod("getPrevalentAbundance", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = assay_name, assay_name = "counts", ...){
        # check assay
        .check_assay_present(assay.type, x)
        #
        getPrevalentAbundance(assay(x,assay.type))
    }
)


############################# agglomerateByPrevalence ##########################

#' Agglomerate data based on population prevalence
#' 
#' @name agglomerateByPrevalence
#'  
#' @inheritParams agglomerateByRank
#' 
#' @param other.name \code{Character scalar}. Used as the label for the
#'   summary of non-prevalent taxa. (default: \code{"Other"})
#' 
#' @param other_label Deprecated. use \code{other.name} instead.
#'
#' @details
#' \code{agglomerateByPrevalence} sums up the values of assays at the taxonomic
#' level specified by \code{rank} (by default the highest taxonomic level
#' available) and selects the summed results that exceed the given population
#' prevalence at the given detection level. The other summed values (below the
#' threshold) are agglomerated in an additional row taking the name indicated by
#' \code{other.name} (by default "Other").
#'
#' @return
#' \code{agglomerateByPrevalence} returns a taxonomically-agglomerated object
#' of the same class as x and based on prevalent taxonomic results.
#'
#' @examples
#' ## Data can be aggregated based on prevalent taxonomic results
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' tse <- transformAssay(tse, method = "relabundance")
#' tse <- agglomerateByPrevalence(
#'     tse,
#'     rank = "Phylum",
#'     assay.type = "relabundance",
#'     detection = 1/100,
#'     prevalence = 50/100)
#'
#' tse
#'
#' # Here data is aggregated at the taxonomic level "Phylum". The five phyla
#' # that exceed the population prevalence threshold of 50/100 represent the
#' # five first rows of the assay in the aggregated data. The sixth and last row
#' # named by default "Other" takes the summed up values of all the other phyla
#' # that are below the prevalence threshold.
#'
#' assay(tse)[,1:5]
#'
#' @export
NULL

#' @rdname agglomerateByPrevalence
#' @export
setMethod("agglomerateByPrevalence", signature = c(x = "SummarizedExperiment"),
    function(x, rank = NULL, other.name = other_label, other_label = "Other",
        ...){
        # input check
        if(!.is_a_string(other.name)){
            stop("'other.name' must be a single character value.",
                call. = FALSE)
        }
        #
        # Check assays that they can be merged safely
        temp <- mapply(.check_assays_for_merge, assayNames(x), assays(x))
        #
        x <- .merge_features(x, rank, check.assays = FALSE, ...)
        pr <- getPrevalent(x, rank = NULL, ...)
        f <- rownames(x) %in% pr
        if(any(!f)){
            other_x <- agglomerateByVariable(x[!f,], by = "rows",
                                            factor(rep(1L,sum(!f))),
                                            check_assays = FALSE)
            rowData(other_x)[,colnames(rowData(other_x))] <- NA
            # set the other label
            rownames(other_x) <- other.name
            if(!is.null(rank)){
                rowData(other_x)[,rank] <- other.name
            }
            x <- rbind(x[f,], other_x)
        }
        x
    }
)

#' @rdname agglomerateByPrevalence
#' @export
setMethod("agglomerateByPrevalence", 
    signature = c(x = "TreeSummarizedExperiment"),
    function(x, rank = NULL, other.name = other_label, other_label = "Other",
            update.tree = FALSE, ...){
        # input check
        if(!.is_a_bool(update.tree)){
            stop("'update.tree' must be TRUE or FALSE.", call. = FALSE)
        }
        # update.refseq is a hidden parameter as for all other agglomeration
        # methods from the agglomerate-methods man page.
        # Here 'list(...)[["update.refseq"]]' is used to access it.
        merge_refseq <- list(...)[["update.refseq"]]
        if( is.null(merge_refseq) ){
            merge_refseq <- FALSE
        }
        if( !.is_a_bool(merge_refseq) ){
            stop("'update.refseq' must be TRUE or FALSE.", call. = FALSE)
        }
        # Agglomerate based on prevalence with SE method
        res <- callNextMethod()
        # If user wants to agglomerate reference sequences. At this point,
        # sequences are only subsetted without finding consensus sequences.
        if( merge_refseq && !is.null(referenceSeq(x))  ){
            # If user wants to agglomerate based on rank
            x <- .merge_features(x, rank, check.assays = FALSE, ...)
            # Find groups that will be used to agglomerate the data
            f <- rownames(x)[ match(rownames(x), rownames(res)) ]
            f[ is.na(f) ] <- other.name
            # Find consensus sequences, and add them to result
            ref_seq <- referenceSeq(x)
            ref_seq <- .merge_refseq_list(ref_seq, f, rownames(res), ...)
            referenceSeq(res) <- ref_seq
        }
        # Update tree if user has specified to do so
        if( update.tree ){
            res <- .agglomerate_trees(res, 1)
        }
        return(res)
    }
)

# Get abundance. Determines if relative abundance is calculated or not.
.to_rel_abund <- function(
        mat, as.relative = as_relative, as_relative = FALSE, ...) {
    # input check
    if( !.is_a_bool(as.relative) ){
        stop("'as.relative' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    if( as.relative ){
        mat <- .calc_rel_abund(mat)
    }
    return(mat)
}
