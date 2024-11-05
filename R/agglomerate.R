#' Agglomerate or merge data using taxonomic information
#'
#' Agglomeration functions can be used to sum-up data based on specific criteria
#' such as taxonomic ranks, variables or prevalence.
#'
#' \code{agglomerateByRank} can be used to sum up data based on associations
#' with certain taxonomic ranks, as defined in \code{rowData}. Only available
#' \code{\link{taxonomyRanks}} can be used.
#'
#' \code{agglomerateByVariable} merges data on rows or columns of a
#' \code{SummarizedExperiment} as defined by a \code{factor} alongside the
#' chosen dimension. This function allows agglomeration of data based on other
#' variables than taxonomy ranks.
#' Metadata from the \code{rowData} or \code{colData} are
#' retained as defined by \code{archetype}.
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} are
#' agglomerated, i.e. summed up. If the assay contains values other than counts
#' or absolute values, this can lead to meaningless values being produced.
#'  
#' @inheritParams getPrevalence
#'
#' @param empty.fields \code{Character vector}. Defines which values should be
#'   regarded as empty. (Default: \code{c(NA, "", " ", "\t")}). They will be
#'   removed if \code{na.rm = TRUE} before agglomeration.
#'
#' @param agglomerateTree Deprecated. Use \code{update.tree} instead.
#' 
#' @param agglomerate.tree Deprecated. Use \code{update.tree} instead.
#'
#' @param ... arguments passed to \code{agglomerateByRank} function for
#'   \code{SummarizedExperiment} objects,
#'   to \code{\link[=agglomerate-methods]{agglomerateByVariable}} and
#'   \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}},
#'   to \code{getPrevalence} and \code{getPrevalentTaxa} and used in
#'   \code{agglomeratebyPrevalence}
#'   \itemize{
#'        \item \code{empty.ranks.rm}: \code{Logical scalar}. Determines
#'        whether to remove those columns of rowData that include only NAs after
#'        agglomeration. (Default: \code{FALSE})
#'        \item \code{make.unique}: \code{Logical scalar}. Determines
#'        whether to make rownames unique. (Default: \code{TRUE})
#'        \item \code{detection}: The threshold value for determining presence
#'        or absence. A value in \code{x} must exceed this threshold to be
#'        considered present.
#'        \item \code{assay.type}: \code{Character scalar}. Specifies the assay
#'        used to
#'        calculate prevalence. (Default: \code{"counts"})
#'        \item \code{prevalence}: Prevalence threshold (in 0 to 1). The
#'        required prevalence is strictly greater by default. To include the
#'        limit, set \code{include.lowest} to \code{TRUE}.
#'        \item \code{update.refseq}: \code{Logical scalar}. Should a
#'        consensus sequence be calculated? If set to \code{FALSE}, the result
#'        from \code{archetype} is returned; If set to \code{TRUE} the result
#'        from
#'        \code{\link[DECIPHER:ConsensusSequence]{DECIPHER::ConsensusSequence}}
#'        is returned. (Default: \code{FALSE})
#'        \item \code{archetype} Of each level of \code{group}, which element
#'        should
#'        be regarded as the archetype and metadata in the columns or rows kept,
#'        while merging? This can be single integer value or an integer vector
#'        of the same length as \code{levels(group)}. (Default:
#'        \code{1L}, which means the first element encountered per
#'        factor level will be kept)
#'    }
#'
#' @param altexp \code{Character scalar} or \code{integer scalar}. 
#'   Specifies an alternative experiment containing the input data.
#'
#' @param altexp.rm \code{Logical scalar}. Should alternative
#'   experiments be removed prior to agglomeration? This prevents too many
#'   nested alternative experiments by default. (Default:
#'   \code{TRUE})
#' 
#' @param strip_altexp Deprecated. Use \code{altexp.rm} instead.
#'
#' @param by \code{Character scalar}. Determines if data is merged
#'   row-wise / for features ('rows') or column-wise / for samples ('cols').
#'   Must be \code{'rows'} or \code{'cols'}.
#'
#' @param group \code{Character scalar}, \code{character vector} or
#' \code{factor vector}. A column name from \code{rowData(x)} or
#' \code{colData(x)} or alternatively a vector specifying how the merging is
#' performed. If vector, the value must be the same length as
#' \code{nrow(x)/ncol(x)}. Rows/Cols corresponding to the same level will be
#' merged. If \code{length(levels(group)) == nrow(x)/ncol(x)}, \code{x} will be
#' returned unchanged.
#'  
#' @param f Deprecated. Use \code{group} instead.
#'
#' @param update.tree \code{Logical scalar}. Should
#'   \code{rowTree()} also be merged? (Default: \code{FALSE})
#' 
#' @param mergeTree Deprecated. Use \code{update.tree} instead.
#'
#' @details
#' 
#' Agglomeration sums up the values of assays at the specified taxonomic level.
#' With
#' certain assays, e.g. those that include binary or negative values, this
#' summing
#' can produce meaningless values. In those cases, consider performing
#' agglomeration
#' first, and then applying the transformation afterwards.
#'
#' \code{agglomerateByVariable} works similarly to
#' \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}}.
#' However, additional support for \code{TreeSummarizedExperiment} was added and
#' science field agnostic names were used. In addition the \code{archetype}
#' argument lets the user select how to preserve row or column data.
#'
#' For merge data of assays the function from \code{scuttle} are used.
#'
#' @return
#' \code{agglomerateByRank} returns a taxonomically-agglomerated,
#' optionally-pruned object of the same class as \code{x}.
#' \code{agglomerateByVariable} returns an object of the same class as \code{x}
#' with the specified entries merged into one entry in all relevant components.
#' \code{agglomerateByRank} returns a taxonomically-agglomerated,
#' optionally-pruned object of the same class as \code{x}.
#'
#' @name agglomerate-methods
#'
#' @seealso
#' \code{\link[=splitOn]{splitOn}}
#' \code{\link[=unsplitOn]{unsplitOn}}
#' \code{\link[=agglomerate-methods]{agglomerateByVariable}},
#' \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}},
#' \code{\link[=agglomerate-methods]{agglomerateByRank}},
#' \code{\link[SingleCellExperiment:altExps]{altExps}},
#' \code{\link[SingleCellExperiment:splitAltExps]{splitAltExps}}
#'
#' @examples
#'
#' ### Agglomerate data based on taxonomic information
#'
#' data(GlobalPatterns)
#' # print the available taxonomic ranks
#' colnames(rowData(GlobalPatterns))
#' taxonomyRanks(GlobalPatterns)
#'
#' # agglomerate at the Family taxonomic rank
#' x1 <- agglomerateByRank(GlobalPatterns, rank="Family")
#' ## How many taxa before/after agglomeration?
#' nrow(GlobalPatterns)
#' nrow(x1)
#'
#' # agglomerate the tree as well
#' x2 <- agglomerateByRank(GlobalPatterns, rank="Family",
#'                        update.tree = TRUE)
#' nrow(x2) # same number of rows, but
#' rowTree(x1) # ... different
#' rowTree(x2) # ... tree
#'
#' # If assay contains binary or negative values, summing might lead to
#' # meaningless values, and you will get a warning. In these cases, you might
#' # want to do agglomeration again at chosen taxonomic level.
#' tse <- transformAssay(GlobalPatterns, method = "pa")
#' tse <- agglomerateByRank(tse, rank = "Genus")
#' tse <- transformAssay(tse, method = "pa")
#'
#' # removing empty labels by setting na.rm = TRUE
#' sum(is.na(rowData(GlobalPatterns)$Family))
#' x3 <- agglomerateByRank(GlobalPatterns, rank="Family", na.rm = TRUE)
#' nrow(x3) # different from x2
#'
#' # Because all the rownames are from the same rank, rownames do not include
#' # prefixes, in this case "Family:".
#' print(rownames(x3[1:3,]))
#'
#' # To add them, use getTaxonomyLabels function.
#' rownames(x3) <- getTaxonomyLabels(x3, with.rank = TRUE)
#' print(rownames(x3[1:3,]))
#'
#' # use 'empty.ranks.rm' to remove columns that include only NAs
#' x4 <- agglomerateByRank(GlobalPatterns, rank="Phylum",
#'                         empty.ranks.rm = TRUE)
#' head(rowData(x4))
#'
#' # If the assay contains NAs, you might want to consider replacing them,
#' # since summing-up NAs lead to NA
#' x5 <- GlobalPatterns
#' # Replace first value with NA
#' assay(x5)[1,1] <- NA
#' x6 <- agglomerateByRank(x5, "Kingdom")
#' head( assay(x6) )
#' # Replace NAs with 0. This is justified when we are summing-up counts.
#' assay(x5)[ is.na(assay(x5)) ] <- 0
#' x6 <- agglomerateByRank(x5, "Kingdom")
#' head( assay(x6) )
#'
#' ## Look at enterotype dataset...
#' data(enterotype)
#' ## Print the available taxonomic ranks. Shows only 1 available rank,
#' ## not useful for agglomerateByRank
#' taxonomyRanks(enterotype)
#'
#' ### Merge TreeSummarizedExperiments on rows and columns
#'
#' data(esophagus)
#' esophagus
#' plot(rowTree(esophagus))
#' # get a factor for merging
#' f <- factor(regmatches(rownames(esophagus),
#'                        regexpr("^[0-9]*_[0-9]*",rownames(esophagus))))
#' merged <- agglomerateByVariable(esophagus, by = "rows", f,
#'                                 update.tree = TRUE)
#' plot(rowTree(merged))
#' #
#' data(GlobalPatterns)
#' GlobalPatterns
#' merged <- agglomerateByVariable(GlobalPatterns, by = "cols",
#'                                 colData(GlobalPatterns)$SampleType)
#' merged
NULL

#' @rdname agglomerate-methods
#' @export
setGeneric("agglomerateByRank",
            signature = "x",
            function(x, ...)
                standardGeneric("agglomerateByRank"))

#' @rdname agglomerate-methods
#' @aliases agglomerateByVariable
#' @export
setGeneric("agglomerateByVariable",
            signature = "x",
            function(x, ...)
                standardGeneric("agglomerateByVariable"))

#' @rdname agglomerate-methods
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#'
#' @export
setMethod("agglomerateByRank", signature = c(x = "SummarizedExperiment"),
    function(x, rank = taxonomyRanks(x)[1], na.rm = TRUE,
        empty.fields = c(NA, "", " ", "\t", "-", "_"), ...){
        # input check
        if(nrow(x) == 0L){
            stop("No data available in `x` ('x' has nrow(x) == 0L.)",
                call. = FALSE)
        }
        if(!.is_non_empty_string(rank)){
            stop("'rank' must be a non-empty single character value",
                call. = FALSE)
        }
        if(!.is_a_bool(na.rm)){
            stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
        }
        if(ncol(rowData(x)) == 0L){
            stop("taxonomyData needs to be populated.", call. = FALSE)
        }
        .check_taxonomic_rank(rank, x)
        .check_for_taxonomic_data_order(x)
        #

        # Make a vector from the taxonomic data.
        col <- which( taxonomyRanks(x) %in% rank )
        tax_cols <- .get_tax_cols_from_se(x)

        # if na.rm is TRUE, remove the empty, white-space, NA values from
        # tree will be pruned later, if update.tree = TRUE
        if( na.rm ){
            x <- .remove_with_empty_taxonomic_info(x, tax_cols[col],
                                                    empty.fields)
        }
        # If rank is the only rank that is available and this data is unique,
        # then the data is already 'aggregated' and no further operations
        # are needed.
        if (length(taxonomyRanks(x)) == 1L &&
            !anyDuplicated(rowData(x)[,taxonomyRanks(x)])) {
            return(x)
        }

        # get groups of taxonomy entries
        tax_factors <- .get_tax_groups(x, col = col, ...)
        # Convert to factors. Use na.rm so that NA values are not preserved.
        # i.e. they are not converted into character values.
        # NA values are handled earlier in this function.
        tax_factors <- .norm_f(nrow(x), tax_factors, na.rm = TRUE)

        # merge taxa
        x <- agglomerateByVariable(
            x, by = "rows", group = tax_factors, na.rm = TRUE, ...)

        # "Empty" the values to the right of the rank, using NA_character_.
        if( col < length(taxonomyRanks(x)) ){
            badcolumns <- tax_cols[seq_along(tax_cols) > col]
            if(length(badcolumns) > 0L){
                row_data <- rowData(x)
                row_data[, badcolumns] <- NA_character_
                rowData(x) <- row_data
            }
        }
        # adjust rownames
        rownames(x) <- getTaxonomyLabels(x, empty.fields, ...,
                                        with.rank = FALSE,
                                        resolve.loops = FALSE)
        # Remove those columns from rowData that include only NAs
        x <- .remove_NA_cols_from_rowdata(x, ...)
        x <- .add_values_to_metadata(x, "agglomerated_by_rank", rank)

        # Order the data in alphabetical order
        x <- x[ order(rownames(x)), ]
    }
)

#' @rdname agglomerate-methods
#' @aliases agglomerateByVariable
#' @export
setMethod("agglomerateByVariable", signature = c(x = "SummarizedExperiment"),
            function(x, by, group = f, f, ...){
                by <- .check_MARGIN(by)
                FUN <- switch(by, .merge_rows, .merge_cols)
                x <- FUN(x, group, ...)
                return(x)
            }
)

#' @rdname agglomerate-methods
#' @aliases agglomerateByVariable
#' @export
setMethod("agglomerateByVariable",
            signature = c(x = "TreeSummarizedExperiment"),
            function(x, by, group = f, f, update.tree = mergeTree,
                    mergeTree = FALSE, ...){
                # Check by
                by <- .check_MARGIN(by)
                # Get function based on by
                FUN <- switch(by, .merge_rows_TSE, .merge_cols_TSE)
                # Agglomerate
                x <- FUN(x, group, update.tree = update.tree, ...)
                return(x)
            }
)

#' @rdname agglomerate-methods
#' @importFrom SingleCellExperiment altExp altExp<- altExps<-
#' @export
setMethod("agglomerateByRank", signature = c(x = "SingleCellExperiment"),
    function(x, ..., altexp = NULL, altexp.rm = strip_altexp,
            strip_altexp = TRUE){
        # input check
        if(!.is_a_bool(altexp.rm)){
            stop("'altexp.rm' mus be TRUE or FALSE.", call. = FALSE)
        }
        #
        if (!is.null(altexp)) {
            x <- altExp(x, altexp)
        }
        if(altexp.rm && is(x, "SingleCellExperiment")){
            altExps(x) <- NULL
        }
        callNextMethod(x, ...)
    }
)

#' @rdname agglomerate-methods
#' @export
setMethod(
    "agglomerateByRank", signature = c(x = "TreeSummarizedExperiment"),
    function(
        x, ..., update.tree = agglomerateTree,
        agglomerate.tree = agglomerateTree, agglomerateTree = FALSE){
                # input check
                if(!.is_a_bool(update.tree)){
                    stop("'update.tree' must be TRUE or FALSE.",
                        call. = FALSE)
                }
                # If there are multipe rowTrees, it might be that multiple
                # trees are preserved after agglomeration even though the
                # dataset could be presented with one tree.
                # --> order the data so that the taxa are searched from one tree
                # first.
                if( length(rowTreeNames(x)) > 1 ){
                    x <- .order_based_on_trees(x)
                }
                # Agglomerate data
                x <- callNextMethod(x, update.tree = update.tree, ...)
                return(x)
            }
)

################################ HELP FUNCTIONS ################################

.remove_with_empty_taxonomic_info <-
    function(x, column, empty.fields = c(NA,""," ","\t","-","_"))
        {
        tax <- as.character(rowData(x)[,column])
        f <- !(tax %in% empty.fields)
        if(any(!f)){
            x <- x[f, , drop=FALSE]
        }
        x
    }

# This function removes empty rank columns from rowdata. (Those that include
# only NA values)
.remove_NA_cols_from_rowdata <- function(x, empty.ranks.rm = remove_empty_ranks,
    remove_empty_ranks = FALSE, ...){
    # Check empty.ranks.rm
    if( !.is_a_bool(empty.ranks.rm) ){
        stop("'empty.ranks.rm' must be a boolean value.",
            call. = FALSE)
    }
    # If user wants to remove those columns
    if( empty.ranks.rm ){
        # Get columns that include taxonomy information
        rank_cols <- taxonomyRanks(x)
        # Get rowData with only taxonomy
        rd <- rowData(x)[ , rank_cols, drop = FALSE]
        # Remove taxonomy from rowData
        rowData(x) <- rowData(x)[
            , !colnames(rowData(x)) %in% rank_cols, drop = FALSE]
        # Subset data so that it includes only rank columns that include data
        non_empty_ranks <- apply(rd, 2, function(x) !all(is.na(x)))
        rd <- rd[ , non_empty_ranks, drop = FALSE]
        # Adding taxonomy back to SE
        rowData(x) <- cbind(rowData(x), rd)
    }
    return(x)
}

# Order the data so that taxa from tree1 comes first, then taxa
# from tree2...
.order_based_on_trees <- function(x){
    # Get rowlinks and unique trees
    links <- DataFrame(rowLinks(x))
    uniq_trees <- sort(unique(links$whichTree))
    # Get row index to the data
    links$row_i <- seq_len(nrow(x))
    # Calculate, how many rows each tree has, and add it to data
    freq <- as.data.frame(table(links$whichTree))
    links <- merge(links, freq, all.x = TRUE, all.y = FALSE,
                    by.x = "whichTree", by.y = "Var1")
    # Factorize the names of trees
    links$whichTree <- factor(links$whichTree, levels = uniq_trees)
    # Order the data back to its original order based on row indices
    links <- links[order(links$row_i), ]
    # Get the order based on size of tree and name
    order <- order(links$whichTree)
    # Order the data
    x <- x[order, ]
    return(x)
}

# Agglomerate all rowTrees found in TreeSE object. Get tips that represent
# rows and remove all others.
#' @importFrom TreeSummarizedExperiment subsetByLeaf
.agglomerate_trees <- function(x, by = 1, ...){
    # Get right functions based on direction
    tree_names_FUN <- switch(
        by, "1" = rowTreeNames, "2" = colTreeNames, stop("."))
    links_FUN <- switch(by, "1" = rowLinks, "2" = colLinks, stop("."))
    tree_FUN <- switch(by, "1" = rowTree, "2" = colTree, stop("."))
    # Get right argument names for subsetByLeaf call
    args_names <- switch(
        by,
        "1" = c("x", "rowLeaf", "whichRowTree", "updateTree"),
        "2" = c("x", "colLeaf", "whichColTree", "updateTree"),
        stop("."))
    # Get names of trees and links between trees and rows
    tree_names <- tree_names_FUN(x)
    row_links <- links_FUN(x)
    # Loop through tree names
    for( name in tree_names ){
        # Get the tree that is being agglomerated
        tree <- tree_FUN(x, name)
        # Get row links that corresponds this specific tree
        links_temp <- row_links[ row_links[["whichTree"]] == name, ]
        # If the tree represents the data, agglomerate it
        if( nrow(links_temp) > 0 ){
            # Get names of nodes that are preserved
            links_temp <- links_temp[["nodeLab"]]
            # Agglomerate the tree
            args <- list(x, links_temp, name, TRUE)
            names(args) <- args_names
            x <- do.call(subsetByLeaf, args)
        }
    }
    return(x)
}

# This function trims tips until all tips can be found from provided set of
# nodes
#' @importFrom ape drop.tip has.singles collapse.singles
.prune_tree <- function(tree, nodes, collapse.singles = TRUE, ...){
    # Check collapse.singles
    if( !.is_a_bool(collapse.singles) ){
        stop("'collapse.singles' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Get those tips that can not be found from provided nodes
    remove_tips <- .get_tips_to_drop(tree, nodes)
    # As long as there are tips to be dropped, run the loop
    while( length(remove_tips) > 0 ){
        # Drop tips that cannot be found. Drop only one layer at the time. Some
        # dataset might have taxa that are not in tip layer but they are in
        # higher rank. If we delete more than one layer at the time, we might
        # loose the node for those taxa. --> The result of pruning is a tree
        # whose all tips can be found provided nodes i.e., rows of TreeSE. Some
        # taxa might be higher rank meaning that all rows might not be in tips
        # even after pruning; these rows have still child-nodes that represent
        # other rows.
        # Suppress warning: drop all tips of the tree: returning NULL
        tree <- tryCatch({
            drop.tip(
                tree, remove_tips,
                trim.internal = FALSE,
                collapse.singles = FALSE)
        }, warning = function(w) {
            # Do nothing on warning
        }, error = function(e) {
            # Try to prune by also pruning internal nodes. Sometimes that is the
            # case; we need to trim also internal nodes in order to prune
            # leaf.
            drop.tip(
                tree, remove_tips,
                trim.internal = TRUE,
                collapse.singles = FALSE)
        })
        # If all tips were dropped, the result is NULL --> stop loop
        if( is.null(tree) ){
            warning("Pruning resulted to empty tree.", call. = FALSE)
            break
        }
        # Again, get those tips of updated tree that cannot be found from
        # provided nodes
        remove_tips <- .get_tips_to_drop(tree, nodes)
    }
    # Simplify the tree structure. Remove nodes that have only single
    # descendant.
    if( !is.null(tree) && length(tree$tip.label) > 1 && has.singles(tree) &&
            collapse.singles ){
        tree <- collapse.singles(tree)
    }
    return(tree)
}

# This function gets tree and nodes as input. As output, it gives set of tips
# that are not in the set of nodes provided as input.
.get_tips_to_drop <- function(tree, nodes){
    # Get those tips cannot be found from node set
    cannot_be_found <- !tree$tip.label %in% nodes
    # Get those tips that are duplicated. Single node should match with only
    # one row.
    dupl <- duplicated(tree$tip.label)
    # Get indices of those tips that are going to be removed
    tips <- which( cannot_be_found | dupl )
    return(tips)
}
