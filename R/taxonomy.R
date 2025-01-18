#' Functions for accessing taxonomic data stored in \code{rowData}.
#'
#' These function work on data present in \code{rowData} and define a way to
#' represent taxonomic data alongside the features of a
#' \code{SummarizedExperiment}.
#'
#' \code{taxonomyRanks} returns, which columns of \code{rowData(x)} are regarded
#' as columns containing taxonomic information.
#'
#' \code{taxonomyRankEmpty} checks, if a selected rank is empty of information.
#'
#' \code{checkTaxonomy} checks, if taxonomy information is valid and whether
#'   it contains any problems. This is a soft test, which reports some
#'   diagnostic and might mature into a data validator used upon object
#'   creation.
#'
#' \code{getTaxonomyLabels} generates a character vector per row consisting of
#'   the lowest taxonomic information possible. If data from different levels,
#'   is to be mixed, the taxonomic level is prepended by default.
#'
#' \code{IdTaxaToDataFrame} extracts taxonomic results from results of
#'   \code{\link[DECIPHER:IdTaxa]{IdTaxa}}.
#'
#' \code{mapTaxonomy} maps the given features (taxonomic groups; \code{taxa})
#'   to the specified taxonomic level (\code{to} argument) in \code{rowData}
#'   of the \code{SummarizedExperiment} data object
#'   (i.e. \code{rowData(x)[,taxonomyRanks(x)]}). If the argument \code{to} is
#'   not provided, then all matching taxonomy rows in \code{rowData} will be
#'   returned. This function allows handy conversions between different
#    taxonomic levels.
#'
#' @inheritParams agglomerate-methods
#'
#' @param with.rank \code{Logical scalar}. Should the level be add as a
#'   suffix? For example: "Phylum:Crenarchaeota". (Default: \code{FALSE})
#'
#' @param with_rank Deprecated. Use \code{with.rank} instead.
#'
#' @param make.unique \code{Logical scalar}. Should the labels be made
#'   unique, if there are any duplicates? (Default: \code{TRUE})
#'
#' @param make_unique Deprecated. Use \code{make.unique} instead.
#'
#' @param resolve.loops \code{Logical scalar}. Should \code{resolveLoops}
#'   be applied to the taxonomic data? Please note that has only an effect,
#'   if the data is unique. (Default: \code{TRUE})
#'
#' @param resolve_loops Deprecated. Use \code{resolve.loops} instead.
#'
#' @param taxa \code{Character vector}. Used for subsetting the
#'   taxonomic information. If no information is found,\code{NULL} is returned
#'   for the individual element. (Default: \code{NULL})
#'
#' @param from
#' \itemize{
#'   \item For \code{mapTaxonomy}: \code{character scalar}. A value which
#'     must be a valid taxonomic rank. (Default: \code{NULL})
#'   \item otherwise a \code{Taxa} object as returned by
#'     \code{\link[DECIPHER:IdTaxa]{IdTaxa}}
#' }
#'
#' @param to \code{Character Scalar}. Must be a valid
#'   taxonomic rank. (Default: \code{NULL})
#'
#' @param use.grepl \code{Logical}. Should pattern matching via
#'   \code{grepl} be used? Otherwise literal matching is used.
#'   (Default: \code{FALSE})
#'
#' @param use_grepl Deprecated. Use \code{use.grepl} instead.
#'
#' @param ... additional arguments
#' \itemize{
#'   \item \code{lowest.rank}: A lowest taxonomy level to be considered in
#'   \code{getTaxonomyLabels}. Ranks lower than this will be collapsed into rank
#'   specified by \code{lowest.rank}. For example, if genus level is specified,
#'   species will be collapsed into genus. If \code{NULL}, the data is not
#'   collapsed. (Default: \code{NULL})
#' }
#'
#' @param ranks \code{Character vector}. A vector of ranks to be set.
#' @details
#' Taxonomic information from the \code{IdTaxa} function of \code{DECIPHER}
#' package are returned as a special class. With \code{as(taxa,"DataFrame")}
#' the information can be easily converted to a \code{DataFrame} compatible
#' with storing the taxonomic information a \code{rowData}. Please note that the
#' assigned confidence information are returned as \code{metatdata} and can
#' be accessed using \code{metadata(df)$confidence}.
#'
#' @return
#' \itemize{
#'   \item \code{taxonomyRanks}: a \code{character} vector with all the
#'     taxonomic ranks found in \code{colnames(rowData(x))}
#'   \item \code{taxonomyRankEmpty}: a \code{logical} value
#'   \item \code{mapTaxonomy}: a \code{list} per element of taxa. Each
#'     element is either a \code{DataFrame}, a \code{character} or \code{NULL}.
#'     If all \code{character} results have the length of one, a single
#'     \code{character} vector is returned.
#' }
#'
#' @name taxonomy-methods
#'
#' @seealso \code{\link[=agglomerate-methods]{agglomerateByRank}},
#' \code{\link[TreeSummarizedExperiment:toTree]{toTree}},
#' \code{\link[TreeSummarizedExperiment:resolveLoop]{resolveLoop}}
#'
#' @examples
#' data(GlobalPatterns)
#' GlobalPatterns
#' taxonomyRanks(GlobalPatterns)
#'
#' checkTaxonomy(GlobalPatterns)
#'
#' table(taxonomyRankEmpty(GlobalPatterns,"Kingdom"))
#' table(taxonomyRankEmpty(GlobalPatterns,"Species"))
#'
#' getTaxonomyLabels(GlobalPatterns[1:20,])
#' # Taxonomy labels represent the lowest taxonomy name that identifies each
#' # taxa. For instance, they can represent OTUs which does no necessarily
#' # tell much. In this case, you might want to get the labels with higher
#' # taxonomy rank
#' getTaxonomyLabels(GlobalPatterns[1:20,], lowest.rank = "Class")
#'
#' # mapTaxonomy
#' ## returns the unique taxonomic information
#' mapTaxonomy(GlobalPatterns)
#' # returns specific unique taxonomic information
#' mapTaxonomy(GlobalPatterns, taxa = "Escherichia")
#' # returns information on a single output
#' mapTaxonomy(GlobalPatterns, taxa = "Escherichia",to="Family")
#'
#' # setTaxonomyRanks
#' tse <- GlobalPatterns
#' colnames(rowData(tse))[1] <- "TAXA1"
#'
#' setTaxonomyRanks(colnames(rowData(tse)))
#' # Taxonomy ranks set to: taxa1 phylum class order family genus species
#'
#' # getTaxonomyRanks is to get/check if the taxonomic ranks is set to "TAXA1"
#' getTaxonomyRanks()
NULL

# This function returns all supported ranks and their prefixes. These ranks are
# used to detect ranks in imported data.
.taxonomy_rank_prefixes <- c(
    domain = "d",
    superkingdom = "sk",
    kingdom = "k",
    phylum = "p",
    class = "c",
    order = "o",
    family = "f",
    genus = "g",
    species = "s",
    strain = "t"
)

# Function to set taxonomy ranks prefixes (not exported)
#' @importFrom utils assignInMyNamespace
setTaxonomyRankPrefixes <- function(prefixes) {
    # Check if prefixes is a character vector with length >= 1 and it has names
    if( !(is.character(prefixes) && length(prefixes) > 0 &&
            !is.null(names(prefixes))) ){
        stop(
            "'prefixes' must be a non-empty character vector and it must have ",
            "names.", call. = FALSE)
    }
    # Replace default value of mia:::.taxonomy_rank_prefixes
    assignInMyNamespace(".taxonomy_rank_prefixes", prefixes)
}

# Function to get taxonomy ranks prefixes (not exported)
getTaxonomyRankPrefixes <- function() {
    return(.taxonomy_rank_prefixes)
}

#' @format a \code{character} vector of length containing all the taxonomy ranks
#'   recognized. In functions this is used as case insensitive.
TAXONOMY_RANKS <- names(.taxonomy_rank_prefixes)

#' @rdname taxonomy-methods
#'
#' @importFrom SummarizedExperiment rowData
#'
#' @export
setMethod("taxonomyRanks", signature = c(x = "SummarizedExperiment"),
    function(x){
        ranks <- colnames(rowData(x))
        ranks[.get_tax_cols(ranks)]
    }
)

#' @rdname taxonomy-methods
#' @aliases taxonomyRankEmpty
#'
#' @importFrom SummarizedExperiment rowData
#'
#' @export
setMethod("taxonomyRankEmpty", signature = c(x = "SummarizedExperiment"),
    function(x, rank = taxonomyRanks(x)[1],
            empty.fields = c(NA, "", " ", "\t", "-", "_")){
    # input check
    if(ncol(rowData(x)) == 0L){
        stop("rowData needs to be populated.", call. = FALSE)
    }
    if(!.is_non_empty_string(rank)){
        stop("'rank' must be an non empty single character value.",
            call. = FALSE)
    }
    if(!is.character(empty.fields) || length(empty.fields) == 0L){
        stop("'empty.fields' must be a character vector with one or ",
            "more value", call. = FALSE)
    }
    .check_taxonomic_rank(rank, x)
    .check_for_taxonomic_data_order(x)
    #
    rowData(x)[,rank] %in% empty.fields
    }
)

#' @rdname taxonomy-methods
#' @aliases checkTaxonomy
#' @export
setMethod("checkTaxonomy", signature = c(x = "SummarizedExperiment"),
    function(x){
        tmp <- try(.check_for_taxonomic_data_order(x), silent = TRUE)
        ans <- !is(tmp,"try-error")
        if(!ans){
            attr(ans, "msg") <- as.character(ans)
        }
        ans
    }
)

#' @rdname taxonomy-methods
#' @importFrom utils assignInMyNamespace
#' @aliases checkTaxonomy
#' @export
# Function to set taxonomy ranks
setTaxonomyRanks <- function(ranks) {
    ranks <- tolower(ranks)
    # Check if rank is a character vector with length >= 1
    if (!is.character(ranks) || length(ranks) < 1
        || any(ranks == "" | ranks == " " | ranks == "\t" | ranks == "-" |
            ranks == "_")
        || any(grepl("\\s{2,}", ranks))) {
        stop("Input 'rank' should be a character vector with non-empty strings,
            no spaces, tabs, hyphens, underscores, and non-continuous spaces.",
            call. = FALSE)
    }
    #Replace default value of mia::TAXONOMY_RANKS
    assignInMyNamespace("TAXONOMY_RANKS", ranks)
}

#' @rdname taxonomy-methods
#' @export
# Function to get taxonomy ranks
getTaxonomyRanks <- function() {
    return(TAXONOMY_RANKS)
}

.check_taxonomic_rank <- function(rank, x){
    if(length(rank) != 1L){
        stop("'rank' must be a single character value.",call. = FALSE)
    }
    if( !(rank %in% taxonomyRanks(x) ) ){
        stop("'rank' must be a value from 'taxonomyRanks()'",call. = FALSE)
    }
}
.check_taxonomic_ranks <- function(ranks, x){
    if(length(ranks) == 0L){
        stop("'ranks' must contain at least one value.",call. = FALSE)
    }
    if( !all(ranks %in% taxonomyRanks(x) ) ){
        stop("'ranks' must contain values from 'taxonomyRanks()'")
    }
}

#' @importFrom SummarizedExperiment rowData
.check_for_taxonomic_data_order <- function(x){
    ranks <- colnames(rowData(x))
    f <- tolower(ranks) %in% TAXONOMY_RANKS
    if(!any(f)){
        stop("no taxonomic ranks detected in rowData(). Columns with one of ",
            "the following names can be used: '",
            paste(TAXONOMY_RANKS, collapse = "', '"), "'", call. = FALSE)
    }
    m <- match(TAXONOMY_RANKS, tolower(ranks[f]))
    m <- m[!is.na(m)]
    # check that taxonomic ranks are in order. If they are all value in check
    # should be 1 or 0
    check <- unique(c(m[-1], m[length(m)]) - m )
    if(!all(check %in% c(1L,0L))){
        stop("Taxonomic ranks are not in order. Please reorder columns, which ",
            "correspond to taxonomic ranks like this:\n'",
            paste(TAXONOMY_RANKS, collapse = "', '"), "'.",
            call. = FALSE)
    }
}

#' @rdname taxonomy-methods
#' @aliases checkTaxonomy
#' @export
setMethod("getTaxonomyLabels", signature = c(x = "SummarizedExperiment"),
    function(x, empty.fields = c(NA, "", " ", "\t", "-", "_"),
            with.rank = with_rank,
            with_rank = FALSE, make.unique = make_unique, make_unique = TRUE,
            resolve.loops = resolve_loops, resolve_loops = FALSE, ...){
        # input check
        if(nrow(x) == 0L){
            stop("No data available in `x` ('x' has nrow(x) == 0L.)",
                call. = FALSE)
        }
        if(ncol(rowData(x)) == 0L){
            stop("rowData needs to be populated.", call. = FALSE)
        }
        .check_for_taxonomic_data_order(x)
        if(!is.character(empty.fields) || length(empty.fields) == 0L){
            stop("'empty.fields' must be a character vector with one or ",
                "more values.", call. = FALSE)
        }
        if(!.is_a_bool(with.rank)){
            stop("'with.rank' must be TRUE or FALSE.", call. = FALSE)
        }
        if(!.is_a_bool(make.unique)){
            stop("'make.unique' must be TRUE or FALSE.", call. = FALSE)
        }
        if(!.is_a_bool(resolve.loops)){
            stop("'resolve.loops' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        # Collapse taxonomy ranks if user has specified so
        x <- .collapse_lowest_taxonomy_ranks(x, ...)

        dup <- duplicated(rowData(x)[,taxonomyRanks(x), drop = FALSE])
        if(any(dup)){
            td <- apply(
                rowData(x)[,taxonomyRanks(x), drop = FALSE],
                1L,
                paste, collapse = "___")
            td_non_dup <- td[!dup]
            m <- match(td, td_non_dup)
        }
        ans <- .get_taxonomic_label(x[!dup,],
                                    empty.fields = empty.fields,
                                    with.rank = with.rank,
                                    resolve.loops = resolve.loops)
        if(any(dup)){
            ans <- ans[m]
        }
        # last resort - this happens, if annotation data contains ambiguous data
        # sometimes labeled as "circles"
        if(make.unique && anyDuplicated(ans)){
            dup <- which(ans %in% ans[which(duplicated(ans))])
            ans[dup] <- make.unique(ans[dup], sep = "_")
        }
        ans
    }
)

# This function is for collapsing the lowest taxonomy ranks into single value.
# For instance, if user specifies genus rank, genus and species are collapsed
# into one, and species rank is removed from the taxonomy table. If family is
# specified, along with these two, also family is incorporated into this value.
.collapse_lowest_taxonomy_ranks <- function(
        x, lowest.rank = NULL, empty.fields = c(NA, "", " ", "\t", "-", "_"),
        sep = "_", ...){
    if( !.is_a_string(sep) ){
        stop("'sep' must be a single character value.", call. = FALSE)
    }
    # By default, we keep the taxonomy table untouched.
    if( !is.null(lowest.rank) ){
        .check_taxonomic_rank(lowest.rank, x)
        # Get available taxonomy ranks
        available_ranks <- taxonomyRanks(x)
        # Get indices of ranks that we are going to collapse into one
        mod_ranks <- seq(
            which(available_ranks == lowest.rank), length(available_ranks))
        # For each row, collapse ranks into one
        new_rank <- apply(rowData(x)[, mod_ranks, drop = FALSE], 1, function(x){
            # Check if empty, and replace with NA if it is
            x[ x %in% empty.fields ] <- NA
            # Collapse values
            x <- paste(na.omit(x), collapse = sep)
            return(x)
        })
        # Remove the collapsed ranks from the original table
        rowData(x) <- rowData(x)[, -mod_ranks, drop = FALSE]
        # Add the collapsed values to the taxonomy table
        new_rank <- unname(new_rank)
        rowData(x)[[lowest.rank]] <- new_rank
    }
    return(x)
}

#' @importFrom IRanges CharacterList LogicalList
.get_tax_ranks_selected <- function(x,rd, tax_cols, empty.fields){
    # We need DataFrame here to handle cases with a single entry in tax_cols
    charlist <- CharacterList(t(rd[,tax_cols, drop=FALSE]))
    tax_ranks_non_empty <- !is.na(charlist) &
        !LogicalList(lapply(charlist,"%in%",empty.fields))

    tax_ranks_non_empty <- t(as(tax_ranks_non_empty,"matrix"))
    tax_ranks_selected <- apply(tax_ranks_non_empty,1L,which)
    # Check if every row got at least some rank information from taxonomy table
    # i.e. the info was not empty.
    if( any(lengths(tax_ranks_selected) == 0L) || length(
        tax_ranks_selected) == 0L){
        if(!anyDuplicated(rownames(x))){
            return(NULL)
        }
        stop("Only empty taxonomic information detected. Some rows contain ",
            "only entries selected by 'empty.fields'. Cannot generated ",
            "labels. Try option na.rm = TRUE in the function call.",
            call. = FALSE)
    }
    #
    if(is.matrix(tax_ranks_selected)){
        tax_ranks_selected <- apply(tax_ranks_selected,2L,max)
    } else if(is.list(tax_ranks_selected)) {
        tax_ranks_selected <- lapply(tax_ranks_selected,max)
        tax_ranks_selected <- unlist(tax_ranks_selected)
    } else if(is.vector(tax_ranks_selected)){
        tax_ranks_selected <- max(tax_ranks_selected)
    } else {
        stop(".")
    }
    tax_ranks_selected
}

.add_taxonomic_type <- function(rd, ans, tax_cols_selected){
    sep <- rep(":", length(ans))
    tax_cols_selected <- unlist(tax_cols_selected)
    # sep[tax_cols_selected != max(tax_cols_selected)] <- "::"
    types <- colnames(rd)[tax_cols_selected]
    ans <- paste0(types, sep, ans)
    ans
}

.get_taxonomic_label <- function(x,
                                empty.fields = c(NA, "", " ", "\t", "-", "_"),
                                with.rank = FALSE,
                                resolve.loops = FALSE){
    rd <- rowData(x)
    tax_cols <- .get_tax_cols_from_se(x)
    tax_ranks_selected <- .get_tax_ranks_selected(x, rd, tax_cols, empty.fields)
    if(is.null(tax_ranks_selected)){
        return(rownames(x))
    }
    tax_cols_selected <- tax_cols[tax_ranks_selected]
    # resolve loops
    if(resolve.loops){
        td <- as.data.frame(rd[,tax_cols])
        td <- suppressWarnings(resolveLoop(td))
        rd[,tax_cols] <- as(td,"DataFrame")
        rm(td)
    }
    #
    all_same_rank <- length(unique(tax_cols_selected)) == 1L
    ans <- mapply("[",
                    as.data.frame(t(as.data.frame(rd))),
                    tax_cols_selected,
                    SIMPLIFY = FALSE)
    ans <- unlist(ans, use.names = FALSE)
    if(with.rank || !all_same_rank){
        ans <- .add_taxonomic_type(rd, ans, tax_cols_selected)
    }
    ans
}

#' Calculate hierarchy tree
#'
#' These functions generate a hierarchy tree using taxonomic information from a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{SummarizedExperiment}}
#' object and add this hierarchy tree into the \code{rowTree}.
#'
#' @inheritParams taxonomy-methods
#'
#' @param ... optional arguments not used currently.
#'
#' @details
#'
#' \code{addHierarchyTree} calculates a hierarchy tree from the available
#'   taxonomic information and add it to \code{rowTree}.
#'
#' \code{getHierarchyTree} generates a hierarchy tree from the available
#'   taxonomic information. Internally it uses
#'   \code{\link[TreeSummarizedExperiment:toTree]{toTree}} and
#'   \code{\link[TreeSummarizedExperiment:resolveLoop]{resolveLoop}} to sanitize
#'   data if needed.
#'
#' Please note that a hierarchy tree is not an actual phylogenetic tree.
#' A phylogenetic tree represents evolutionary relationships among features.
#' On the other hand, a hierarchy tree organizes species into a hierarchical
#' structure based on their taxonomic ranks.
#'
#' @return
#' \itemize{
#'   \item \code{addHierarchyTree}: a \code{TreeSummarizedExperiment} whose
#'   \code{phylo} tree represents the hierarchy among available taxonomy
#'   information.
#'   \item \code{getHierarchyTree}: a \code{phylo} tree representing the
#'   hierarchy among available taxonomy information.
#' }
#'
#' @name hierarchy-tree
#'
#' @examples
#' # Generate a tree based on taxonomic rank hierarchy (a hierarchy tree).
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' getHierarchyTree(tse)
#'
#' # Add a hierarchy tree to a TreeSummarizedExperiment.
#' # Please note that any tree already stored in rowTree() will be overwritten.
#' tse <- addHierarchyTree(tse)
#' tse
NULL

#' @rdname hierarchy-tree
#' @aliases getHierarchyTree
#' @export
#' @importFrom ape drop.tip
setMethod("getHierarchyTree", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        # Input check
        # If there is no rowData it is not possible to create rowTree
        if( ncol(rowData(x)) == 0L ){
            stop("'x' does not have rowData. Tree cannot be created.",
                call. = FALSE)
        }
        # If there are no taxonomy ranks
        if( length(taxonomyRanks(x)) < 2L ){
            stop(
                "'x' does not contain adequate taxonomy information, and ",
                "hierarchy tree cannot be created. Check rowData and consider ",
                "using setTaxonomyRanks() if ranks differ from defaults..",
                call. = FALSE)
        }
        #
        # Get rowData as data.frame
        td <- rowData(x)[, taxonomyRanks(x), drop = FALSE]
        td <- as.data.frame(td)
        # Get information on empty nodes. It will be used later to polish the
        # created tree.
        td_NA <- .get_empty_nodes(td, ...)
        # Replace empty cells with NA (also "" can be empty value)
        for( i in seq_len(ncol(td_NA)) ){
            td[td_NA[[i]], i] <- NA
        }
        # Remove empty taxonomic levels
        td <- td[ , !vapply(td_NA, all, logical(1)), drop = FALSE]
        # Check if there is no taxonomy information left after removing empty
        # columns
        if( ncol(td) < 2L ){
            stop(
                "'x' does not contain adequate taxonomy information, and ",
                "hierarchy tree cannot be created. Check rowData and consider ",
                "using setTaxonomyRanks() if ranks differ from defaults.",
                call. = FALSE)
        }
        # Make cells unique. Add suffix, if duplicated values are found from
        # certain rank.
        td <- suppressWarnings(resolveLoop(td))
        # Build tree
        tree <- toTree(td)
        tree$tip.label <- paste0(colnames(td)[ncol(td)],":",tree$tip.label)
        # remove empty nodes
        for(i in rev(seq_len(ncol(td)))){
            if(any(td_NA[,i])){
                to_drop <- paste0(colnames(td)[i],":",td[,i][td_NA[,i]])
                tree <- drop.tip(
                    tree,
                    to_drop,
                    trim.internal = FALSE,
                    collapse.singles = FALSE)
            }
        }
        return(tree)
    }
)

#' @rdname hierarchy-tree
#' @export
setMethod("addHierarchyTree", signature = c(x = "SummarizedExperiment"),
    function(x, ...){
        #
        # Get the tree
        tree <- getHierarchyTree(x, ...)
        # Ensure that the object has rowTree slot
        x <- as(x,"TreeSummarizedExperiment")
        # Get node labs: which row represents which node in the tree?
        node_labs <- getTaxonomyLabels(
            x, with.rank = TRUE, resolve.loops = TRUE, make.unique = FALSE)
        # Add tree
        x <- changeTree(x, tree, node_labs)
        return(x)
    }
)

#' @importFrom BiocGenerics %in% grepl
.get_taxa_row_match <- function(taxa, td, from, use.grepl = FALSE){
    if(is.na(taxa)){
        r_f <- is.na(td[[from]])
    } else {
        if(use.grepl){
            r_f <- grepl(taxa, td[[from]], ignore.case = TRUE)
        } else {
            r_f <- td[[from]] %in% taxa
        }
    }
    r_f[is.na(r_f)] <- FALSE
    r_f
}

#' @importFrom BiocGenerics %in% grepl
.get_taxa_any_match <- function(taxa, td, use.grepl = FALSE){
    if(is.na(taxa)){
        r_f <- is.na(td)
    } else {
        if(use.grepl){
            r_f <- vapply(
                td, grepl, logical(nrow(td)), pattern=taxa, ignore.case = TRUE)
        } else {
            r_f <- t(as.matrix(td %in% taxa))
        }
    }
    r_f <- rowSums(r_f, na.rm = TRUE) > 0
    r_f[is.na(r_f)] <- FALSE
    r_f
}

#' @rdname taxonomy-methods
#' @importFrom BiocGenerics %in%
#' @export
setMethod("mapTaxonomy", signature = c(x = "SummarizedExperiment"),
    function(x, taxa = NULL, from = NULL, to = NULL, use.grepl = use_grepl,
            use_grepl = FALSE){
        # input check
        if(!checkTaxonomy(x)){
            stop("Non compatible taxonomic information found. ",
                "checkTaxonomy(x) must be TRUE.",
                call. = FALSE)
        }
        if(!is.null(taxa)){
            if(!is.character(taxa)){
                stop("'taxa' must be a character vector.",
                    call. = FALSE)
            }
        }
        if(!is.null(from)){
            if(!.is_a_string(from)){
                stop("'from' must be a single character value.",
                    call. = FALSE)
            }
            if(!(from %in% taxonomyRanks(x))){
                stop("'from' must be an element of taxonomyRanks(x).",
                    call. = FALSE)
            }
        }
        if(!is.null(to)){
            if(!.is_a_string(to)){
                stop("'to' must be a single character value.",
                    call. = FALSE)
            }
            if(!(to %in% taxonomyRanks(x))){
                stop("'to' must be an element of taxonomyRanks(x).",
                    call. = FALSE)
            }
        }
        if(!is.null(from) && !is.null(to)){
            if(from  == to){
                stop("'from' and 'to' must be different values.", call. = FALSE)
            }
        }
        if(!.is_a_bool(use.grepl)){
            stop("'use.grepl' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        td <- rowData(x)[,taxonomyRanks(x)]
        if(is.null(taxa)){
            return(unique(td))
        }
        #
        r_fs <- NULL
        c_f <- rep(TRUE,ncol(td))
        # Get unique taxa to search
        taxa <- unique(taxa)
        if(!is.null(from)){
            r_fs <- lapply(
                taxa, .get_taxa_row_match, td = td, from = from,
                use.grepl = use.grepl)
            names(r_fs) <- taxa
        } else {
            r_fs <- lapply(
                taxa, .get_taxa_any_match, td = td,
                use.grepl = use.grepl)
            names(r_fs) <- taxa
        }
        if(!is.null(to)) {
            c_f <- colnames(td) == to
            c_f[is.na(c_f)] <- FALSE
        }
        # assemble the result
        ans <- lapply(r_fs, .get_map_result, td = td, c_f = c_f)
        names(ans) <- names(r_fs)
        u_len <- unique(lengths(ans))
        if(length(u_len) == 1L && u_len == 1L){
            ans <- unlist(ans)
        }
        #
        ans
    }
)

.get_map_result <- function(r_f, td, c_f){
    ans <- td[r_f,c_f]
    ans <- unique(ans)
    if(is(ans,"DataFrame") && nrow(ans) == 0L){
        return(NULL)
    }
    ans
}

################################################################################
# helper functions

.get_tax_cols_logical <- function(x){
    tolower(x) %in% TAXONOMY_RANKS
}

.get_tax_cols <- function(x){
    which(.get_tax_cols_logical(x))
}

#' @importFrom SummarizedExperiment rowData
.get_tax_cols_from_se <- function(x){
    .get_tax_cols(colnames(rowData(x)))
}

#' @importFrom SummarizedExperiment rowData
.get_tax_groups <- function(x, col, ignore.taxonomy = onRankOnly,
    onRankOnly = FALSE, ...){
    # input check
    if(!.is_a_bool(ignore.taxonomy)){
        stop("'ignore.taxonomy' must be TRUE or FALSE.", call. = FALSE)
    }

    tax_cols <- .get_tax_cols_from_se(x)
    tax_col_n <- seq_along(tax_cols)
    if(length(tax_col_n) < col){
        stop(".")
    }
    if(ignore.taxonomy){
        groups <- rowData(x)[,tax_cols[tax_col_n == col],drop=TRUE]
    } else {
        groups <- rowData(x)[,tax_cols[tax_col_n <= col],drop=FALSE]
        groups <- apply(groups,1L,paste,collapse="_")
    }
    factor(groups, unique(groups))
}

################################################################################
# IDTAXA to DataFrame conversion

#' @importFrom IRanges CharacterList NumericList
#' @importFrom S4Vectors pc DataFrame
.idtaxa_to_DataFrame <- function(from){
    ranks <- CharacterList(lapply(from,"[[","rank"))
    conf <- NumericList(lapply(from,"[[","confidence"))
    taxa <- CharacterList(lapply(from,"[[","taxon"))
    # even out the lengths
    l <- lengths(ranks)
    ml <- max(l)
    diff <- ml - l
    add <- CharacterList(lapply(diff,rep,x=NA))
    ranks <- pc(ranks,add)
    conf <- pc(conf,as(add,"NumericList"))
    taxa <- pc(taxa,add)
    # convert to DataFrame
    names <- unique(unlist(ranks))
    names <- names[!is.na(names)]
    taxa <- DataFrame(as.matrix(taxa))
    colnames(taxa) <- names
    conf <- DataFrame(as.matrix(conf))
    colnames(conf) <- names
    # subset to valid taxonomic information
    f <- tolower(names) %in% TAXONOMY_RANKS
    taxa <- taxa[,f]
    conf <- conf[,f]
    # combine with confidence data
    metadata(taxa)$confidence <- conf
    #
    taxa
}

#' @rdname taxonomy-methods
#' @export
IdTaxaToDataFrame <- .idtaxa_to_DataFrame

# This function gives information on if cell in taxonomy table is empty or not.
.get_empty_nodes <- function(
        td, empty.fields = c(NA, "", " ", "\t", "-", "_"), ...){
    # Check empty.fields
    if(!is.character(empty.fields) || length(empty.fields) == 0L){
        stop(
            "'empty.fields' must be a character vector with one or ",
            "more value", call. = FALSE)
    }
    #
    # Loop over columns. For each cell, get info if the cell is empty or not.
    is_empty <- lapply(td, function(x){
        temp <- x %in% empty.fields
        return(temp)
    })
    # Convert to data.frame
    is_empty <- as.data.frame(is_empty)
    return(is_empty)
}
