#' Merge SE objects into single SE object.
#' 
#' @inheritParams rarefyAssay
#' @inheritParams getDominant
#' 
#' @param y a \code{\link{SummarizedExperiment}} object when \code{x} is a
#' \code{\link{SummarizedExperiment}} object. Disabled when \code{x} is a list.
#' 
#' @param join \code{Character scalar}. A value for selecting the joining
#' method. Must be 'full', 'inner', 'left', or 'right'. 'left' and 'right' are
#' disabled when more than two objects are being merged.
#' (Default: \code{"full"})
#' 
#' @param missing.values \code{NA}, \code{0} or \code{Character scalar}.
#' Specifies the notation of missing values. (By default: \code{NA})
#' 
#' @param missing_values Deprecated. Use \code{missing.values} instead.
#' 
#' @param collapse.cols \code{Logical scalar}. Determines whether to collapse
#' identically named samples to one. (Default: \code{FALSE})
#' 
#' @param collapse_samples Deprecated. Use \code{collapse.cols} instead.
#' 
#' @param collapse.rows \code{Logical scalar}. Selects whether to collapse
#' identically named features to one. Since all taxonomy information is
#' taken into account,
#' this concerns rownames-level (usually strain level) comparison. Often
#' OTU or ASV level is just an arbitrary number series from sequencing machine
#' meaning that the OTU information is not comparable between studies. With this
#' option, it is possible to specify whether these strains are combined if their
#' taxonomy information along with OTU number matches.
#' (Default: \code{TRUE})
#' 
#' @param collapse_features Deprecated. Use \code{collapse.rows} instead.
#'
#' @param ... optional arguments (not used).
#'
#' @return A single \code{SummarizedExperiment} object.
#'
#' @details
#' This function merges multiple \code{SummarizedExperiment} objects. It
#' combines
#' \code{rowData}, \code{assays}, and \code{colData} so that the output includes
#' each unique row and column ones. The merging is done based on
#' \code{rownames} and
#' \code{colnames}. \code{rowTree} and \code{colTree} are preserved if linkage
#' between rows/cols and the tree is found.
#' 
#' Equally named rows are interpreted as equal. Further
#' matching based on \code{rowData} is not done. For samples, collapsing 
#' is disabled by default meaning that equally named samples that are stored 
#' in different objects are interpreted as unique. Collapsing can be enabled 
#' with \code{collapse.cols = TRUE} when equally named samples describe the same
#' sample. 
#' 
#' If, for example, all rows are not shared with
#' individual objects, there are missing values in \code{assays}.
#' The notation of missing
#' can be specified with the \code{missing.values} argument. If input consists
#' of
#' \code{TreeSummarizedExperiment} objects, also \code{rowTree}, \code{colTree},
#' and
#' \code{referenceSeq} are preserved if possible. The data is preserved if 
#' all the rows or columns can be found from it.
#' 
#' Compared to \code{cbind} and \code{rbind} \code{mergeSEs} 
#' allows more freely merging since \code{cbind} and \code{rbind} expect 
#' that rows and columns are matching, respectively.
#' 
#' You can choose joining methods from \code{'full'}, \code{'inner'},
#'  \code{'left'}, and  \code{'right'}. In all the methods, all the samples are 
#'  included in the result object. However, with different methods, it is
#'  possible 
#'  to choose which rows are included.
#' 
#' \itemize{
#'   \item{\code{full} -- all unique features}
#'   \item{\code{inner} -- all shared features}
#'   \item{\code{left} -- all the features of the first object}
#'   \item{\code{right} -- all the features of the second object}
#' }
#' 
#' The output depends on the input. If the input contains
#' \code{SummarizedExperiment}
#' object, then the output will be \code{SummarizedExperiment}. When all the
#' input
#' objects belong to \code{TreeSummarizedExperiment}, the output will be 
#' \code{TreeSummarizedExperiment}.
#'
#' @seealso
#' \itemize{
#'   \item{\code{TreeSummarizedExperiment::cbind}}
#'   \item{\code{TreeSummarizedExperiment::rbind}}
#'   \item{\code{\link[dplyr:full_join]{full_join}}}
#'   \item{\code{\link[dplyr:inner_join]{inner_join}}}
#'   \item{\code{\link[dplyr:left_join]{left_join}}}
#'   \item{\code{\link[dplyr:right_join]{right_join}}}
#' }
#'
#' @name mergeSEs
#' @export
#' 
#' @examples
#' data(GlobalPatterns)
#' data(esophagus)
#' data(enterotype)
#' 
#' # Take only subsets so that it wont take so long
#' tse1 <- GlobalPatterns[1:100, ]
#' tse2 <- esophagus
#' tse3 <- enterotype[1:100, ]
#' 
#' # Merge two TreeSEs
#' tse <- mergeSEs(tse1, tse2)
#' 
#' # Merge a list of TreeSEs
#' list <- SimpleList(tse1, tse2, tse3)
#' tse <- mergeSEs(list, assay.type = "counts", missing.values = 0)
#' tse
#' 
#' # With 'join', it is possible to specify the merging method. Subsets are used
#' # here just to show the functionality
#' tse_temp <- mergeSEs(tse[1:10, 1:10], tse[5:100, 11:20], join = "left")
#' tse_temp
#' 
#' # If your objects contain samples that describe one and same sample,
#' # you can collapse equally named samples to one by specifying 'collapse.cols'
#' tse_temp <- mergeSEs(list(tse[1:10, 1], tse[1:20, 1], tse[1:5, 1]), 
#'                        collapse.cols = TRUE,
#'                        join = "inner")
#' tse_temp
#' 
#' # Merge all available assays
#' tse <- transformAssay(tse, method="relabundance")
#' ts1 <- transformAssay(tse1, method="relabundance")
#' tse_temp <- mergeSEs(tse, tse1, assay.type = assayNames(tse))
#' 
NULL

###################### Function for SimpleList of TreeSEs ######################

#' @rdname mergeSEs
#' @export
setMethod("mergeSEs", signature = c(x = "SimpleList"),
        function(x, assay.type="counts", assay_name = NULL, join = "full",
                missing.values = missing_values, missing_values = NA, 
                collapse.cols = collapse_samples, collapse_samples = FALSE, 
                collapse.rows = collapse_features, collapse_features = TRUE,
                verbose = TRUE, ... ){
            ################## Input check ##################
            # Check the objects 
            class <- .check_objects_and_give_class(x)
            if (!is.null(assay_name) & is.null(assay.type)) {
                .Deprecated(new="assay.type", old="assay_name", msg=paste0(
                    "The argument assay_name is deprecated and replace ",
                "with assay.type"))
            assay.type <- assay_name
            } else if (!is.null(assay_name) & !is.null(assay.type)) {
                warning("The assay.type argument is used and assay_name is ",
                        "ignored")
            }
            # CHeck which assays can be found, and if any --> FALSE
            assay.type <- .assays_cannot_be_found(assay.type = assay.type, x)
            if( .is_a_bool(assay.type) && assay.type == FALSE ){
                stop("'assay.type' must specify an assay from assays. ",
                    "'assay.type' cannot be found at least in one SE object.",
                    call. = FALSE)
            }

            # Check join
            if( !(.is_a_string(join) &&
                join %in% c("full", "inner", "left", "right") ) ){
                stop("'join' must be 'full', 'inner', 'left', or 'right'.",
                    call. = FALSE)
            }
            # Check if join is not available
            if( length(x) > 2 &&
                join %in% c("left", "right") ){
                stop("Joining method 'left' and 'right' are not available ",
                    "when more than two objects are being merged.",
                    call. = FALSE)
            }
            # Is missing.values one of the allowed ones
            missing_values_bool <- length(missing.values) == 1L &&
                (is.numeric(missing.values) && missing.values == 0) ||
                .is_a_string(missing.values) || is.na(missing.values)
            # If not then give error
            if(  !missing_values_bool ){
                stop("'missing.values' must be 0, NA, or a single character ",
                    "value.", call. = FALSE)
            }
            # Check collapse.cols
            if( !.is_a_bool(collapse.cols) ){
                stop("'collapse.cols' must be TRUE or FALSE.", call. = FALSE)
            }
            # Check collapse.cols
            if( !.is_a_bool(collapse.rows) ){
                stop("'collapse.rows' must be TRUE or FALSE.", call. = FALSE)
            }
            # Check verbose
            if( !.is_a_bool(verbose) ){
                stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
            }

            ################ Input check end ################
            # Give message if TRUE
            if( verbose ){
                message("Merging with ", join, " join...")
                message("1/", length(x), appendLF = FALSE)
            }
            # Merge objects
            tse <- .merge_SEs(
                x, class, join, assay.type, missing.values, collapse.cols,
                collapse.rows, verbose, ...)
            return(tse)
        }
)

########################### Function for two TreeSEs ###########################

#' @rdname mergeSEs
#' @export
setMethod("mergeSEs", signature = c(x = "SummarizedExperiment"),
        function(x, y = NULL, ...){
            ################## Input check ##################
            # Check y
            if( !(is(y, "SummarizedExperiment")) ){
                stop("'y' must be a 'SummarizedExperiment' object.",
                    call. = FALSE)
            } 
            ################ Input check end ################
            # Create a list based on TreeSEs
            list <- SimpleList(x, y)
            # Call the function for list
            mergeSEs(list, ...)
        }
)

########################### Function for list TreeSEs ##########################

#' @rdname mergeSEs
#' @export
setMethod("mergeSEs", signature = c(x = "list"),
    function(x, ...){
        # Convert into a list
        x <- SimpleList(x)
        # Call the function for list
        mergeSEs(x, ...)
    }
)

################################ HELP FUNCTIONS ################################

################################## .merge_SEs ##################################
# This function merges SE objects into one SE

# Input: A list of SEs
# Output: SE

#' @importFrom SingleCellExperiment SingleCellExperiment
.merge_SEs <- function(
        x, class, join, assay.type, missing.values, collapse.cols,
        collapse.rows, verbose, ...){

    # Take first element and remove it from the list
    tse <- x[[1]]
    x[[1]] <- NULL
    # Add rowData info to rownames
    rownames_name <- "rownames_that_will_be_used_to_adjust_names"
    tse <- .add_rowdata_to_rownames(tse, rownames_name = rownames_name)

    # Initialize a list for TreeSE-specific slots
    tse_args <- list(
        rowTrees = NULL,
        colTrees = NULL,
        refSeqs = NULL
    )
    # If the class is TreeSE, get TreeSE-specific slots
    if( class == "TreeSummarizedExperiment" ){
        tse_args <- .get_TreeSE_args(tse, tse_args)
    }

    # Get the data in a list
    args <- .get_SummarizedExperiment_data(tse = tse, assay.type = assay.type)

    # Get the function based on class
    FUN_constructor <- switch(
        class,
        TreeSummarizedExperiment = TreeSummarizedExperiment,
        SingleCellExperiment = SingleCellExperiment,
        SummarizedExperiment = SummarizedExperiment
    )
    # Create an object
    tse <- do.call(FUN_constructor, args = args)

    # Loop through individual TreeSEs and add them to tse
    if( length(x) > 0 ){
        for( i in seq_len(length(x)) ){
            # Give message if TRUE
            if( verbose ){
                message("\r", i+1, "/", length(x)+1, appendLF = FALSE)
            }
            
            # Get the ith object
            temp <- x[[i]]
            # Add rownames to rowData so that full matches are found
            temp <- .add_rowdata_to_rownames(
                temp, rownames_name = rownames_name)

            # Modify names if specified
            if( !collapse.cols ){
                temp <- .get_unique_names(tse, temp, "col")
            }
            if( !collapse.rows ){
                temp <- .get_unique_names(tse, temp, "row")
            }
            # Merge data
            args <- .merge_SummarizedExperiments(
                tse1 = tse,
                tse2 = temp,
                join = join,
                assay.type = assay.type,
                missing.values = missing.values
                )
            # If class is TreeSE, get trees and links, and reference sequences
            if( class == "TreeSummarizedExperiment" ){
                tse_args <- .get_TreeSE_args(temp, tse_args)
            }
            
            # Create an object
            tse <- do.call(FUN_constructor, args = args)
        }
    }
    # Add new line to, so that possible warning or  message has new line
    if( verbose ){
        message("")
    }
    # Get the data
    rowTrees <- tse_args$rowTrees
    colTrees <- tse_args$colTrees
    refSeqs <- tse_args$refSeqs
    # If data includes rowTrees, add them
    if( !is.null(rowTrees) ){
        tse <- .check_and_add_trees(tse, rowTrees, "row", verbose, ...)
    }
    # If data includes colTrees, add them
    if( !is.null(colTrees) ){
        tse <- .check_and_add_trees(tse, colTrees, "col", verbose, ...)
    }
    # If data includes reference sequences, add them
    if( !is.null(refSeqs) ){
        tse <- .check_and_add_refSeqs(tse, refSeqs, verbose)
    }
    # Adjust rownames
    rownames(tse) <- rowData(tse)[[rownames_name]]
    rowData(tse)[[rownames_name]] <- NULL
    # Ensure that rownames and colnames are unique. TreeSE allows equal
    # rownames. However, if in this points rownames are equal, it means that
    # rows are not really equal. They might have different taxonomy data, but
    # their ID might be same. This can cause problems.
    rownames(tse) <- make.unique(rownames(tse))
    colnames(tse) <- make.unique(colnames(tse))
    return(tse)
}

########################### .add_rowdata_to_rownames ###########################
# This function adds taxonomy information to rownames to enable more specific
# match between rows

# Input: (Tree)SE, name of the column that is being added to rowData
# Output: (Tree)SE with rownames that include all taxonomy information
.add_rowdata_to_rownames <- function(
        x, rownames_name, only.taxonomy = TRUE, ...){
    # Input check
    if( !.is_a_bool(only.taxonomy) ){
        stop("'only.taxonomy' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Add rownames to rowData
    rowData(x)[[rownames_name]] <- rownames(x)
    # Get rowData
    rd <- rowData(x)
    # If specified, get specify column only based on detected taxonomy columns.
    # Otherwise, use all the data to specify rows.
    if( only.taxonomy ){
        tax_cols <- match(
            c(tolower(TAXONOMY_RANKS), rownames_name),
            tolower(colnames(rd)), nomatch = 0 )
        tax_cols <- colnames(rd)[tax_cols]
    } else{
        tax_cols <- colnames(rd)
    }
    # Get taxonomy_info
    taxonomy_info <- rd[ , tax_cols, drop = FALSE]
    # Combine taxonomy info
    rownames <- apply(
        taxonomy_info, 1, function(x) paste0(x[!is.na(x)], collapse = "_"))
    # Add new rownames
    rownames(x) <- rownames
    return(x)
}

############################ .check_and_add_refSeqs ############################
# This function check if reference sequences can be added, and adds them if it
# is possible

# Input: reference sequences and TreeSE
# Output: TreeSE
.check_and_add_refSeqs <- function(tse, refSeqs, verbose){
    # Give message if wanted
    if( verbose ){
        message("Adding referenceSeqs...")
    }
    
    # Get the rownames that are included in reference sequences
    rows_that_have_seqs <- lapply(refSeqs, FUN = function(x){
        names(x[[1]])
    })
    rows_that_have_seqs <- unlist(rows_that_have_seqs)
    # Check that all the rownames are included
    if( !all(rownames(tse) %in% rows_that_have_seqs) ||
            is.null(rownames(tse)) ){
        warning("referenceSeqs do not match with the data so they are ",
                "discarded.", call. = FALSE)
        return(tse)
    }
    # Get the maximum number of DNA sets that individual TreeSE had / max
    # number of sets that individual rownames set had.
    max_numrow <- max(lengths(refSeqs))
    
    # Initialize a list
    result_list <- list()
    # Loop from 1 to max number of DNA sets
    for(i in seq_len(max_numrow) ){
        # Loop over DNA set list. Each element is found from unique TreeSE
        temp_seqs <- lapply(refSeqs, FUN = function(x){
            # If the ith element cannot be found, give the last
            if( i > length(x) ){
                return(x[[length(x)]])
            } else{
                # Otherwise give the ith element
                return(x[[i]])
            }
        })
        # Combine the list that includes DNA sets from unique TreeSEs.
        temp_seqs <- do.call(c, temp_seqs)
        # Get only those taxa that are included in TreeSE
        temp_seqs <- temp_seqs[ match(rownames(tse), names(temp_seqs)), ]
        # Add combined sequences into a list
        result_list <- c(result_list, temp_seqs)
    }
    # Create a DNAStrinSetList if there are more than one element
    if(length(result_list) > 1){
        result <- do.call(DNAStringSetList, result_list)
    } else{
        # Otherwise, give the only DNA set as it is
        result <- result_list[[1]]
    }
    # Add it to the correct slot
    referenceSeq(tse) <- result
    return(tse)
}

############################# .check_and_add_trees #############################
# This function check if tree can be added, and adds it if it can

# Input: tree data and TreeSE
# Output: TreeSE
#' @importFrom TreeSummarizedExperiment changeTree
.check_and_add_trees <- function(
        tse, trees_and_links, MARGIN = "row", verbose = FALSE, ...){
    # Give a message if verbose is specified
    if( verbose ){
        message("Merging ", MARGIN, "Tree...")
    }
    # Get trees
    trees <- trees_and_links$trees
    # Get links
    links <- trees_and_links$links
    # Based on margin, get rownames or colnames of the TreeSE object; to check
    # if the data matches with trees
    if(MARGIN == "row"){
        names <- rownames(tse)
    } else{
        names <- colnames(tse)
    }
    # All rownames/colnames should be included in trees/links
    if( !all(names %in% links[["names"]]) ||
        is.null(names) || length(names) == 0 ){
        warning(MARGIN, "Tree(s) does not match with the data so it ", 
                "is discarded.", call. = FALSE)
        return(tse)
    }
    # Combine trees into single tree.
    tree <- .merge_trees(trees, links, ...)
    # Order links so that the order matches with TreeSE
    links <- links[match(rownames(tse), links[["names"]]), ]
    # Add the data in correct slot based on MARGIN
    args <- list(tse, tree, links[["nodeLab"]])
    arg_names <- switch(
        MARGIN,
        "row" = c("x", "rowTree", "rowNodeLab"),
        "col" = c("x", "colTree", "colNodeLab"))
    names(args) <- arg_names
    tse <- do.call(changeTree, args)
    return(tse)
}

################################# .merge_trees #################################
# This function merges list of trees into single tree.

# Input: list of trees and link DataFrame
# Output: single tree
#' @importFrom ape bind.tree as.phylo
#' @importFrom dplyr as_tibble
.merge_trees <- function(trees, links, ...){
    # Bind trees to combine one large tree
    # Take first tree
    tree <- trees[[1]]
    trees[[1]] <- NULL
    # Loop through trees and merge them. The merging is done incrementally to
    # avoid creating huge trees that exceed memory.
    for( t in trees ){
        # Bind from root node if available. If not, then bind from node 0.
        tree <- bind.tree(tree, t)
        # Prune the tree so that it includes rows in tips. This step removes
        # additional tips, i.e., only tips that are in rows are preserved. Also
        # it simplifies the structure preserving the necessary information on
        # the dataset. Moreover, it ensures that there are no duplicated
        # tips which might be the case if the merged trees had shared taxa in
        # addition to unique taxa.
        tree <- .prune_tree(tree, links[["nodeLab"]], ...)
        # The trees are bind without merging. This means that we can have
        # duplicated nodes and branches. For instance, there can be a node
        # "family x" which is present in two trees that were merged. Moreover,
        # descendant nodes of these "family x" nodes can differ, which means
        # that we cannot just remove duplicated nodes. Instead, we have to
        # relink nodes so that each node label is present only one time and
        # all its child nodes are preserved.
        if( any(duplicated( c(tree$tip.label, tree$node.label) )) ){
            # Convert to table so that we can modify the data
            old_tree <- tree <- as_tibble(tree)
            # Remove duplicated nodes
            tree <- tree[ !duplicated(tree[["label"]]), ]
            # Reindex nodes
            tree[["node"]] <- seq_len(nrow(tree))
            # Reorder the old tree to match new trees parent node order
            old_tree <- old_tree[ match(tree[["parent"]], old_tree[["node"]]), ]
            # Reindex parent nodes of new tree
            parent <- tree[ match(old_tree[["label"]], tree[["label"]]), ]
            parent <- parent[["node"]]
            tree[["parent"]] <- parent
            # Convert back to phylo object
            tree <- as.phylo(tree)
        }
    }
    return(tree)
}

############################### .get_TreeSE_args ###############################
# This function fetches TreeSummarizedExperiment specific data: rowTree,
# colTree, and referenceSeq

# Input: TreeSE and argument list
# Output: An argument list
.get_TreeSE_args <- function(tse, tse_args){
    # If rowTree slot is not NULL
    if( !is.null(tse@rowTree) ){
        # Get trees that will be added
        trees_add <- tse@rowTree
        # Get rowLinks, convert them to basic DataFrame, 
        # so that additional column can be added
        links <- DataFrame(rowLinks(tse))
        # Add rownames as one of the columns
        links$names <- rownames(tse)
        
        # If there is no data yet / if rowTree arguments are NULL
        if( is.null(tse_args$rowTrees) ){
            # Get the tree data as a list. Tree is as a list, and links as DF
            rowTrees <- list(
                trees = trees_add,
                links = links
            )
            # Replace NULL with tree data
            tse_args$rowTrees <- rowTrees
        } else{
            # If tree data already exist
            # How many trees there already are
            tree_num_before <- length(tse_args$rowTrees$tree)
            # Get unique names
            unique_names <- make.unique( 
                names( c(tse_args$rowTrees$tree, trees_add) )
            )
            # Update the names of current data
            names(tse_args$rowTrees$tree) <- unique_names[ tree_num_before ]
            # Get unique names of trees that will be added
            unique_names_add <- unique_names[ -seq_len(tree_num_before) ]
            # Get corresponding current names
            names_add <- names(trees_add)
            # Update tree names from links
            links[ , "whichTree" ] <- 
                unique_names_add[ match( links[ , "whichTree" ], names_add ) ]
            # Update tree names
            names(trees_add) <- unique_names_add
            # Add data to a list
            tse_args$rowTrees <- list( 
                trees = c(tse_args$rowTrees$trees, trees_add),
                links = rbind(tse_args$rowTrees$links, links)
                )
        }
    }
    # If colTree slot is not NULL
    if( !is.null(tse@colTree) ){
        # Get trees that will be added
        trees_add <- tse@rowTree
        # Get colLinks, convert them to basic DataFrame, 
        # so that additional column can be added
        links <- DataFrame(colLinks(tse))
        # Add colnames as one of the columns
        links$names <- colnames(tse)
        
        # If there is no data yet / if colTree arguments are NULL
        if( is.null(tse_args$colTrees) ){
            # Get the tree data as a list. Tree is as a list, and links as DF
            colTrees <- list(
                trees = trees_add,
                links = links
            )
            # Replace NULL with tree data
            tse_args$colTrees <- colTrees
        } else{
            # If tree data already exist
            # How many trees there already are
            tree_num_before <- length(tse_args$colTrees$tree)
            # Get unique names
            unique_names <- make.unique( 
                names( c(tse_args$colTrees$tree, trees_add) )
            )
            # Update the names of current data
            names(tse_args$colTrees$tree) <- unique_names[ tree_num_before ]
            # Get unique names of trees that will be added
            unique_names_add <- unique_names[ -seq_len(tree_num_before) ]
            # Get corresponding current names
            names_add <- names(trees_add)
            # Update tree names from links
            links[ , "whichTree" ] <- 
                unique_names_add[ match( links[ , "whichTree" ], names_add ) ]
            # Update tree names
            names(trees_add) <- unique_names_add
            # Add data to a list
            tse_args$rowTrees <- list( 
                trees = c(tse_args$colTrees$trees, trees_add),
                links = rbind(tse_args$colTrees$links, links)
            )
        }
    }
    # If reference sequences exist
    if( !is.null(referenceSeq(tse)) ){
        # Get the data
        refSeq <- referenceSeq(tse)
        # Check if it is a individual set
        if( is(refSeq, "DNAStringSet") ){
            # Convert individual set to a list, so that all refseqs are in same
            # format
            refSeq <- DNAStringSetList(refSeq)
        }
        # Add data to a list
        refSeqs <- list(
            refSeq
        )
        # If there is no data yet, replace the NULL
        if( is.null(tse_args$refSeqs) ){
            tse_args$refSeqs <- refSeqs
        } else{
            # otherwise add data to a list
            tse_args$refSeqs <- c( tse_args$refSeqs, refSeqs ) 
        }
    }
    return(tse_args)
}

######################## .get_SummarizedExperiment_data ########################
# This function gets the desired data from one SE object and creates a list of 
# arguments containing the data
# Arguments of SCE and TreeSE are also fetched with this function.
# TreeSE-specific slots are collected with different function so that they are
# merged at the end.

# Input: SE
# Output: A list of arguments
.get_SummarizedExperiment_data <- function(tse, assay.type){
    # Remove all information but rowData, colData, metadata and assay
    row_data <- rowData(tse)
    col_data <- colData(tse)
    assays <- assays(tse)[ assay.type ]
    metadata <- metadata(tse)
    # Create a list of arguments
    args <- list(assays = assays,
                rowData = row_data,
                colData = col_data,
                metadata = metadata
    )
    return(args)
    
}

######################## .check_objects_and_give_class #########################
# This function checks that the object are in correct format

# Input: a list of objects
# Output: A shared class of objects
.check_objects_and_give_class <- function(x){
    # Allowed classes
    allowed_classes <- c("TreeSummarizedExperiment", "SingleCellExperiment",
        "SummarizedExperiment")
    
    # Get the class based on hierarchy TreeSE --> SCE --> SE
    # and check that objects are in correct format
    classes <- lapply(x, .check_object_for_merge)
    classes <- unlist(classes)
    # Get the shared class that is highest in hierarchy
    if( all( classes %in% allowed_classes[1] ) ){
        class <- allowed_classes[1]
    } else if( all( classes %in% allowed_classes[c(1,2)] ) ){
        class <- allowed_classes[2]
    } else {
        class <- allowed_classes[3]
    }
    
    # If there are multiple classes, give a warning
    if( length(unique( classes )) > 1 ){
        warning("The Input consist of multiple classes. ",
                "The output is '", class, "'.",
                call. = FALSE)
    }
    return(class)
}

########################### .check_object_for_merge ############################
# This function checks an object that it is in correct format. Additionally, it
# returns its class

# Input: (Tree)SE
# Output: Class of (Tree)SE
.check_object_for_merge <- function(x){
    # Check that the class matches with supported ones
    if( !is(x, "SummarizedExperiment") ){
        stop("Input includes an object that is not 'SummarizedExperiment'.",
            call. = FALSE)
    }
    # Check that there are no object with no dimensions
    if( ncol(x) == 0 || nrow(x) == 0 ){
        stop("Input includes an object that has either no columns or/and ",
            "no rows.", call. = FALSE)
    }
    # Check that object has row/colnames
    if( is.null(rownames(x)) || is.null(colnames(x)) ){
        stop("Input includes object(s) whose rownames and/or colnames is ",
            "NULL. Please add them.", call. = FALSE)
    }
    # Check if the col/rownames are duplicated
    if( any(duplicated(rownames(x))) || any(duplicated(colnames(x))) ){
        stop("Input includes object(s) whose rownames and/or colnames include ",
            "duplicates. Please make them unique.", call. = FALSE)
    }
    # Get class
    class <- class(x)
    return(class)
    
}
########################### .assays_cannot_be_found ############################
# This function checks that the assay(s) can be found from TreeSE objects of
# a list.

# Input: the name of the assay and a list of TreeSE objects
# Output: A list of assay.types that can be found or FALSE if any
.assays_cannot_be_found <- function(assay.type, x){
    # Loop through objects
    assays <- lapply(x, FUN = function(tse){
        # Check if the assay.types can be found. If yes, then TRUE. If not,
        # then FALSE
        temp <- lapply(assay.type, .assay_cannot_be_found, tse = tse)
        # Unlist and return
        return( unlist(temp) )
    })
    # Create a data.frame from the result
    assays <- as.data.frame(assays, row.names = assay.type)
    colnames(assays) <- paste0("tse", seq_len(length(assays)))
    # Which assays can be found from all the objects?
    assays <- rownames(assays)[ rowSums(assays) == ncol(assays) ]
    # If none of assays were found, return FALSE
    if( length(assays) == 0 ){
        assays <- FALSE
    }
    # Give warning if assays were dropped
    if( length(assays) < length(assay.type) ){
        warning("The following assay(s) was not found from all the objects ",
                "so it is dropped from the output: ",
                paste0("'", setdiff(assay.type, assays), sep = "'",
                    collapse = ", "),
                call. = FALSE)
    }
    return(assays)
}

############################ .assay_cannot_be_found ############################
# This function checks that the assay can be found from TreeSE. If it can be
# found --> TRUE, if it cannot be found --> FALSE

# Input: the name of the assay and TreSE object
# Output: TRUE or FALSE
.assay_cannot_be_found <- function(assay.type, tse){
    # Check if the assay.type can be found. If yes, then TRUE. If not, then
    # FALSE
    tryCatch(
        {
            .check_assay_present(assay.type, tse)
            return(TRUE)
            
        },
        error = function(cond) {
            return(FALSE)
        }
    )
}

########################### ..get_unique_names ###########################
# This function convert colnames unique

# Input: TreeSEs and MARGIN
# Output: One TreeSE with unique sample names compared to other TreeSE
.get_unique_names <- function(tse1, tse2, MARGIN, suffix=2){
    # Based on MARGIN, get right names
    if( MARGIN == "row" ){
        names1 <- rownames(tse1)
        names2 <- rownames(tse2)
    } else{
        names1 <- colnames(tse1)
        names2 <- colnames(tse2)
    }
    # If there are duplicated names
    if( any(names2 %in% names1) ){
        # Get duplicated names
        ind <- names2 %in% names1
        temp_names2 <- names2[ind]
        # Get unique suffix
        while( any(paste0(names2, ".", suffix) %in% names1) ){
            suffix <- suffix + 1
        }
        temp_names2 <- paste0(temp_names2, ".", suffix)
        # Assign names back
        if( MARGIN == "row" ){
            rownames(tse2)[ind] <- temp_names2
        } else{
            colnames(tse2)[ind] <- temp_names2
        }
    }
    return(tse2)
}

######################## .merge_SummarizedExperiments ##########################
# This function merges the data of two SE objects into one set of arguments that
# can be feed to create a single object.
# TreeSE and SCE are all merged with this function since SCE or TreeSE-specific
# slots are not merged at this point. TreeSE-specific slots are collected and
# merged at the end.

# Input: Two SEs
# Output: A list of arguments
.merge_SummarizedExperiments <- function(
        tse1, tse2, join, assay.type, missing.values){
    # Merge rowData
    rowdata <- .merge_rowdata(tse1, tse2, join)
    # Merge colData
    coldata <- .merge_coldata(tse1, tse2, join)
    # Merge assays
    assays <- lapply(assay.type, .merge_assay,
                    tse1 = tse1, tse2 = tse2,
                    join = join, missing.values = missing.values,
                    rd = rowdata, cd = coldata)
    assays <- SimpleList(assays)
    names(assays) <- assay.type
    # Combine metadata
    metadata <- c( metadata(tse1), metadata(tse2) )
    
    # Create a list of data
    args <- list(assays = assays,
                rowData = rowdata,
                colData = coldata,
                metadata = metadata)
    return(args)
}

################################ .merge_assay ##################################
# This function merges assays.

# Input: Two TreeSEs, the name of the assay, joining method, value to denote
# missing values, merged rowData, and merged colData
# Output: Merged assay
.merge_assay <- function(
        tse1, tse2, assay.type, join, missing.values, rd, cd){
    # Take assays
    assay1 <- assay(tse1, assay.type)
    assay2 <- assay(tse2, assay.type)
    
    # Merge two assays into one
    assay <- .join_two_tables(assay1, assay2, join)
    
    # Convert into matrix
    assay <- as.matrix(assay)
    
    # Fill missing values
    assay[ is.na(assay) ] <- missing.values
    
    # Order the assay based on rowData and colData
    assay <- assay[ match(rownames(rd), rownames(assay)), , drop = FALSE ]
    assay <- assay[ , match(rownames(cd), colnames(assay)), drop = FALSE]
    return(assay)
}

############################### .merge_rowdata #################################
# This function merges rowDatas,

# Input: Two TreeSEs and joining method
# Output: Merged rowData
.merge_rowdata <- function(tse1, tse2, join){
    # Take rowDatas
    rd1 <- rowData(tse1)
    rd2 <- rowData(tse2)
    
    # Convert column names to lower
    if( length(colnames(rd1)) > 0 ){
        colnames(rd1) <- tolower(colnames(rd1))
    }
    if( length(colnames(rd2)) > 0 ){
        colnames(rd2) <- tolower(colnames(rd2))
    }
    
    # Merge rowdata
    rd <- .join_two_tables(rd1, rd2, join)
    
    # There might be duplicated rownames. This might occur when there are
    # features with equal taxonomy data but merged datasets have some
    # additional info that do not match with each other. --> collapse
    # duplicated rows/features into one row.
    dupl_rows <- rownames(rd)[ duplicated(rownames(rd)) ]
    if( length(dupl_rows) > 0 ){
        for( r in dupl_rows ){
            # Get duplicated rows
            temp <- rd[rownames(rd) %in% r, , drop = FALSE]
            # Remove columns with no info
            temp <- temp[, vapply(temp, function(x)
                !all(is.na(x)), logical(1)), drop = FALSE]
            # 1st row is kept if the value is not NA
            keep1 <- vapply(temp, function(x) !is.na(x[1]),  logical(1))
            # 2nd row is kept if it is not NA and the value differs from 1st row
            keep2 <- vapply(temp, function(x)
                !is.na(x[2]) && (is.na(x[1]) ||
                (!is.na(x[1]) && x[2] != x[[1]])), logical(1))
            # Get the rows
            keep1 <- temp[1, keep1, drop = FALSE]
            keep2 <- temp[2, keep2, drop = FALSE]
            # If the two different, previously merged TreeSEs had same row,
            # but the row had variales with unequal values, there are variable
            # and variable.1 columns. Check if the value can be found from those
            # columns.
            for( col in colnames(keep2) ){
                # Get those column that have been added with suffix
                col_add <- grep(paste0(col, "."), colnames(keep1))
                if( length(col_add) > 0 ){
                    col_add <- keep1[, col_add, drop = FALSE]
                    # If the value can be found already, remove it, i.e., do
                    # not add it again
                    if( keep2[, col] %in% unlist(col_add) ){
                        keep2[, col] <- NULL
                    }
                }
            }
            # Combine rows together
            temp <- cbind(keep1, keep2)
            colnames(temp) <- make.unique(colnames(temp))
            # Remove the rows from the original rowData
            rd <- rd[!rownames(rd) %in% r, , drop = FALSE]
            # Add the data back
            rd[["rownames_merge_ID"]] <- rownames(rd)
            temp[["rownames_merge_ID"]] <- r
            rd <- merge(rd, temp, all = TRUE)
            # Ensure that the rowData is DF (it can handle duplicated rownames)
            rd <- DataFrame(rd)
            # Add original rownames
            rownames(rd) <- rd[["rownames_merge_ID"]]
            rd[["rownames_merge_ID"]] <- NULL
        }
        
    }
    
    # Get column indices that match with taxonomy ranks
    ranks_ind <- match( TAXONOMY_RANKS, colnames(rd) )
    # Remove NAs
    ranks_ind <- ranks_ind[ !is.na(ranks_ind) ]
    # If ranks were found
    if( length(ranks_ind) != 0 ){
        # Get the data in correct order, take only column that have ranks
        rd_rank <- rd[ , ranks_ind, drop = FALSE]
        # Take other columns
        rd_other <- rd[ , -ranks_ind, drop = FALSE]
        # Get rank names
        rank_names <- colnames(rd_rank)
        # Convert names that they have capital letters
        new_rank_names <- paste(
            toupper(substr(rank_names, 1, 1)),
            substr(rank_names, 2, nchar(rank_names)), sep = "")
        # Add new names to colnames of rd_rank
        colnames(rd_rank) <- new_rank_names
        
        # Combine columns
        rd <- cbind(rd_rank, rd_other)
    }
    return(rd)
}

############################### .merge_coldata #################################
# This function merges colDatas,

# Input: Two TreeSEs and joining method
# Output: Merged colData
.merge_coldata <- function(tse1, tse2, join){
    # Take colDatas
    cd1 <- colData(tse1)
    cd2 <- colData(tse2)
    
    # Merge coldata
    cd <- .join_two_tables(cd1, cd2, join = "full")
    # Convert into DataFrame
    cd <- DataFrame(cd)
    
    return(cd)
}

############################## .join_two_tables ################################
# This general function is used to merge rowDatas, colDatas, and assays.

# Input: Two tables and joining method
# Output: One merged table

#' @importFrom dplyr coalesce
.join_two_tables <- function(df1, df2, join){
    # Get parameter based on join
    all.x <- switch(join,
                    full = TRUE,
                    inner = FALSE,
                    left = TRUE,
                    right = FALSE
    )
    all.y <- switch(join,
                    full = TRUE,
                    inner = FALSE,
                    left = FALSE,
                    right = TRUE
    )
    
    # Ensure that the data is always data.frame
    df1 <- as.data.frame(df1)
    df2 <- as.data.frame(df2)
    
    # STEP 1: Check whether variables can be merged; if their classes are equal.
    # Adjust colnames if their classes differ
    if( ncol(df1) > 0 && ncol(df2) > 0 ){
        # Get classes of variables of df1
        c1 <- lapply(df1, FUN = function(x){c(class(x), !all(is.na(x)))})
        colnames1 <- names(c1)
        c1 <- data.frame(t(data.frame(c1)))
        # Get classes of variables of df2
        c2 <- lapply(df2, FUN = function(x){c(class(x), !all(is.na(x)))})
        colnames2 <- names(c2)
        c2 <- data.frame(t(data.frame(c2)))
        # Rename columns and convert columns to correct class
        colnames(c1) <- colnames(c2) <- c("class", "not_na")
        class(c1$class) <- class(c2$class) <- "character"
        class(c1$not_na) <- class(c2$not_na) <- "logical"
        # Add original colnames to one column
        c1$rownames <- colnames1
        c2$rownames <- colnames2
        # Merge class information into one df
        classes <- merge(c1, c2, by="rownames", all=TRUE)
        # Add info whether certain variable was found from df1/df2
        classes$found1 <- !is.na(classes$class.x)
        classes$found2 <- !is.na(classes$class.y)
        classes$found_both <- classes$found1 & classes$found2
        # Add information whether the classes of equally named variables match
        # Take into account if other df has only NA values in column --> class
        # of that variable can be wrong --> use known class
        classes$no_match <- FALSE
        classes[classes$found_both, "no_match"] <- 
            classes[classes$found_both, "class.x"] != classes[
                classes$found_both, "class.y"] & classes$not_na.x[
                    classes$found_both] & classes$not_na.y[classes$found_both]
        # Add new colnames to columns. If equally named variables' classes
        # differ add also class information to colnames
        classes$colnames1 <- classes$rownames
        classes$colnames2 <- classes$rownames
        classes[classes$no_match, "colnames1"] <-
            paste0(classes[classes$no_match, "colnames1"], "_",
                classes[classes$no_match, "class.x"])
        classes[classes$no_match, "colnames2"] <-
            paste0(classes[classes$no_match, "colnames2"], "_",
                classes[classes$no_match, "class.y"])
        # Give warning if there were missmatch between equally named variables
        # and their classes
        if( any(classes$no_match) ){
            warning("Datasets include equally named variables called '",
                "'but their class differ. In the output, variables are not ",
                "combined and they are renamed based on their class.",
                "Please check the following columns:\n",
                paste0("'", paste(
                    classes[classes$no_match, "rownames"], collapse = "', '"),
                    "'"),
                call. = FALSE)
        }
        # Add new column names to df1
        colnames <- classes[classes$found1, "rownames"]
        colnames <- classes[classes$found1, "colnames1"][
            match(colnames(df1), colnames)]
        colnames(df1) <- colnames
        # Add new column names to df2
        colnames <- classes[classes$found2, "rownames"]
        colnames <- classes[classes$found2, "colnames2"][
            match(colnames(df2), colnames)]
        colnames(df2) <- colnames
    }
    # Add rownames to one of the columns
    df1$rownames_merge_ID <- rownames(df1)
    df2$rownames_merge_ID <- rownames(df2)
    
    # STEP 2: merge
    # Finally merge data frames into one data frame
    df <- merge(df1, df2, all.x = all.x, all.y = all.y)
    
    # STEP 3: polish
    # Get column names because they can be changed when class is changed
    colnames <- colnames(df)
    # Convert to DF because there can be equally named columns (DataFrame does
    # allow those)
    df <- DataFrame(df)
    colnames(df) <- colnames
    # Add original rownames and remove additional column
    rownames(df) <- df$rownames_merge_ID
    df$rownames_merge_ID <- NULL
    return(df)
}
