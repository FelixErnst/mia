#' Split \code{TreeSummarizedExperiment} column-wise or row-wise based on
#' grouping variable
#'
#' @inheritParams agglomerate-methods
#'
#' @param group \code{Character scalar}, \code{character vector} or
#' \code{factor vector}. A column name from \code{rowData(x)} or
#' \code{colData(x)} or alternatively a vector specifying how the merging is
#' performed. If vector, the value must be the same length as
#' \code{nrow(x)/ncol(x)}. Rows/Cols corresponding to the same level will be
#' merged. If \code{length(levels(group)) == nrow(x)/ncol(x)}, \code{x} will be
#' returned unchanged. If \code{group} matches with both dimensions,
#' \code{by} must be specified. (Default: \code{NULL})
#'   
#' @param f Deprecated. Use \code{group} instead.
#' 
#' @param update_rowTree Deprecated. Use \code{update.tree} instead.
#'   
#' @param altexp \code{Character vector}. Specify the alternative experiments
#'   to be unsplit. (Default: \code{names(altExps(x))})
#' 
#' @param altExpNames Deprecated. Use \code{altexp} instead.
#'   
#' @param ... Arguments passed to \code{agglomerateByVariable} function for
#'   \code{SummarizedExperiment} objects and other functions.
#'   See \code{\link[=agglomerate-methods]{agglomerateByVariable}} for more 
#'   details.
#'   \itemize{
#'     \item \code{use.names}: \code{Logical scalar}. Specifies whether to name
#'     elements of
#'     list by their group names. (Default: \code{TRUE})
#'   }
#'
#'
#' @details
#' \code{splitOn} split data based on grouping variable. Splitting can be done
#' column-wise or row-wise. The returned value is a list of
#' \code{SummarizedExperiment} objects; each element containing members of each
#' group.
#'
#' @return
#' For \code{splitOn}: \code{SummarizedExperiment} objects in a
#' \code{SimpleList}.
#'
#' For \code{unsplitOn}: \code{x}, with \code{rowData} and \code{assay}
#' data replaced by the unsplit data. \code{colData} of x is kept as well
#' and any existing \code{rowTree} is dropped as well, since existing
#' \code{rowLinks} are not valid anymore.
#' 
#' @name splitOn
#' @seealso
#' \code{\link[=agglomerate-methods]{agglomerateByRanks}}
#' \code{\link[=agglomerate-methods]{agglomerateByVariable}},
#' \code{\link[scuttle:sumCountsAcrossFeatures]{sumCountsAcrossFeatures}},
#' \code{\link[=agglomerate-methods]{agglomerateByRank}},
#' \code{\link[SingleCellExperiment:altExps]{altExps}},
#' \code{\link[SingleCellExperiment:splitAltExps]{splitAltExps}}
#'
#' @export
#' 
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' # Split data based on SampleType. 
#' se_list <- splitOn(tse, group = "SampleType")
#' 
#' # List of SE objects is returned. 
#' se_list
#' 
#' # Create arbitrary groups
#' rowData(tse)$group <- sample(1:3, nrow(tse), replace = TRUE)
#' colData(tse)$group <- sample(1:3, ncol(tse), replace = TRUE)
#' 
#' # Split based on rows
#' # Each element is named based on their group name. If you don't want to name
#' # elements, use use_name = FALSE. Since "group" can be found from rowdata and
#' # colData you must use `by`.
#' se_list <- splitOn(tse, group = "group", use.names = FALSE, by = 1)
#' 
#' # When column names are shared between elements, you can store the list to
#' # altExps
#' altExps(tse) <- se_list
#' 
#' altExps(tse)
#' 
#' # If you want to split on columns and update rowTree, you can do
#' se_list <- splitOn(tse, group = colData(tse)$group, update.tree = TRUE)
#' 
#' # If you want to combine groups back together, you can use unsplitBy
#' unsplitOn(se_list)
#' 
NULL

# This function collects group (grouping variable), by, and 
# use.names and returns them as a list.
.norm_args_for_split_by <- function(
        x, group = f, f, by = MARGIN, MARGIN = NULL, use.names = use_names,
        use_names = TRUE, ...){
    # input check
    # Check group
    if(is.null(group)){
        stop("'group' must either be a single non-empty character value or",
            " vector coercible to factor alongside the one of the dimensions ",
            "of 'x'", call. = FALSE)
    }
    # Check by
    if( !is.null(by) ){
        by <- .check_MARGIN(by)
    }
    # If group is a vector containing levels
    if( !.is_non_empty_string(group) ){
        # Convert into factors
        group <- factor(group, unique(group))
        # Check if the length of group matches with one of the dimensions
        if(!length(group) %in% dim(x)){
            stop("'group' must either be a single non-empty character value or",
                " vector coercible to factor alongside the on of the ",
                "dimensions of 'x'.",
                call. = FALSE)
        # If it matches with both dimensions, give error if by is not specified
        } else if( is.null(by) && all(length(group) == dim(x)) ){
            stop("The length of 'group' matches with nrow and ncol. ",
                "Please specify 'by'.", call. = FALSE)
        # If by is specified but it does not match with length of group
        } else if( !is.null(by) && (length(group) !=  dim(x)[[by]]) ){
            stop("'group' does not match with ", 
                ifelse(by==1, "nrow", "ncol"), ". Please check 'by'.",
                call. = FALSE)
        # IF group matches with nrow
        } else if(length(group) == dim(x)[[1]] && is.null(by)  ){
            by <- 1L
        # If group matches with ncol
        } else if( is.null(by) ){
            by <- 2L
        }
    # Else if group is a character specifying column from rowData or colData  
    } else {
        # If by is specified
        if( !is.null(by) ){
            # Search from rowData or colData based on by
            dim_name <- switch(by,
                                "1" = "rowData",
                                "2" = "colData")
            # Specify right function
            dim_FUN <- switch(by,
                                "1" = retrieveFeatureInfo,
                                "2" = retrieveCellInfo)
            # Try to get information
            tmp <- try({dim_FUN(x, group, search = dim_name)},
                        silent = TRUE)
            # Give error if it cannot be found
            if(is(tmp,"try-error")){
                stop("'group' is not found. ",
                    "Please check that 'group' specifies a column from ",
                    dim_name, ".", call. = FALSE)
            }
            # Get values
            group <- tmp$value
        # Else if by is not specified
        } else{
            # Try to get information from rowData
            tmp_row <- try({retrieveFeatureInfo(x, group, search = "rowData")},
                            silent = TRUE)
            # Try to get information from colData
            tmp_col <- try({retrieveCellInfo(x, group, search = "colData")}, 
                            silent = TRUE)
            
            # If it was not found 
            if( is(tmp_row, "try-error") && is(tmp_col, "try-error") ){
                stop("'group' is not found. ",
                    "Please check that 'group' specifies a column from ",
                    "rowData or colData.", 
                    call. = FALSE)
                # If group was found from both
            } else if( !is(tmp_row, "try-error") && !is(tmp_col, "try-error") ){
                stop("'group' can be found from both rowData and colData. ",
                    "Please specify 'by'.",
                    call. = FALSE)
                # If it was found from rowData
            } else if( !is(tmp_row, "try-error") ){
                by <- 1L
                # Get values
                group <- tmp_row$value
                # Otherwise, it was found from colData
            } else{
                by <- 2L
                # Get values
                group <- tmp_col$value
            }
        }
        # Convert values into factors
        group <- factor(group, unique(group))
        
        # If there are NAs, add NA as level
        if( any(is.na(group)) ){
            group <- addNA(group)
        }
    }
    # Check use.names
    if( !.is_a_bool(use.names) ){
        stop("'use.names' must be TRUE or FALSE.",
            call. = FALSE)
    }
    # Create a list from arguments
    list(f = group,
        by = by,
        use.names = use.names)
}

# PErform the split
.split_on <- function(x, args, ...){
    # Get grouping variable and its values
    f <- args[["f"]]
    # Choose nrow or ncol based on by
    dim_FUN <- switch(args[["by"]],
                        "1" = nrow,
                        "2" = ncol)
    # Get indices from 1 to nrow/ncol
    idx <- seq_len(dim_FUN(x))
    # Split indices into groups based on grouping variable
    idxs <- split(idx, f)
    # Subset function takes SE and list of groups which have indices
    # It divides the data into groups
    subset_FUN <- function(x, i = TRUE, j = TRUE){
        x[i, j]
    }
    # Based on by, divide data in row-wise or column-wise
    if(args[["by"]] == 1){
        ans <- SimpleList(lapply(idxs, subset_FUN, x = x))
    } else {
        ans <- SimpleList(lapply(idxs, subset_FUN, x = x, i = TRUE))
    }
    # If user do not want to use names, unname
    if(!args[["use.names"]]){
        ans <- unname(ans)
    # Otherwise convert NAs to "NA", if there is a level that do not have name
    } else{
        names(ans)[ is.na(names(ans)) ] <- "NA"
    }
    ans
}

#' @rdname splitOn
#' @export
setMethod("splitOn", signature = c(x = "SummarizedExperiment"),
    function(x, group = f, f = NULL,  ...){
        # Get arguments
        args <- .norm_args_for_split_by(x, group = group, ...)
        # Split data
        .split_on(x, args, ...)
    }
)

#' @rdname splitOn
#' @export
setMethod("splitOn", signature = c(x = "SingleCellExperiment"),
    function(x, group = f, f = NULL, ...){
        # Get arguments
        args <- .norm_args_for_split_by(x, group = group, ...)
        # Should alternative experiment be removed? --> yes
        args[["altexp.rm"]] <- TRUE
        # Split data
        .split_on(x, args, ...)
    }
)

#' @rdname splitOn
#' @export
setMethod("splitOn", signature = c(x = "TreeSummarizedExperiment"),
    function(x, group = f, f = NULL, update.tree = update_rowTree,
            update_rowTree = TRUE, ...){
        # Input check
        # Check update.tree
        if( !.is_a_bool(update.tree) ){
            stop("'update.tree' must be TRUE or FALSE.",
                call. = FALSE)
        }
        # Input check end
        # Split data
        x <- callNextMethod()
        # Manipulate rowTree or not?
        if( update.tree ){
            # Update both colTree and rowTree
            for( direction in c(1, 2) ){
                # If the returned value is a list, go through all of them
                if( is(x, "SimpleList") ){
                    x <- lapply(x, function(y){
                        .agglomerate_trees(y, MARGIN = direction, ...)})
                    x <- SimpleList(x)
                } else {
                    # Otherwise, the returned value is TreeSE
                    x <- .agglomerate_trees(x, MARGIN = direction, ...)
                }
            }
        }
        x
    }
)

################################################################################
# unsplitOn

# Perform the unsplit
.list_unsplit_on <- function(
        ses, update.tree = FALSE, by = MARGIN, MARGIN = NULL, ...){
    # Input check
    is_check <- vapply(ses,is,logical(1L),"SummarizedExperiment")
    if(!all(is_check)){
        stop("Input must be a list of SummarizedExperiment or derived objects ",
            "only.",
            call. = FALSE)
    }
    # Check update.tree
    if( !.is_a_bool(update.tree) ){
        stop("'update.tree' must be TRUE or FALSE.",
            call. = FALSE)
    }
    # Check by
    if( !is.null(by) ){
        by <- .check_MARGIN(by)
    }
    # Input check end
    # If list contains only one element, return it
    if( length(ses) == 1 ){
        return(ses[[1]])
    }
    # Get dimensions of each SE in the list
    dims <- vapply(ses, dim, integer(2L))
    # Based on which dimension SE objects share, select `by`.
    # If they share rows, then by is col, and vice versa
    if( is.null(by) ){
        if( length(unique(dims[1L,])) == 1 && length(unique(dims[2L,])) == 1 ){
            stop("The dimensions match with row and column-wise. ",
                "Please specify 'by'.", call. = FALSE)
        } else if(length(unique(dims[1L,])) == 1L){
            by <- 2L
        } else if(length(unique(dims[2L,])) == 1L) {
            by <- 1L
        } else {
            stop("The dimensions are not equal across all elements. ", 
                "Please check that either number of rows or columns match.", 
                call. = FALSE)
        }
    } else{
        # Get correct dimension, it is opposite of by
        dim <- ifelse(by == 1, 2, 1)
        if( length(unique(dims[dim,])) != 1L ){
            stop("The dimensions are not equal across all elements.",
                call. = FALSE)
        }
    }
    
    # Get the class of objects SCE, SE or TreeSE
    class_x <- class(ses[[1L]])
    # Combine assays
    args <- list(assays = .unsplit_assays(ses, MARGIN = by))
    # Combine rowData if data share columns
    if(by == 1L){
        rd <- .combine_rowData(ses)
        # Add rownames since they are missing after using combining
        rownames(rd) <- unlist(unname(lapply(ses, rownames)))
        rr <- .combine_rowRanges(ses)
        args$rowRanges <- rr
        args$colData <- colData(ses[[1L]])
    # Combine colData if data share rows
    } else {
        args$colData <- .combine_colData(ses)
        args$rowRanges <- rowRanges(ses[[1L]])
        rd <- rowData(ses[[1L]])
    }
    # Create a object specified by class_x from the data
    ans <- do.call(class_x, args)
    # Add rowData
    rowData(ans) <- rd
    # Update rownames
    rownames(ans) <- rownames(rd)
    
    # IF the object is TreeSE. add rowTree
    if( class_x == "TreeSummarizedExperiment" ){
        # Add both colTree and rowTree
        for( direction in c("row", "col") ){
            ans <- .add_trees(ans, ses, direction, update.tree = update.tree)
        }
    }
    return(ans)
}

# This function adds rowTrees and colTrees to result.
.add_trees <- function(tse, tses, MARGIN, update.tree){
    # Get functions based on MARGIN
    tree_FUN <- switch(MARGIN, "row" = rowTree, "col" = colTree)
    tree_name_FUN <- switch(MARGIN, "row" = rowTreeNames, "col" = colTreeNames)
    tree_assign_FUN <- switch(MARGIN, "row" = `rowTree<-`, "col" = `colTree<-`)
    link_FUN <- switch(MARGIN, "row" = rowLinks, "col" = colLinks)
    name_FUN <- switch(MARGIN, "row" = rownames, "col" = colnames)
    # Update or add old tree from the first element of list
    if( update.tree ){
        # Get trees as list.
        trees <- list()
        for( x in tses ){
            temp <- tree_FUN(x, tree_name_FUN(x))
            if( !is(temp, "list") ){
                temp <- list(temp)
                names(temp) <- tree_name_FUN(x)
            }
            trees <- c(trees, temp)
        }
        # If there are trees available
        if( any(lengths(trees) > 0) ){
            # Get links
            links <- lapply(tses, link_FUN)
            # Links must be single DF and it must have "name" column denoting
            # row/colnames
            links <- do.call(rbind, links)
            links <- DataFrame(links)
            nams <- lapply(tses, name_FUN)
            nams <- unlist( nams[lengths(trees) > 0] )
            links[["names"]] <- nams
            # Combine trees into single tree
            args <- list(trees = trees, links = links)
            tse <- .check_and_add_trees(tse, args, MARGIN = MARGIN)
        }
    } else{
        tree <- tree_FUN(tses[[1L]])
        tse <- tree_assign_FUN(tse, value = tree)
    }
    return(tse)
}

#' @importFrom SummarizedExperiment colData
#' @importFrom BiocGenerics rbind
.combine_colData <- function(ses) {
    # Get colDatas of objects
    cds <- lapply(ses, colData)
    # Bind them together row-wise
    cd <- do.call(rbind,unname(cds))
    # Add sample names
    rownames(cd) <- unlist(unname(lapply(ses, colnames)))
    cd
}


#' @rdname splitOn
#' @importFrom SingleCellExperiment altExpNames altExp altExps
#' @export
setMethod("unsplitOn", signature = c(x = "list"),
    function(x, update.tree = update_rowTree, update_rowTree = TRUE, ...){
        # Unsplit list and create SCE, SE, or TreeSE from it
        .list_unsplit_on(x, update.tree, ...)
    }
)
#' @rdname splitOn
#' @importFrom SingleCellExperiment altExpNames altExp altExps
#' @export
setMethod("unsplitOn", signature = c(x = "SimpleList"),
    function(x, update.tree = update_rowTree, update_rowTree = TRUE, ...){
        unsplitOn(as.list(x), update.tree, ...)
    }
)

#' @rdname splitOn
#' @importFrom SingleCellExperiment altExpNames altExp altExps reducedDims<-
#' @export
setMethod("unsplitOn", signature = c(x = "SingleCellExperiment"),
    function(x, altexp = altExpNames, altExpNames = names(altExps(x)),
            keep.dimred = keep_reducedDims,
            keep_reducedDims = FALSE, ...){
        # input check
        if(!.is_a_bool(keep.dimred)){
            stop("'keep.dimred' must be TRUE or FALSE.", call. = FALSE)
        }
        # Get alternative experiment names since data is located there
        ae_names <- names(altExps(x))
        # Get only those experiments that user has specified
        ae_names <- ae_names[ae_names %in% altexp]
        if(length(ae_names) == 0L){
            stop("No altExp matching 'altexp' in name.", call. = FALSE)
        }
        # Get alternative experiments as a list
        ses <- altExps(x)[ae_names]
        # And unsplit the data
        se <- .list_unsplit_on(ses, ...)
        # Add reducedDims if specified
        if( keep.dimred ){
            reducedDims(se) <- reducedDims(x)
        }
        return(se)
    }
)
