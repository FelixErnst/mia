#' @importFrom ape drop.tip
#' @importFrom rbiom unifrac
.get_unifrac <- function(
        x, tree, weighted = FALSE, node.label = nodeLab, nodeLab = NULL, ...){
    # Transpose the matrix so that the orientation is the same as in other
    # dissimilatity methods
    x <- t(x)
    # Check x
    if( !is.matrix(as.matrix(x)) ){
        stop("'x' must be a matrix", call. = FALSE)
    }
    # Check tree
    if( is.null(tree) ){
        stop("Please provide 'tree'.", call. = FALSE)
    }
    # input check
    if(!.is_a_bool(weighted)){
        stop("'weighted' must be TRUE or FALSE.", call. = FALSE)
    }
    if(is.null(colnames(x)) || is.null(rownames(x))){
        stop("colnames and rownames must not be NULL", call. = FALSE)
    }
    # node.label should be NULL or character vector specifying links between 
    # rows and tree labels
    if( !(is.null(node.label) ||
            (is.character(node.label) && length(node.label) == nrow(x) &&
            all(node.label[ !is.na(node.label) ] %in% tree$tip.label))) ){
        stop(
            "'node.label' must be NULL or character specifying links between ",
            "abundance table and tree labels.", call. = FALSE)
    }
    # check that matrix and tree are compatible
    if( is.null(node.label) && !all(rownames(x) %in% c(tree$tip.label)) ) {
        stop(
            "Incompatible tree and abundance table! Please try to provide ",
            "'node.label'.", call. = FALSE)
    }
    # Merge rows, so that rows that are assigned to same tree node are
    # agglomerated
    # together. If nodeLabs were provided, merge based on those. Otherwise merge
    # based on rownames
    if( is.null(node.label) ){
        node.label <- rownames(x)
    }
    # Prune tree if there are nodes that cannot be found from tips or if there
    # are tips that cannot be found from abundance matrix. It might be
    # that certain row is linked to internal node or that the tree has extra
    # tips that do not match with rows (e.g. after subsetting).
    if( any( !node.label %in% tree$tip.label ) ||
            any( !tree$tip.label %in% node.label) ){
        tree <- .prune_tree(tree, node.label, ...)
        warning("Pruning tree...", call. = FALSE)
    }
    # If node labels cannot be found from tips even after pruning, give error.
    # This kind of tree cannot be used in unifrac since it expects that every
    # row is linked to tips.
    if( any( !node.label %in% tree$tip.label ) ){
        stop(
            "Unifrac cannot be calculated since tree is not compatible. ",
            "Each row must be linked to tip of the tree.", call. = FALSE)
    }
    # Merge assay so that each row represent single tip. It might be that
    # multiple rows are linked to single tip.
    x <- .merge_assay_by_rows(x, node.label, ...)

    # Calculate unifrac. Use implementation from rbiom package
    res <- unifrac(x, tree = tree, weighted = weighted)
    return(res)
}

# Aggregate matrix based on nodeLabs. At the same time, rename rows based on
# node.label
# --> each row represent specific node of tree
#' @importFrom scuttle sumCountsAcrossFeatures
.merge_assay_by_rows <- function(x, node.label, average = FALSE, ...){
    if( !.is_a_bool(average) ){
        stop("'average' must be TRUE or FALSE.", call. = FALSE)
    }
    # Merge assay based on nodeLabs
    x <- sumCountsAcrossFeatures(
        x, ids = node.label, subset.row = NULL, subset.col = NULL,
        average = average)
    # Remove NAs from node.label
    node.label <- node.label[ !is.na(node.label) ]
    # Get the original order back
    x <- x[ node.label, ]
    return(x)
}
