
.estimate_faith <- function(
        x, mat, tree = NULL, tree.name = "phylo", node.label = node_lab,
        node_lab = NULL, ...){
    # Input check
    # If tree is NULL, then we must have TreeSE object that has rowTree
    if( is.null(tree) &&
        !(is(x, "TreeSummarizedExperiment") && !is.null(rowTree(x))) ){
        stop("'tree' must be provided.", call. = FALSE)
    }
    # If tree is not specified, then we get rowTree
    if( is.null(tree) ){
        .check_rowTree_present(tree.name, x)
        tree <- rowTree(x, tree.name)
        # When we get rowTree, we know the linking between rows and nodes.
        node.label <- rowLinks(x)[ , "nodeLab" ]
        node.label[ rowLinks(x)[, "whichTree"] != tree.name ] <- NA
    }
    # Check that tree is correct
    if( is.null(tree) || is.null(tree$edge.length) ){
        stop(
            "'tree' is NULL or it does not have any branches. ",
            "The Faith's alpha diversity index is not possible to ",
            "calculate.", call. = FALSE)
    }

    # Check that node.label is NULL or it specifies links between rownames and
    # node labs
    if( !( is.null(node.label) ||
            is.character(node.label) && length(node.label) == nrow(x) ) ){
        stop(
            "'node.label' must be NULL or a vector specifying links between ",
            "rownames and node labs of 'tree'.",
            call. = FALSE)
    }

    # Subset rows of the assay to correspond node_labs (if there are any NAs
    # in node labels)
    if( !is.null(node.label) && any(is.na(node.label)) ){
        warning(
            "The tree named does not include all the ",
            "rows. 'x' is subsetted.", call. = FALSE)
        mat <- mat[ !is.na(node.label), ]
        node.label <- node.label[ !is.na(node.label) ]
    }
    # If there are node labels (any TreeSE should have because they are rowLinks
    # by default), rename the features in matrix to match with labels found in
    # tree.
    if( !is.null(node.label) ){
        rownames(mat) <- node.label
    }

    # To calculate faith, the assay must have rownames. TreeSE has always
    # rownames at this point, but if the object is SE, it might be that it is
    # missing rownames.
    if( is.null(rownames(mat)) ){
        stop("'x' must have rownames.", call. = FALSE)
    }
    # Calculates Faith index
    res <- .calc_faith(mat, tree, ...)
    return(res)
}

.calc_shannon <- function(mat, ...){
    vegan::diversity(t(mat), index="shannon")
}

# NOTE: vegan::diversity(x, index = "simpson")
# gives Simpson diversity, also called Gini-Simpson
# index: 1-lambda, where lambda is the Simpson index
# (lambda). This may cause confusion if your familiarity
# with diversity indices is limited.
# Moreover, Simpson's lambda is simply the
# squared sum of relative abundances so we can
# just use that for clarity and simplicity.
#.get_simpson <- function(x, ...){
.simpson_lambda <- function(mat, ...){

    # Convert table to relative values
    rel <- .calc_rel_abund(mat)

    # Squared sum of relative abundances
    colSums2(rel^2)
}

.calc_gini_simpson <- function(mat, ...){
    1 - .simpson_lambda(mat, ...)
}

.calc_inverse_simpson <- function(mat, ...){
    1 / .simpson_lambda(mat, ...)
}

.calc_coverage <- function(mat, threshold = 0.9, ...){

    # Threshold must be a numeric value between 0-1
    if( !( is.numeric(threshold) && (threshold >= 0 && threshold <= 1) ) ){
        stop("'threshold' must be a numeric value between 0-1.",
            call. = FALSE)
    }

    # Convert table to relative values
    rel <- .calc_rel_abund(mat)

    # Number of groups needed to have threshold (e.g. 50 %) of the
    # ecosystem occupied
    coverage <- apply(rel, 2, function(x) {
        min(which(cumsum(rev(sort(x/sum(x)))) >= threshold))
    })
    names(coverage) <- colnames(rel)
    coverage
}

.calc_fisher <- function(mat, ...){
    vegan::fisher.alpha(t(mat))
}

.calc_faith <- function(mat, tree, only.tips = FALSE, ...){
    # Input check
    if( !.is_a_bool(only.tips) ){
        stop("'only.tips' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Remove internal nodes if specified
    if( only.tips ){
        mat <- mat[ rownames(mat) %in% tree$tip.label, ]
    }
    # To ensure that the function works with NA also, convert NAs to 0.
    # Zero means that the taxon is not present --> same as NA (no information)
    mat[ is.na(mat) ] <- 0

    # Gets vector where number represent nth sample
    samples <- seq_len(ncol(mat))

    # Repeats taxa as many times there are samples, i.e. get all the
    # taxa that are analyzed in each sample.
    taxa <- rep(rownames(mat), length(samples))

    # Gets those taxa that are present/absent in each sample.
    # Gets one big list that combines
    # taxa from all the samples.
    present_combined <- taxa[ mat[, samples] > 0 ]

    # Gets how many taxa there are in each sample.
    # After that, determines indices of samples' first taxa with cumsum.
    split_present <- as.vector(cumsum(colSums(mat > 0)))

    # Determines which taxa belongs to which sample by first determining
    # the splitting points,
    # and after that giving every taxa number which tells their sample.
    split_present <- as.factor(cumsum((seq_along(present_combined)-1) %in%
                        split_present))

    # Assigns taxa to right samples based on their number that they got from
    # previous step, and deletes unnecessary names.
    present <- unname(split(present_combined, split_present))

    # If there were samples without any taxa present/absent, the length of the
    # list is not the number of samples since these empty samples are missing.
    # Add empty samples as NULL.
    names(present) <- names(which(colSums2(mat) > 0))
    present[names(which(colSums2(mat) == 0))] <- list(NULL)
    present <- present[colnames(mat)]

    # Assign NA to all samples
    faiths <- rep(NA,length(samples))

    # If there are no taxa present, then faith is 0
    ind <- lengths(present) == 0
    faiths[ind] <- 0

    # If there are taxa present
    ind <- lengths(present) > 0
    # Loop through taxa that were found from each sample
    faiths_for_taxa_present <- lapply(present[ind], function(x){
        # Trim the tree
        temp <- .prune_tree(tree, x, ...)
        # Sum up all the lengths of edges
        temp <- sum(temp$edge.length)
        return(temp)
    })
    faiths_for_taxa_present <- unlist(faiths_for_taxa_present)
    faiths[ind] <- faiths_for_taxa_present
    return(faiths)
}

.calc_log_modulo_skewness <- function(mat, quantile = 0.5,
    nclasses = num_of_classes, num_of_classes = 50, ...){
    # quantile must be a numeric value between 0-1
    if( !( is.numeric(quantile) && (quantile >= 0 && quantile <= 1) ) ){
        stop("'quantile' must be a numeric value between 0-1.",
            call. = FALSE)
    }
    # nclasses must be a positive numeric value
    if( !( is.numeric(nclasses) && nclasses > 0 ) ){
        stop("'nclasses' must be a positive numeric value.",
            call. = FALSE)
    }
    # Determine the quantile point.
    quantile_point <- quantile(max(mat), quantile)
    # Tabulate the arithmetic abundance classes. Use the same classes
    # for all samples for consistency
    cutpoints <- c(seq(0, quantile_point, length=nclasses), Inf)
    # Calculates sample-wise frequencies. How many taxa in each interval?
    freq_table <- table(cut(mat, cutpoints), col(mat))
    # Calculates the skewness of frequency table. Returns skewness for each
    # sample
    r <- .calc_skewness(freq_table)
    # Return log-modulo
    log(1 + r)
}

#' @importFrom DelayedMatrixStats rowSums2 rowMeans2
.calc_skewness <- function(x) {
    # Transposes the table
    x <- t(x)
    # Each value is subtracted by sample-wise mean, which is raised to the
    # power of 3.
    # Then the sample-wise sum is taken from these values.
    numerator <- rowSums2((x - rowMeans2(x))^3)
    # Sample-wise sum is divided by number of taxa that are not NA.
    numerator <- numerator/rowSums2(!is.na(x))
    # Each value is subtracted by sample-wise mean, which is raises to the
    # power of 2.
    # Then the sample-wise sum is taken from these values.
    denominator <- rowSums2((x - rowMeans2(x))^2)
    # Sample-wise sum is divided by number of taxa that are not NA. Then
    # these values
    # are raised to the power of 3/2.
    denominator <- (denominator/rowSums2(!is.na(x)))^(3/2)
    # Result
    result <- numerator/denominator
    return(result)
}

#' @importFrom SummarizedExperiment assay assays
.get_diversity_values <- function(index, x, mat, ...){
    FUN <- switch(index,
        shannon = .calc_shannon,
        gini_simpson = .calc_gini_simpson,
        inverse_simpson = .calc_inverse_simpson,
        coverage = .calc_coverage,
        fisher = .calc_fisher,
        faith = .estimate_faith,
        log_modulo_skewness = .calc_log_modulo_skewness
        )
    res <- FUN(x = x, mat = mat, ...)
    res <- unname(res)
    return(res)
}
