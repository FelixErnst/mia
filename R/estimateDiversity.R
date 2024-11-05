#' Estimate (alpha) diversity measures
#'
#' Several functions for calculating (alpha) diversity indices, including 
#' the \code{vegan} package options and some others.
#'
#' The available indices include the \sQuote{Coverage}, 
#' \sQuote{Faith's phylogenetic diversity}, \sQuote{Fisher alpha},
#' \sQuote{Gini-Simpson}, 
#' \sQuote{Inverse Simpson}, \sQuote{log-modulo skewness}, and \sQuote{Shannon} 
#' indices. See details for more information and references.
#'
#' @param x a \code{\link{SummarizedExperiment}} object or
#'   \code{\link{TreeSummarizedExperiment}}. The latter is recommended for
#'   microbiome data sets and tree-based alpha diversity indices.
#' 
#' @param tree A phylogenetic tree that is used to calculate 'faith' index.
#'   If \code{x} is a \code{TreeSummarizedExperiment}, \code{rowTree(x)} is 
#'   used by default.
#'
#' @param assay.type the name of the assay used for
#'   calculation of the sample-wise estimates.
#'   
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'
#' @param index a \code{character} vector, specifying the diversity measures
#'   to be calculated.
#'
#' @param name a name for the column(s) of the colData the results should be
#'   stored in. By default this will use the original names of the calculated
#'   indices.
#'   
#' @param tree.name a single \code{character} value for specifying which
#'   rowTree will be used to calculate faith index. 
#'   (By default: \code{tree.name = "phylo"})
#' 
#' @param tree_name Deprecated. Use \code{tree.name} instead.
#'   
#' @param node.label NULL or a character vector specifying the links between rows and 
#'   node labels of \code{tree}. If a certain row is not linked with the tree, missing 
#'   instance should be noted as NA. When NULL, all the rownames should be found from
#'   the tree. (By default: \code{node.label = NULL})
#' 
#' @param node_lab Deprecated. Use \code{node.label} instead.
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'
#' @param ... optional arguments:
#' \itemize{
#'   \item threshold: A numeric value in the unit interval,
#'   determining the threshold for coverage index. By default,
#'   \code{threshold} is 0.9.
#'   \item quantile: Arithmetic abundance classes are evenly cut up to to
#'   this quantile of the data. The assumption is that abundances higher than
#'   this are not common, and they are classified in their own group.
#'   By default, \code{quantile} is 0.5.
#'   \item nclasses: The number of arithmetic abundance classes
#'   from zero to the quantile cutoff indicated by \code{quantile}. 
#'   By default, \code{nclasses} is 50.
#'   \item num_of_classes Deprecated. Use \code{nclasses} instead.
#'   \item only.tips: A boolean value specifying whether to remove internal
#'   nodes when Faith's index is calculated. When \code{only.tips=TRUE}, those
#'   rows that are not tips of tree are removed.
#'   (By default: \code{only.tips=FALSE})
#' }
#'
#' @return \code{x} with additional \code{\link{colData}} named \code{*name*}
#'
#' @details
#'
#' Alpha diversity is a joint quantity that combines elements or community
#' richness and evenness. Diversity increases, in general, when species richness
#' or evenness increase.
#'
#' By default, this function returns all indices.
#'
#' \itemize{
#' 
#' \item 'coverage': Number of species needed to cover a given fraction of
#' the ecosystem (50 percent by default). Tune this with the threshold
#' argument.
#' 
#' \item 'faith': Faith's phylogenetic alpha diversity index measures how
#' long the taxonomic distance is between taxa that are present in the sample.
#' Larger values represent higher diversity. Using this index requires
#' rowTree. (Faith 1992)
#' 
#' If the data includes features that are not in tree's tips but in
#' internal nodes, there are two options. First, you can keep those features,
#' and prune the tree to match features so that each tip can be found from
#' the features. Other option is to remove all features that are not tips.
#' (See \code{only.tips} parameter)
#' 
#' \item 'fisher': Fisher's alpha; as implemented in
#' \code{\link[vegan:diversity]{vegan::fisher.alpha}}. (Fisher et al. 1943)
#' 
#' \item 'gini_simpson': Gini-Simpson diversity i.e. \eqn{1 - lambda},
#' where \eqn{lambda} is the
#' Simpson index, calculated as the sum of squared relative abundances.
#' This corresponds to the diversity index
#' 'simpson' in \code{\link[vegan:diversity]{vegan::diversity}}.
#' This is also called Gibbs–Martin, or Blau index in sociology,
#' psychology and management studies. The Gini-Simpson index (1-lambda)
#' should not be
#' confused with Simpson's dominance (lambda), Gini index, or
#' inverse Simpson index (1/lambda).
#' 
#' \item 'inverse_simpson': Inverse Simpson diversity:
#' \eqn{1/lambda} where \eqn{lambda=sum(p^2)} and p refers to relative
#' abundances.
#' This corresponds to the diversity index
#' 'invsimpson' in vegan::diversity. Don't confuse this with the
#' closely related Gini-Simpson index
#'
#' \item 'log_modulo_skewness': The rarity index characterizes the
#' concentration of species at low abundance. Here, we use the skewness of
#' the frequency 
#' distribution of arithmetic abundance classes (see Magurran & McGill 2011).
#' These are typically right-skewed; to avoid taking log of occasional
#' negative skews, we follow Locey & Lennon (2016) and use the log-modulo
#' transformation that adds a value of one to each measure of skewness to
#' allow logarithmization.
#'
#' \item 'shannon': Shannon diversity (entropy).
#' 
#' }
#'
#' @references
#'
#' Beisel J-N. et al. (2003)
#' A Comparative Analysis of Diversity Index Sensitivity.
#' _Internal Rev. Hydrobiol._ 88(1):3-15.
#' \url{https://portais.ufg.br/up/202/o/2003-comparative_evennes_index.pdf}
#'
#' Bulla L. (1994)
#' An  index of diversity and its associated diversity measure.
#' _Oikos_ 70:167--171
#'
#' Faith D.P. (1992)
#' Conservation evaluation and phylogenetic diversity.
#' _Biological Conservation_ 61(1):1-10.
#'
#' Fisher R.A., Corbet, A.S. & Williams, C.B. (1943)
#' The relation between the number of species and the number of individuals in
#' a random sample of animal population.
#' _Journal of Animal Ecology_ *12*, 42-58.
#' 
#' Locey K.J. & Lennon J.T. (2016)
#' Scaling laws predict global microbial diversity.
#' _PNAS_ 113(21):5970-5975.
#'
#' Magurran A.E., McGill BJ, eds (2011)
#' Biological Diversity: Frontiers in Measurement and Assessment.
#' (Oxford Univ Press, Oxford), Vol 12.
#'
#' Smith B. & Wilson JB. (1996)
#' A Consumer's Guide to Diversity Indices.
#' _Oikos_ 76(1):70-82.
#'
#' @seealso
#' \code{\link[scater:plotColData]{plotColData}}
#' \itemize{
#'   \item \code{\link[=estimateRichness]{estimateRichness}}
#'   \item \code{\link[=estimateEvenness]{estimateEvenness}}
#'   \item \code{\link[=estimateDominance]{estimateDominance}}
#'   \item \code{\link[vegan:diversity]{diversity}}
#'   \item \code{\link[vegan:specpool]{estimateR}}
#' }
#'
#' @name .estimate_diversity
#' @noRd
#' 
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # All index names as known by the function
#' index <- c("shannon","gini_simpson","inverse_simpson", "coverage", "fisher", 
#' "faith",  "log_modulo_skewness")
#'
#' # Corresponding polished names
#' name <- c("Shannon","GiniSimpson","InverseSimpson", "Coverage", "Fisher", 
#' "Faith",  "LogModSkewness")
#'
#' # Calculate diversities
#' tse <- .estimate_diversity(tse, index = index)
#'
#' # The colData contains the indices with their code names by default
#' colData(tse)[, index]
#'
#' # Removing indices
#' colData(tse)[, index] <- NULL
#' 
#' # 'threshold' can be used to determine threshold for 'coverage' index
#' tse <- .estimate_diversity(tse, index = "coverage", threshold = 0.75)
#' # 'quantile' and 'nclasses' can be used when
#' # 'log_modulo_skewness' is calculated
#' tse <- .estimate_diversity(tse, index = "log_modulo_skewness",
#'        quantile = 0.75, nclasses = 100)
#'
#' # It is recommended to specify also the final names used in the output.
#' tse <- .estimate_diversity(tse,
#'         index = c("shannon", "gini_simpson", "inverse_simpson", "coverage",
#'                    "fisher", "faith", "log_modulo_skewness"),
#'         name = c("Shannon", "GiniSimpson",  "InverseSimpson",  "Coverage",
#'                    "Fisher", "Faith", "LogModSkewness"))
#' # The colData contains the indices by their new names provided by the user
#' colData(tse)[, name]
#'
#' # Compare the indices visually
#' pairs(colData(tse)[, name])
#'
#' # Plotting the diversities - use the selected names
#' library(scater)
#' plotColData(tse, "Shannon")
#' # ... by sample type
#' plotColData(tse, "Shannon", "SampleType")
#' \donttest{
#' # combining different plots
#' library(patchwork)
#' plot_index <- c("Shannon","GiniSimpson")
#' plots <- lapply(plot_index,
#'                plotColData,
#'                object = tse,
#'                x = "SampleType",
#'                colour_by = "SampleType")
#' plots <- lapply(plots,"+",
#'    theme(axis.text.x = element_text(angle=45,hjust=1)))
#' names(plots) <- plot_index
#' plots$Shannon + plots$GiniSimpson + plot_layout(guides = "collect")
#' }
NULL

.estimate_diversity <- function(
        x, assay.type = assay_name, assay_name = "counts",
        index = c(
            "coverage", "faith", "fisher", "gini_simpson", "inverse_simpson",
            "log_modulo_skewness", "shannon"),
        name = index, BPPARAM = SerialParam(), ...){
    # Check index
    supported_index <- c(
        "coverage", "faith", "fisher", "gini_simpson", "inverse_simpson",
        "log_modulo_skewness", "shannon")
    index_string <- paste0(
        "'", paste0(supported_index, collapse = "', '"), "'")
    if ( !all(index %in% supported_index) || !(length(index) > 0)) {
        stop("'", paste0(
            "'index' must be from the following options: '", index_string),
            "'", call. = FALSE)
    }
    # Check name
    if(!.is_non_empty_character(name) || length(name) != length(index)){
        stop("'name' must be a non-empty character value and have the ",
            "same length as 'index'.",
            call. = FALSE)
    }
    # Check assay and check that vegan package is available since it is
    # required in some of the diversity calculations.
    .check_assay_present(assay.type, x)
    #
    # Calculate specified diversity indices
    dvrsts <- BiocParallel::bplapply(
        index,
        .get_diversity_values,
        x = x,
        mat = assay(x, assay.type),
        BPPARAM = BPPARAM,
        ...)
    # Add them to colData
    x <- .add_values_to_colData(x, dvrsts, name)
    return(x)
}

################################################################################

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
