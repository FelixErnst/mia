#' Calculate dissimilarities
#'
#' These functions are designed to calculate dissimilarities on data stored 
#' within a 
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object. For overlap, Unifrac, and Jensen-Shannon Divergence (JSD) 
#' dissimilarities, the functions use mia internal functions, while for other 
#' types of dissimilarities, they rely on \code{\link[vegan:vegdist]{vegdist}}
#' by default.
#'
#' @param x
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' or \code{matrix}.
#'
#' @param method \code{Character scalar}. Specifies which dissimilarity to 
#' calculate. (Default: \code{"bray"})
#' 
#' @param name \code{Character scalar}. The name to be used to store the result 
#' in metadata of the output. (Default: \code{method})
#'
#' @param assay.type \code{Character scalar}. Specifies which assay to use for 
#' calculation. (Default: \code{"counts"})
#' 
#' @param niter The number of iterations performed. If \code{NULL},
#' rarefaction is disabled. (Default: \code{NULL})
#'   
#' @param transposed \code{Logical scalar}. Specifies if x is transposed with
#' cells in rows. (Default: \code{FALSE})
#'
#' @param ... other arguments passed into \code{\link[vegan:avgdist]{avgdist}},
#' \code{\link[vegan:vegdist]{vegdist}}, or into mia internal functions:
#' 
#' \itemize{
#'   \item \code{sample}: The sampling depth in rarefaction.
#'   (Default: \code{min(rowSums2(x))})
#'   
#'   \item \code{dis.fun}: \code{Character scalar}. Specifies the dissimilarity
#'   function to be used.
#'   
#'   \item \code{tree.name}: (Unifrac)  \code{Character scalar}. Specifies the
#'   name of the tree from \code{rowTree(x)} that is used in calculation.
#'   Disabled when \code{tree} is specified. (Default: \code{"phylo"})
#'   
#'   \item \code{tree}: (Unifrac) \code{phylo}. A phylogenetic tree used in
#'   calculation. (Default: \code{NULL})
#'   
#'   \item \code{weighted}: (Unifrac) \code{Logical scalar}. Should use
#'   weighted-Unifrac calculation? 
#'   Weighted-Unifrac takes into account the relative abundance of
#'   species/taxa shared between samples, whereas unweighted-Unifrac only
#'   considers presence/absence. Default is \code{FALSE}, meaning the
#'   unweighted-Unifrac dissimilarity is calculated for all pairs of samples.
#'   (Default: \code{FALSE})
#'   
#'   \item \code{node.label} (Unifrac) \code{character vector}. Used only if
#'   \code{x} is a matrix. Specifies links between rows/columns and tips of 
#'   \code{tree}. The length must equal the number of rows/columns of \code{x}. 
#'   Furthermore, all the node labs must be present in \code{tree}.
#'   
#'   \item \code{chunkSize}: (JSD) \code{Integer scalar}. Defines the size of
#'   data  send to the individual worker. Only has an effect, if \code{BPPARAM} 
#'   defines more than one worker. (Default: \code{nrow(x)})
#'   
#'   \item \code{BPPARAM}: (JSD)
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}.
#'   Specifies whether the calculation should be parallelized.
#'   
#'   \item \code{detection}: (Overlap) \code{Numeric scalar}.
#'   Defines detection threshold for absence/presence of features. Feature that 
#'   has abundance under threshold in either of samples, will be discarded when 
#'   evaluating overlap between samples. (Default: \code{0}) 
#' }
#'
#' @return 
#' \code{getDissimilarity} returns a sample-by-sample dissimilarity matrix.
#' 
#' \code{addDissimilarity} returns \code{x} that includes dissimilarity matrix 
#' in its metadata. 
#'   
#' @details 
#' Overlap reflects similarity between sample-pairs. When overlap is 
#' calculated using relative abundances, the higher the value the higher the 
#' similarity is. When using relative abundances, overlap value 1 means that 
#' all the abundances of features are equal between two samples, and 0 means 
#' that samples have completely different relative abundances. 
#'   
#' Unifrac is calculated with \code{\link[rbiom:unifrac]{rbiom:unifrac()}}.
#' 
#' If rarefaction is enabled, \code{\link[vegan:avgdist]{vegan:avgdist()}} is
#' utilized.
#' 
#' For JSD implementation:
#' Susan Holmes \email{susan@@stat.stanford.edu}.
#' Adapted for phyloseq by Paul J. McMurdie.
#' Adapted for mia by Felix G.M. Ernst
#'   
#' @name getDissimilarity
#' 
#' @seealso
#' \url{http://en.wikipedia.org/wiki/Jensen-Shannon_divergence}
#' 
#' @references
#' For unifrac dissimilarity: \url{http://bmf.colorado.edu/unifrac/}
#' 
#' See also additional descriptions of Unifrac in the following articles:
#'
#' Lozupone, Hamady and Knight, ``Unifrac - An Online Tool for Comparing
#' Microbial Community Diversity in a Phylogenetic Context.'', BMC
#' Bioinformatics 2006, 7:371
#'
#' Lozupone, Hamady, Kelley and Knight, ``Quantitative and qualitative (beta)
#' diversity measures lead to different insights into factors that structure
#' microbial communities.'' Appl Environ Microbiol. 2007
#'
#' Lozupone C, Knight R. ``Unifrac: a new phylogenetic method for comparing
#' microbial communities.'' Appl Environ Microbiol. 2005 71 (12):8228-35.
#' 
#' For JSD dissimilarity: 
#' Jensen-Shannon Divergence and Hilbert space embedding.
#' Bent Fuglede and Flemming Topsoe University of Copenhagen,
#' Department of Mathematics
#' \url{http://www.math.ku.dk/~topsoe/ISIT2004JSD.pdf}
#'
#' @export
#'
#' @examples
#' library(mia)
#' library(scater)
#' 
#' # load dataset
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' 
#' ### Overlap dissimilarity
#' 
#' tse <- addDissimilarity(tse, method = "overlap", detection = 0.25)
#' metadata(tse)[["overlap"]][1:6, 1:6]
#' 
#' ### JSD dissimilarity
#' 
#' tse <- addDissimilarity(tse, method = "jsd")
#' metadata(tse)[["jsd"]][1:6, 1:6]
#' 
#' # Multi Dimensional Scaling applied to JSD dissimilarity matrix
#' tse <- runMDS(tse, FUN = getDissimilarity, method = "overlap", 
#'               assay.type = "counts")
#' metadata(tse)[["MDS"]][1:6, ]
#'               
#' ### Unifrac dissimilarity
#' 
#' res <- getDissimilarity(tse, method = "unifrac", weighted = FALSE)
#' dim(as.matrix((res)))
#' 
#' tse <- addDissimilarity(tse, method = "unifrac", weighted = TRUE)
#' metadata(tse)[["unifrac"]][1:6, 1:6]
#' 
NULL

#' @rdname getDissimilarity
#' @export
setMethod(
    "addDissimilarity", signature = c(x = "SummarizedExperiment"),
    function(x, method = "bray", name = method, ...){
    #
    res <- getDissimilarity(x, method = method, ...)
    # Add matrix to original SE
    x <- .add_values_to_metadata(x, names = name, values = as.matrix(res))
    return(x)
    }
)

#' @rdname getDissimilarity
#' @export
setMethod(
    "getDissimilarity", signature = c(x = "SummarizedExperiment"),
    function(
        x, method = "bray", assay.type = "counts", niter = NULL,
        transposed = FALSE, ...){
    # Input checks
    .check_assay_present(assay.type, x)
    if( !.is_non_empty_string(method) ){
        stop("'method' must be a non-empty single character value",
            call. = FALSE)
    }
    if( !.is_a_bool(transposed) ){
        stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Get arguments
    mat <- assay(x, assay.type)
    if( !transposed ){
        mat <- t(mat)
    }
    args <- c(
        list(x = mat, method = method, niter = niter), list(...))
    # Calculate dissimilarity based on matrix
    mat <- do.call(getDissimilarity, args)
    return(mat)
    }
)

#' @rdname getDissimilarity
#' @export
setMethod(
    "getDissimilarity", signature = c(x = "TreeSummarizedExperiment"),
    function(
        x, method = "bray", assay.type = "counts", niter = NULL,
        transposed = FALSE, ...){
    # Input checks
    .check_assay_present(assay.type, x)
    if( !.is_non_empty_string(method) ){
        stop("'method' must be a non-empty single character value",
            call. = FALSE)
    }
    if( !.is_a_bool(transposed) ){
        stop("'transposed' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Retrieve tree arguments from TreeSE object, if method is unifrac
    if( method %in% c("unifrac") ){
        args <- .get_tree_args(
            x,  method = method, assay.type = assay.type,
            transposed = transposed, ...)
    } else{
        # For other cases, do not fetch tree data from TreeSE
        mat <- assay(x, assay.type)
        if( !transposed ){
            mat <- t(mat)
        }
        args <- c(
            list(x = mat, method = method, niter = niter), list(...))
    }
    # Calculate dissimilarity
    mat <- do.call(getDissimilarity, args)
    return(mat)
    }
)

#' @rdname getDissimilarity
#' @export
setMethod(
    "getDissimilarity", signature = c(x = "ANY"), function(
        x, method = "bray", niter = NULL, ...){
    # Input check
    if( !.is_a_string(method) ){
        stop("'method' must be a single character value.", call. = FALSE)
    }
    #
    # Calculate dissimilarity
    mat <- .calculate_dissimilarity(
        mat = x, method = method, niter = niter, ...)
    return(mat)
    }
)

# This function chooses right method and calculates dissimilarity matrix.
#' @importFrom MatrixGenerics rowSums2
#' @importFrom vegan vegdist avgdist
.calculate_dissimilarity <- function(
        mat, method, niter, dis.fun = distfun, distfun = FUN, FUN = NULL,
        sample = min(rowSums2(mat)), ...){
    # input check
    if( !(is.null(dis.fun) || is.function(dis.fun)) ){
        stop("'dis.fun' must be NULL or a function.", call. = FALSE)
    }
    if( !(is.null(niter) || .is_an_integer(niter)) ){
        stop("'niter' must be NULL or an integer.", call. = FALSE)
    }
    # sample is only used when niter is specified
    if( !is.null(niter) && !.is_an_integer(sample) ){
        stop("'sample' must be an integer.", call. = FALSE)
    } 
    #
    # If the dissimilarity function is not specified, get default choice
    if( is.null(dis.fun) ){
        if( method %in% c("overlap") ){
            dis.fun <- .get_overlap
        } else if( method %in% c("unifrac")  ){
            dis.fun <- .get_unifrac
        } else if( method %in% c("jsd")  ){
            dis.fun <- .get_jsd
        } else{
            dis.fun <- vegdist
        }
    }
    # Initialize an argument list
    args <- c(list(mat), list(...))
    # If rarefaction is specified, calculate dissimilarity with vegan::avgdist
    # function that utilizes the specified dissimilarity function. Otherwise,
    # call the specified function directly.
    if( !is.null(niter) ){
        # Remove arguments that will overlap with arguments added below
        args <- args[ !names(args) %in% c("dmethod", "iterations") ]
        # Add arguments specific for avgdist
        args <- c(args, list(
            dmethod = method, iterations = niter, sample = sample,
            distfun = dis.fun))
        # Calculate dissimilarities
        res <- do.call(avgdist, args)
    } else{
        args <- c(args, list(method = method))
        res <- do.call(dis.fun, args)
    }
    return(res)
}

# If user want to calculate unifrac dissimilarity, this function gathers tree
# data.
.get_tree_args <- function(
        x, method, assay.type = "counts", transposed = FALSE, tree = NULL, ...){
    # Check tree.
    if( !(is.null(tree) || is(tree, "phylo")) ){
        stop("'tree' must be NULL or phylo.", call. = FALSE)
    }
    #
    # Create an argument list that includes matrix, and tree-related parameters.
    args <- list(method = method)
    args <- c(args, list(...))
    # Either add tree that was provided by user, or get tree from TreeSE
    if( !is.null(tree) ){
        mat <- assay(x, assay.type)
        if( !transposed ){
            mat <- t(mat)
        }
        tree_args <- list(x = mat, tree = tree)
    } else{
        tree_args <- .get_tree_args_from_TreeSE(x, transposed = transposed, 
            assay.type = assay.type, ...)
    }
    args <- c(args, tree_args)
    return(args)
}

# This function fetches tree arguments fro TreeSE slots.
.get_tree_args_from_TreeSE <- function(
        x, transposed, tree.name = "phylo", assay.type, ...){
    # Get functions and parameters based on direction
    tree_present_FUN <- if (transposed) .check_colTree_present
        else .check_rowTree_present
    tree_FUN <- if (transposed) colTree else rowTree
    links_FUN <- if (transposed) colLinks else rowLinks
    margin_name <- if (transposed) "col" else "row"
    # Check tree.name
    tree_present_FUN(tree.name, x)
    #
    # Select only those features/samples that are in the tree
    links <- links_FUN(x)
    present_in_tree <- links[, "whichTree"] == tree.name
    if( any(!present_in_tree) ){
        warning(
            "Not all ", margin_name, "s were present in the ", margin_name,
            "Tree specified by 'tree.name'. 'x' is subsetted.",
            call. = FALSE)
        # Subset the data
        if( transposed ){
            x <- x[, present_in_tree]
        } else{
            x <- x[present_in_tree, ]
        }
    }
    # Get tree
    tree <- tree_FUN(x, tree.name)
    # Get links and take only nodeLabs
    links <- links_FUN(x)
    links <- links[ , "nodeLab"]
    node.label <- links
    
    # Get assay. By default, dissimilarity between samples is calculated. In
    # dissimilarity functions, features must be in columns and samples in rows
    # in this case.
    mat <- assay(x, assay.type)
    if( !transposed ){
        mat <- t(mat)
    }
    # Return a list of tree arguments
    args <- list(x = mat, tree = tree, node.label = node.label)
    return(args)
}
