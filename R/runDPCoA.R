#' Calculation of Double Principal Correspondance analysis
#'
#' Double Principal Correspondance analysis is made available via the
#' \code{ade4} package in typical fashion. Results are stored in the
#' \code{reducedDims} and are available for all the expected functions.
#'
#' @inheritParams getDominant
#' @inheritParams getDissimilarity
#'
#' @details
#'   For \code{addDPCoA} a \linkS4class{TreeSummarizedExperiment} containing the
#'   expression values as well as a \code{rowTree} to calculate \code{y} using
#'   \code{\link[ape:cophenetic.phylo]{cophenetic.phylo}}.
#'
#' @param y a \code{dist} or a symmetric \code{matrix} compatible with
#'   \code{ade4:dpcoa}
#'
#' @param ncomponents \code{Numeric scalar}. Indicates the number of DPCoA
#' dimensions to obtain. (Default: \code{2})
#'
#' @param ntop \code{Numeric scalar}. Specifies the number of features with the
#' highest variances to use for dimensionality reduction. Alternatively
#' \code{NULL}, if all features should be used. (Default: \code{NULL})
#'
#' @param subset.row \code{Character Vector}. Specifies the subset of features
#' to use for dimensionality reduction. This can be a character vector of row
#' names, an integer vector of row indices or a logical vector.
#' (Default: \code{NULL})
#' 
#' @param subset_row Deprecated. Use \code{subset.row} instead.
#'
#' @param scale \code{Logical scalar}. Should the expression values be
#' standardized? (Default: \code{FALSE})
#' 
#' @param name \code{Character scalar}. A name for the column of the 
#'   \code{colData} where results will be stored. (Default: \code{"DPCoA"})
#' 
#' @param altexp \code{Character scalar} or \code{integer scalar}. Specifies an
#'   alternative experiment containing the input data. (Default: \code{NULL})
#'   
#' @param exprs_values Deprecated. Use \code{assay.type} instead.
#' 
#' @param tree.name \code{Character scalar}. Specifies the name of the
#'   tree to be included in the phyloseq object that is created, 
#'   (Default: \code{"phylo"})
#'   
#' @param tree_name Deprecated. Use \code{tree.name} instead.
#'
#' @param ... Currently not used.
#'
#' @details
#' In addition to the reduced dimension on the features, the reduced dimension
#' for samples are returned as well as \code{sample_red} attribute.
#' \code{eig}, \code{feature_weights} and \code{sample_weights} are
#' returned as attributes as well.
#'
#' @returns
#' For \code{getDPCoA} a matrix with samples as rows and CCA dimensions as
#' columns
#'
#' For \code{addDPCoA} a modified \code{x} with the results stored in
#' \code{reducedDim} as the given \code{name}
#'
#'
#' @name runDPCoA
#' @seealso
#' \code{\link[scater:plotReducedDim]{plotReducedDim}}
#' \code{\link[SingleCellExperiment:reducedDims]{reducedDims}}
#'
#' @examples
#' data(esophagus)
#' dpcoa <- getDPCoA(esophagus)
#' head(dpcoa)
#'
#' esophagus <- addDPCoA(esophagus)
#' reducedDims(esophagus)
#'
#' library(scater)
#' plotReducedDim(esophagus, "DPCoA")
NULL

.calculate_dpcoa <- function(
        x, y, ncomponents = 2, ntop = NULL,
        subset.row = subset_row, subset_row = NULL, scale = FALSE,
        transposed = FALSE, ...)
{
    .require_package("ade4")
    # input check
    # Check ncomponents
    if( !(.is_an_integer(ncomponents) && ncomponents > 0) ){
        stop("'ncomponents' must be a single integer value specifying the ",
            "number of DPCoA dimensions.", call. = FALSE)
    }
    # Check ntop
    if( !(is.null(ntop) || (.is_an_integer(ntop) && ntop > 0))  ){
        stop("'ntop' must be NULL or a single integer value specifying the ",
            "number of features with the highest variance.", call. = FALSE)
    }
    y <- as.matrix(y)
    if(length(unique(dim(y))) != 1L){
        stop("'y' must be symmetric.", call. = FALSE)
    }
    #
    # Get NAs. ade4:dpcoa lead to an error if there are any NAs
    if( any( is.na(x) ) ){
        stop("'x' includes NAs. Please try to convert them into numeric ",
            "values.", call. = FALSE)
    }
    if(!transposed) {
        if(is.null(ntop)){
            ntop <- nrow(x)
        }
        x <- .get_mat_for_reddim(
            x, subset_row = subset.row, ntop = ntop, scale = scale)
    }
    y <- y[colnames(x), colnames(x), drop = FALSE]
    if(nrow(y) != ncol(x)){
        stop("x and y must have corresponding dimensions.", call. = FALSE)
    }
    y <- sqrt(y)
    y <- as.dist(y)
    #
    dpcoa <- ade4::dpcoa(data.frame(x), y, scannf = FALSE, nf = ncomponents)
    ans <- as.matrix(dpcoa$li)
    rownames(ans) <- rownames(x)
    colnames(ans) <- NULL
    attr(ans,"eig") <- dpcoa$eig
    tmp <- as.matrix(dpcoa$dls)
    rownames(tmp) <- colnames(x)
    colnames(tmp) <- NULL
    attr(ans,"sample_red") <- tmp
    attr(ans,"feature_weights") <- unname(dpcoa$dw)
    attr(ans,"sample_weights") <- unname(dpcoa$lw)
    ans
}

#' @export
#' @rdname runDPCoA
setMethod("getDPCoA", c("ANY","ANY"), .calculate_dpcoa)

#' @export
#' @importFrom ape cophenetic.phylo
#' @rdname runDPCoA
setMethod("getDPCoA", signature = c("TreeSummarizedExperiment","missing"),
    function(x, ..., assay.type = assay_name, assay_name = exprs_values, 
        exprs_values = "counts", tree.name = tree_name, tree_name = "phylo")
    {
        .require_package("ade4")
        # Check assay.type
        .check_assay_present(assay.type, x)
        # Check tree.name
        .check_rowTree_present(tree.name, x)
        #
        # Get tree
        tree <- rowTree(x, tree.name)
        # Select only those features that are in the rowTree
        whichTree <- rowLinks(x)[ , "whichTree"] == tree.name
        if( any(!whichTree) ){
            warning("Not all rows were present in the rowTree specified by ",
                    "'tree.name'. 'x' is subsetted.", call. = FALSE)
            # Subset the data
            x <- x[ whichTree, ]
        }
        dist <- cophenetic.phylo(tree)
        # Get assay
        mat <- assay(x, assay.type)
        getDPCoA(mat, dist, ...)
    }
)

#' @export
#' @rdname runDPCoA
#' @aliases getDPCoA
calculateDPCoA <- function(x,...){
    getDPCoA(x,...)
}

#' @export
#' @rdname runDPCoA
#' @importFrom SingleCellExperiment reducedDim<-
addDPCoA <- function(x, ..., altexp = NULL, name = "DPCoA"){
    # Input check
    # Check and get altExp if altexp is not NULL
    if( !is.null(altexp) ){
        # Check altexp
        .check_altExp_present(altexp, x)
        # Get altExp
        y <- altExp(x, altexp)
    } else {
        y <- x
    }
    # Check name
    if( !.is_a_string(name) ){
        stop("'name' must be a single character value specifying a name of ",
            "reducedDim where the result will be stored.", call. = FALSE)
    }
    reducedDim(x, name) <- getDPCoA(y, ...)
    x
}

#' @export
#' @rdname runDPCoA
#' @aliases addDPCoA
runDPCoA <- function(x,...){
    addDPCoA(x,...)
}
