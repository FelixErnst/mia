#' Clustering wrapper
#'
#' This function returns a \code{SummarizedExperiment} with clustering 
#'   information in its colData or rowData
#'
#' @param x A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#' 
#' @param by \code{Character scalar}. Determines if association is calculated
#'   row-wise / for features ('rows') or column-wise / for samples ('cols').
#'   Must be \code{'rows'} or \code{'cols'}.
#'   
#' @param MARGIN Deprecated. Use \code{by} instead.
#'   
#' @param clust.col \code{Character scalar}. Indicates the name of the 
#'   \code{rowData} (or \code{colData}) where the data will be stored.
#'   (Default: \code{"clusters"})
#'   
#' @param ... Additional parameters to use altExps for example
#' @inheritParams bluster::clusterRows
#' @inheritParams runDMN
#' @inheritParams transformAssay
#'   
#'   
#' @details
#' This is a wrapper for the \code{clusterRows} function from the 
#' \link[https://bioconductor.org/packages/release/bioc/html/bluster.html]{bluster} 
#' package.
#'
#' When setting \code{full = TRUE}, the clustering information will be stored in
#' the metadata of the object.
#' 
#' By default, clustering is done on the features.
#'
#' @return
#' \code{addCluster} returns an object of the same type as the \code{x}
#' parameter with clustering information named \code{clusters} stored in
#' \code{colData} or \code{rowData}. 
#'
#' @name addCluster
#' @export
#'
#' @examples
#' library(bluster)
#' data(GlobalPatterns, package = "mia")
#' tse <- GlobalPatterns
#'
#' # Cluster on rows using Kmeans
#' tse <- addCluster(tse, KmeansParam(centers = 3))
#' 
#' # Clustering done on the samples using Hclust
#' tse <- addCluster(tse, 
#'                by = "samples", 
#'                HclustParam(metric = "bray", dist.fun = vegan::vegdist))
#' 
#' # Getting the clusters
#' colData(tse)$clusters
#' 
NULL

#' @rdname addCluster
#' @export
#' @importFrom bluster clusterRows
setMethod("addCluster", signature = c(x = "SummarizedExperiment"),
    function(
            x, BLUSPARAM, assay.type = assay_name, 
            assay_name = "counts", by = MARGIN, MARGIN = "rows", full = FALSE, 
            name = "clusters", clust.col = "clusters", ...) {
        .require_package("bluster")
        # Checking parameters
        by <- .check_MARGIN(by)
        se <- .check_and_get_altExp(x, ...)
        .check_assay_present(assay.type, se)
        if( !.is_a_string(name) ){
            stop("'name' must be a non-empty single character value.",
                call. = FALSE)
        }
        if( !.is_a_string(clust.col) ){
            stop("'clust.col' must be a non-empty single character value.",
                call. = FALSE)
        }
        if( !.is_a_bool(full) ){
            stop("'full' must be TRUE or FALSE.", call. = FALSE)
        }
        #
        # Get assay
        mat <- assay(se, assay.type)
        # Transpose if clustering on the columns
        if(by == 2){
            mat <- t(mat)
        }
        # Get clusters
        result <- clusterRows(mat, BLUSPARAM, full)
        # If user has specified full=TRUE, result includes additional info
        # that will be stored to metadata.
        if (full) {
            clusters <- result$clusters
            x <- .add_values_to_metadata(x, name, result$objects, ...)
        } else {
            clusters <- result
        }
        # Setting clusters in the object. The adding function requires data as
        # list
        clusters <- list(clusters)
        x <- .add_values_to_colData(
            x, clusters, clust.col, MARGIN = by, colname = "clust.col", ...)
        return(x)
    }
)
