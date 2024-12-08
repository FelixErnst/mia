#' Create a \code{TreeSummarizedExperiment} object from a phyloseq object
#' 
#' @inheritParams convertFromBIOM
#' 
#' @details 
#' \code{convertFromPhyloseq} converts \code{phyloseq}
#' objects into 
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}} objects.
#' All data stored in a \code{phyloseq} object is transferred.
#'
#' @return 
#' \code{convertFromPhyloseq} returns an object of class 
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @importFrom S4Vectors SimpleList DataFrame make_zero_col_DFrame
#' @importFrom SummarizedExperiment colData colData<-
#'
#' @export
#'
#' @rdname convertFromPhyloseq
#'
#' @examples
#' 
#' ### Coerce a phyloseq object to a TreeSE object
#' if (requireNamespace("phyloseq")) {
#'     data(GlobalPatterns, package="phyloseq")
#'     convertFromPhyloseq(GlobalPatterns)
#'     data(enterotype, package="phyloseq")
#'     convertFromPhyloseq(enterotype)
#'     data(esophagus, package="phyloseq")
#'     convertFromPhyloseq(esophagus)
#' }
convertFromPhyloseq <- function(x) {
    # input check
    .require_package("phyloseq")
    if(!is(x,"phyloseq")){
        stop("'x' must be a 'phyloseq' object")
    }
    #
    # Get the assay
    counts <- phyloseq::otu_table(x)
    # Check the orientation, and transpose if necessary
    if( !phyloseq::taxa_are_rows(x) ){
        counts <- t(counts)
    }
    # Create a list of assays
    assays <- SimpleList(counts = counts)
    
    rowData <- tryCatch(phyloseq::tax_table(x), error = function(e) NULL) |>
        data.frame() |> DataFrame()
    if( nrow(rowData) == 0L ){
        rowData <- S4Vectors::make_zero_col_DFrame(nrow(assays$counts))
        rownames(rowData) <- rownames(assays$counts)
    }
    colData <- tryCatch(phyloseq::sample_data(x), error = function(e) NULL) |>
        data.frame() |> DataFrame()
    if( nrow(colData) == 0L ){
        colData <- S4Vectors::make_zero_col_DFrame(ncol(assays$counts))
        rownames(colData) <- colnames(assays$counts)
    }
    rowTree <- tryCatch(phyloseq::phy_tree(x), error = function(e) NULL)
    referenceSeq <- tryCatch(phyloseq::refseq(x), error = function(e) NULL)
    TreeSummarizedExperiment(assays = assays,
                            rowData = rowData,
                            colData = colData,
                            rowTree = rowTree,
                            referenceSeq = referenceSeq)
}
