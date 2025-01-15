#' Create a phyloseq object from a \code{TreeSummarizedExperiment} object
#'  
#' @inheritParams convertFromBIOM
#' 
#' @param x a 
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'   object
#'
#' @param assay.type \code{Character scalar}. Specifies the name of assay 
#'   used. (Default: \code{"counts"})
#'   
#' @param assay_name Deprecated. Use \code{assay.type} instead.
#'   
#' @param tree.name \code{Character scalar}. Specifies the name of the
#'   tree to be included in the phyloseq object that is created, 
#'   (Default: \code{"phylo"})
#'   
#' @param tree_name Deprecated. Use \code{tree.name} instead.
#'
#' @details 
#' \code{convertToPhyloseq} creates a phyloseq object from a 
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object. By using \code{assay.type}, it is possible to specify which table
#' from \code{assay} is added to the phyloseq object.
#'
#' @return
#' \code{convertToPhyloseq} returns an object of class 
#' \code{\link[phyloseq:phyloseq-class]{phyloseq}}
#'
#' @name convertFromPhyloseq
#' @export
#'
#' @examples
#' 
#' ### Coerce a TreeSE object to a phyloseq object
#' # Get tse object
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#'
#' # Create a phyloseq object from it
#' phy <- convertToPhyloseq(tse)
#' phy
#'
#' # By default the chosen table is counts, but if there are other tables,
#' # they can be chosen with assay.type.
#'
#' # Counts relative abundances table
#' tse <- transformAssay(tse, method = "relabundance")
#' phy2 <- convertToPhyloseq(tse, assay.type = "relabundance")
#' phy2
#' 
#' @export
NULL

#' @rdname convertFromPhyloseq
#' @export
setMethod("convertToPhyloseq", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts", assay_name = NULL, ...){
        # Input check
        .require_package("phyloseq")
        # Check that tse do not have zero rows
        if(!all(dim(x) > 0)){
            stop("'x' contains zero rows. 'x' can not be converted
                to a phyloseq object.", call. = FALSE)
        }

        if (!is.null(assay_name)) {
            .Deprecated(old="assay_name", new="assay.type",
                "Now assay_name is deprecated. Use assay.type instead.")
        }
        
        # Check assay.type
        .check_assay_present(assay.type, x)
        
        # phyloseq object requires nonduplicated rownames. If there are 
        # duplicated rownames, they are converted so that they are unique
        if( any(duplicated(rownames(x))) ){
            rownames(x) <- getTaxonomyLabels(x)
        }
        # List of arguments
        args <- list()
        # Gets the abundance data from assay, and converts it to otu_table
        otu_table <- as.matrix(assay(x, assay.type))
        otu_table <- phyloseq::otu_table(otu_table, taxa_are_rows = TRUE)
        # Adds to the list
        args[["otu_table"]] <- otu_table

        # If rowData includes information
        if(!( length(rowData(x)[,taxonomyRanks(x)]) == 0 ||
                is.null((rowData(x)[,taxonomyRanks(x)])) )){
            # Converts taxonomy table to characters if it's not already
            rowData(x) <- DataFrame(lapply(rowData(x), as.character))
            # Gets the taxonomic data from rowData, and converts it to tax_table
            tax_table <- as.matrix(rowData(x)[,taxonomyRanks(x),drop=FALSE])
            tax_table <- phyloseq::tax_table(tax_table)
            # Adds to the list
            args[["tax_table"]] <- tax_table
        }
        
        # If colData includes information
        if(!( length(colData(x)) == 0 || is.null(ncol(colData(x))) )){
            # Gets the feature_data from colData and converts it to sample_data
            sample_data <- as.data.frame(colData(x))
            sample_data <- phyloseq::sample_data(sample_data)
            # Adds to the list
            args[["sample_data"]] <- sample_data
        }
        
        # Creates a phyloseq object
        phyloseq <- do.call(phyloseq::phyloseq, args)
        return(phyloseq)
    }
)

#' @rdname convertFromPhyloseq
#' @export
setMethod("convertToPhyloseq", signature = c(x = "TreeSummarizedExperiment"),
    function(x, tree.name = tree_name, tree_name = "phylo", ...){
        # If rowTrees exist, check tree.name
        if( length(rowTreeNames(x)) > 0 ){
            .check_rowTree_present(tree.name, x)
            # Subset the data based on the tree
            x <- x[ rowLinks(x)$whichTree == tree.name, ]
            add_phy_tree <- TRUE
        } else{
            add_phy_tree <- FALSE
        }
        #
        
        # phyloseq and tree objects require nonduplicated rownames. If there are
        # duplicated rownames, they are converted so that they are unique
        if( any(duplicated(rownames(x))) ){
            rownames(x) <- getTaxonomyLabels(x)
        }
        # Gets otu_table object from the function above this, if tse contains
        # only abundance table.
        # Otherwise, gets a phyloseq object with otu_table, and tax_table
        # and/or sample_data
        obj <- callNextMethod()
        # List of arguments
        args <- list()
        # Adds to the list of arguments, if 'obj' is not a phyloseq object
        # i.e. is an otu_table
        if(!is(obj,"phyloseq")){
            # Adds otu_table to the list
            args[["otu_table"]] <- obj
        }
        
        # Add phylogenetic tree
        if( add_phy_tree ){
            phy_tree <- .get_rowTree_for_phyloseq(x, tree.name)
            # If the object is a phyloseq object, adds phy_tree to it
            if(is(obj,"phyloseq")){
                phyloseq::phy_tree(obj) <- phy_tree
            } else{
                # Adds to the list
                args[["phy_tree"]] <- phy_tree
            }
        }
        
        # If referenceSeq has information, stores it to refseq and converts is
        # to phyloseq's refseq.
        if( !is.null(referenceSeq(x)) ){
            # Get referenceSeqs
            refseq <- .get_referenceSeq_for_phyloseq(x, ...)
            # IF refSeq passed the test, add it
            if( !is.null(refseq) ){
                # Convert it to phyloseq object
                refseq <- phyloseq::refseq(refseq)
                # If the object is a phyloseq object, adds refseq to it
                if(is(obj,"phyloseq")){
                    obj <- phyloseq::merge_phyloseq(obj, refseq)
                } else{
                    # Adds to the list
                    args[["refseq"]] <- refseq
                }
            }
        }
        
        # If 'obj' is not a phyloseq object, creates one.
        if(!is(obj,"phyloseq")){
            # Creates a phyloseq object
            phyloseq <- do.call(phyloseq::phyloseq, args)
        } else{
            phyloseq <- obj
        }
        phyloseq
    }
)

################################ HELP FUNCTIONS ################################
# If tips do not match with rownames, prune the tree
.get_x_with_pruned_tree <- function(x, tree.name){
    # Get rowLinks
    row_links <- rowLinks(x)
    # Gets node labels
    node_labs <- row_links[ , "nodeLab"]
    # Prunes the tree
    tree_pruned <- ape::keep.tip(rowTree(x), node_labs)
    # Replace tip labels with corresponding rownames
    tree_pruned$tip.label <- rownames(x)
    # Assigns the pruned tree back to TSE object
    rowTree(x) <- tree_pruned
    warning("Tips of rowTree are renamed to match rownames.", call. = FALSE)
    return(x)
}

# In phyloseq, tips and rownames must match
.get_rowTree_for_phyloseq <- function(x, tree.name){
    # Check if the rowTree's tips match with rownames:
    # tips labels are found from rownames
    if( any(!( rowTree(x, tree.name)$tip.label) %in% rownames(x)) ){
        # If rowtree do not match, tree is pruned
        x <- .get_x_with_pruned_tree(x, tree.name)
    }
    # Get rowTree
    phy_tree <- rowTree(x, tree.name)
    # Convert rowTree to phyloseq object
    phy_tree <- phyloseq::phy_tree(phy_tree)
        
    return(phy_tree)
}

.get_referenceSeq_for_phyloseq <- function(x, referenceSeq = 1, ...){
    # Get reference seqs
    refSeqs <- referenceSeq(x)
    # Is referenceSeq a list / does it contain multiple DNA sets
    is_list <- is(refSeqs, "DNAStringSetList")
    
    # Take only one set, if it is a list
    if( is_list ){
        # Check referenceSeq
        if( !( (.is_non_empty_string(referenceSeq) &&
                referenceSeq %in% names(refSeqs)) ||
            (.is_an_integer(referenceSeq) && (referenceSeq>0 &&
                referenceSeq<=length(refSeqs))) )
            ){
            stop("'referenceSeq' must be a non-empty single character value ",
                "or an integer specifying the DNAStringSet from ",
                "DNAStringSetList.", call. = FALSE)
        }
        # Get specified referenceSeq
        refSeqs <- refSeqs[[referenceSeq]]
        warning("Use 'referenceSeq' to specify DNA set from ",
                "DNAStringSetList. Current choice is '", referenceSeq, "'.",
                call. = FALSE)
    }
    # Check if all rownames have referenceSeqs
    if( !(all(rownames(x) %in% names(refSeqs)) &&
        all(names(refSeqs) %in% rownames(x) )) ){
        warning("referenceSeq does not match with rownames so they are ",
                "discarded.", call. = FALSE)
        refSeqs <- NULL
    }
    return(refSeqs)
}