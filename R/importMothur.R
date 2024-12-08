#' Import Mothur results as a \code{TreeSummarizedExperiment}
#' 
#' This method creates a \code{TreeSummarizedExperiment} object from
#' \code{Mothur} files provided as input. 
#'
#' @param assay.file \code{Character scalar}. Defines the file
#'   path of the feature table to be imported. The File has to be in 
#'   \code{shared file} format as defined in Mothur documentation.
#' 
#' @param sharedFile Deprecated. Use \code{assay.file} instead.
#'
#' @param row.file \code{Character scalar}. Defines the file path of
#'   the taxonomy table to be imported. The File has to be in \code{taxonomy
#'   file} or \code{constaxonomy file} format  as defined in Mothur
#'   documentation. (Default: \code{NULL}).
#' 
#' @param taxonomyFile Deprecated. Use \code{row.file} instead.
#'
#' @param col.file \code{Character scalar}. Defines file path of
#'   the sample metadata to be imported. The File has to be in \code{design
#'   file} format as defined in Mothur documentation. (Default: \code{NULL}).
#' 
#' @param designFile Deprecated. Use \code{col.file} instead.
#'
#' @details
#' Results exported from Mothur can be imported as a
#' \code{SummarizedExperiment} using \code{importMothur}. Except for the
#' \code{assay.file}, the other data types, \code{row.file}, and
#' \code{col.file}, are optional, but are highly encouraged to be provided.
#'
#' @return  A
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object
#'
#' @name importMothur
#' @seealso
#' \code{\link[=convertFromPhyloseq]{convertFromPhyloseq}}
#' \code{\link[=convertFromBIOM]{convertFromBIOM}}
#' \code{\link[=convertFromDADA2]{convertFromDADA2}}
#' \code{\link[=importQIIME2]{importQIIME2}}
#'
#' @references
#' \url{https://mothur.org/}
#' \url{https://mothur.org/wiki/shared_file/}
#' \url{https://mothur.org/wiki/taxonomy_file/}
#' \url{https://mothur.org/wiki/constaxonomy_file/}
#' \url{https://mothur.org/wiki/design_file/}
#' 
#' @examples
#' # Abundance table
#' counts <- system.file("extdata", "mothur_example.shared", package = "mia")
#' # Taxa table (in "cons.taxonomy" or "taxonomy" format)
#' taxa <- system.file(
#'     "extdata", "mothur_example.cons.taxonomy", package = "mia")
#' #taxa <- system.file("extdata", "mothur_example.taxonomy", package = "mia")
#' # Sample meta data
#' meta <- system.file("extdata", "mothur_example.design", package = "mia")
#' 
#' # Creates se object from files
#' se <- importMothur(assay.file = counts, row.file = taxa, col.file = meta)
#' # Convert SE to TreeSE
#' tse <- as(se, "TreeSummarizedExperiment")
#' tse
NULL

#' @rdname importMothur
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors make_zero_col_DFrame
#' @export
importMothur <- function(assay.file = sharedFile,
                        sharedFile, 
                        taxonomyFile = NULL,
                        row.file = taxonomyFile,
                        designFile = NULL,
                        col.file = designFile) {

    # input check
    if(!.is_non_empty_string(assay.file)){
        stop("'assay.file' must be a single character value.", 
            call. = FALSE)
    }
    if(!is.null(row.file) && !.is_non_empty_string(row.file)){
        stop("'row.file' must be a single character value or NULL.", 
            call. = FALSE)
    }
    if(!is.null(col.file) && !.is_non_empty_string(col.file)){
        stop("'col.file' must be a single character value or NULL.",
            call. = FALSE)
    }
    
    # Reads the assay.file 
    feature_tab_and_data_to_colData <- .read_mothur_feature(assay.file)
    # Extracts feature_tab
    feature_tab <- feature_tab_and_data_to_colData$assay
    # Extracts data that goes to colData
    data_to_colData <- feature_tab_and_data_to_colData$colData
    
    # If rowData information exists, gets that. Otherwise, tax_tab is just 
    # data.frame without information
    if (!is.null(row.file)) {
        taxa_tab <- .read_mothur_taxonomy(row.file, feature_tab)
    } else {
        taxa_tab <- S4Vectors::make_zero_col_DFrame(nrow(feature_tab))
        rownames(taxa_tab) <- rownames(feature_tab)
    }
    
    # If colData informationor data_to_colData exists, gets that. Otherwise, 
    # sample_tab is just data frame without information
    if (!is.null(col.file) && !is.null(data_to_colData)) {
        sample_meta <- .read_mothur_sample_meta(col.file, data_to_colData)
    } else {
        sample_meta <- S4Vectors::make_zero_col_DFrame(ncol(feature_tab))
        rownames(sample_meta) <- colnames(feature_tab)
    }

    TreeSummarizedExperiment(
        assays = S4Vectors::SimpleList(counts = feature_tab),
        rowData = taxa_tab,
        colData = sample_meta)
}

# These extra information must be added to colData. Return list of assay and 
# extra info
.read_mothur_feature <- function(assay.file){

    if (!.is_mothur_shared_file(assay.file)) {
        stop("The input '", assay.file, "' must be in `shared` format.",
            call. = FALSE)
    }

    # Stores name of columns will be included in colData not in assays
    MOTHUR_NON_ASSAY_COLS <- c("label","numOtus","Group")
    data <- read.table(assay.file, check.names=FALSE, header=TRUE,
                        sep="\t", stringsAsFactors=FALSE)
    # Checks that colnames contain information and it is not NULL
    if ( !(length(colnames(data)) > 0) || is.null(colnames(data)) ){
        stop("'shared' does not include names of taxa.",
            call. = FALSE)
    }
    # Takes all columns but not those that goes to colData, 
    # and transforms the data frame to matrix
    assay <- as.matrix(data[,!(colnames(data) %in% MOTHUR_NON_ASSAY_COLS )])
    # Transposes the matrix --> taxa to rows
    assay <- t(assay)
    # Gets those data that goes to colData, and creates a data frame it
    colData <- DataFrame(data[,MOTHUR_NON_ASSAY_COLS])
    return(list(assay = assay,
                colData = colData))
}

.read_mothur_taxonomy <- function(row.file, feature_tab){
    
    # If the file is in "cons.taxonomy" format
    if (.is_mothur_constaxonomy_file(row.file, feature_tab)) {
        data <- read.table(row.file, check.names=FALSE,
                            header=TRUE, sep="\t", stringsAsFactors=FALSE)
    } 
    # If the file is in "taxonomy" format, adds column names
    else if (.is_mothur_taxonomy_file(row.file, feature_tab)){
        data <- read.table(row.file, check.names=FALSE,
                            header=FALSE, sep="\t", 
                            stringsAsFactors=FALSE, 
                            col.names = c("OTU", "Taxonomy"))
    }
    # Else the file is not either gives an error
    else{
        stop("The input '", row.file, "' must be provided in the ",
            "`taxonomy` or `cons.taxonomy` format. In addition, it must ",
            "match the data of the 'assay.file'",
            call. = FALSE)
    }

    # Column that includes taxonomical information
    MOTHUR_TAX_COL <- "Taxonomy"

    # Checks that colnames contain information, it is not NULL, and taxonomical 
    # information is present 
    if ( !(length(colnames(data)) > 0) || 
        is.null(colnames(data)) || 
        is.null(data[[MOTHUR_TAX_COL]]) ){
        stop("'taxonomy' does not include taxonomical information.",
            call. = FALSE)
    }
    
    # Removes additional characters between taxa
    data [,MOTHUR_TAX_COL] <- gsub("[\"]", "", data [,MOTHUR_TAX_COL])
    data [,MOTHUR_TAX_COL] <- gsub("[(1-100)]", "", data [,MOTHUR_TAX_COL])
    # Splits taxa level into separate columns
    into <- c("Kingdom", "Phylum", "Order", "Class", "Family", "Genus")
    tax <- tidyr::separate(data,
                            MOTHUR_TAX_COL,
                            into=into,
                            sep=";",
                            extra="merge")
    # Removes ";" from the end of genus level names
    tax$Genus <- gsub(";", "", tax$Genus)
    rowData <- tax
    # Adds rownames
    rownames(rowData) <- rowData$OTU
    return(rowData)
}

.read_mothur_sample_meta <- function(col.file, data_to_colData){
    # Checks if file is in "design" format. data_to_colData$Group includes 
    # sample names that were extracted from assay, i.e. assay.file
    if (!.is_mothur_design_file(col.file, data_to_colData$Group)) {
        stop("The input '", col.file, "' must be in `design` format, 
            and it must inlude same sample names as 'assay.file'.",
            call. = FALSE)
    }

    # Reads the file
    colData <- read.table(col.file, check.names=FALSE,
                            header=TRUE, sep="\t",
                            stringsAsFactors=FALSE)
    
    # Combines the extracted colData and data from the assay
    colData <- cbind(colData, data_to_colData)
    
    # Adds sample names to rownames of the data 
    rownames(colData) <- colData$group
    return(colData)
}

#' extract file extension
#' @noRd
.is_mothur_shared_file <- function(file){
    result <- FALSE
    # Columns that every shared file should have
    columns_that_must_be_found <- c("label","Group", "numOtus")
    
    # Reads the data
    data <- read.table(file, check.names=FALSE, header=TRUE,
                        sep="\t", stringsAsFactors=FALSE)
    
    # If data contains column names, then it is shared file
    if( identical(colnames(data)[seq_len(3)], columns_that_must_be_found) ){
        result <- TRUE
    }
    return(result)
}

#' extract file extension
#' @noRd
.is_mothur_taxonomy_file <- function(file, feature_tab){
    result <- FALSE
    
    # Reads the data
    data <- read.table(file, check.names=FALSE, header=FALSE,
                        sep="\t", stringsAsFactors=FALSE)
    
    # If data contains 2 columns and the first column includes same taxa as 
    # feature_tab in rownames, then it is taxonomy file
    if( ncol(data) == 2 && identical(data[,1], rownames(feature_tab)) ){
        result <- TRUE
    }
    return(result)
}

#' extract file extension
#' @noRd
.is_mothur_constaxonomy_file <- function(file, feature_tab){
    result <- FALSE
    # Columns that every constaxonomy file should have
    columns_that_must_be_found <- c("OTU", "Size", "Taxonomy")
    
    # Reads the data
    data <- read.table(file, check.names=FALSE, header=TRUE,
                        sep="\t", stringsAsFactors=FALSE)
    
    # If data contains column names, and "OTU" column that includes same taxa as
    # feature_tab, 
    # then it is constaxonomy file
    if( identical(colnames(data), columns_that_must_be_found) && 
        identical(data$OTU, rownames(feature_tab)) ){
        result <- TRUE
    }
    return(result)
}

.is_mothur_design_file <- function(file, sample_names){
    result <- FALSE
    # Reads the data
    data <- read.table(file, check.names=FALSE, header=TRUE,
                        sep="\t", stringsAsFactors=FALSE)
    
    # If data contains "group" column that include sample names, then it is
    # design file
    if( identical(data$group, sample_names) ){
        result <- TRUE
    }
    return(result)
}
