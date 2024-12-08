#' Import HUMAnN results to \code{TreeSummarizedExperiment}
#
#' @param file \code{Character scalar}. Defines the file
#' path of the HUMAnN file. The file must be in merged HUMAnN format.
#'
#' @param col.data a DataFrame-like object that includes sample names in
#' rownames, or a single \code{character} value defining the file
#' path of the sample metadata file. The file must be in \code{tsv} format
#' (Default: \code{NULL}).
#' 
#' @param colData Deprecated. Use \code{col.data} instead.
#' 
#' @param ... additional arguments:
#' \itemize{
#'   \item \code{assay.type}:  \code{Character scalar}. Specifies the name of
#'   the assay used in calculation. (Default: \code{"counts"})
#'   \item \code{prefix.rm}:  \code{Logical scalar}. Should
#'     taxonomic prefixes be removed? (Default: \code{FALSE})
#'   \item \code{remove.suffix}: \code{Logical scalar}. Should
#'     suffixes of sample names be removed? HUMAnN pipeline adds suffixes
#'     to sample names. Suffixes are formed from file names. By selecting
#'     \code{remove.suffix = TRUE}, you can remove pattern from end of sample
#'     names that is shared by all. (Default: \code{FALSE})
#' }
#'
#' @details
#' Import HUMAnN (currently version 3.0 supported) results of functional
#' predictions based on metagenome composition (e.g. pathways or gene families).
#' The input must be in merged HUMAnN format. (See
#' \href{https://github.com/biobakery/humann#humann_join_tables}{the HUMAnN
#' documentation and \code{humann_join_tables} method.})
#' 
#' The function parses gene/pathway information along with taxonomy information
#' from the input file. This information is stored to \code{rowData}. Abundances
#' are stored to \code{assays}.
#' 
#' Usually the workflow includes also taxonomy data from Metaphlan. See
#' \link[=importMetaPhlAn]{importMetaPhlAn} to load the data to \code{TreeSE}.
#'
#' @return  A
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object
#'
#' @name importHUMAnN
#' 
#' @seealso
#' \code{\link[=importMetaPhlAn]{importMetaPhlAn}}
#' \code{\link[=convertFromPhyloseq]{convertFromPhyloseq}}
#' \code{\link[=convertFromBIOM]{convertFromBIOM}}
#' \code{\link[=convertFromDADA2]{convertFromDADA2}}
#' \code{\link[=importQIIME2]{importQIIME2}}
#' \code{\link[=importMothur]{importMothur}}
#'
#' @export
#' 
#' @references
#' Beghini F, McIver LJ, Blanco-Míguez A, Dubois L, Asnicar F, Maharjan S,
#' Mailyan A, Manghi P, Scholz M, Thomas AM, Valles-Colomer M, Weingart G,
#' Zhang Y, Zolfo M, Huttenhower C, Franzosa EA, & Segata N (2021)
#' Integrating taxonomic, functional, and strain-level profiling of diverse
#' microbial communities with bioBakery 3. \emph{eLife}. 10:e65088.
#' 
#' @examples
#' # File path
#' file_path <- system.file("extdata", "humann_output.tsv", package = "mia")
#' # Import data
#' tse <- importHUMAnN(file_path)
#' tse
#' 
NULL

importHUMAnN <- function(file, col.data = colData, colData = NULL, ...){
    ################################ Input check ###############################
    if(!.is_non_empty_string(file)){
        stop("'file' must be a single character value.",
            call. = FALSE)
    }
    if (!file.exists(file)) {
        stop(file, " does not exist", call. = FALSE)
    }
    if(!is.null(col.data) &&
        !(.is_non_empty_string(col.data) || is.data.frame(col.data) ||
            is.matrix(col.data) || is(col.data, "DataFrame")) ){
        stop("'col.data' must be a single character value, DataFrame or NULL.",
            call. = FALSE)
    }
    ############################## Input check end #############################
    # Humann files has these columns that goes to rowData
    rowdata_col <- c("Pathway", "Gene_Family")
    # Read metaphlan data
    data <- .read_humann(file, rowdata_col, ...)
    # Create TreeSE from the data
    tse <- .create_tse_from_humann(data, rowdata_col, ...)
    # Add col.data if provided
    if( !is.null(col.data) ){
        tse <- .add_coldata(tse, col.data)
    }
    return(tse)
}

################################ HELP FUNCTIONS ################################

# Read Humann file, catch error if it occurs
#' @importFrom utils read.delim
.read_humann <- function(file, rowdata_col, remove.suffix = FALSE, ...){
    ################################ Input check ###############################
    if(!.is_a_bool(remove.suffix)){
        stop("'remove.suffix' must be TRUE or FALSE.", call. = FALSE)
    }
    ############################## Input check end #############################
    # Read the table. Catch error and give more informative message
    table <- tryCatch(
        {
            read.delim(file, check.names = FALSE)
        },
        error = function(condition){
            stop("Cannot read the file: ", file,
                "\nPlease check that the file is in merged HUMAnN file ",
                "format.", call. = FALSE)
        }
    )
    # In the first column name, there is "# " prefix. Remove it
    colnames(table)[1] <- gsub("# ", "", colnames(table)[1])
    # Replace spaces with underscore
    colnames(table)[1] <- gsub(" ", "_", colnames(table)[1])
    # Remove possible suffix from the colnames if user has specified
    if( remove.suffix ){
        table <- .remove_suffix(table, rowdata_col)
    }
    # Add rownames
    rownames(table) <- table[, 1] 
    # Check that file is in right format
    if( .check_metaphlan(table, rowdata_col) ){
        stop("Cannot read the file: ", file,
            "\nPlease check that the file is in merged HUMAnN file format.",
            call. = FALSE)
    }
    return(table)
}

# This function parses humann file and creates tse from it.
.create_tse_from_humann <- function(
        data, rowdata_col, assay.type = "counts", ...){
    # Get rowdata columns
    rowdata_id <- unlist(lapply(rowdata_col, grep, colnames(data)))
    rowdata <- data[ , rowdata_id, drop = FALSE]
    # Get columns that go to assay
    assay <- data[ , -rowdata_id, drop = FALSE]
    
    # Parse rowdata. The data includes gene/pathway info along with taxonomy
    # info. Separate the information to own columns.
    rowdata_id <- unlist(lapply(rowdata_col, grep, colnames(rowdata)))
    # Get the column and split gene/pathway and taxonomy info
    rowdata_temp <- strsplit(rowdata[[1]], "\\|",)
    # Are some rows missing taxonomy info? Get their indices.
    missing_taxa <- lengths(rowdata_temp) == 1
    # Create a df from the list of splitted data.
    rowdata_temp <- data.frame(t(data.frame(rowdata_temp)))
    # Replace missing taxonomy info with NA
    rowdata_temp[missing_taxa, 2] <- NA
    # Replace column names
    colnames(rowdata_temp) <- c(colnames(rowdata), "Taxonomy")
    
    # Now we have rowdata that includes gene/pathway info in one column and
    # taxonomy info in other. Let's parse the taxonomy info so that species
    # genus etc levels are in unique columns.
    taxonomy <- .parse_taxonomy(
        rowdata_temp, column_name = "Taxonomy", sep = "\\.", ...)
    
    # Convert all data to DataFrame
    rowdata <- DataFrame(rowdata)
    rowdata_temp <- DataFrame(rowdata_temp)
    
    # Create rowdata from the information that is parsed
    colnames(rowdata) <- paste0(colnames(rowdata), "_long")
    rowdata <- cbind(rowdata, rowdata_temp[, 1, drop = FALSE])
    rowdata <- cbind(rowdata, taxonomy)
    
    # Create assays
    assay <- as.matrix(assay)
    assays <- SimpleList(counts = assay)
    names(assays) <- assay.type

    # Create TreeSE
    tse <- TreeSummarizedExperiment(assays = assays, rowData = rowdata)
    return(tse)
}
