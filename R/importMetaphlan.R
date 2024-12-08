#' Import Metaphlan results to \code{TreeSummarizedExperiment}
#
#' @param file \code{NULL} or \code{DataFrame}-like object. Defines the file
#'   path of the Metaphlan file. The file must be in merged Metaphlan format.
#'
#' @param col.data a DataFrame-like object that includes sample names in
#'   rownames, or a single \code{character} value defining the file
#'   path of the sample metadata file. The file must be in \code{tsv} format
#'   (Default: \code{NULL}).
#' 
#' @param colData Deprecated. use \code{col.data} instead.
#'   
#' @param sample_meta Deprecated. Use \code{col.data} instead.
#' 
#' @param tree.file \code{Character scalar}. Defines the file
#'   path of the phylogenetic tree.
#'   (Default: \code{NULL}).
#' 
#' @param phy_tree Deprecated. Use \code{tree.file} instead.
#'   
#' @param ... additional arguments:
#' \itemize{
#'   \item \code{assay.type}: \code{Character scalar}. Specifies the name of
#'   the assay used in calculation. (Default: \code{"counts"})
#'   \item \code{prefix.rm}: \code{Logical scalar}. Should
#'     taxonomic prefixes be removed? (Default: \code{FALSE})
#'   \item \code{remove.suffix}: \code{Logical scalar}. Should
#'     suffixes of sample names be removed? Metaphlan pipeline adds suffixes
#'     to sample names. Suffixes are formed from file names. By selecting
#'     \code{remove.suffix = TRUE}, you can remove pattern from end of sample
#'     names that is shared by all. (Default: \code{FALSE})
#'   \item \code{set.ranks}: \code{Logical scalar}. Should the columns
#'     in the rowData that are treated as taxonomy ranks be updated according to
#'     the ranks found in the imported data?
#'     (Default: \code{FALSE})
#' }
#'
#' @details
#' Import Metaphlan (versions 2, 3 and 4 supported) results.
#' Input must be in merged Metaphlan format.
#' (See
#' \href{https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4#merging-tables}{
#' the Metaphlan documentation and \code{merge_metaphlan_tables} method.})
#' Data is imported so that data at the lowest rank is imported as a 
#' \code{TreeSummarizedExperiment} object. Data at higher rank is imported as a
#' \code{SummarizedExperiment} objects which are stored to \code{altExp} of
#' \code{TreeSummarizedExperiment} object.
#' 
#' Currently Metaphlan versions 2, 3, and 4 are supported.
#'
#' @return  A
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object
#'
#' @name importMetaPhlAn
#' @seealso
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
#' Integrating taxonomic, functional, 
#' and strain-level profiling of diverse microbial communities with bioBakery 3.
#' \emph{eLife}. 10:e65088. doi: 10.7554/eLife.65088
#'
#' @examples
#' # (Data is from tutorial
#' # https://github.com/biobakery/biobakery/wiki/metaphlan3#merge-outputs)
#' 
#' # File path
#' file_path <- system.file(
#'     "extdata", "merged_abundance_table.txt", package = "mia")
#' # Import data
#' tse <- importMetaPhlAn(file_path)
#' # Data at the lowest rank
#' tse
#' # Data at higher rank is stored in altExp
#' altExps(tse)
#' # Higher rank data is in SE format, for example, Phylum rank
#' altExp(tse, "phylum")
#' 
NULL

importMetaPhlAn <- function(
        file, col.data = colData, colData = sample_meta,
        sample_meta = NULL, tree.file = phy_tree, phy_tree = NULL, ...){
    
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
    if(!is.null(tree.file) && !.is_non_empty_string(tree.file)){
        stop("'tree.file' must be a single character value or NULL.",
            call. = FALSE)
    }
    ############################## Input check end #############################
    # Get rowdata columns. metaphlan v2 has ID column. Metaphlan > v2 has
    # clade_name for taxonomy names. Some has taxonomy.
    rowdata_col <- c("clade_name", "ID", "_id", "taxonomy")
    # Read metaphlan data
    data <- .read_metaphlan(file, rowdata_col, ...)
    # Parse data into separate tables, which include data at certain taxonomy
    # rank
    tables <- .parse_metaphlan(data, ...)

    # Create multiple SE objects at different rank from the data
    se_objects <- lapply(tables, function(x){
        .create_se_from_metaphlan(x, rowdata_col, ...)
        })
    
    # Get the object with lowest rank
    tse <- se_objects[[ length(se_objects) ]]
    # Convert it to TreeSE so that it has altExp
    tse <- as(tse, "TreeSummarizedExperiment")
    # Remove it, so that it is not added multiple times
    se_objects[[ length(se_objects) ]] <- NULL
    # Add rest of the objects to altExp
    if( length(se_objects) > 0 ){
        for( rank in names(se_objects) ){
            # Add SEs to altExp of tse, give name according to rank
            altExp(tse, rank) <- se_objects[[rank]]
        }
    }
    # Set taxonomy ranks using .set_taxonomy_ranks
    .set_ranks_based_on_rowdata(tse,...)
    
    # Load sample meta data if it is provided
    if( !is.null(col.data) ) {
        tse <- .add_coldata(tse, col.data)
    }
    
    
    # Load tree if it is provided
    if (!is.null(tree.file)) {
        tree <- ape::read.tree(tree.file)
        rowTree(tse) <- tree
    } 
    
    return(tse)
}

################################ HELP FUNCTIONS ################################

# Read Metaphlan file, catch error if it occurs
.read_metaphlan <- function(file, rowdata_col, remove.suffix = FALSE, ...){
    ################################ Input check ###############################
    if(!.is_a_bool(remove.suffix)){
        stop("'remove.suffix' must be TRUE or FALSE.", call. = FALSE)
    }
    ############################## Input check end #############################
    # Read the table. Catch error and give more informative message
    table <- tryCatch(
        {
            read.table(
                file, header = TRUE, comment.char = "#", check.names = FALSE)
        },
        error = function(condition){
            stop("Cannot read the file: ", file,
                "\nPlease check that the file is in merged Metaphlan file ",
                "format.", call. = FALSE)
        }
    )
    # Check that file is in right format
    if( .check_metaphlan(table, rowdata_col) ){
        stop("Cannot read the file: ", file,
            "\nPlease check that the file is in merged Metaphlan file format.",
            call. = FALSE)
    }
    # Remove possible suffix from the colnames if user has specified
    if( remove.suffix ){
        table <- .remove_suffix(table, rowdata_col)
    }
    return(table)
}

# Check that metaphlan file contains correct information
.check_metaphlan <- function(data, rowdata_col){
    # Get rowData columns
    rowdata_id <- unlist(lapply(rowdata_col, grep, colnames(data)))
    rowdata_columns <- data[ , rowdata_id, drop = FALSE]
    # Get columns that go to assay
    assay_columns <- data[ , -rowdata_id, drop = FALSE]
    # Initialize result 
    result <- TRUE
    
    # Check rowdata column names that they contain right information, and check
    # that rest of the columns represents abundances in samples.
    # If these requirements are met, give FALSE. Otherwise, give TRUE.
    if( length(rowdata_id) > 0 &&
            all(unlist(lapply(assay_columns, is.numeric))) ){
        result <- FALSE
    }
    return(result)
}

# Get metaphlan table as input and return multiple tables which each include
# data at certain taxonomy rank
.parse_metaphlan <- function(table, ...){
    # ID in Metaphlan v2, > 2 clade_name
    col <- colnames(table) %in% c("clade_name", "ID")
    if( sum(col) != 1 ){
        stop("Cannot parse Metaphlan file.", call. = FALSE)
    }
    # Get the lowest level of each row
    
    levels <- lapply(table[ , col], FUN = .get_lowest_taxonomic_level)
    # Convert list to vector
    levels <- unlist(levels)
    # Split table so that each individual table contains information only
    # at specific rank
    tables <- split(table, levels)
    # Get the order
    metaphlan_tax <- names(.taxonomy_rank_prefixes)
    indices <- match(tolower(metaphlan_tax), tolower(names(tables)))
    # Remove NAs which occurs if rank is not included
    indices <- indices[!is.na(indices)]
    # Order tables 
    tables <- tables[indices]
    return(tables)
}

# Get the lowest level of the string that contains multiple taxonomic levels
# with prefixes Output is single character that specifies the rank, e.g,
# "s" == "Species"
.get_lowest_taxonomic_level <- function(string){
    # List all ranks and what prefix they correspond
    ranks <- .taxonomy_rank_prefixes
    # Get indices that specify location of rank prefixes 
    ranks_pattern <- paste(ranks, collapse = "|")
    ranks_pattern <- paste0("(", ranks_pattern, ")__")
    levels <- gregexpr(ranks_pattern, string)[[1]]
    # Get the location of lowest rank
    lowest_level_ind <- levels[length(levels)]
    # Get the lowest rank that was found
    lowest_level <- substr(
        string, start = lowest_level_ind, stop = lowest_level_ind)
    
    # Convert prefix into full rank name
    lowest_level <- names(ranks[ match(lowest_level, ranks) ])
    return(lowest_level)
}

# Create SE object that include rowdata and assay, from the metaphlan table
.create_se_from_metaphlan <- function(
        table, rowdata_col, assay.type = assay_name, assay_name = "counts",
        ...){
    # Check assay.type
    if( !.is_non_empty_character(assay.type) ){
        stop("'assay.type' must be a non-empty character value.",
            call. = FALSE)
    }
    # Get rowdata columns
    rowdata_id <- unlist(lapply(rowdata_col, grep, colnames(table)))
    # Get those columns that belong to rowData
    rowdata <- table[, rowdata_id, drop = FALSE]
    # Get those columns that belong to assay
    assay <- table[, -rowdata_id, drop = FALSE]
    # Parse taxonomic levels
    taxonomy <- .parse_taxonomy(
        rowdata[ , 1, drop = FALSE], sep = "\\|",
        column_name = colnames(rowdata)[[1]], ...)
    # Add parsed taxonomy level information to rowdata
    rowdata <- cbind(taxonomy, rowdata)
    # Ensure that rowData is DataFrame
    rowdata <- DataFrame(rowdata, check.names = FALSE)
    
    # Ensure that assay is matrix
    assay <- as.matrix(assay)
    # Create assays list and add assay with specific name
    assays <- S4Vectors::SimpleList()
    assays[[assay.type]] <- assay
    
    # Create TreeSE
    se <- TreeSummarizedExperiment(
        assays = assays, rowData = rowdata)
    # Add taxonomy labels
    rownames(se) <- getTaxonomyLabels(se)
    return(se)
}

# This function can be used to add col.data to TreeSE. It checks that sample
# names match (full or partial) and adds the metadata to altExps also.
.add_coldata <- function(tse, coldata){
    # If the coldata is character specifying the path
    if( .is_non_empty_character(coldata) ){
        coldata <- read.table(file = coldata, header = TRUE, sep = "\t")
        rownames(coldata) <- coldata[, 1]
    }
    # Ensure that the coldata is DF
    coldata <- DataFrame(coldata)

    # Usually when metaphlan sample names are created, the workflow adds file
    # name as a suffix. Because of that sample names do not match anymore
    # necessarily. Try partial match if full match is not found.
    if( all(rownames(coldata) %in% colnames(tse)) ){
        sample_names <- rownames(coldata)
        names(sample_names) <- sample_names
    } else{
        sample_names <- lapply(rownames(coldata), function(x){
            x <- colnames(tse)[grep(x, colnames(tse))]
            if( length(x) != 1 ){
                x <- NULL
            }
            return(x)
        })
        names(sample_names) <- rownames(coldata)
        sample_names <- unlist(sample_names)
    }

    # Check if all samples were found. (In partial match certain sample name
    # is missing if one matching name was not found.). In this part, all
    # colnames should be found if data sets are matching. (More samples in
    # metadata is allowed.)
    if( !(all(colnames(tse) %in% sample_names) &&
            length(sample_names) == ncol(tse)) ){
        warning(
            "The sample names in 'col.data' do not match with the data. ",
            "The sample metadata is not added.", call. = FALSE
            )
        return(tse)
    }

    # Reorder sample metadata based on the data
    sample_names <- sample_names[match(colnames(tse), sample_names)]
    coldata <- coldata[names(sample_names), ]

    # Give warning if partial match was used
    if( !all(rownames(coldata) %in% colnames(tse)) ){
        warning("Partial match was used to match sample names between ",
                "'col.data' and the data. Please check that they are correct.",
                call. = FALSE
        )
        # Replace colnames with names from sample metadata. They are without
        # additional suffix.
        coldata[["colnames_data"]] <- colnames(tse)
        coldata[["colnames_metadata"]] <- rownames(coldata)
        colnames(tse) <- rownames(coldata)
    }

    # Add sample metadata
    colData(tse) <- coldata
    # If there are altExps, add the metadata also there
    if( length(altExps(tse)) ){
        altexps <- altExps(tse)
        altexps <- lapply(altexps, function(x){
            colnames(x) <- rownames(coldata)
            colData(x) <- coldata
            return(x)
            })
        altexps <- SimpleList(altexps)
        altExps(tse) <- altexps
        
    }
    return(tse)
}

# In humann/metaphlan pipeline, suffix is added to colnames. The suffix comes
# from file names. This can cause problems when, e.g., taxonomy and pathway
# information is combined. Because their suffixes differ, the sample names
# differ. The purpose of the function is to remove those file names.
.remove_suffix <- function(
        data, rowdata_col = c("clade_name", "ID", "_id", "taxonomy")){
    # Get rowdata columns
    rowdata_id <- unlist(lapply(rowdata_col, grep, colnames(data)))
    # Get sample names
    sample_names <- colnames(data)[-rowdata_id]
    # Get the possible suffixes, i.e., words that are divided by underscores
    # based on first sample name.
    suffix <- strsplit(sample_names[[1]], "_")[[1]]
    # Loop over these suffixes. From the end, add words one by one. Store those
    # suffixes that match with all sample names. The result is suffix that is
    # the longest.
    shared <- NULL
    for( i in length(suffix):1 ){
        pattern <- paste0("_", paste(suffix[i:length(suffix)], collapse="_"))
        if( all(grepl(pattern, sample_names)) ){
            shared <- pattern
        }
    }
    # Remove suffix from sample names
    if( !is.null(shared) ){
        sample_names <- gsub(shared, "", sample_names)
        colnames(data)[-rowdata_id] <- sample_names
    }
    return(data)
}
