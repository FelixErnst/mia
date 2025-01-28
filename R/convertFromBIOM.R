#' Convert a \code{TreeSummarizedExperiment} object to/from \sQuote{BIOM}
#' results
#'
#' For convenience, a few functions are available to convert BIOM, DADA2 and
#' phyloseq objects to
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' objects, and
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' objects to phyloseq objects.
#'
#' @param prefix.rm \code{Logical scalar}. Should
#' taxonomic prefixes be removed? The prefixes is removed only from detected
#' taxa columns meaning that \code{rank.from.prefix} should be enabled in the
#' most cases. (Default: \code{FALSE})
#'
#' @param removeTaxaPrefixes Deprecated. Use \code{prefix.rm} instead.
#'
#' @param rank.from.prefix \code{Logical scalar}. If file does not have
#' taxonomic ranks on feature table, should they be scraped from prefixes?
#' (Default: \code{FALSE})
#'
#' @param rankFromPrefix Deprecated.Use \code{rank.from.prefix} instead.
#'
#' @param artifact.rm \code{Logical scalar}. If file have
#' some taxonomic character naming artifacts, should they be removed.
#' (default (Default: \code{FALSE})
#'
#' @param remove.artifacts Deprecated. Use \code{artifact.rm} instead.
#'
#' @details
#' \code{convertFromBIOM} coerces a \code{\link[biomformat:biom-class]{biom}}
#' object to a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
#'
#' \code{convertToBIOM} coerces a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object to a \code{\link[biomformat:biom-class]{biom}} object.
#'
#' @return
#' \code{convertFromBIOM} returns an object of class
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @name importBIOM
#'
#' @seealso
#' \code{\link[=importQIIME2]{importQIIME2}}
#' \code{\link[=importMothur]{importMothur}}
#'
#' @examples
#'
#' # Convert BIOM results to a TreeSE
#' # Load biom file
#' library(biomformat)
#' biom_file <- system.file("extdata", "rich_dense_otu_table.biom",
#'                          package = "biomformat")
#'
#' # Make TreeSE from BIOM object
#' biom_object <- biomformat::read_biom(biom_file)
#' tse <- convertFromBIOM(biom_object)
#'
#' # Convert TreeSE object to BIOM
#' biom <- convertToBIOM(tse)
#'
NULL

#' Import BIOM results to \code{TreeSummarizedExperiment}
#'
#' @param file BIOM file location
#'
#' @param ... additional arguments to be passed to \code{convertFromBIOM}
#'
#' @details
#' \code{importBIOM} loads a BIOM file and creates a
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object from the BIOM object contained in the loaded file.
#'
#' @return
#' \code{importBIOM} returns an object of class
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @name importBIOM
#'
#' @seealso
#' \code{\link[=importMetaPhlAn]{importMetaPhlAn}}
#' \code{\link[=convertFromPhyloseq]{convertFromPhyloseq}}
#' \code{\link[=convertFromBIOM]{convertFromBIOM}}
#' \code{\link[=convertFromDADA2]{convertFromDADA2}}
#' \code{\link[=importQIIME2]{importQIIME2}}
#' \code{\link[=importMothur]{importMothur}}
#' \code{\link[=importHUMAnN]{importHUMAnN}}
#'
#' @examples
#' # Load biom file
#' library(biomformat)
#' biom_file <- system.file(
#'     "extdata", "rich_dense_otu_table.biom", package = "biomformat")
#'
#' # Make TreeSE from biom file
#' tse <- importBIOM(biom_file)
#'
#' # Get taxonomyRanks from prefixes and remove prefixes
#' tse <- importBIOM(
#'     biom_file, rank.from.prefix = TRUE, prefix.rm = TRUE)
#'
#' # Load another biom file
#' biom_file <- system.file(
#'    "extdata", "Aggregated_humanization2.biom", package = "mia")
#'
#' # Clean artifacts from taxonomic data
#' tse <- importBIOM(biom_file, artifact.rm = TRUE)
#'
#' @export
importBIOM <- function(file, ...) {
    .require_package("biomformat")
    biom <- biomformat::read_biom(file)
    convertFromBIOM(biom, ...)
}

#' @rdname importBIOM
#'
#' @param x object of type \code{\link[biomformat:biom-class]{biom}}
#'
#' @export
#' @importFrom S4Vectors make_zero_col_DFrame DataFrame
#' @importFrom dplyr %>% bind_rows
convertFromBIOM <- function(
        x, prefix.rm = removeTaxaPrefixes,
        removeTaxaPrefixes = FALSE, rank.from.prefix = rankFromPrefix,
        rankFromPrefix = FALSE,
        artifact.rm = remove.artifacts, remove.artifacts = FALSE, ...){
    # input check
    .require_package("biomformat")
    if(!is(x,"biom")){
        stop("'x' must be a 'biom' object", call. = FALSE)
    }
    if( !.is_a_bool(prefix.rm) ){
        stop("'prefix.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(rank.from.prefix) ){
        stop("'rank.from.prefix' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(artifact.rm) ){
        stop("'artifact.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    counts <- as(biomformat::biom_data(x), "matrix")
    sample_data <- biomformat::sample_metadata(x)
    feature_data <- biomformat::observation_metadata(x)

    # colData is initialized with empty tables with rownames if it is NULL
    if( is.null(sample_data) ){
        sample_data <- S4Vectors::make_zero_col_DFrame(ncol(counts))
        rownames(sample_data) <- colnames(counts)
    # Otherwise convert it into correct format if it is a list
    } else if( is(sample_data, "list") ){
        # Merge list of data.frames into one
        sample_data <- bind_rows(sample_data)
        sample_data < as.data.frame(sample_data)
    }
    # rowData is initialized with empty tables with rownames if it is NULL
    if( is.null(feature_data) ){
        feature_data <- S4Vectors::make_zero_col_DFrame(nrow(counts))
        rownames(feature_data) <- rownames(counts)
    # Otherwise convert it into correct format if it is a list
    } else if( is(feature_data, "list") ){
        # Feature data is a list of taxa info. Dfs are merged together
        # differently than sample metadata since the column names are only
        # "Taxonomy". If there is only one taxonomy level, the column name does
        # not get a suffix.
        # --> bind rows based on the index of column.

        # Get the maximum length of list
        max_length <- max( lengths(feature_data) )
        # Get the column names from the taxa info that has all the levels that
        # occurs in the data
        colnames <- names( head(
            feature_data[ lengths(feature_data) == max_length ], 1)[[1]])
        # Convert the list so that all individual taxa info have the max length
        # of the list objects. All vectors are appended with NAs, if they do not
        # have all the levels. E.g., if only Kingdom level is found, all lower
        # ranks are now NA
        feature_data <- lapply(feature_data, function(x){
            length(x) <- max_length
            return(x)
        })
        # Create a data.frame from the list
        feature_data <- do.call(rbind, feature_data)
        # Transposing feature_data and make it df object
        feature_data <- as.data.frame(feature_data)
        # Add column that includes all the data
        merged_col <- apply(feature_data, 1, paste0, collapse = ";")
        feature_data <- cbind(feature_data, merged_col)
        # Add correct colnames
        colnames(feature_data) <- c(colnames, "taxonomy_unparsed")
    }
    # If there is only one column in the feature data, it is the most probable
    # that the taxonomy is not parsed. Try to parse it.
    if( ncol(feature_data) == 1 ){
        colnames(feature_data) <- "taxonomy_unparsed"
        tax_tab <- .parse_taxonomy(
            feature_data, column_name = colnames(feature_data))
        feature_data <- cbind(tax_tab, feature_data)
        feature_data <- as.data.frame(feature_data)
    }

    # Clean feature_data from possible character artifacts if specified.
    if( artifact.rm ){
        feature_data <- .detect_taxa_artifacts_and_clean(feature_data, ...)
    }
    # Replace taxonomy ranks with ranks found based on prefixes. Test first
    # if columns are already ranks. If they are, do not change them.
    cols_not_already_ranks <- lapply(
        colnames(feature_data), function(x) !x %in% TAXONOMY_RANKS) |>
        unlist() |> all()
    if( rank.from.prefix && cols_not_already_ranks ){
        # Find ranks
        ranks <- lapply(
            colnames(feature_data), .replace_colnames_based_on_prefix,
            x = feature_data)
        # Replace old ranks with found ranks
        colnames(feature_data) <- unlist(ranks)
    }
    # Adjust row and colnames
    rownames(counts) <- rownames(feature_data) <- biomformat::rownames(x)
    colnames(counts) <- rownames(sample_data) <- biomformat::colnames(x)

    # Convert into DataFrame
    sample_data <- DataFrame(sample_data)
    feature_data <- DataFrame(feature_data)
    # Convert into list
    assays <- SimpleList(counts = counts)
    # Create TreeSE
    tse <- TreeSummarizedExperiment(
        assays = assays,
        colData = sample_data,
        rowData = feature_data)
    # Set ranks based on rowData columns if user has specified to do so
    temp <- .set_ranks_based_on_rowdata(tse, ...)
    # Remove prefixes if specified and rowData includes info
    if(prefix.rm && ncol(rowData(tse)) > 0){
        rowData(tse) <- .remove_prefixes_from_taxa(rowData(tse), ...)
    }
    return(tse)
}

#' @rdname importBIOM
#'
#' @param x
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @param assay.type \code{Character scaler}. The name of assay.
#' (Default: \code{"counts"})
#'
#' @param ... Additional arguments. Not used currently.
#'
#' @export
setMethod(
    "convertToBIOM", signature = c(x = "SummarizedExperiment"),
        function(x, assay.type = "counts", ...){
    ############################## Input check #############################
    .require_package("biomformat")
    .check_assay_present(assay.type, x)
    ############################ Input check end ###########################
    x <- .convert_se_to_biom(x, assay.type)
    return(x)
    }
)

################################ HELP FUNCTIONS ################################
# This function removes prefixes from taxonomy names
.remove_prefixes_from_taxa <- function(
        feature_tab,
        prefixes = paste0(
            "(", paste0(.taxonomy_rank_prefixes, collapse = "|"), ")__"),
        only.taxa.col = TRUE, ignore.col = "taxonomy_unparsed", ...){
    #
    if( !.is_a_bool(only.taxa.col) ){
        stop("'only.taxa.col' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    if( !.is_a_string(prefixes) ){
        stop("'prefixes' must be a single character value.", call. = FALSE)
    }
    #
    if( !(is.character(ignore.col) || is.null(ignore.col)) ){
        stop("'ignore.col' must be a character value or NULL.", call. = FALSE)
    }
    #
    # Subset by taking only taxonomy info if user want to remove the pattern
    # only from those. (Might be too restricting, e.g., if taxonomy columns are
    # not detected in previous steps. That is way the default is FALSE)
    ind <- rep(TRUE, ncol(feature_tab))
    if( only.taxa.col ){
        ind <- tolower(colnames(feature_tab)) %in% TAXONOMY_RANKS
    }
    # Also remove those columns that are specified to ignore (by default the
    # column that has unparsed data)
    ind <- ind & !colnames(feature_tab) %in% ignore.col
    temp <- feature_tab[, ind, drop = FALSE]
    # If there are columns left for removing the pattern
    if( ncol(temp) > 0 ){
        # Remove patterns
        temp <- lapply(
            temp, gsub, pattern = prefixes, replacement = "")
        temp <- as.data.frame(temp)
        # If cell had only prefix, it is now empty string. Convert to NA
        temp[ temp == "" ] <- NA
        # Combine table
        feature_tab[, ind] <- temp
    }
    # Ensure that returned value is DF
    feature_tab <- DataFrame(feature_tab)
    return(feature_tab)
}

# Find taxonomy rank based on prefixes. If found, return
# corresponding rank. Otherwise, return the original
# rank that is fed to function.
.replace_colnames_based_on_prefix <- function(colname, x){
    # Get column
    col <- x[ , colname]
    # List prefixes
    all_ranks <- .taxonomy_rank_prefixes
    prefixes <- paste0("^", all_ranks, "__")
    names(prefixes) <- names(all_ranks)
    # Find which prefix is found from each column value, if none.
    found_rank <- lapply(prefixes, FUN = function(pref){
        all(grepl(pattern = pref, col) | is.na(col)) && !all(is.na(col))
        })
    found_rank <- unlist(found_rank)
    # If only one prefix was found (like it should be), get the corresponding
    # rank name.
    if( sum(found_rank) == 1 ){
        colname <- names(prefixes)[found_rank]
    }
    return(colname)
}

# Detect and clean non wanted characters from Taxonomy data if needed.
#' @importFrom stringr str_extract_all
.detect_taxa_artifacts_and_clean <- function(
        x, pattern = "auto", ignore.col = "taxonomy_unparsed", ...) {
    #
    if( !.is_non_empty_character(pattern) ){
        stop("'pattern' must be a single character value.", call. = FALSE)
    }
    #
    if( !(is.character(ignore.col) || is.null(ignore.col)) ){
        stop("'ignore.col' must be a character value or NULL.", call. = FALSE)
    }
    #
    # Store rownames because they are lost during modifications
    row_names <- rownames(x)
    # Get those columns that are not modified
    not_mod <- x[ , colnames(x) %in% ignore.col, drop = FALSE]
    # Subset the data being modified so that it includes only those that user
    # wants to modify
    x <- x[ , !colnames(x) %in% ignore.col, drop = FALSE]
    # If there are still columns to clean after subsetting
    if( ncol(x) > 0 ){
        # Remove artifacts
        if( pattern == "auto" ){
            # Remove all but these characters
            pattern <- "[[:alnum:]]|-|_|\\[|\\]|,|;\\||[[:space:]]"
            x <- lapply(x, function(col){
                # Take all specified characters as a matrix where each column
                # is a character
                temp <- str_extract_all(
                    col, pattern = pattern, simplify = TRUE)
                # Collapse matrix to strings
                temp <- apply(temp, 1, paste, collapse = "")
                # Now NAs are converted into characters. Convert them back
                temp[ temp == "NA" ] <- NA
                # Convert also empty strings to NA
                temp[ temp == "" ] <- NA
                return(temp)
            })
        } else{
            # Remove pattern specified by user
            x <- lapply(x, gsub, pattern = pattern, replacement = "")
        }
    }
    # Ensure that the data is data.frame
    x <- as.data.frame(x)
    # Combine modified and non-modified data
    x <- cbind(x, not_mod)
    # Add rownames because they are dropped while removing artifacts
    rownames(x) <- row_names
    return(x)
}

# This function get's SE object as input and returns BIOM object.
.convert_se_to_biom <- function(x, assay.type){
    # Get assay, rowData and colData
    assay <- assay(x, assay.type)
    rd <- rowData(x)
    cd <- colData(x)
    # Convert rowData and colData to data.frame since biomformat cannot handle
    # DF.
    rd <- as.data.frame(rd)
    cd <- as.data.frame(cd)
    # Check if assay contains integers or floats. biom constructor
    # requires that information since the default value is "int".
    mat_type <- ifelse(all(assay %% 1 == 0), "int", "float")

    # Create argument list
    args <- list(data = assay, matrix_element_type = mat_type)
    # Add rowData and colData only if they contain information
    if( ncol(rd) != 0 ){
        args[["observation_metadata"]] <- rd
    }
    if( ncol(cd) != 0 ){
        args[["sample_metadata"]] <- cd
    }
    # Create biom
    biom <- do.call(biomformat::make_biom, args)
    return(biom)
}
