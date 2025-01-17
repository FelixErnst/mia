################################################################################
# This function gives user a message when the package is loaded into the session
.onAttach <- function(libname, pkgname) {
    pkg_version <- utils::packageDescription(pkgname, fields = "Version")
    msg <- paste0(
        "This is ", pkgname, " version ", pkg_version, "\n",
        "- Online documentation and vignettes: https://microbiome.github.io/",
        pkgname, "/",
        "\n",
        "- Online book 'Orchestrating Microbiome Analysis (OMA)': ",
        "https://microbiome.github.io/OMA/docs/devel/"
    )
    packageStartupMessage(msg)
}

################################################################################
# internal methods loaded from other packages

.get_mat_from_sce <- scater:::.get_mat_from_sce
.get_mat_for_reddim <- scater:::.get_mat_for_reddim

################################################################################
# integration with other packages

.require_package <- function(pkg){
    if(!requireNamespace(pkg, quietly = TRUE)){
    stop("'",pkg,"' package not found. Please install the '",pkg,"' package ",
        "to use this function.", call. = FALSE)
    }
}

################################################################################
# testing

.is_a_bool <- function(x){
    is.logical(x) && length(x) == 1L && !is.na(x)
}

.is_non_empty_character <- function(x){
    is.character(x) && all(nzchar(x))
}

.is_non_empty_string <- function(x){
    .is_non_empty_character(x) && length(x) == 1L
}

.is_a_string <- function(x){
    is.character(x) && length(x) == 1L
}

.is_integer <- function(x){
    is.numeric(x) && all(x%%1==0)
}

.is_an_integer <- function(x){
    .is_integer(x) && x%%1==0
}

.are_whole_numbers <- function(x){
    tol <- 100 * .Machine$double.eps
    abs(x - round(x)) <= tol && !is.infinite(x)
}

.is_a_numeric <- function(x){
    is.numeric(x) && length(x) == 1L
}

.is_numeric_string <- function(x){
    x <- as.character(x)
    suppressWarnings({x <- as.numeric(x)})
    !is.na(x)
}

.is_function <- function(x){
    typeof(x) == "closure" && is(x, "function")
}

.all_are_existing_files <- function(x){
    all(file.exists(x))
}

.get_name_in_parent <- function(x) {
    .safe_deparse(do.call(substitute, list(substitute(x), parent.frame())))
}

.safe_deparse <- function (expr, ...) {
    paste0(deparse(expr, width.cutoff = 500L, ...), collapse = "")
}

################################################################################
# checks

#' @importFrom SummarizedExperiment assays
.check_assay_present <- function(
        assay.type, x, name = .get_name_in_parent(assay.type)){
    if(!.is_non_empty_string(assay.type)){
        stop("'",name,"' must be a single non-empty character value.",
            call. = FALSE)
    }
    if(!(assay.type %in% names(assays(x)))){
        stop("'",name,"' must be a valid name of assays(x)", call. = FALSE)
    }
}

.check_rowTree_present <- function(
        tree.name, x, name = .get_name_in_parent(tree.name) ){
    if( !.is_non_empty_string(tree.name) ){
        stop("'", name, "' must be a single non-empty character value.",
            call. = FALSE)
    }
    if( !(tree.name %in% rowTreeNames(x)) ){
        stop("'", name, "' must specify a tree from 'rowTreeNames(x)'.",
            call. = FALSE)
    }
}

.check_colTree_present <- function(
        tree.name, x, name = .get_name_in_parent(tree.name) ){
    if( !.is_non_empty_string(tree.name) ){
        stop("'", name, "' must be a single non-empty character value.",
            call. = FALSE)
    }
    if( !(tree.name %in% colTreeNames(x)) ){
        stop("'", name, "' must specify a tree from 'colTreeNames(x)'.",
            call. = FALSE)
    }
}

# Check if alternative experiment can be found from altExp slot.
.check_altExp_present <- function(
        altexp, tse, altExpName = .get_name_in_parent(altexp),
        tse_name = .get_name_in_parent(tse), .disable.altexp = FALSE, ...){
    # Disable altExp if specified
    if( !.is_a_bool(.disable.altexp) ){
        stop("'.disable.altexp' must be TRUE or FALSE.", call. = FALSE)
    }
    if( .disable.altexp ){
        altexp <- NULL
    }
    # Check that altexp.name must be an integer or name
    if( !(.is_a_string(altexp) || .is_an_integer(altexp) || is.null(altexp)) ){
        stop(
            "'", altExpName, "' must be a string or an integer.", call. = FALSE)
    }
    # If is not NULL, but the object does not have altExp slot
    if( !is.null(altexp) && !is(tse, "SingleCellExperiment") ){
        stop(
            "'", altExpName, "', is specified but '", tse_name, "' does not ",
            "have altExp slot.", call. = FALSE)
    }
    # Then check that altExp can be found; name or index.
    if( !is.null(altexp) && !altexp %in% c(
            altExpNames(tse), seq_len(length(altExps(tse)))) ){
        stop(
            "'", altExpName, "', does not specify an experiment from altExp ",
            "slot of '", tse_name, "'.", call. = FALSE)
    }
}

# Check whether dimred is present in tse
.check_dimred_present <- function(dimred, x){
    specifies_index <- .is_integer(dimred) && dimred > 0 &&
        dimred <= length(reducedDims(x))
    specifies_name <- .is_a_string(dimred) && dimred %in% reducedDimNames(x)
    if( !specifies_index && !specifies_name ){
        stop("'dimred' must specify name or index from reducedDims(x).",
             call. = FALSE)
    }
    return(NULL)
}

# Check MARGIN parameters. Should be defining rows or columns.
.check_MARGIN <- function(MARGIN, name = .get_name_in_parent(MARGIN)) {
    # MARGIN must be one of the following options
    if( !(length(MARGIN) == 1L && tolower(MARGIN) %in% c(
            1, 2, "1", "2", "features", "samples", "columns", "col", "row",
            "rows", "cols")) ) {
        stop("'", name,"' must be 'rows' or 'cols'.", call. = FALSE)
    }
    # Convert MARGIN to numeric if it is not.
    MARGIN <- ifelse(tolower(MARGIN) %in% c(
        "samples", "columns", "col", 2, "cols"), 2, 1)
    return(MARGIN)
}

################################################################################
# Internal wrappers for getters

# Input: (Tree)SE
# Output: (Tree)SE
.check_and_get_altExp <- function(
        x, altexp = NULL, ...){
    # If altexp is specified, check and get it.
    # Otherwise return the original object
    if( !is.null(altexp) ){
        # Check altexp
        .check_altExp_present(altexp, x, ...)
        # Get altExp and return it
        x <- altExp(x, altexp)
    }
    return(x)
}

################################################################################
# Internal wrappers for setters

# This function adds values to colData (or rowData). The data must be in a list.
# Each element of list represent a column to be added to col/rowData.
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<-
#' @importFrom S4Vectors DataFrame
.add_values_to_colData <- function(
        x, values, name, altexp = NULL, MARGIN = default.MARGIN,
        default.MARGIN = 2, transpose.MARGIN = FALSE, colname = "name",
        ...){
    #
    if( !.is_a_string(colname) ){
        stop("'colname' must be a string.", call. = FALSE)
    }
    #
    # Check if altExp can be found
    .check_altExp_present(altexp, x)
    # Check that MARGIN is correct
    MARGIN <- .check_MARGIN(MARGIN)
    #
    # If transpose.MARGIN is TRUE, transpose MARGIN, i.e. 1 --> 2, and 2 --> 1.
    # In certain functions, values calculated by rows (MARGIN=1) are stored to
    # colData (MARGIN=2) and vice versa.
    if( transpose.MARGIN ){
        MARGIN <- ifelse(MARGIN == 1, 2, 1)
    }
    # converts each value:name pair into a DataFrame
    values <- mapply(
        function(value, n){
            value <- DataFrame(value)
            colnames(value)[1L] <- n
            if(ncol(value) > 1L){
                i <- seq.int(2,ncol(value))
                colnames(value)[i] <- paste0(n,"_",colnames(value)[i])
            }
            value
        },
        values,
        name)
    values <- do.call(cbind, values)
    
    # Based on MARGIN, get rowDatra or colData
    FUN <- switch(MARGIN, rowData, colData)
    # If altexp.name was not NULL, then we know that it specifies correctly
    # altExp from the slot. Take the colData/rowData from experiment..
    if( !is.null(altexp) ){
        cd <- FUN( altExp(x, altexp) )
    } else{
        cd <- FUN(x)
    }
    
    # check for duplicated values
    f <- colnames(cd) %in% colnames(values)
    FUN_name <- switch(MARGIN, "rowData", "colData")
    if(any(f)) {
        warning(
            "The following values are already present in `", FUN_name,
            "` and will be overwritten: '",
            paste(colnames(cd)[f], collapse = "', '"),
            "'. Consider using the '", colname,
            "' argument to specify alternative names.",
            call. = FALSE)
    }
    # Keep only unique values
    cd <- cbind( (cd)[!f], values )
    
    # Replace colData with new one
    x <- .add_to_coldata(x, cd, altexp = altexp, MARGIN = MARGIN)
    return(x)
}

# Get feature or sample metadata. Allow hidden usage of MARGIN and altExp.
#' @importFrom SummarizedExperiment rowData colData
.add_to_coldata <- function(
        x, cd, altexp = NULL, .disable.altexp = FALSE,
        MARGIN = default.MARGIN, default.MARGIN = 1, ...){
    #
    if( !.is_a_bool(.disable.altexp) ){
        stop("'.disable.altexp' must be TRUE or FALSE.", call. = FALSE)
    }
    # Check if altExp can be found
    .check_altExp_present(altexp, x, ...)
    # Check that MARGIN is correct
    MARGIN <- .check_MARGIN(MARGIN)
    # Based on MARGIN, add result to rowData or colData
    FUN <- switch(MARGIN, `rowData<-`, `colData<-`)
    # If altexp was specified, add result to altExp. Otherwise add it directly
    # to x.
    if( !is.null(altexp) && !.disable.altexp ){
        altExp(x, altexp) <- FUN( altExp(x, altexp), value = cd )
    } else{
        x <- FUN(x, value = cd)
    }
    return(x)
}

#' @importFrom S4Vectors metadata metadata<-
.add_values_to_metadata <- function(
        x, names, values, altexp = NULL, metadata.name = "name", ...){
    #
    if( !.is_a_string(metadata.name) ){
        stop("'metadata.name' must be a string.", call. = FALSE)
    }
    # Check if altExp can be found
    .check_altExp_present(altexp, x)
    #
    # Create a list and name elements
    add_metadata <- list(values)
    names(add_metadata) <- names
    # Get old metadata
    if( !is.null(altexp) ){
        old_metadata <- metadata( altExp(x, altexp) )
    } else{
        old_metadata <- metadata(x)
    }
    # Check if names match with elements that are already present
    f <- names(old_metadata) %in% names(add_metadata)
    if( any(f) ){
        warning(
            "The following values are already present in `metadata` and will ",
            "be overwritten: '",
            paste(names(old_metadata)[f], collapse = "', '"),
            "'. Consider using the '", metadata.name,
            "' argument to specify alternative ", "names.", call. = FALSE)
    }
    # keep only unique values
    add_metadata <- c( old_metadata[!f], add_metadata )
    # Add metadata to altExp or directly to x
    if( !is.null(altexp) ){
        metadata( altExp(x, altexp) ) <- add_metadata
    } else{
        metadata(x) <- add_metadata
    }
    return(x)
}

# This function can be used to add values to altExp
.add_to_altExps <- function(x, values, name = names(values), ...){
    # Check values
    if( !((is(values, "list") || is(values, "SimpleList")) &&
            length(values) > 0) ){
        stop("'values' must be non-empty list.", call. = FALSE)
    }
    # Check names
    if( !is.character(name) && length(name) > 1L ){
        stop("'name' must be a character value.", call. = FALSE)
    }
    # Names must match with list
    if( length(values) != length(name) ){
        stop("Lenght of 'name' must match with 'values'.", call. = FALSE)
    }
    #
    # If the object is SE, convert it to TreeSE
    if( !is(x, "SingleCellExperiment") ){
        x <- as(x, "TreeSummarizedExperiment")
        warning(
            "SummarizedExperiment does not have altExps slot. ",
            "Therefore, it is converted to TreeSummarizedExperiment.",
            call. = FALSE)
    }
    #
    # Add names to values
    names(values) <- name
    # Get altExps
    old_altexp <- altExps(x)
    # Check if names match with elements that are already present
    f <- names(old_altexp) %in% names(values)
    if( any(f) ){
        warning(
            "The following values are already present in `altExps` and will ",
            "be overwritten: '",
            paste(names(old_altexp)[f], collapse = "', '"),
            "'. Consider using the 'name' argument to specify alternative ",
            "names.", call. = FALSE)
    }
    # Keep only unique values
    values <- c( old_altexp[!f], values )
    # Add to altExps
    altExps(x) <- values
    return(x)
}

# This function can be used to add values to reducedDims
.add_values_to_reducedDims <- function(x, values, name, ...){
    # Check values
    if( !((is(values, "matrix") || is(values, "dist")) && length(values) > 0) ){
        stop("'values' must be a matrix.", call. = FALSE)
    }
    # Check names
    if( !.is_a_string(name) ){
        stop("'name' must be a character value.", call. = FALSE)
    }
    #
    # If the object is SE, convert it to TreeSE
    if( !is(x, "SingleCellExperiment") ){
        x <- as(x, "TreeSummarizedExperiment")
        warning(
            "SummarizedExperiment does not have reducedDims slot. ",
            "Therefore, it is converted to TreeSummarizedExperiment.",
            call. = FALSE)
    }
    if( !identical(rownames(as.matrix(values)), colnames(x)) ){
        stop("Rownames of the matrix should match with colnames(x).",
            " The result is not added to reducedDims.", call. = FALSE)
    }
    # Throw warning if values of reducedDim are overwritten
    if( name %in% names(reducedDims(x)) ){
        warning(
            "The following values are already present in `reducedDims` and", 
            " will be overwritten: '", name,
            "'. Consider using the 'name' argument to specify alternative ",
            "names.", call. = FALSE)
    }
    reducedDim(x, name) <- values
    return(x)
}

################################################################################
# Other common functions

# keep dimnames of feature table (assay) consistent with the meta data 
# of sample (colData) and feature (rowData)
.set_feature_tab_dimnames <- function(
        feature_tab, sample_meta, feature_meta, refseq = NULL,
        include.all = FALSE, ...) {
    #
    if( !.is_a_bool(include.all) ){
        stop("'include.all' must be TRUE or FALSE.", call. = FALSE)
    }
    # Sample and feature names must be present
    if( is.null(colnames(feature_tab)) || is.null(rownames(sample_meta)) ){
        stop("Sample ids must be present.", call. = FALSE)
    }
    if( is.null(rownames(feature_tab)) || is.null(rownames(feature_meta)) ){
        stop("Feature ids must be present.", call. = FALSE)
    }
    #
    # Metadata can include features that are not present in abundance table.
    # With 'include.all', it is possible to include them in the dataset.
    if( include.all ){
        samples <- union(colnames(feature_tab), rownames(sample_meta))
        features <- union(rownames(feature_tab), rownames(feature_meta))
        # Add additional columns and sort data
        feature_tab <- feature_tab[
            match(features, rownames(feature_tab)),
            match(samples, colnames(feature_tab)), drop = FALSE]
        sample_meta <- sample_meta[
            match(samples, rownames(sample_meta)), , drop = FALSE]
        feature_meta <- feature_meta[
            match(features, rownames(feature_meta)), , drop = FALSE]
        # Add correct names
        colnames(feature_tab) <- rownames(sample_meta) <- samples
        rownames(feature_tab) <- rownames(feature_meta) <- features
    }
    
    # The abundance table takes the precedence, and metadata is modified
    # accordingly. Give warning if there are samples or features that are not
    # included in the metadata.
    not_found <- sum(!colnames(feature_tab) %in% rownames(sample_meta))
    if( not_found != 0 ){
        warning("The dataset includes ", not_found, " samples that do not ",
                "have metadata. Please check for errors.", call. = FALSE)
    }
    not_found <- sum(!rownames(feature_tab) %in% rownames(feature_meta))
    if( not_found != 0 ){
        warning("The dataset includes ", not_found, " features that are not ",
                "included in the taxonomy table. Please check for errors.",
                call. = FALSE)
    }
    # If the metadata includes samples or features that will be removed give
    # warning for the user.
    not_found <- sum(!rownames(sample_meta) %in% colnames(feature_tab))
    if( not_found != 0 ){
        warning("The sample metadata includes ", not_found, " samples that ",
                "are not included in the abundance table and thus being ",
                "removed. Please check for errors.", call. = FALSE)
    }
    not_found <- sum(!rownames(feature_meta) %in% rownames(feature_tab))
    if( not_found != 0 ){
        warning("The taxonomy table includes ", not_found, " features that ",
                "are not present in the abundance table and thus being ",
                "removed. Please check for errors.", call. = FALSE)
    }
    
    # We order the metadata based on abundance table. Moreover, we subset
    # the metadata to match with abundance table if there are additional data.
    ind <- match(colnames(feature_tab), rownames(sample_meta))
    sample_meta <- sample_meta[ind, , drop = FALSE]
    rownames(sample_meta) <- colnames(feature_tab)
    ind <- match(rownames(feature_tab), rownames(feature_meta))
    feature_meta <- feature_meta[ind, , drop = FALSE]
    rownames(feature_meta) <- rownames(feature_tab)
    # Reference sequences are optional and all the features must have sequences.
    # This is because DNAStringSet object cannot have empty element for
    # features that were not included.
    if( !is.null(refseq) && all(rownames(feature_tab) %in% names(refseq)) ){
        refseq <- refseq[ rownames(feature_tab) ]
    } else if( !is.null(refseq) ){
        warning("Reference sequences are incompatible with the data.",
                call. = FALSE)
        refseq <- NULL
    } else{
        refseq <- NULL
    }
    # Return a list of these tables
    data_list <- list(
        assay = feature_tab,
        rowData = feature_meta,
        colData = sample_meta,
        referenceSeq = refseq
    )
    return(data_list)
}

#' Parse taxa in different taxonomic levels
#' @param taxa_tab `data.frame` object.
#' 
#' @param sep character string containing a regular expression, separator
#'  between different taxonomic levels, defaults to one compatible with both
#'  GreenGenes and SILVA `; |;"`.
#'  
#' @param col.name a single \code{character} value defining the column of
#' taxa_tab that includes taxonomical information.
#'  
#' @param prefix.rm {\code{TRUE} or \code{FALSE}: Should 
#'  taxonomic prefixes be removed? (default: \code{prefix.rm = FALSE})}
#'  
#' @return  a `data.frame`.
#' @keywords internal
#' @importFrom IRanges CharacterList IntegerList
#' @importFrom S4Vectors DataFrame
#' @noRd
.parse_taxonomy <- function(
    taxa_tab, sep = "; |;", col.name = column_name, column_name = "Taxon",
    remove.prefix = prefix.rm, prefix.rm = removeTaxaPrefixes,
    removeTaxaPrefixes = FALSE, ...) {
    ############################### Input check ################################
    # Check sep
    if(!.is_non_empty_string(sep)){
        stop("'sep' must be a single character value.", call. = FALSE)
    }
    # Check col.name
    if( !(.is_non_empty_string(col.name) && col.name %in% colnames(taxa_tab)) ){
        stop("'col.name' must be a single character value defining column ",
            "that includes information about taxonomic levels.", call. = FALSE)
    }
    # Check remove.prefix
    if(!.is_a_bool(remove.prefix)){
        stop("'remove.prefix' must be TRUE or FALSE.", call. = FALSE)
    }
    ############################## Input check end #############################
    
    #  work with any combination of taxonomic ranks available
    all_ranks <- .taxonomy_rank_prefixes
    all_prefixes <- paste0(all_ranks, "__")
    names(all_prefixes) <- names(all_ranks)
    
    # split the taxa strings
    taxa_split <- CharacterList(strsplit(taxa_tab[, col.name],sep))
    # extract present prefixes
    taxa_prefixes <- lapply(
        taxa_split, gsub, pattern = "(^[a-zA-Z]+__).*", replacement = "\\1")
    # match them to the order given by present_prefixes
    taxa_prefixes_match <- lapply(taxa_prefixes, match, x = all_prefixes)
    taxa_prefixes_match <- IntegerList(taxa_prefixes_match)
    # get the taxa values without prefixes
    if(remove.prefix){
        pattern <- paste0("(", paste0(all_ranks, collapse = "|"), ")__")
        taxa_split <- lapply(
            taxa_split, gsub, pattern = pattern, replacement = "")
        taxa_split <- CharacterList(taxa_split)
    }
    # extract by order matches
    taxa_split <- taxa_split[taxa_prefixes_match]
    #
    if(length(unique(lengths(taxa_split))) != 1L){
        stop("Something went wrong while splitting taxonomic ",
            "levels. Please check that 'sep' is correct.", call. = FALSE)
    }
    taxa_tab <- DataFrame(as.matrix(taxa_split))
    colnames(taxa_tab) <- names(all_ranks)
    
    # Subset columns so that include only those columns that have some
    # information
    non_empty <- colSums(is.na(taxa_tab)) != nrow(taxa_tab)
    taxa_tab <- taxa_tab[ , non_empty, drop = FALSE]
    
    return(taxa_tab)
}

################################################################################
# internal wrappers for agglomerateByRank/agglomerateByVariable
.merge_features <- function(x, merge.by = rank, rank = NULL, ...) {
    # Check if merge.by parameter belongs to taxonomyRanks
    if( .is_a_string(merge.by) && merge.by %in% taxonomyRanks(x) ){
        # Merge using agglomerateByRank
        x <- agglomerateByRank(x, rank = merge.by, ...)
    } else if( !is.null(merge.by) ){
        # Merge using agglomerateByVariable
        x <- agglomerateByVariable(x, by = "rows", group = merge.by, ...)
    }
    return(x)
}

################################################################################
# This function sets taxonomy ranks based on rowData of TreeSE. With this,
# user can automatically set ranks based on imported data.
.set_ranks_based_on_rowdata <- function(
        tse, set.ranks = FALSE, verbose = TRUE,
        ignore.col = "taxonomy_unparsed", ...){
    #
    if( !.is_a_bool(set.ranks) ){
        stop("'set.ranks' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    if( !.is_a_bool(verbose) ){
        stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    if( !(is.character(ignore.col) || is.null(ignore.col)) ){
        stop("'ignore.col' must be a character value or NULL.", call. = FALSE)
    }
    #
    # Get ranks from rowData
    rd <- rowData(tse)
    # Remove those columns that are ignored. By default, the column
    # containing unparsed taxonomy.
    rd <- rd[ , !colnames(rd) %in% ignore.col, drop = FALSE]
    # Ranks must be character columns
    is_char <- lapply(
        rd, function(x) is.character(x) || is.factor(x))
    is_char <- unlist(is_char)
    rd <- rd[ , is_char, drop = FALSE]
    # If user wants to set ranks and there are ranks after filtering out
    # those columns that are not characters.
    if( set.ranks && ncol(rd) > 0L ){
        # Finally, set ranks and give message
        ranks <- colnames(rd)
        temp <- setTaxonomyRanks(ranks)
        if( verbose ){
            message(
                "TAXONOMY_RANKS set to: '",
                paste0(ranks, collapse = "', '"), "'")
        }
    }
    # If user wanted to set ranks but there were no suitable columns in rowData,
    # give warning
    if( set.ranks && ncol(rd) == 0L ){
        warning(
            "Ranks cannot be set. rowData(x) does not include columns ",
            "specifying character values.", call. = FALSE)
    }
    return(NULL)
}

################################################################################
# This function converts vector of character values to capitalized.
.capitalize <- function(x){
    paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}
