# This function can be used to unify the group id vector. It can be any
# kind of vector, but this converts it to factor.
.norm_f <- function(
        i, f, dim.type = c("rows","columns"), empty.rm = FALSE, ...){
    if(!.is_a_bool(empty.rm)){
        stop("'empty.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    dim.type <- match.arg(dim.type)
    if(!is.character(f) && !is.factor(f)){
        stop("'f' must be a factor or character vector coercible to a ",
            "meaningful factor.",
            call. = FALSE)
    }
    if(i != length(f)){
        stop("'f' must have the same number of ",dim.type," as 'x'",
            call. = FALSE)
    }
    # This is done otherwise we lose NA values
    if( !empty.rm && any(is.na(f)) ){
        f <- as.character(f)
        f[ is.na(f) ] <- "NA"
    }
    if(is.character(f)){
        f <- factor(f)
    }
    f
}

# When we merge rows or columns, first member of group is kept as default
# (in colData or rowData). This function controls this and allows user to
# specify some other element than the first one.
.norm_archetype <- function(f, archetype){
    if(length(archetype) > 1L){
        if(length(levels(f)) != length(archetype)){
            stop("length of 'archetype' must have the same length as ",
                "levels('f')",
                call. = FALSE)
        }
    }
    f_table <- table(f)
    if(!is.null(names(archetype))){
        if(anyNA(names(archetype)) || anyDuplicated(names(archetype))){
            stop("If 'archetype' is named, names must be non-NA and unqiue.",
                call. = FALSE)
        }
        archetype <- archetype[names(f_table)]
    }
    if(any(f_table < archetype)){
        stop("'archetype' out of bounds for some levels of 'f'. The maximum of",
            " 'archetype' is defined as table('f')", call. = FALSE)
    }
    if(length(archetype) == 1L){
        archetype <- rep(archetype,length(levels(f)))
    }
    archetype
}

# This function returns the index/position of rows/columns that are kept
# after merging.
#' @importFrom S4Vectors splitAsList
.get_element_pos <- function(f, archetype){
    archetype <- as.list(archetype)
    f_pos <- seq.int(1L, length(f))
    f_pos_split <- S4Vectors::splitAsList(f_pos, f)
    f_pos <- unlist(f_pos_split[archetype])
    f_pos
}

# This function merges assays and row/colData.
#' @importFrom S4Vectors SimpleList
#' @importFrom scuttle sumCountsAcrossFeatures
.merge_rows_or_cols <- function(
        x, f, by, archetype = 1L, average = FALSE, BPPARAM = SerialParam(),
        check.assays = TRUE, na.rm = FALSE, ...){
    # input check
    if( !.is_a_bool(average) ){
        stop("'average' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(check.assays) ){
        stop("'check.assays' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_bool(na.rm) ){
        stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Get correct functions based on whether we agglomerate rows or cols
    rowData_FUN <- switch(by, rowData, colData)
    nrow_FUN <- switch(by, nrow, ncol)
    rownames_FUN <- switch(by, rownames, colnames)
    rownames_ass_FUN <- switch(by, `rownames<-`, `colnames<-`)
    # If user specified column name from row/colData, get the values
    if( .is_a_string(f) && f %in% colnames(rowData_FUN(x)) ){
        f <- rowData_FUN(x)[[ f ]]
    }
    # Check that the group ID vector is specifying groups for each element
    f <- .norm_f(nrow_FUN(x), f, ...)
    # If the data is already agglomerated at each group
    if(length(levels(f)) == nrow_FUN(x)){
        return(x)
    }
    # In merging, first element of certain group is kept by default. archetype,
    # can control this behavior; it can specify the preserved rows for every
    # group or index.
    archetype <- .norm_archetype(f, archetype)
    
    # Get assays
    assays <- assays(x)
    # We check whether the assays include values that cannot be summed. For
    # instance, summing negative values do not make sense.
    if( check.assays ){
        temp <- lapply(seq_len(length(assays)), function(i)
            .check_assays_for_merge(names(assays)[[i]], assays[[i]]))
    }
    
    # Transpose if we are merging columns
    if( by == 2L ){
        assays <- lapply(assays, function(mat) t(mat))
    }
    # Get the aggregation function based on whether user wants to exclude NAs
    # and if there are any NAs. scuttle::sumCountsAcrossFeatures cannot handle
    # NAs so if user wants to exclude them, we use own implementation.
    FUN <- if( na.rm && anyNA(assays[[1]])) .sum_counts_accross_features_na else
        sumCountsAcrossFeatures
    # Agglomerate assays
    assays <- lapply(assays, FUN, average = average, ids = f, BPPARAM = BPPARAM)
    # Transpose back to original orientation
    if( by == 2L ){
        assays <- lapply(assays, function(mat) t(mat))
    }
    # Convert to SimpleList
    assays <- assays |> SimpleList()
    
    # Now we have agglomerated assays, but TreeSE has still the original form.
    # We take specified rows/columns from the TreeSE.
    idx <- .get_element_pos(f, archetype = archetype)
    if( by == 1L ){
        x <- x[idx, ]
    } else{
        x <- x[ , idx]
    }
    
    # Add assays back to TreeSE
    assays(x, withDimnames = FALSE) <- assays
    # Change row/colnames. Currently, they have same names as in original data
    # but just certain rows. Change them to represent groups
    x <- rownames_ass_FUN(x, rownames_FUN(assays[[1]]))
    return(x)
}

# This function works similarly to scuttle::sumCountsAcrossFeatures but this
# excludes NAs from the data. The scuttle function cannot handle NAs.
#' @importFrom DelayedArray DelayedArray type rowsum
.sum_counts_accross_features_na <- function(x, average, ids, ...){
    # Which cell is not NA?
    is_not_na <- !is.na(x)
    type(is_not_na) <- "integer"
    # Aggregate data to certain groups
    x <- rowsum(x, ids, na.rm = TRUE)
    # Calculate average if specified
    if( average ){
        x <- x/rowsum(is_not_na, ids)
    }
    return(x)
}

# This functions checks if assay has negative or binary values. It does not
# make sense to sum them, so we give warning to user.
.check_assays_for_merge <- function(assay.type, assay){
    # Check if assays include binary or negative values
    if( all(assay == 0 | assay == 1) ){
        warning("'", assay.type, "'", " includes binary values.",
                "\nAgglomeration of it might lead to meaningless values.",
                "\nCheck the assay, and consider doing transformation again",
                "manually with agglomerated data.",
                call. = FALSE)
    }
    if( !all( assay >= 0 | is.na(assay) ) ){
        warning("'", assay.type, "'", " includes negative values.",
                "\nAgglomeration of it might lead to meaningless values.",
                "\nCheck the assay, and consider doing transformation again",
                "manually with agglomerated data.",
                call. = FALSE)
    }
    return(assay)
}

#' @importFrom Biostrings DNAStringSetList
.merge_refseq_list <- function(sequences_list, f, names, ...){
    threshold <- list(...)[["threshold"]]
    if(is.null(threshold)){
        threshold <- 0.05
    }
    if(!is(sequences_list,"DNAStringSetList")){
        return(.merge_refseq(sequences_list, f, names, threshold))
    }
    names <- names(sequences_list)
    seqs <- DNAStringSetList(lapply(sequences_list, .merge_refseq, f, names,
                                    threshold))
    names(seqs) <- names
    seqs
}

#' @importFrom Biostrings DNAStringSetList
#' @importFrom DECIPHER ConsensusSequence
.merge_refseq <- function(sequences, f, names, threshold){
    sequences <- split(sequences,f)
    seq <- unlist(DNAStringSetList(lapply(sequences, ConsensusSequence,
                                            threshold = threshold)))
    seq
}

.merge_rows_TSE <- function(x, f, archetype = 1L, update.tree = FALSE,
    update.refseq = mergeRefSeq, mergeRefSeq = FALSE, ...){
    # input check
    if(!.is_a_bool(update.tree)){
        stop("'update.tree' must be TRUE or FALSE.", call. = FALSE)
    }
    if(!.is_a_bool(update.refseq)){
        stop("'update.refseq' must be TRUE or FALSE.", call. = FALSE)
    }
    # for optionally merging referenceSeq
    refSeq <- NULL
    if(update.refseq){
        refSeq <- referenceSeq(x)
    }
    #
    x <- .merge_rows_or_cols(x, f, by = 1L, archetype = 1L, ...)
    # optionally merge rowTree
    if( update.tree ){
        x <- .agglomerate_trees(x, 1, ...)
    }
    # optionally merge referenceSeq
    if(!is.null(refSeq)){
        referenceSeq(x) <- .merge_refseq_list(refSeq, f, rownames(x), ...)
    }
    x
}

.merge_cols_TSE <- function(x, f, archetype = 1L, update.tree = FALSE, ...){
    # input check
    if(!.is_a_bool(update.tree)){
        stop("'update.tree' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    x <- .merge_rows_or_cols(x, f, by = 2L, archetype = 1L, ...)
    # optionally merge colTree
    if( update.tree ){
        x <- .agglomerate_trees(x, 2, ...)
    }
    return(x)
}
