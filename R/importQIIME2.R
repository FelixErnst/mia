#' Import QIIME2 results to \code{TreeSummarizedExperiment}
#'
#' Results exported from QIMME2 can be imported as a
#' \code{TreeSummarizedExperiment} using \code{importQIIME2}. Except for the
#' \code{assay.file}, the other data types, \code{row.file},
#' \code{refseq.file} and \code{tree.file}, are optional, but are highly
#' encouraged to be provided.
#'
#' @param assay.file \code{Character scalar}. Defines the file
#'   path of the feature table to be imported.
#' 
#' @param featureTableFile Deprecated. use \code{assay.file} instead.
#'
#' @param row.file \code{Character scalar} or \code{NULL}. Defines the file
#'   path of the taxonomy table to be imported. (default:
#'   \code{NULL}).
#' 
#' @param taxonomyTableFile Deprecated. use \code{row.file} instead.
#'
#' @param col.file \code{Character scalar} or \code{NULL}. Defines the file path
#'   of the sample metadata to be imported. The file has to be in tsv format.
#'   (Default: \code{NULL}).
#' 
#' @param sampleMetaFile Deprecated. Use \code{col.file} instead.
#'
#' @param as.refseq \code{Logical scalar} or \code{NULL}. Should the feature
#'   names of the feature table be regarded as reference sequences? This setting
#'   will be disregarded, if \code{refseq.file} is not \code{NULL}. If the
#'   feature names do not contain valid DNA characters only, the reference
#'   sequences will not be set.
#' 
#' @param featureNamesAsRefSeq Deprecated. Use \code{as.refseq} instead.
#'
#' @param refseq.file \code{Character scalar} or \code{NULL}. Defines the file
#' path of the reference sequences for each feature. (Default: \code{NULL}).
#' 
#' @param refSeqFile Deprecated. Use \code{refseq.file} instead.
#'
#' @param tree.file \code{Character scalar}. Defines the file path of
#'   the phylogenetic tree. (Default: \code{NULL}).
#' 
#' @param phyTreeFile Deprecated. Use \code{tree.file} instead.
#'
#' @param ... additional arguments:
#' \itemize{
#'   \item \code{temp.dir}: the temporary directory used for decompressing the
#'     data. (default: \code{tempdir()})
#'   \item \code{prefix.rm}: \code{TRUE} or \code{FALSE}: Should
#'     taxonomic prefixes be removed? (default:
#'     \code{prefix.rm = FALSE})
#' }
#'
#' @details
#' Both arguments \code{as.refseq} and \code{refseq.file} can be used
#' to define reference sequences of features. \code{as.refseq} is
#' only taken into account, if \code{refseq.file} is \code{NULL}. No reference
#' sequences are tried to be created, if \code{featureNameAsRefSeq} is
#' \code{FALSE} and \code{refseq.file} is \code{NULL}.
#'
#' @return  A
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object
#'
#' @name importQIIME2
#' @seealso
#' \code{\link[=convertFromPhyloseq]{convertFromPhyloseq}}
#' \code{\link[=convertFromBIOM]{convertFromBIOM}}
#' \code{\link[=convertFromDADA2]{convertFromDADA2}}
#' \code{\link[=importMothur]{importMothur}}
#'
#' @export
#' 
#' @author Yang Cao
#'
#' @references
#' Bolyen E et al. 2019: Reproducible, interactive, scalable and extensible
#' microbiome data science using QIIME 2. Nature Biotechnology 37: 852–857.
#' \url{https://doi.org/10.1038/s41587-019-0209-9}
#'
#' \url{https://qiime2.org}
#'
#' @examples
#' assay.file <- system.file("extdata", "table.qza", package = "mia")
#' row.file <- system.file("extdata", "taxonomy.qza", package = "mia")
#' col.file <- system.file("extdata", "sample-metadata.tsv", package = "mia")
#' tree.file <- system.file("extdata", "tree.qza", package = "mia")
#' refseq.file <- system.file("extdata", "refseq.qza", package = "mia")
#' tse <- importQIIME2(
#'   assay.file = assay.file,
#'   row.file = row.file,
#'   col.file = col.file,
#'   refseq.file = refseq.file,
#'   tree.file = tree.file
#' )
#'
#' tse

#' @importFrom S4Vectors make_zero_col_DFrame
importQIIME2 <- function(assay.file = featureTableFile,
                            featureTableFile,
                            row.file = taxonomyTableFile,
                            taxonomyTableFile = NULL,
                            col.file = sampleMetaFile,
                            sampleMetaFile = NULL,
                            as.refseq = featureNamesAsRefSeq,
                            featureNamesAsRefSeq = TRUE,
                            refseq.file = refSeqFile,
                            refSeqFile = NULL,
                            tree.file = phyTreeFile,
                            phyTreeFile = NULL,
                            ...) {
    .require_package("yaml")
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
    if(!.is_a_bool(as.refseq)){
        stop("'as.refseq' must be TRUE or FALSE.", call. = FALSE)
    }
    if(!is.null(refseq.file) && !.is_non_empty_string(refseq.file)){
        stop("'refseq.file' must be a single character value or NULL.",
            call. = FALSE)
    }
    if(!is.null(tree.file) && !.is_non_empty_string(tree.file)){
        stop("'tree.file' must be a single character value or NULL.",
            call. = FALSE)
    }
    #

    feature_tab <- importQZA(assay.file, ...)

    if (!is.null(row.file)) {
        taxa_tab <- importQZA(row.file, ...)
        taxa_tab <- .subset_taxa_in_feature(taxa_tab, feature_tab)
    } else {
        taxa_tab <- S4Vectors::make_zero_col_DFrame(nrow(feature_tab))
        rownames(taxa_tab) <- rownames(feature_tab)
    }

    if (!is.null(col.file)) {
        sample_meta <- .read_q2sample_meta(col.file)
    } else {
        sample_meta <- S4Vectors::make_zero_col_DFrame(ncol(feature_tab))
        rownames(sample_meta) <- colnames(feature_tab)
    }

    if (!is.null(tree.file)) {
        tree <- importQZA(tree.file, ...)
    } else {
        tree <- NULL
    }

    # if row.names(feature_tab) is a DNA sequence,  set it as refseq
    if (!is.null(refseq.file)){
        refseq <- importQZA(refseq.file, ...)
    } else if (as.refseq) {
        refseq <- .rownames_as_dna_seq(rownames(feature_tab))
    } else {
        refseq <- NULL
    }
    
    feature_tab <- .set_feature_tab_dimnames(feature_tab, sample_meta, taxa_tab)
    TreeSummarizedExperiment(
        assays = S4Vectors::SimpleList(counts = feature_tab),
        rowData = taxa_tab,
        colData = sample_meta,
        rowTree = tree,
        referenceSeq = refseq
    )
}

#' Read the qza file output from QIIME2
#'
#' Import the QIIME2 artifacts to R.
#'
#' @param file character, path of the input qza file. Only files in format of
#'   `BIOMV210DirFmt` (feature table), `TSVTaxonomyDirectoryFormat` (taxonomic
#'   table), `NewickDirectoryFormat` (phylogenetic tree ) and
#'   `DNASequencesDirectoryFormat` (representative sequences) are supported
#'    right now.
#' @param temp.dir character, a temporary directory in which the qza file will
#' be decompressed to, default `tempdir()`.
#' 
#' @param temp Deprecated. Use \code{temp.dir} instead.
#' 
#' @return `matrix` object for feature table, `DataFrame` for taxonomic table,
#'   [`ape::phylo`] object for phylogenetic tree,
#'   [`Biostrings::DNAStringSet-class`] for representative sequences of taxa.
#'   
#' @name importQIIME2
#' @export
#'
#' @examples 
#' # Read individual files
#' assay.file <- system.file("extdata", "table.qza", package = "mia")
#' row.file <- system.file("extdata", "taxonomy.qza", package = "mia")
#' col.file <- system.file("extdata", "sample-metadata.tsv", package = "mia")
#' 
#' assay <- importQZA(assay.file)
#' rowdata <- importQZA(row.file, prefix.rm = TRUE)
#' coldata <- read.table(col.file, header = TRUE, sep = "\t", comment.char = "")
#' 
#' # Assign rownames 
#' rownames(coldata) <- coldata[, 1]
#' coldata[, 1] <- NULL
#' 
#' # Order coldata based on assay
#' coldata <- coldata[match(colnames(assay), rownames(coldata)), ]
#' 
#' # Create SE from individual files
#' se <- SummarizedExperiment(
#'     assays = list(assay), rowData = rowdata, colData = coldata)
#' se
#' 
#' @importFrom utils unzip
#' @importFrom ape read.tree
#' @importFrom Biostrings readDNAStringSet
importQZA <- function(file, temp.dir = temp, temp = tempdir(), ...) {
    if (!file.exists(file)) {
        stop(file, " does not exist", call. = FALSE)
    }
    if (.get_ext(file) != "qza") {
        stop("The input '", file, "' must be in `qza` format ",
            "(QIIME2 Artifact)", call. = FALSE)
    }

    unzipped_file <- unzip(file, exdir = temp.dir)
    on.exit(unlink(c(unzipped_file,unique(dirname(unzipped_file))),
                    recursive = TRUE))
    meta_file <- grep("metadata.yaml", unzipped_file, value = TRUE)
    metadata <- yaml::read_yaml(meta_file[1])
    uuid <- metadata$uuid

    format <- metadata$format
    # support for multiple BIOM formats:  V100, V210
    if (grepl("BIOMV", format)) {
        format <- "BIOMV"
    }

    format_files <- c(
        "feature-table.biom", "taxonomy.tsv",
        "tree.nwk", "dna-sequences.fasta"
    )
    formats <- c(
        "BIOMV", "TSVTaxonomyDirectoryFormat",
        "NewickDirectoryFormat", "DNASequencesDirectoryFormat"
    )
    file <- file.path(
        temp.dir, uuid, "data", format_files[match(format, formats)])

    res <- switch (
        format,
        BIOMV = .read_q2biom(file),
        TSVTaxonomyDirectoryFormat = .read_q2taxa(file, ...),
        NewickDirectoryFormat = read.tree(file),
        DNASequencesDirectoryFormat = readDNAStringSet(file),
        stop(
            "Only files in format of 'BIOMV210DirFmt', ",
            "'TSVTaxonomyDirectoryFormat', NewickDirectoryFormat' and ",
            "'DNASequencesDirectoryFormat' are supported.",
            call. = FALSE
        )
    )
    res
}

#' Read QIIME2 feature table
#' @param file character, file name of the biom file.
#' @noRd
.read_q2biom <- function(file) {
    .require_package("biomformat")
    biomobj <- suppressWarnings(biomformat::read_biom(file))
    feature_tab <- as(biomformat::biom_data(biomobj),"matrix")

    feature_tab
}

#' Read QIIME2 taxa file
#' @keywords internal
#' @importFrom utils read.table
#' @noRd
.read_q2taxa <- function(file, ...) {
    taxa_tab <- utils::read.table(file, sep = '\t', header = TRUE)
    
    confidence <- NULL
    featureID <- NULL
    # make sure confidence in numeric
    if ("Confidence" %in% colnames(taxa_tab)) {
        confidence <- as.numeric(taxa_tab[,"Confidence"])
    }
    if("Feature.ID" %in% names(taxa_tab)) {
        featureID <- taxa_tab[,"Feature.ID"]
    }
    
    taxa_tab <- .parse_taxonomy(
        taxa_tab, sep = "; |;", column_name = "Taxon", ...)
    
    rownames(taxa_tab) <- featureID
    taxa_tab$Confidence <- confidence
    
    taxa_tab
}

#' Read QIIME2 sample meta data file
#' @keywords internal
#' @importFrom utils read.table
#' @noRd
.read_q2sample_meta <- function(file) {
    sam <- read.table(file = file, header = TRUE, sep = "\t", comment.char = "")
    rownames(sam) <- as.character(sam[, 1])

    # Find if there is #q2:types row, and store its index
    idx <- which(sam == "#q2:types", arr.ind = TRUE)

    # If the length is over zero, "#q2:types" row was found. Then it is removed.
    if(!(length(idx)==0)){
        sam <- sam[-idx[, "row"], ]
    }

    S4Vectors::DataFrame(sam)
}


#' Subset taxa according to the taxa in feature table
#' @keywords internal
#' @noRd
.subset_taxa_in_feature <- function(taxa_tab, feature_tab) {
    idx <- match( rownames(feature_tab), rownames(taxa_tab) )
    taxa_tab <- taxa_tab[idx, , drop = FALSE]

    taxa_tab
}

#' check the row.names of feature table is DNA sequence or not
#' @keywords internal
#' @importFrom Biostrings DNAStringSet
#' @noRd
.rownames_as_dna_seq <- function(seq){
    names(seq) <- paste0("seq_", seq_along(seq))
    seq <- try({DNAStringSet(seq)}, silent = TRUE)
    if (is(seq, "try-error")) {
        return(NULL)
    }

    seq
}

#' extract file extension
#' @noRd
.get_ext <- function(file) {
    ex <- strsplit(basename(file), split = ".", fixed = TRUE)[[1]]
    ex[length(ex)]
}
