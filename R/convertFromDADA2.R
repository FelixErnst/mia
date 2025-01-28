#' Create a \code{TreeSummarizedExperiment} object from \sQuote{DADA2} results
#'
#' @param ... Additional arguments. For \code{convertFromDADA2}, see
#' \code{mergePairs} function for more details.
#'
#' @details
#'
#' \code{convertFromDADA2} is a wrapper for the
#' \code{mergePairs} function from the \code{dada2} package.
#' A count matrix is constructed via
#' \code{makeSequenceTable(mergePairs(...))} and rownames are dynamically
#' created as \code{ASV(N)} with \code{N} from 1 to \code{nrow} of the count
#' tables. The colnames and rownames from the output of \code{makeSequenceTable}
#' are stored as \code{colnames} and in the \code{referenceSeq} slot of the
#' \code{TreeSummarizedExperiment}, respectively.
#'
#' @return \code{convertFromDADA2} returns an object of class
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @importFrom S4Vectors SimpleList
#' @importFrom Biostrings DNAStringSet
#'
#' @rdname convertFromDADA2
#'
#' @export
#'
#' @examples
#'
#' ### Coerce DADA2 results to a TreeSE object
#' if(requireNamespace("dada2")) {
#'   fnF <- system.file("extdata", "sam1F.fastq.gz", package="dada2")
#'   fnR = system.file("extdata", "sam1R.fastq.gz", package="dada2")
#'   dadaF <- dada2::dada(fnF, selfConsist=TRUE)
#'   dadaR <- dada2::dada(fnR, selfConsist=TRUE)
#'
#'   tse <- convertFromDADA2(dadaF, fnF, dadaR, fnR)
#'   tse
#' }
#' @importFrom stringr str_pad
convertFromDADA2 <- function(...) {
    # input checks
    .require_package("dada2")
    #
    mergers <- dada2::mergePairs(...)
    seqtab <- dada2::makeSequenceTable(mergers)
    seqtab <- t(seqtab)
    # generate row and col names
    rName <- paste0("ASV",
        str_pad(seq.int(1L,nrow(seqtab)), nchar(nrow(seqtab)) + 1L, pad="0"))
    cName <- colnames(seqtab)
    # retrieve count data and reference sequence
    assays <- S4Vectors::SimpleList(counts = unname(seqtab))
    refseq <- Biostrings::DNAStringSet(rownames(seqtab))
    # construct ME an name rows and cols
    output <- TreeSummarizedExperiment(assays = assays,
                                        referenceSeq = refseq)
    colnames(output) <- cName
    rownames(output) <- rName
    output
}
