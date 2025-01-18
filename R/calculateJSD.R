# written by Susan Holmes \email{susan@@stat.stanford.edu}.
# Adapted for phyloseq by Paul J. McMurdie.
# Adapted for mia by Felix G.M. Ernst
#' @importFrom DelayedArray rowSums
.JSD <- function(x, y){
    # Function to compute Shannon-Jensen Divergence
    # x and y are the frequencies for the same p categories
    # Assumes relative abundance transformation already happened
    # (for efficiency)

    # Define the mean point
    m <- (x+y)/2
    # Define each samples component
    P1 <- x*log(x/m)
    P2 <- y*log(y/m)
    # In the case of zeroes entries log is undefined, JSD is defined as zero
    P1[!is.finite(P1)] <- 0
    P2[!is.finite(P2)] <- 0
    d <- (P1+P2)/2
    return(rowSums(d, na.rm = TRUE))
}

#' @importFrom utils combn
#' @importFrom stats as.dist
#' @importFrom BiocParallel SerialParam register bplapply bpisup bpstart bpstop
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
#' 
.get_jsd <- function(x, BPPARAM = SerialParam(), chunkSize = nrow(x), ...){
    # input check
    if(is.null(rownames(x))){
        rownames(x) <- seq_len(nrow(x))
    }
    if(missing(chunkSize) || is.na(chunkSize) || is.null(chunkSize) ||
        !is.integer(chunkSize)){
        chunkSize <- nrow(x)
    } else if(length(chunkSize) != 1L) {
        chunkSize <- chunkSize[1L]
    }
    #
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }
    # Coerce to relative abundance by sample (row)
    x <- sweep(x, 1L, rowSums(x), "/")
    # create N x 2 matrix of all pairwise combinations of samples.
    spn <- utils::combn(rownames(x), 2, simplify = TRUE)
    #
    N <- ncol(spn)
    f <- ceiling(seq_len(N)/chunkSize)
    A <- split(spn[1L,], f)
    B <- split(spn[2L,], f)
    FUN <- function(X, a, b){
        .JSD(X[a,,drop=FALSE], X[b,,drop=FALSE])
    }
    distlist <- BiocParallel::bpmapply(
        FUN, A, B,
        MoreArgs = list(X = x),
        BPPARAM = BPPARAM,
        SIMPLIFY = FALSE)
    distlist <- do.call(c, unname(distlist))
    # reformat
    # initialize distmat with NAs
    distmat <- matrix(NA_real_, nrow(x), nrow(x))
    rownames(distmat) <- colnames(distmat) <- rownames(x)
    matIndices <- matrix(c(unlist(B), unlist(A)), ncol = 2)
    distmat[matIndices] <- distlist
    #
    stats::as.dist(distmat)
}
