
.calc_observed <- function(mat, detection = 0, ...){
    if( !.is_a_numeric(detection) ){
        stop("'detection' must be numeric value.", call. = FALSE)
    }
    colSums(mat > detection)
}

.calc_chao1 <- function(mat, ...){
    # Required to work with DelayedArray
    if(is(mat, "DelayedArray")) {
        mat <- matrix(mat, nrow = nrow(mat))
    }

    ans <- t(vegan::estimateR(t(mat))[c("S.chao1","se.chao1"),])
    colnames(ans) <- c("", "_se")
    ans
}

.calc_ace <- function(mat, ...){
    # Required to work with DelayedArray
    if(is(mat, "DelayedArray")) {
        mat <- matrix(mat, nrow = nrow(mat))
    }

    ans <- t(vegan::estimateR(t(mat))[c("S.ACE","se.ACE"),])
    colnames(ans) <- c("", "_se")
    ans
}

.calc_hill <- function(mat, ...){
    # Exponent of Shannon diversity
    exp(vegan::diversity(t(mat), index="shannon"))
}

.get_richness_values <- function(index, mat, detection, ...) {
    FUN <- switch(index,
        observed = .calc_observed,
        chao1 = .calc_chao1,
        ace = .calc_ace,
        hill = .calc_hill
    )
    res <- FUN(mat = mat, detection = detection, ...)
    if( is.matrix(res) ){
        rownames(res) <- NULL
    } else{
        res <- unname(res)
    }
    return(res)
}
