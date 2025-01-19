
.gini_dominance <- function(x, w=rep(1, length(x))) {
    # See also reldist::gini for an independent implementation
    x <- as.vector(x)
    o <- order(x)
    x <- x[o]
    w <- w[o]/sum(w)
    p <- cumsum(w)
    nu <- cumsum(w * x)
    n <- length(nu)
    nu <- nu/nu[[n]]
    sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}

.calc_gini_dominance <- function(mat, ...){
    apply(mat, 2L, .gini_dominance)
}

.calc_core_dominance <- function(mat, ...){
    getPrevalentAbundance(mat, detection = 0, as.relative = TRUE)
}

.calc_dominance <- function(mat, index, ntaxa = 1L, aggregate = TRUE, ...){

    # Check ntaxa
    if(!(.is_an_integer(ntaxa) && ntaxa>0 && ntaxa<3)){
        stop("'ntaxa' must be a numerical value 1 or 2.", call. = FALSE)
    }
    # Check aggregate
    if(!.is_a_bool(aggregate)){
        stop("'aggregate' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    if (index == "absolute") {
        # ntaxa=1 by default but can be tuned
        as_relative <- FALSE
    } else if (index == "relative") {
        # ntaxa=1 by default but can be tuned
        as_relative <- TRUE
    } else if (index == "dbp") {
        # Berger-Parker: if selected fix the following values
        ntaxa <- 1
        as_relative <- TRUE
    } else if (index == "dmn") {
        # McNaughton's dominance: if selected fix the following values
        ntaxa <- 2
        aggregate <- TRUE
        as_relative <- TRUE
    }

    if (as_relative) {
        # Calculates the relative abundance per sample
        mat <- .calc_rel_abund(mat)
    }

    # Aggregate or not
    if (!aggregate) {
        idx <- apply(mat, 2L,
                    function(mc) {
                        order(as.vector(mc), decreasing = TRUE)[[ntaxa]]
                    })
    } else {
        idx <- apply(mat, 2L,
                    function(mc) {
                        order(as.vector(mc), decreasing = TRUE)[seq_len(ntaxa)]
                    })
        idx <- split(as.vector(idx),
                    unlist(lapply(seq_len(length(idx) / ntaxa),rep.int,ntaxa)))
    }

    ans <- lapply(mapply(function(i,j,x){x[i,j]},
                        i = idx,
                        j = seq_len(ncol(mat)),
                        MoreArgs = list(x = mat),
                        SIMPLIFY = FALSE),
                    sum)
    ans <- unlist(ans)

    # Adds sample names to the table
    names(ans) <- colnames(mat)
    ans
}

.get_dominance_values <- function(
        index, mat, ntaxa = 1, aggregate = TRUE, ...) {
    FUN <- switch(index,
        simpson_lambda = .simpson_lambda,
        core_abundance = .calc_core_dominance,
        gini = .calc_gini_dominance,
        absolute = .calc_dominance,
        relative = .calc_dominance,
        dbp = .calc_dominance,
        dmn = .calc_dominance
        )
    res <- FUN(index, mat = mat, ntaxa = ntaxa, aggregate = aggregate, ...)
    res <- unname(res)
    return(res)
}
