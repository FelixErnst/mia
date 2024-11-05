#' Estimate dominance measures
#'
#' This function calculates community dominance indices.
#' This includes the \sQuote{Absolute}, \sQuote{Berger-Parker},
#' \sQuote{Core abundance},
#' \sQuote{Gini}, \sQuote{McNaughton’s}, \sQuote{Relative}, and
#' \sQuote{Simpson's} indices.
#'
#' @param x a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{
#'   SummarizedExperiment}} object
#'
#' @param assay.type A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}}
#'   to calculate the sample-wise estimates.
#'
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'   
#' @param index a \code{character} vector, specifying the indices to be
#'   calculated.
#'
#' @param ntaxa Optional and only used for the \code{Absolute} and
#'   \code{Relative} dominance indices: The n-th position of the dominant taxa
#'   to consider (default: \code{ntaxa = 1}). Disregarded for the indices
#'   \dQuote{dbp},
#'   \dQuote{core_abundance}, \dQuote{Gini}, \dQuote{dmn}, and \dQuote{Simpson}.
#'
#' @param aggregate Optional and only used for the \code{Absolute}, \code{dbp},
#'   \code{Relative}, and \code{dmn} dominance indices:
#'   Aggregate the values for top members selected by \code{ntaxa} or not. If
#'   \code{TRUE}, then the sum of relative abundances is returned. Otherwise the
#'   relative abundance is returned for the single taxa with the indicated rank
#'   (default: \code{aggregate = TRUE}). Disregarded for the indices
#'   \dQuote{core_abundance}, \dQuote{gini}, \dQuote{dmn}, and \dQuote{simpson}.
#'
#' @param name A name for the column(s) of the colData where the calculated
#'   Dominance indices should be stored in.
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'   (Currently not used)
#'
#' @param ... additional arguments currently not used.
#'
#' @details
#'
#' A dominance index quantifies the dominance of one or few species in a
#' community. Greater values indicate higher dominance.
#'
#' Dominance indices are in general negatively correlated with alpha diversity
#' indices (species richness, evenness, diversity, rarity). More dominant
#' communities are less diverse.
#'
#' \code{.estimate_dominance} calculates the following community dominance
#' indices:
#'
#' \itemize{
#' 
#' \item 'absolute': Absolute index equals to the absolute abundance of the
#' most dominant n species of the sample (specify the number with the argument
#' \code{ntaxa}). Index gives positive integer values.
#' 
#' \item 'dbp': Berger-Parker index (See Berger & Parker 1970) calculation
#' is a special case of the 'relative' index. dbp is the relative abundance of
#' the most
#' abundant species of the sample. Index gives values in interval 0 to 1,
#' where bigger value represent greater dominance.
#'
#' \deqn{dbp = \frac{N_1}{N_{tot}}}{%
#' dbp = N_1/N_tot} where \eqn{N_1} is the absolute abundance of the most
#' dominant species and \eqn{N_{tot}} is the sum of absolute abundances of all
#' species.
#' 
#' \item 'core_abundance': Core abundance index is related to core species.
#' Core species are species that are most abundant in all samples, i.e., in
#' whole data set. Core species are defined as those species that have
#' prevalence over 50\%. It means that in order to belong to core species,
#' species must be prevalent in 50\% of samples. Core species are used to
#' calculate the core abundance index. Core abundance index is sum of relative
#' abundances of core species in the sample. Index gives values in interval
#' 0 to 1, where bigger value represent greater dominance.
#'
#' \deqn{core_abundance = \frac{N_{core}}{N_{tot}}}{%
#' core_abundance = N_core/N_tot} where \eqn{N_{core}} is the sum of absolute
#' abundance of the core species and \eqn{N_{tot}} is the sum of absolute
#' abundances of all species.
#' 
#' \item 'gini':  Gini index is probably best-known from socioeconomic
#' contexts (Gini 1921). In economics, it is used to measure, for example, how
#' unevenly income is distributed among population. Here, Gini index is used
#' similarly, but income is replaced with abundance. 
#' 
#' If there is small group of species
#' that represent large portion of total abundance of microbes, the inequality
#' is large and Gini index closer to 1. If all species has equally large
#' abundances, the equality is perfect and Gini index equals 0. This index
#' should not be confused with Gini-Simpson index, which quantifies diversity.
#'
#' \item 'dmn': McNaughton’s index is the sum of relative abundances of the two
#' most abundant species of the sample (McNaughton & Wolf, 1970). Index gives
#' values in the unit interval:
#'
#' \deqn{dmn = (N_1 + N_2)/N_tot}
#'
#' where \eqn{N_1} and \eqn{N_2} are the absolute
#' abundances of the two most dominant species and \eqn{N_{tot}} is the sum of
#' absolute abundances of all species.
#'
#' \item 'relative': Relative index equals to the relative abundance of the
#' most dominant n species of the sample (specify the number with the
#' argument \code{ntaxa}).
#' This index gives values in interval 0 to 1.
#'
#' \deqn{relative = N_1/N_tot}
#'
#' where \eqn{N_1} is the absolute abundance of the most
#' dominant species and \eqn{N_{tot}} is the sum of absolute abundances of all
#' species.
#'
#' \item 'simpson_lambda': Simpson's (dominance) index or Simpson's lambda is
#' the sum of squared relative abundances. This index gives values in the unit interval.
#' This value equals the probability that two randomly chosen individuals
#' belongs to the
#' same species. The higher the probability, the greater the dominance (See
#' e.g. Simpson 1949).
#'
#' \deqn{lambda = \sum(p^2)}
#'
#' where p refers to relative abundances.
#'
#' There is also a more advanced Simpson dominance index (Simpson 1949).
#' However, this is not provided and the simpler squared sum of relative
#' abundances is used instead as the alternative index is not in the unit
#' interval and it is highly
#' correlated with the simpler variant implemented here.
#' 
#' }
#'
#' @references
#'
#' Berger WH & Parker FL (1970)
#' Diversity of Planktonic Foraminifera in Deep-Sea Sediments.
#' _Science_ 168(3937):1345-1347. doi: 10.1126/science.168.3937.1345
#'
#' Gini C (1921)
#' Measurement of Inequality of Incomes.
#' _The Economic Journal_ 31(121): 124-126. doi: 10.2307/2223319
#'
#' McNaughton, SJ and Wolf LL. (1970).
#' Dominance and the niche in ecological systems.
#' _Science_ 167:13, 1--139
#'
#' Simpson EH (1949)
#' Measurement of Diversity.
#' _Nature_ 163(688). doi: 10.1038/163688a0
#'
#' @return \code{x} with additional \code{\link{colData}} named
#'   \code{*name*}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[=estimateRichness]{estimateRichness}}
#'   \item \code{\link[=estimateEvenness]{estimateEvenness}}
#'   \item \code{\link[=estimateDiversity]{estimateDiversity}}
#' }
#'
#' @name .estimate_dominance
#' @noRd
#'
#' @examples
#' data(esophagus)
#'
#' # Calculates Simpson's lambda (can be used as a dominance index)
#' esophagus <- .estimate_dominance(esophagus, index="simpson_lambda")
#'
#' # Shows all indices
#' colData(esophagus)
#'
#' # Indices must be written correctly (e.g. dbp, not dbp), otherwise an error
#' # gets thrown
#' \donttest{esophagus <- .estimate_dominance(esophagus, index="dbp")}
#' # Calculates dbp and Core Abundance indices
#' esophagus <- .estimate_dominance(esophagus, index=c("dbp", "core_abundance"))
#' # Shows all indices
#' colData(esophagus)
#' # Shows dbp index
#' colData(esophagus)$dbp
#' # Deletes dbp index
#' colData(esophagus)$dbp <- NULL
#' # Shows all indices, dbp is deleted
#' colData(esophagus)
#' # Deletes all indices
#' colData(esophagus) <- NULL
#'
#' # Calculates all indices
#' esophagus <- .estimate_dominance(esophagus)
#' # Shows all indices
#' colData(esophagus)
#' # Deletes all indices
#' colData(esophagus) <- NULL
#'
#' # Calculates all indices with explicitly specified names
#' esophagus <- .estimate_dominance(esophagus,
#'     index = c("dbp", "dmn", "absolute", "relative",
#'               "simpson_lambda", "core_abundance", "gini"),
#'     name  = c("BergerParker", "McNaughton", "Absolute", "Relative",
#'               "SimpsonLambda", "CoreAbundance", "Gini")
#' )
#' # Shows all indices
#' colData(esophagus)
#'
NULL

setGeneric(
    ".estimate_dominance", signature = c("x"),
    function(x, ...) standardGeneric(".estimate_dominance"))

setGeneric(
    ".estimate_dominance",signature = c("x"),
    function(
        x, assay.type = assay_name, assay_name = "counts",
        index = c(
            "absolute", "dbp", "core_abundance", "gini", "dmn", "relative",
            "simpson_lambda"),
        ntaxa = 1, aggregate = TRUE, name = index, BPPARAM = SerialParam(), ...)
    standardGeneric(".estimate_dominance"))

.estimate_dominance <- function(
        x, assay.type = assay_name, assay_name = "counts",
        index = c(
            "absolute", "dbp", "core_abundance", "gini", "dmn", "relative",
            "simpson_lambda"),
        ntaxa = 1, aggregate = TRUE, name = index, BPPARAM = SerialParam(),
        ...){
    # Input check
    # Check assay.type
    .check_assay_present(assay.type, x)
    # Check indices
    index <- match.arg(index, several.ok = TRUE)
    if(!.is_non_empty_character(name) || length(name) != length(index)){
        stop("'name' must be a non-empty character value and have the 
             same length as 'index'",
             call. = FALSE)
    }

    # Check aggregate
    if(!.is_a_bool(aggregate)){
        stop("'aggregate' must be TRUE or FALSE.", call. = FALSE)
    }

    # Calculates dominance indices
    dominances <- BiocParallel::bplapply(
        index,
        FUN = .get_dominance_values,
        mat = assay(x,assay.type),
        ntaxa = ntaxa,
        aggregate = aggregate,
        BPPARAM = BPPARAM)

    # Add dominance indices to colData
    x <- .add_values_to_colData(x, dominances, name)
    return(x)
}

#---------------------------Help functions--------------------------------------

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

.calc_dominance <- function(mat, ntaxa, aggregate, index){

    # Check ntaxa
    if(!(ntaxa>0 && ntaxa<3)){
        stop("'ntaxa' must be a numerical value 1 or 2.", call. = FALSE)
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
