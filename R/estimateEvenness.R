#' Estimate Evenness measures
#'
#' This function calculates community evenness indices.
#' These include the \sQuote{Camargo}, \sQuote{Pielou}, \sQuote{Simpson},
#' \sQuote{Evar} and \sQuote{Bulla} evenness measures.
#' See details for more information and references.
#'
#' @param x a \code{\link{SummarizedExperiment}} object
#'
#' @param assay.type A single character value for selecting the
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}} used
#'   for calculation of the sample-wise estimates.
#'   
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'
#' @param index a \code{character} vector, specifying the evenness measures to
#'   be calculated.
#'
#' @param name a name for the column(s) of the colData the results should be
#'   stored in.
#'
#' @param BPPARAM A
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   object specifying whether calculation of estimates should be parallelized.
#'
#' @param ... optional arguments:
#' \itemize{
#'   \item threshold:  a numeric threshold. assay values below or equal
#'     to this threshold will be set to zero.
#' }
#'
#' @return \code{x} with additional \code{\link{colData}} named \code{*name*}
#'
#' @details
#' Evenness is a standard index in community ecology, and it quantifies how
#' evenly the abundances of different species are distributed. The following
#' evenness indices are provided:
#'
#' By default, this function returns all indices.
#'
#' The available evenness indices include the following (all in lowercase):
#' \itemize{
#'   \item 'camargo': Camargo's evenness (Camargo 1992)
#'   \item 'simpson_evenness': Simpson’s evenness is calculated as inverse
#'   Simpson diversity (1/lambda) divided by
#'   observed species richness S: (1/lambda)/S.
#'   \item 'pielou': Pielou's evenness (Pielou, 1966), also known as Shannon or
#'   Shannon-Weaver/Wiener/Weiner
#'   evenness; H/ln(S). The Shannon-Weaver is the preferred term; see
#'   Spellerberg and Fedor (2003).
#'   \item 'evar': Smith and Wilson’s Evar index (Smith & Wilson 1996).
#'   \item 'bulla': Bulla’s index (O) (Bulla 1994).
#' }
#'   
#' Desirable statistical evenness metrics avoid strong bias towards very
#' large or very small abundances; are independent of richness; and range
#' within the unit interval with increasing evenness (Smith & Wilson 1996).
#' Evenness metrics that fulfill these criteria include at least camargo,
#' simpson, smith-wilson, and bulla. Also see Magurran & McGill (2011)
#' and Beisel et al. (2003) for further details.
#'
#' @references
#'
#' Beisel J-N. et al. (2003)
#' A Comparative Analysis of Evenness Index Sensitivity.
#' _Internal Rev. Hydrobiol._ 88(1):3-15.
#' URL: \url{https://portais.ufg.br/up/202/o/2003-comparative_evennes_index.pdf}
#'
#' Bulla L. (1994)
#' An  index  of  evenness  and  its  associated  diversity  measure.
#' _Oikos_ 70:167--171.
#'
#' Camargo, JA. (1992)
#' New diversity index for assessing structural alterations in aquatic
#' communities.
#' _Bull. Environ. Contam. Toxicol._ 48:428--434.
#'
#' Locey KJ and Lennon JT. (2016)
#' Scaling laws predict global microbial diversity.
#' _PNAS_ 113(21):5970-5975; doi:10.1073/pnas.1521291113.
#'
#' Magurran AE, McGill BJ, eds (2011)
#' Biological Diversity: Frontiers in Measurement and Assessment
#' (Oxford Univ Press, Oxford), Vol 12.
#'
#' Pielou, EC. (1966)
#' The measurement of diversity in different types of
#' biological collections. _J Theoretical Biology_ 13:131--144.
#'
#' Smith B and Wilson JB. (1996)
#' A Consumer's Guide to Evenness Indices.
#' _Oikos_ 76(1):70-82.
#'
#' Spellerberg and Fedor (2003).
#' A tribute to Claude Shannon (1916 –2001) and a plea for more rigorous use of
#' species richness, species diversity and the ‘Shannon–Wiener’ Index.
#' _Alpha Ecology & Biogeography_ 12, 177–197.
#'
#' @seealso
#' \code{\link[scater:plotColData]{plotColData}}
#' \itemize{
#'   \item{\code{\link[=estimateRichness]{estimateRichness}}}
#'   \item{\code{\link[=estimateDominance]{estimateDominance}}}
#'   \item{\code{\link[=estimateDiversity]{estimateDiversity}}}
#' }
#'
#' @name .estimate_evenness
#' @noRd
#' @examples
#' data(esophagus)
#' tse <- esophagus
#'
#' # Specify index and their output names
#' index <- c("pielou", "camargo", "simpson_evenness", "evar", "bulla")
#' name  <- c("Pielou", "Camargo", "SimpsonEvenness",  "Evar", "Bulla")
#'
#' # Estimate evenness and give polished names to be used in the output
#' tse <- estimate_evenness(tse, index = index, name = name)
#' # Check the output
#' head(colData(tse))
#'
NULL

.estimate_evenness <- function(
        x, assay.type = assay_name, assay_name = "counts",
        index = c("camargo", "pielou", "simpson_evenness", "evar", "bulla"),
        name = index, ..., BPPARAM = SerialParam()){
    # input check
    index <- match.arg(index, several.ok = TRUE)
    if(!.is_non_empty_character(name) || length(name) != length(index)){
        stop("'name' must be a non-empty character value and have the ",
            "same length as 'index'.",
            call. = FALSE)
    }
    .check_assay_present(assay.type, x)
    #
    vnss <- BiocParallel::bplapply(
        index,
        .get_evenness_values,
        mat = assay(x, assay.type),
        BPPARAM = BPPARAM, ...)
    x <- .add_values_to_colData(x, vnss, name)
    return(x)
}

.calc_bulla_evenness <- function(mat) {
    # Species richness (number of species)
    S <- colSums2(mat > 0, na.rm = TRUE)

    # Relative abundances
    p <- t(mat)/colSums2(mat, na.rm = TRUE)

    i <- seq_len(nrow(p))
    O <- vapply(i,function(i){sum(pmin(p[i,], 1/S[i]))},numeric(1))

    # Bulla's Evenness
    (O - 1/S)/(1 - 1/S)
}

# Camargo's eveness x: species counts zeroes: include zeros Inspired
# by code from Pepijn de Vries and Zhou Xiang at
# researchgate.net/post/How_can_we_calculate_the_Camargo_evenness_index_in_R
# but rewritten here
.calc_camargo_evenness <- function(mat) {
    N <- colSums2(mat > 0, na.rm = TRUE)

    seq <- IntegerList(lapply(N - 1,seq_len))

    x <- mapply(
        function(i, n, s){
            xx <- 0
            for (j in s) {
                xx <- xx + sum(abs(mat[(j + 1):n,i] - mat[j,i]))
            }
            xx
        },
        seq_along(N),
        N,
        seq)
    # Return
    1 - x/(colSums2(mat, na.rm = TRUE) * N)
}

# x: Species count vector
.calc_simpson_evenness <- function(mat) {

    # Species richness (number of detected species)
    S <- colSums2(mat > 0, na.rm = TRUE)

    # Simpson evenness (Simpson diversity per richness)
    .calc_inverse_simpson(mat)/S
}

# x: Species count vector
.calc_pielou_evenness <- function(mat) {
    # Remove zeroes
    mat[mat == 0] <- NA

    # Species richness (number of detected species)
    S <- colSums2(mat > 0, na.rm = TRUE)

    # Relative abundances
    p <- t(mat)/colSums2(mat, na.rm = TRUE)

    # Shannon index
    H <- (-rowSums2(p * log(p), na.rm = TRUE))

    # Simpson evenness
    H/log(S)
}

# Smith and Wilson’s Evar index
.calc_evar_evenness <- function(mat) {
    N <- colSums2(mat, na.rm = TRUE)

    # Log abundance
    a <- log(mat)
    a[is.na(a) | is.infinite(a)] <- 0

    # Richness
    S <- colSums2(mat > 0, na.rm = TRUE)

    c <- colSums2(a, na.rm = TRUE)/S
    d <- t((t(a) - c)^2/S)
    d[mat == 0] <- 0

    f <- colSums2(d, na.rm = TRUE)

    (1 - 2/pi * atan(f))
}

.get_evenness_values <- function(index, mat, threshold = 0, ...){
    if(!is.numeric(threshold) || length(threshold) != 1L){
        stop("'threshold' must be a single numeric value.", call. = FALSE)
    }
    if(threshold > 0){
        mat[mat <= threshold] <- 0
    }
    FUN <- switch(index,
        camargo = .calc_camargo_evenness,
        pielou = .calc_pielou_evenness,
        simpson_evenness = .calc_simpson_evenness,
        evar = .calc_evar_evenness,
        bulla = .calc_bulla_evenness
        )
    res <- FUN(mat = mat, ...)
    res <- unname(res)
    return(res)
}

