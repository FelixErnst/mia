#' @name
#' getPERMANOVA
#' 
#' @title
#' Calculate PERMANOVA (Permutational Multivariate Analysis of Variance)
#' 
#' @description
#' These functions perform PERMANOVA to assess the significance of group
#' differences based on a specified dissimilarity matrix. The results can be
#' returned directly or added to metadata in an object of class
#' \code{TreeSummarizedExperiment}.
#' 
#' @details
#' PERMANOVA is a non-parametric method used to test whether the centroids of 
#' different groups (as defined by the formula or covariates) are significantly 
#' different in terms of multivariate space.
#' 
#' PERMANOVA relies on the assumption of group homogeneity, meaning the groups
#' should be distinct and have similar variances within each group. This
#' assumption is essential as PERMANOVA is sensitive to differences in
#' within-group dispersion, which can otherwise confound results.
#' This is why the functions return homogeneity test results by default.
#' 
#' The functions utilize \code{\link[vegan:adonis2]{vegan::adonis2}} to compute
#' PERMANOVA. For homogeneity testing,
#' \code{\link[vegan:betadisper]{vegan::betadisper}} along
#' with \code{\link[vegan:permutest]{vegan::permutest}} are utilized by default,
#' which allow testing for equal dispersion across groups and validate the
#' homogeneity assumption.
#' 
#' PERMANOVA and distance-based redundancy analysis (dbRDA) are closely related
#' methods for analyzing multivariate data. PERMANOVA is non-parametric, making
#' fewer assumptions about the data. In contrast, dbRDA assumes a linear
#' relationship when constructing the ordination space, although it also
#' employs a PERMANOVA-like approach to test the significance of predictors
#' within this space. dbRDA offers a broader scope overall, as it provides
#' visualization of the constrained ordination, which can reveal patterns and
#' relationships. However, when the underlying data structure is non-linear,
#' the results from the two methods can differ significantly due to dbRDA's
#' reliance on linear assumptions.
#' 
#' @return
#' \code{getPERMANOVA} returns the PERMANOVA results or a list containing the
#' PERMANOVA results and homogeneity test results if
#' \code{test.homogeneity = TRUE}. \code{addPERMANOVA} adds these results to
#' metadata of \code{x}.
#' 
#' @inheritParams runCCA
#' 
#' @param name \code{Character scalar}. A name for the results that will be
#' stored to metadata. (Default: \code{"permanova"})
#' 
#' @param method \code{Character scalar}. A dissimilarity metric used in
#' PERMANOVA and group dispersion calculation. (Default: \code{"bray"})
#' 
#' @param test.homogeneity \code{Logical scalar}. Should the homogeneity of
#' group dispersions be evaluated? (Default: \code{TRUE})
#'
#' @param ... additional arguments passed to \code{vegan::adonis2}.
#' \itemize{
#'   \item \code{by}: \code{Character scalar}. Specifies how significance is
#'   calculated. (Default: \code{"margin"})
#'   
#'   \item \code{na.action}: \code{function}. Action to take when missing
#'   values for any of the variables in \code{formula} are encountered.
#'   (Default: \code{na.fail})
#'   
#'   \item \code{full} \code{Logical scalar}. should all the results from the
#'   homogeneity calculations be returned. When \code{FALSE}, only
#'   summary tables are returned. (Default: \code{FALSE})
#'   
#'   \item \code{homogeneity.test}: \code{Character scalar}. Specifies
#'   the significance test used to analyse
#'   \code{\link[vegan:betadisper]{vegan::betadisper}} results.
#'   Options include 'permanova'
#'   (\code{\link[vegan:permutest]{vegan::permutest}}), 'anova'
#'   (\code{\link[stats:anova]{stats::anova}}) and 'tukeyhsd'
#'   (\code{\link[stats:TukeyHSD]{stats::TukeyHSD}}).
#'   (Default: \code{"permanova"})
#' }
#'
#' @examples
#' 
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' 
#' # Apply relative transformation
#' tse <- transformAssay(tse, method = "relabundance")
#' # Perform PERMANOVA
#' tse <- addPERMANOVA(
#'     tse,
#'     assay.type = "relabundance",
#'     method = "bray",
#'     formula = x ~ SampleType,
#'     permutations = 99
#'     )
#' # The results are stored to metadata
#' metadata(tse)[["permanova"]]
#' 
#' # Calculate dbRDA
#' rda_res <- getRDA(
#'     tse, assay.type = "relabundance", method = "bray",
#'     formula = x ~ SampleType, permutations = 99)
#' # Significance results are similar to PERMANOVA
#' attr(rda_res, "significance")
#' 
#' @seealso
#' For more details on the actual implementation see
#' \code{\link[vegan:adonis2]{vegan::adonis2}},
#' \code{\link[vegan:betadisper]{vegan::betadisper}}, and
#' \code{\link[vegan:permutest]{vegan::permutest}}. See also
#' \code{\link[=runCCA]{addCCA}} and \code{\link[=runCCA]{addRDA}}
#'
NULL

#' @rdname getPERMANOVA
setGeneric("getPERMANOVA", signature = c("x"), function(x, ...)
    standardGeneric("getPERMANOVA"))

#' @rdname getPERMANOVA
setGeneric("addPERMANOVA", signature = c("x"), function(x, ...)
    standardGeneric("addPERMANOVA"))

#' @export
#' @rdname getPERMANOVA
setMethod("getPERMANOVA", "SingleCellExperiment", function(x,  ...){
    # Get altexp if specified
    x <- .check_and_get_altExp(x, ...)
    res <- callNextMethod(x, ...)
    return(res)
    }
)

#' @export
#' @rdname getPERMANOVA
setMethod("getPERMANOVA", "SummarizedExperiment",
    function(
        x, assay.type = "counts", formula = NULL, col.var = NULL, ...){
        ############################# Input check ##############################
        # Assay must be present
        .check_assay_present(assay.type, x)
        # Formula must be either correctly specified or not specified
        if( !(is.null(formula) || is(formula, "formula")) ){
            stop("'formula' must be formula or NULL.", call. = FALSE)
        }
        # User can also specify covariates with col.var
        if( !(is.null(col.var) || (is.character(col.var) &&
                all(col.var %in% colnames(colData(x))))) ){
            stop("'col.var' must specify column from colData(x) or be NULL.",
                call. = FALSE)
        }
        ########################### Input check end ############################
        # Get abundance table, formula and related sample metadata as DF
        mat <- assay(x, assay.type)
        temp <- .get_formula_and_covariates(x, formula, col.var)
        formula <- temp[["formula"]]
        covariates <- temp[["variables"]]
        # Calculate PERMANOVA with matrix method
        res <- getPERMANOVA(mat, formula = formula, data = covariates, ...)
        return(res)
    }
)

#' @export
#' @rdname getPERMANOVA
setMethod("getPERMANOVA", "ANY", function(
        x, formula, data, method = "bray", test.homogeneity = TRUE, ...){
    if( !is.matrix(x) ){
        stop("'x' must be matrix.", call. = FALSE)
    }
    if( !is(formula, "formula") ){
        stop("'formula' must be formula or NULL.", call. = FALSE)
    }
    if( !(is.data.frame(data) || is.matrix(data) || is(data, "DFrame")) ){
        stop("'data' must be data.frame or coarcible to one.", call. = FALSE)
    }
    if( ncol(x) != nrow(data) ){
        stop("Number of columns in 'x' should match with number of rows in ",
            "'data'.", call. = FALSE)
    }
    if( !.is_a_string(method) ){
        stop("'method' must be a single character value.", call. = FALSE)
    }
    if( !.is_a_bool(test.homogeneity) ){
        stop("'test.homogeneity' must be TRUE or FALSE.", call. = FALSE)
    }
    # Calculate PERMANOVA
    res <- .calculate_permanova(
        x, formula = formula, data = data, method = method, ...)
    # Test homogeneity
    if( test.homogeneity ){
        homogeneity <- .calculate_homogeneity(x, data, method = method, ...)
        res <- list(permanova = res, homogeneity = homogeneity)
    }
    return(res)
})

#' @export
#' @rdname getPERMANOVA
setMethod("addPERMANOVA", "SummarizedExperiment",
    function(x, name = "permanova", ...){
    if( !.is_a_string(name) ){
        stop("'name' must be a single character value.", call. = FALSE)
    }
    # Calculate permanova
    res <- getPERMANOVA(x, ...)
    # Addresults to metadata
    x <- .add_values_to_metadata(x, name, res)
    return(x)
    }
)

################################ HELP FUNCTIONS ################################

# This function is internal function to perform PERMANOVA from abundance
# matrix, formula and sample metadata table.
#' @importFrom vegan adonis2
.calculate_permanova <- function(
        x, formula, data, by = "margin", na.action = na.fail, ...){
    #
    if( !.is_a_string(by) ){
        stop("'by' must be a single character value.", call. = FALSE)
    }
    #
    # Get abundance data into correct orientation. Samples must be in rows and
    # features in columns. Also ensure that the abundance table is matrix.
    x <- as.matrix(t(x))
    # Covariate table must be data.frame
    data <- data.frame(data, check.names = FALSE)
    # Instead of letting na.action pass through, give informative error
    # about missing values.
    if( any(is.na(data)) && isTRUE(all.equal(na.action, na.fail)) ){
        stop("Variables contain missing values. Set na.action to na.exclude ",
            "to remove samples with missing values.", call. = FALSE)
    }
    # This next step ensures that the formula points correctly to abundance
    # table
    formula <- as.character(formula)
    formula[[2]] <- "x"
    formula <- as.formula(
        paste(as.character(formula)[c(2,1,3)], collapse = " "))
    # Calculate PERMANOVA
    res <- adonis2(
        formula = formula, data = data, by = by, na.action = na.action, ...)
    return(res)
}
