#' Canonical Correspondence Analysis and Redundancy Analysis
#'
#' These functions perform Canonical Correspondence Analysis on data stored
#' in a \code{SummarizedExperiment}.
#'
#' @inheritParams getDominant
#' @inheritParams getDissimilarity
#'
#' @param formula \code{formula}. If \code{x} is a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' a formula can be supplied. Based on the right-hand side of the given formula
#' \code{colData} is subset to \code{col.var}.
#'   
#' \code{col.var} and \code{formula} can be missing, which turns the CCA
#' analysis into a CA analysis and dbRDA into PCoA/MDS.
#'
#' @param data \code{data.frame} or coarcible to one. The covariance table
#' including covariates defined by \code{formula}.
#'
#' @param col.var \code{Character scalar}. When \code{x} is a
#' \code{SummarizedExperiment},\code{col.var} can be used to specify variables
#' from \code{colData}.
#'   
#' @param variables Deprecated. Use \code{col.var} instead.
#' 
#' @param test.signif \code{Logical scalar}. Should the PERMANOVA and analysis
#' of multivariate homogeneity of group dispersions be performed.
#' (Default: \code{TRUE})
#'   
#' @param altexp \code{Character scalar} or \code{integer scalar}. Specifies an
#' alternative experiment containing the input data.
#' 
#' @param name \code{Character scalar}. A name for the \code{reducedDim()}  
#' where results will be stored. (Default: \code{"CCA"})
#' 
#' @param exprs_values Deprecated. Use \code{assay.type} instead.
#'
#' @param ... additional arguments passed to vegan::cca or vegan::dbrda and
#' other internal functions.
#' \itemize{
#'   \item \code{method} a dissimilarity measure to be applied in dbRDA and
#'   possible following homogeneity test. (Default: \code{"euclidean"})
#'   
#'   \item \code{scale}: \code{Logical scalar}. Should the expression values be
#'   standardized? \code{scale} is disabled when using \code{*RDA} functions.
#'   Please scale before performing RDA. (Default: \code{TRUE})
#'   
#'   \item \code{na.action}: \code{function}. Action to take when missing
#'   values for any of the variables in \code{formula} are encountered.
#'   (Default: \code{na.fail})
#'   
#'   \item \code{full} \code{Logical scalar}. Should all the results from the
#'   significance calculations be returned. When \code{FALSE}, only
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
#'   
#'   \item \code{permutations} a numeric value specifying the number of
#'   permutations for significance testing in \code{vegan::anova.cca}.
#'   (Default: \code{999})
#' }
#' 
#' @details
#' *CCA functions utilize \code{vegan:cca} and *RDA functions
#' \code{vegan:dbRDA}. By default, dbRDA is done with euclidean distances, which
#' is equivalent to RDA. \code{col.var} and \code{formula} can be missing,
#' which turns the CCA analysis into a CA analysis and dbRDA into PCoA/MDS.
#'   
#' Significance tests are done with \code{vegan:anova.cca} (PERMANOVA). Group
#' dispersion, i.e., homogeneity within groups is analyzed with 
#' \code{\link[vegan:betadisper]{vegan::betadisper}}
#' (multivariate homogeneity of groups dispersions
#' (variances)) and statistical significance of homogeneity is tested with a
#' test specified by \code{homogeneity.test} parameter.
#'
#' @return
#' For \code{getCCA} a matrix with samples as rows and CCA dimensions as
#' columns. Attributes include output from \code{\link[vegan:scores]{scores}},
#' eigenvalues, the \code{cca}/\code{rda} object and significance analysis
#' results.
#'
#' For \code{addCCA} a modified \code{x} with the results stored in
#' \code{reducedDim} as the given \code{name}.
#'
#' @name runCCA
#' @seealso
#' For more details on the actual implementation see
#' \code{\link[vegan:cca]{cca}} and \code{\link[vegan:dbrda]{dbrda}}
#'
#' @examples
#' library(miaViz)
#' data("enterotype", package = "mia")
#' tse <- enterotype
#'
#' # Perform CCA and exclude any sample with missing ClinicalStatus
#' tse <- addCCA(
#'     tse,
#'     formula = data ~ ClinicalStatus,
#'     na.action = na.exclude
#'     )
#'
#' # Plot CCA
#' plotCCA(tse, "CCA", colour_by = "ClinicalStatus")
#'
#' # Fetch significance results
#' attr(reducedDim(tse, "CCA"), "significance")
#'
#' tse <- transformAssay(tse, method = "relabundance")
#'
#' # Specify dissimilarity measure
#' tse <- addRDA(
#'     tse,
#'     formula = data ~ ClinicalStatus,
#'     assay.type = "relabundance",
#'     method = "bray",
#'     name = "RDA_bray",
#'     na.action = na.exclude
#'     )
#'
#' # To scale values when using *RDA functions, use
#' # transformAssay(MARGIN = "features", ...) 
#' tse <- transformAssay(tse, method = "standardize", MARGIN = "features")
#'
#' # Data might include taxa that do not vary. Remove those because after
#' # z-transform their value is NA
#' tse <- tse[rowSums(is.na(assay(tse, "standardize"))) == 0, ]
#'
#' # Calculate RDA
#' tse <- addRDA(
#'    tse,
#'    formula = data ~ ClinicalStatus,
#'    assay.type = "standardize",
#'    name = "rda_scaled",
#'    na.action = na.omit
#'    )
#'
#' # Plot RDA
#' plotRDA(tse, "rda_scaled", colour_by = "ClinicalStatus")
#'
#' # A common choice along with PERMANOVA is ANOVA when statistical significance
#' # of homogeneity of groups is analysed. Moreover, full significance test
#' # results can be returned.
#' tse <- addRDA(
#'     tse,
#'     formula = data ~ ClinicalStatus,
#'     homogeneity.test = "anova",
#'     full = TRUE
#'     )
#'
#' # Example showing how to pass extra parameters, such as 'permutations',
#' # to anova.cca
#' tse <- addRDA(
#'     tse,
#'     formula = data ~ ClinicalStatus,
#'     permutations = 500
#'     )
#'
NULL

#' @rdname runCCA
setGeneric("getCCA", signature = c("x"), function(x, ...)
    standardGeneric("getCCA"))

#' @rdname runCCA
setGeneric("addCCA", signature = c("x"), function(x, ...)
    standardGeneric("addCCA"))

#' @rdname runCCA
setGeneric("getRDA", signature = c("x"), function(x, ...)
    standardGeneric("getRDA"))

#' @rdname runCCA
setGeneric("addRDA", signature = c("x"), function(x, ...)
    standardGeneric("addRDA"))


############################# Exported CCA methods #############################

#' @export
#' @rdname runCCA
#' @aliases getCCA
calculateCCA <- function(x, ...){
    getCCA(x, ...)
}

#' @export
#' @rdname runCCA
#' @aliases addCCA
runCCA <- function(x, ...){
    addCCA(x, ...)
}

#' @export
#' @rdname runCCA
setMethod("getCCA", "ANY", function(x, formula, data, ...){
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
    #
    res <- .calculate_rda(
        x, formula = formula, data = data, ord.method = "CCA", ...)
    return(res)
    })

#' @export
#' @rdname runCCA
setMethod("getCCA", "SummarizedExperiment",
    function(
        x, formula = NULL, col.var = variables, variables = NULL,
        test.signif = TRUE, assay.type = assay_name, assay_name = exprs_values,
        exprs_values = "counts", ...){
        ############################# Input check ##############################
        if( !(is.null(formula) || is(formula, "formula")) ){
            stop("'formula' must be formula or NULL.", call. = FALSE)
        }
        if( !(is.null(col.var) ||
                (is.character(col.var) &&
                    all(col.var %in% colnames(colData(x))))) ){
            stop("'col.var' must specify column from colData(x) or be NULL.",
                call. = FALSE)
        }
        if( !is.null(formula) && !is.null(col.var) ){
            stop("Specify either 'formula' or 'col.var'.", call. = FALSE)
        }
        .check_assay_present(assay.type, x)
        if( !.is_a_bool(test.signif) ){
            stop("'test.signif' must be TRUE or FALSE.", call. = FALSE)
        }
        ########################### Input check end ############################
        # Get assay
        mat <- assay(x, assay.type)
        # Get formula and variables
        temp <- .get_formula_and_covariates(x, formula, col.var)
        formula <- temp[["formula"]]
        covariates <- temp[["variables"]]
        # Calculate CCA
        cca <- getCCA(mat, formula = formula, data = covariates, ...)
        # Test significance if specified
        if( test.signif ){
            res <- .test_rda(mat, attr(cca, "obj"), covariates, ...)
            attr(cca, "significance") <- res
        }
        return(cca)
        }
    )

#' @export
#' @rdname runCCA
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("addCCA", "SingleCellExperiment",
    function(x, altexp = NULL, name = "CCA", ...){
        ############################# Input check ##############################
        if( !(is.null(altexp) ||
            (length(altexp) == 1L && is.character(altexp) &&
                all(altexp %in% altExpNames(x))) ||
            (.is_an_integer(altexp) && altexp > 0 &&
                altexp <= length(altExps(x))) ) ){
            stop("'altexp' must specify an alternative experiment from ",
                "altExp(x).", call. = FALSE)
        }
        ########################### Input check end ############################
        # Get TreeSE from altexp if specified.
        if( !is.null(altexp) ){
            y <- altExp(x, altexp)
        } else {
            y <- x
        }
        # Calculate CCA
        cca <- getCCA(y, ...)
        # Add object to reducedDim
        x <- .add_object_to_reduceddim(x, cca, name = name, ...)
        return(x)
        }
    )

############################# Exported RDA methods #############################

#' @export
#' @rdname runCCA
#' @aliases getRDA
calculateRDA <- function(x, ...){
    getRDA(x, ...)
}

#' @export
#' @rdname runCCA
#' @aliases addRDA
runRDA <- function(x, ...){
    addRDA(x, ...)
}

#' @export
#' @rdname runCCA
setMethod("getRDA", "ANY", function(x, formula, data, ...){
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
    #
    res <- .calculate_rda(
        x, formula = formula, data = data, ord.method = "RDA", ...)
    return(res)
    })

#' @export
#' @rdname runCCA
setMethod("getRDA", "SummarizedExperiment",
    function(
        x, formula = NULL, col.var = variables, variables = NULL,
        test.signif = TRUE, assay.type = assay_name, assay_name = exprs_values,
        exprs_values = "counts", ...){
        ############################# Input check ##############################
        if( !(is.null(formula) || is(formula, "formula")) ){
            stop("'formula' must be formula or NULL.", call. = FALSE)
        }
        if( !(is.null(col.var) ||
                (is.character(col.var) &&
                    all(col.var %in% colnames(colData(x))))) ){
            stop("'col.var' must specify column from colData(x) or be NULL.",
                call. = FALSE)
        }
        if( !is.null(formula) && !is.null(col.var) ){
            stop("Specify either 'formula' or 'col.var'.", call. = FALSE)
        }
        .check_assay_present(assay.type, x)
        if( !.is_a_bool(test.signif) ){
            stop("'test.signif' must be TRUE or FALSE.", call. = FALSE)
        }
        ########################### Input check end ############################
        # Get assay
        mat <- assay(x, assay.type)
        # Get formula and variables
        temp <- .get_formula_and_covariates(x, formula, col.var)
        formula <- temp[["formula"]]
        covariates <- temp[["variables"]]
        # Calculate CCA
        rda <- getRDA(mat, formula = formula, data = covariates, ...)
        # Test significance if specified
        if( test.signif ){
            res <- .test_rda(mat, attr(rda, "obj"), covariates, ...)
            attr(rda, "significance") <- res
        }
        return(rda)
        }
    )

#' @export
#' @rdname runCCA
#' @importFrom SingleCellExperiment reducedDim<-
setMethod("addRDA", "SingleCellExperiment",
    function(x, altexp = NULL, name = "RDA", ...){
        ############################# Input check ##############################
        if( !(is.null(altexp) ||
            (length(altexp) == 1L && is.character(altexp) &&
                all(altexp %in% altExpNames(x))) ||
            (.is_an_integer(altexp) && altexp > 0 &&
                altexp <= length(altExps(x))) ) ){
            stop("'altexp' must specify an alternative experiment from ",
                "altExp(x).", call. = FALSE)
        }
        ########################### Input check end ############################
        # Get TreeSE from altexp if specified.
        if( !is.null(altexp) ){
            y <- altExp(x, altexp)
        } else {
            y <- x
        }
        # Calculate RDA
        rda <- getRDA(y, ...)
        # Add object to reducedDim
        x <- .add_object_to_reduceddim(x, rda, name = name, ...)
        return(x)
        }
    )

################################ HELP FUNCTIONS ################################

# THis function creates a formula and covariate table that can be then
# forwarded to subsequent functions. The formula/variables are created based
# on usr-specified formual or variables. If neither is specified, the function
# returns empty table with corresponding formula.
.get_formula_and_covariates <- function(x, formula, col.var){
    # If formula is specified
    if( !is.null(formula) ){
        # Otherwise if formula is provided, get variables based on formula
        variables <- .get_variables_based_on_formula(x, formula)
    } else if( !is.null(col.var) ){
        # If column variables are specified, create a formula based on them and
        # get variables
        formula <- as.formula(
            paste0("mat ~ `", paste(col.var, collapse = "` + `"), "`"))
        # Get the data from colData
        variables <- colData(x)[ , col.var, drop = FALSE]
    } else{
        # Else either formula nor column variables were specified. Get empty
        # variable table and formula that does not incorporate covariates.
        variables <- colData(x)[, 0]
        formula <- mat ~ 1
    }
    # Return a list that holds both formula and covariates
    res <- list(formula = formula, variables = variables)
    return(res)
}

# This function fetch variables from colData based on formula. If formula
# was not specified, the function returns an empty table.
#' @importFrom stats terms
#' @importFrom SummarizedExperiment colData
.get_variables_based_on_formula <- function(x, formula){
    # If user specified "all" covariates, give error. User should use col.var
    # to do that.
    if( as.character(formula)[[3]] == "." ){
        stop("To incorporate all variables from colData(x) with 'x~.', use ",
            "'col.var' and leave 'formula' unspecified.", call. = FALSE)
    }
    # Get variables from formula
    terms <- rownames(attr(terms(formula), "factors"))
    terms <- terms[terms != as.character(formula)[2L]]
    terms <- .remove_special_functions_from_terms(terms)
    # Check that all variables specify a column from colData
    if( !all(terms %in% colnames(colData(x))) ){
        stop("All variables on the right hand side of 'formula' must be ",
            "present in colData(x).", call. = FALSE)
    }
    # Get the variables from colData
    df <- colData(x)[, terms, drop = FALSE]
    return(df)
}

# This function parses right-hand side formula so that it now includes only
# the covariates.
.remove_special_functions_from_terms <- function(terms){
    names(terms) <- terms
    m <- regexec("^Condition\\(([^\\(\\)]*)\\)$|^([^\\(\\)]*)$", terms)
    m <- regmatches(terms, m)
    terms <- vapply(m, function(n){
        n <- n[seq.int(2L,length(n))]
        n <- n[n != ""]
        return(n)
        }, character(1))
    return(terms)
}

# This function performs dbRDA or CCA. It returns side scores with other
# information scores in attributes.
#' @importFrom stats as.formula na.fail
#' @importFrom vegan cca dbrda sppscores<- eigenvals scores
.calculate_rda <- function(
        x, formula, data, scores, scale = TRUE, na.action = na.fail,
        method = distance, distance = "euclidean", ord.method = "CCA", ...){
    # input check
    if(!.is_a_bool(scale)){
        stop("'scale' must be TRUE or FALSE.", call. = FALSE)
    }
    if( !.is_a_string(method) ){
        stop("'method' must be a single string value.", call. = FALSE)
    }
    if( !(.is_a_string(ord.method) && ord.method %in% c("CCA", "RDA")) ){
        stop("'ord.method' must be 'CCA' or 'RDA'.", call. = FALSE)
    }
    #
    # Get data in correct orientation. Samples should be in rows in abundance
    # table.
    x <- as.matrix(t(x))
    data <- data.frame(data, check.names = FALSE)
    # Instead of letting na.action pass through, give informative error
    # about missing values.
    if( any(is.na(data)) && isTRUE(all.equal(na.action, na.fail)) ){
        stop("Variables contain missing values. Set na.action to na.exclude",
            " to remove samples with missing values.", call. = FALSE)
    }
    # Ensure that the left hand side of the formula points to abundance matrix
    formula <- as.character(formula)
    formula[[2]] <- "x"
    # Create a formula from string
    formula <- as.formula(
        paste(as.character(formula)[c(2,1,3)], collapse = " "))
    
    # Initialize an argument list with common arguments
    args <- c(
        list(formula = formula, data = data, na.action = na.action), list(...))
    # vegan::cca function has a scale parameter that is not in dbrda
    if( ord.method == "CCA" ){
        args <- c(args, list(scale = scale))
    }
    # vegan::dbRDA has additional distance parameter that specifies
    # dissimilarity metric
    if( ord.method == "RDA" ){
        args <- c(args, list(distance = method))
    }
    # Perform CCA or RDA
    ord_FUN <- if (ord.method == "CCA") cca else dbrda
    res_obj <- do.call(ord_FUN, args)
    # The function stores the call. However, as we use do.call(), the function
    # sores the call incorrectly (e.g., it prints whole abundance table). That
    # is why we modify the call to have only necessary information.
    na_action_name <- deparse(substitute(na.action))
    na_action_name <- paste(na_action_name, collapse = " ")
    na_action_name <- sub('.*UseMethod\\("([^"]+)"\\).*', '\\1', na_action_name)
    res_obj$call <- paste0(
        res_obj$method, "(formula = ", deparse(formula),
        ", data = data",
        ifelse(ord.method == "CCA", paste0(", scale = ", scale), ""),
        ifelse(ord.method == "RDA", paste0(", distance = ", method), ""),
        ", na.action = ", na_action_name,  ", ...)")
    # Add species scores to rda object since they are missing
    # (in cca they are included). If we used na.action, some samples might be
    # removed. That is why we have to subset the abundance table first by
    # filtering missing values.
    if( ord.method == "RDA" ){
        if( !is.null(res_obj$na.action) ){
            x <- x[-res_obj$na.action, ]
        }
        sppscores(res_obj) <- x
    }
    
    # Get eigenvalues from the object
    eig <- eigenvals(res_obj)
    # Get total number of coordinates. There might be imaginary axes, exclude
    # them because otherwise error occurs in scores() call.
    tot_num_coord <- sum(eig >= 0)
    # Now we know how many coordinates there are. We use that info to get
    # scores for all coordinates.
    res_values <- vegan::scores(res_obj, seq_len(tot_num_coord))
    # Get the site scores. They represent samples' coordinates in ordinated
    # space.
    res_mat <- res_values[["sites"]]
    # Add other values as attributes
    res_values <- c(res_values, list(eig = eig))
    attributes(res_mat) <- c(
        attributes(res_mat), list(obj = res_obj), res_values)
    return(res_mat)
}

# Perform PERMANOVA and homogeneity analysis to RDA/CCA object
#' @importFrom vegan anova.cca betadisper
#' @importFrom stats anova
.test_rda <- function(mat, rda, variables, ...){
    # Perform permanova for whole model and for variables
    permanova_model <- anova.cca(rda, by = NULL)
    if( !is.null(variables) ){
        res <- .test_rda_vars(
            mat, rda, variables, permanova_model, ...)
    } else{
        res <- list(permanova = permanova_model)
    }
    return(res)
}

# Test association of variables to ordination
#' @importFrom vegan anova.cca
.test_rda_vars <- function(
        mat, rda, variables, permanova_model, by = "margin", full = FALSE, ...){
    # Check full parameter
    if( !.is_a_bool(full) ){
        stop("'full' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # Perform PERMANOVA
    permanova <- anova.cca(rda, by = by, ...)
    # Create a table from the results
    # PERMANOVAs
    table_model <- as.data.frame(permanova_model)
    permanova_tab <- as.data.frame(permanova)
    # Take only model results
    permanova_tab <- rbind(table_model[1, ], permanova_tab)
    # Add info about total variance
    permanova_tab[ , "Total variance"] <- permanova_tab["Model", 2] +
        permanova_tab["Residual", 2]
    # Add info about explained variance
    permanova_tab[ , "Explained variance"] <- permanova_tab[ , 2] /
        permanova_tab[ , "Total variance"]
    
    # Perform homogeneity analysis
    homogeneity <- .calculate_homogeneity(mat, variables, full = full, ...)
    
    # Return whole data or just a tables
    permanova_res <- permanova_tab
    if( full ){
        permanova_res <- list(
            summary = permanova_tab,
            model = permanova_model,
            variables = permanova)
    }
    res <- list(permanova = permanova_res, homogeneity = homogeneity)
    return(res)
}

#' @importFrom vegan vegdist betadisper
.calculate_homogeneity <- function(
        mat, variables, method = distance, distance = "euclidean",
        homogeneity.test = "permanova", full = FALSE, ...){
    # Check homogeneity.test
    if( !(.is_a_string(homogeneity.test) &&
            homogeneity.test %in% c("permanova", "anova", "tukeyhsd")) ){
        stop("'homogeneity.test' must be 'permanova', 'anova', or 'tukeyhsd'.",
            call. = FALSE)
    }
    # Check full parameter
    if( !.is_a_bool(full) ){
        stop("'full' must be TRUE or FALSE.", call. = FALSE)
    }
    # Check that there are more than one group per variable. If there are only
    # one group, we cannot calculate homogeneity for those variables.
    cols <- vapply(
        variables, function(x) length(unique(x))>1, FUN.VALUE = logical(1))
    if( any(!cols) ){
        warning("The following variables contain only one group. Cannot ",
            "calculate homogeneity for them: '",
            paste0(names(cols)[!cols], collapse = "','"), "'", call. = FALSE)
        variables <- variables[, cols, drop = FALSE]
    }
    #
    # Calculate dissimilarity matrix
    mat <- t(mat)
    diss_mat <- vegdist(mat, method = method, ...)
    # For all variables run the analysis
    homogeneity <- lapply(colnames(variables), function(x){
        # Get variable values
        var <- variables[[x]]
        # Run betadisper.
        # Suppress possible warnings: "some squared distances are negative and
        # changed to zero" Suppress possible messages:
        # "missing observations due to 'group' removed"
        suppressWarnings(
            suppressMessages(
                betadisper_res <- betadisper(diss_mat, group = var)
            )
        )
        # Run significance test
        significance <- .homogeneity_significance(
            betadisper_res, homogeneity.test, ...)
        # Return the results as a list
        models <- list(betadisper_res, significance[["obj"]])
        names(models) <- c("betadisper", homogeneity.test)
        res <- list(
            models = models,
            table = significance[["table"]]
        )
        return(res)
    })
    names(homogeneity) <- colnames(variables)
    # Get models
    homogeneity_model <- lapply(homogeneity, function(x) x[["models"]])
    names(homogeneity_model) <- colnames(variables)
    # Get tables
    homogeneity_tab <- lapply(homogeneity, function(x) x[["table"]])
    # Combine tables
    res <- do.call(rbind, homogeneity_tab)
    # Return either only summary table or whole result including fitterd models
    if( full ){
        res <- list(summary = res, variables = homogeneity_model)    
    }
    return(res)
}

# Perform statistical test for group homogeneity results
#' @importFrom stats anova TukeyHSD
#' @importFrom vegan permutest
.homogeneity_significance <- function(betadisper_res, homogeneity.test, ...){
    # Run specified significance test
    FUN <- switch(homogeneity.test,
        "permanova" = permutest,
        "anova" = anova,
        "tukeyhsd" = TukeyHSD
    )
    res <- do.call(FUN, c(list(betadisper_res), ...))
    # Get summary table from the results
    if( homogeneity.test == "anova" ){
        tab <- as.data.frame(res)
    } else if( homogeneity.test == "tukeyhsd" ){
        tab <- res[["group"]]
    } else{
        tab <- res[["tab"]]
    }
    # Modify permanova/anova table
    if( homogeneity.test != "tukeyhsd" ){
        # Add info about total variance
        tab[ , "Total variance"] <- tab["Groups", "Sum Sq"] +
            tab["Residuals", "Sum Sq"]
        # Add info about explained variance
        tab[ , "Explained variance"] <- tab[ , "Sum Sq"] /
            tab[ , "Total variance"]
        # Get only groups row (drop residuals row)
        tab <- tab[1, ]
    }
    # Create list from the object and results
    res <- list(
        obj = res,
        table = tab
    )
    return(res)
}

# Add RDA/CCA to reducedDim
.add_object_to_reduceddim <- function(
        tse, rda, name, subset.result = TRUE, ...){
    # Test subset
    if( !.is_a_bool(subset.result) ){
        stop("'subset.result' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # If samples do not match / there were samples without appropriate metadata
    # and they are now removed
    if( !all(colnames(tse) %in% rownames(rda)) && subset.result ){
        # Take a subset
        tse <- tse[ , rownames(rda) ]
        # Give a message
        warning(
            "Certain samples are removed from the result because they did ",
            "not include sufficient metadata.", call. = FALSE)
    } else if( !all(colnames(tse) %in% rownames(rda)) && !subset.result ){
        # If user do not want to subset the data
        # Save attributes from the object
        attr <- attributes(rda)
        attr <- attr[ !names(attr) %in% c("dim", "dimnames")]
        # Find samples that are removed
        samples_not_found <- setdiff(colnames(tse), rownames(rda))
        # Create an empty matrix
        mat <- matrix(nrow = length(samples_not_found), ncol=ncol(rda))
        rownames(mat) <- samples_not_found
        # Combine the data and order it in correct order
        rda <- rbind(rda, mat)
        rda <- rda[colnames(tse), ]
        # Add attributes
        attr <- c(attributes(rda), attr)
        attributes(rda) <- attr
    }
    # Add object to reducedDIm
    reducedDim(tse, name) <- rda
    return(tse)
}
