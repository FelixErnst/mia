#' Perform mediation analysis
#'
#' \code{getMediation} and \code{addMediation} provide a wrapper of
#' \code{\link[mediation:mediate]{mediate}} for
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
#'
#' @param x a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
#'
#' @param outcome \code{Character scalar}. Indicates the colData variable used
#'   as outcome in the model.
#'
#' @param treatment \code{Character scalar}. Indicates the colData variable
#'   used as treatment in the model.
#'
#' @param mediator \code{Character scalar}. Indicates the colData variable used
#'   as mediator in the model. (Default: \code{NULL})
#'
#' @param assay.type \code{Character scalar}. Specifies the assay used for
#'   feature-wise mediation analysis. (Default: \code{NULL})
#'
#' @param dimred \code{Character scalar}. Indicates the reduced dimension
#'   result in \code{reducedDims(object)} for component-wise mediation analysis.
#'   (Default: \code{NULL})
#'
#' @param family \code{Character scalar}. A specification for the outcome model
#' link function. (Default: \code{gaussian("identity")})
#'
#' @param covariates \code{Character scalar} or \code{character vector}.
#' Indicates the colData variables used as covariates in the model.
#' (Default: \code{NULL})
#'
#' @param p.adj.method \code{Character scalar}. Selects adjustment method
#'   of p-values. Passed to `p.adjust` function.
#'   (Default: \code{"holm"})
#'
#' @param add.metadata \code{Logical scalar}. Should the model metadata be
#'   returned. (Default: \code{TRUE})
#'
#' @param sort \code{Logical scalar}. Should the results be sorted by decreasing
#'   significance in terms of ACME_pval. (Default: \code{FALSE})
#'
#' @param verbose \code{Logical scalar}. Should execution messages be printed.
#'   (Default: \code{TRUE})
#'
#' @param name \code{Character scalar}. A name for the column of the
#'   \code{colData} where results will be stored. (Default: \code{"mediation"})
#'
#' @param ... additional parameters that can be passed to
#'   \code{\link[mediation:mediate]{mediate}}.
#'
#' @details
#' This wrapper of \code{\link[mediation:mediate]{mediate}} for
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' provides a simple method to analyse the effect of a treatment variable on an
#' outcome variable found in \code{colData(se)} through the mediation of either
#' another variable in colData (argument \code{mediator}) or an assay
#' (argument \code{assay.type}) or a reducedDim (argument \code{dimred}).
#' Importantly, those three arguments are mutually exclusive.
#'
#' @return
#' \code{getMediation} returns a \code{data.frame} of adjusted p-values and
#' effect sizes for the ACMEs and ADEs of every mediator given as input, whereas
#' \code{addMediation} returns an updated
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' instance with the same \code{data.frame} stored in the metadata as
#' "mediation" or as specified in the \code{name} argument. Its columns include:
#'
#' \describe{
#'   \item{mediator}{the mediator variable}
#'   \item{acme}{the Average Causal Mediation Effect (ACME) estimate}
#'   \item{acme_pval}{the original p-value for the ACME estimate}
#'   \item{acme_lower}{the lower bound of the CI for the ACME estimate}
#'   \item{acme_upper}{the upper bound of the CI for the ACME estimate}
#'   \item{ade}{the Average Direct Effect (ADE) estimate}
#'   \item{ade_pval}{the original p-value for the ADE estimate}
#'   \item{ade_lower}{the lower bound of the CI for the ADE estimate}
#'   \item{ade_upper}{the upper bound of the CI for the ADE estimate}
#'   \item{total}{the Total Effect estimate}
#'   \item{total_lower}{the lower bound of the CI for the Total Effect estimate}
#'   \item{total_upper}{the upper bound of the CI for the Total Effect estimate}
#'   \item{total_pval}{the original p-value for the Total Effect estimate}
#'   \item{acme_padj}{the adjusted p-value for the ACME estimate}
#'   \item{ade_padj}{the adjusted p-value for the ADE estimate}
#'   \item{total_padj}{the adjusted p-value for the Total Effect estimate}
#' }
#'
#' The original output of \code{\link[mediation:mediate]{mediate}} for each
#' analysed mediator is stored as the "model_metadata" attribute of the
#' resulting \code{data.frame} and can be accessed via the \code{attr} function.
#'
#' @name getMediation
#'
#' @examples
#' \dontrun{
#' # Import libraries
#' library(mia)
#' library(miaViz)
#' library(scater)
#'
#' # Load dataset
#' data(hitchip1006, package = "miaTime")
#' tse <- hitchip1006
#'
#' # Agglomerate features by family (merely to speed up execution)
#' tse <- agglomerateByRank(tse, rank = "Phylum")
#' # Convert BMI variable to numeric
#' tse$bmi_group <- as.numeric(tse$bmi_group)
#'
#' # Analyse mediated effect of nationality on BMI via alpha diversity
#' # 100 permutations were done to speed up execution, but ~1000 are recommended
#' med_df <- getMediation(
#'     tse,
#'     outcome = "bmi_group",
#'     treatment = "nationality",
#'     mediator = "diversity",
#'     covariates = c("sex", "age"),
#'     treat.value = "Scandinavia",
#'     control.value = "CentralEurope",
#'     boot = TRUE, sims = 100)
#'
#' # Visualise model statistics
#' plotMediation(med_df)
#'
#' # Apply clr transformation to counts assay
#' tse <- transformAssay(tse, method = "clr", pseudocount = 1)
#'
#' # Analyse mediated effect of nationality on BMI via clr-transformed features
#' # 100 permutations were done to speed up execution, but ~1000 are recommended
#' tse <- addMediation(
#'     tse, name = "assay_mediation",
#'     outcome = "bmi_group",
#'     treatment = "nationality",
#'     assay.type = "clr",
#'     covariates = c("sex", "age"),
#'     treat.value = "Scandinavia",
#'     control.value = "CentralEurope",
#'     boot = TRUE, sims = 100,
#'     p.adj.method = "fdr")
#'
#' # Show results for first 5 mediators
#' head(metadata(tse)$assay_mediation, 5)
#'
#' # Perform ordination
#' tse <- runMDS(
#'     tse, name = "MDS", method = "euclidean",
#'     assay.type = "clr", ncomponents = 3)
#'
#' # Analyse mediated effect of nationality on BMI via NMDS components
#' # 100 permutations were done to speed up execution, but ~1000 are recommended
#' tse <- addMediation(
#'     tse, name = "reddim_mediation",
#'     outcome = "bmi_group",
#'     treatment = "nationality",
#'     dimred = "MDS",
#'     covariates = c("sex", "age"),
#'     treat.value = "Scandinavia",
#'     control.value = "CentralEurope",
#'     boot = TRUE, sims = 100,
#'     p.adj.method = "fdr")
#'
#' # Show results for first 5 mediators
#' head(metadata(tse)$reddim_mediation, 5)
#'
#' # Access model metadata
#' attr(metadata(tse)$reddim_mediation, "model_metadata")
#' }
#'
NULL

#' @rdname getMediation
#' @export
#' @importFrom stats gaussian
setMethod("addMediation", signature = c(x = "SummarizedExperiment"),
        function(x, outcome, treatment, name = "mediation",
            mediator = NULL, assay.type = NULL, dimred = NULL,
            family = gaussian(), covariates = NULL, p.adj.method = "holm",
            add.metadata = TRUE, verbose = TRUE, ...) {

            med_df <- getMediation(
                x, outcome, treatment,
                mediator, assay.type, dimred,
                family, covariates, p.adj.method,
                add.metadata, verbose, ...
            )

            x <- .add_values_to_metadata(x, name, med_df)
            return(x)
        }
)

#' @rdname getMediation
#' @export
#' @importFrom stats gaussian
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
setMethod("getMediation", signature = c(x = "SummarizedExperiment"),
        function(x, outcome, treatment,
            mediator = NULL, assay.type = NULL, dimred = NULL,
            family = gaussian(), covariates = NULL, p.adj.method = "holm",
            add.metadata = TRUE, sort = FALSE, verbose = TRUE, ...) {
        ###################### Input check ########################
        if( !outcome %in% names(colData(x)) ){
            stop(outcome, " not found in colData(x).", call. = FALSE)
        }
        if( !treatment %in% names(colData(x)) ){
            stop(treatment, " not found in colData(x).", call. = FALSE)
        }
        if( !is.null(covariates) && !all(covariates %in% names(colData(x))) ){
            stop("covariates not found in colData(x).", call. = FALSE)
        }
        if( !.is_a_bool(add.metadata) ){
            stop("add.metadata must be TRUE or FALSE.", call. = FALSE)
        }
        if( !.is_a_bool(verbose) ){
            stop("verbose must be TRUE or FALSE.", call. = FALSE)
        }

        # Check that arguments can be passed to mediate and filter samples
        x <- .check_mediate_args(
            x, outcome, treatment, mediator, covariates, verbose, ...
        )

        # Check which mediator was provided (colData, assay or reducedDim)
        med_opts <- unlist(lapply(
            list(mediator = mediator, assay.type = assay.type, dimred = dimred),
            function(x) !is.null(x)
        ))

        if( sum(med_opts) != 1 ){
            # Throw error if none or multiple mediation options are specified
            stop(
                "The arguments mediator, assay.type and dimred are mutually ",
                "exclusive, but ", sum(med_opts), " were provided.",
                call. = FALSE
            )
        }

        if ( med_opts[[1]] ){
            # Check that mediator is in colData
            if( !mediator %in% names(colData(x)) ) {
                stop(mediator, " not found in colData(x).", call. = FALSE)
            }
            # Use mediator for analysis
            mediators <- mediator
            mat <- NULL
        } else if( med_opts[[2]] ){
            # Check that assay is in assays
            .check_assay_present(assay.type, x)
            # Define matrix for analysis
            mat <- assay(x, assay.type)
            # Use assay for analysis
            mediators <- rownames(mat)
        } else if( med_opts[[3]] ){
            # Check that reducedDim is in reducedDims
            if(!dimred %in% reducedDimNames(x)){
                stop(dimred, " not found in reducedDims(x).", call. = FALSE)
            }
            # Define matrix for analysis
            mat <- t(reducedDim(x, dimred))
            # Give component names to matrix rows
            rownames(mat) <- paste0(dimred, seq(1, nrow(mat)))
            # Use reducedDim for analysis
            mediators <- rownames(mat)
        }

        # Create template list of models
        models <- list()
        # Set initial index
        i <- 0

        for( mediator in mediators ){
            # Update index
            i <- i + 1
            if( verbose ){
                message("\rMediator ", i, " out of ",
                        length(mediators), ": ", mediator
                )
            }
            # Run mediation analysis for current mediator
            med_out <- .run_mediate(
                x, outcome, treatment, mediator,
                family = family, mat = mat,
                covariates = covariates, ...
            )
            # Update list of models
            models <- c(models, list(med_out))
        }

        # Name models by mediators
        names(models) <- mediators
        # Combine results into dataframe
        med_df <- .make_output(models, p.adj.method, add.metadata, sort)

        return(med_df)
    }
)

# Check that arguments can be passed to mediate and remove unused samples
#' @importFrom stats na.omit
.check_mediate_args <- function(
    x, outcome, treatment, mediator, covariates, verbose = TRUE, ...) {

    # Create dataframe from selected columns of colData
    df <- as.data.frame(colData(x)[ , names(colData(x)) %in% c(
        outcome, treatment, mediator, covariates)])
    # Store kwargs into variable
    kwargs <- list(...)

    # Remove missing data from df
    df <- na.omit(df)
    diff <- ncol(x) - nrow(df)

    if( diff != 0 ){
        # Remove missing data from se
        x <- x[ , rownames(df)]

        if( verbose ){
            message(diff, " samples removed because of missing data.")
        }
    }

    # If boot is TRUE and treatment variable is discrete and has 3+ levels
    if( !is.null(kwargs[["boot"]]) && !is.numeric(df[[treatment]]) &&
        length(unique((df[[treatment]]))) > 2 ) {

        ## if control and treat value are not specified
        if( any(vapply(kwargs[c("control.value", "treat.value")],
                is.null, logical(1))) ){
            stop(
                "Too many treatment levels. Consider specifing a treat.value ",
                "and a control.value", call. = FALSE
            )
        }
        ## but if they are specified
        # if they appear in the treatment variable
        if( !all(kwargs[c("control.value", "treat.value")] %in%
                unique(df[[treatment]])) ){
            stop(
                "treat.value and/or control.value not found in the levels of ",
                "the treatment variable.", call. = FALSE
            )
        }

        # Find indices of samples that belong to either control or treatment
        keep <- df[[treatment]] %in% kwargs[c("control.value", "treat.value")]

        # Remove samples different from control and treatment from df
        df <- df[keep, ]
        diff <- ncol(x) - nrow(df)

        # Remove samples different from control and treatment from se
        x <- x[ , rownames(df)]

        if( verbose ){
            message(
                diff, " samples removed because different ",
                "from control and treatment."
            )
        }
    }
    return(x)
}

# Run mediation analysis
#' @importFrom stats lm formula glm
.run_mediate <- function(x, outcome, treatment, mediator = NULL, mat = NULL,
                        family = gaussian(), covariates = NULL, ...) {
    .require_package("mediation")
    # Create initial dataframe with outcome and treatment variables
    df <- data.frame(
        Outcome = colData(x)[[outcome]], Treatment = colData(x)[[treatment]])

    if( is.null(mat) ){
        # If matrix not given, fetch mediator from colData
        df[["Mediator"]] <- colData(x)[[mediator]]
    } else {
        # If matrix given, use it as mediators
        df[["Mediator"]] <- mat[mediator, ]
    }

    # Define basic formula mediation model
    relation_m <- "Mediator ~ Treatment"
    # Define basic formula outcome model
    relation_dv <- "Outcome ~ Treatment + Mediator"

    if( !is.null(covariates) ){
        # Fetch covariates from colData and store them in dataframe
        df <- cbind(df, colData(x)[covariates])
        # Add covariate to formula of mediation model
        relation_m <- paste(
            relation_m, "+",
            paste(covariates,collapse = " + ")
        )
        # Add covariate to formula of outcome model
        relation_dv <- paste(
            relation_dv, "+",
            paste(covariates, collapse = " + ")
        )
    }

    # Fit mediation model
    fit_m <- do.call(lm, list(formula = formula(relation_m), data = df))
    # Fit outcome model
    fit_dv <- do.call(
        glm,
        list(formula = formula(relation_dv), family = family, data = df)
    )
    # Run mediation analysis
    med_out <- mediation::mediate(
        fit_m, fit_dv,
        treat = "Treatment", mediator = "Mediator",
        covariates = covariates, ...
    )

    return(med_out)
}

# Combine results into output dataframe
#' @importFrom dplyr select
#' @importFrom tidyr starts_with ends_with unnest_wider
#' @importFrom stringr str_extract str_replace_all
#' @importFrom stats p.adjust
.make_output <- function(models, p.adj.method, add.metadata, sort) {
    # Combine results
    res <- do.call(rbind, models) |> as.data.frame()
    res[["mediator"]] <- names(models)
    # Select certain data types
    res <- res |>
        select(mediator, starts_with(c("d.avg", "z.avg", "tau"))) |>
        select(-ends_with(c("sims")))
    # Get columns with scalar values and turn them into vector instead of list
    cols <- vapply(res, function(col) all(lengths(col) == 1L), logical(1L))
    res[, cols] <- lapply(res[, cols], unlist)
    # Unnest confidence interval columns
    res <- res |> unnest_wider(ends_with(".ci"), names_sep = "_")
    # Replace confidence interval values with "lower" and "upper"
    limits <- str_extract(colnames(res), "([0-9.]+)%") |> unique()
    limits <- limits[!is.na(limits)]
    limits <- limits[order(as.numeric(gsub("%", "", limits)))]
    lookup <- c("lower", "upper")
    names(lookup) <- limits
    colnames(res) <- str_replace_all(colnames(res), lookup)
    # Tidy other names also
    lookup <- c("d.avg" = "acme", "z.avg" = "ade", "tau" = "total",
        "\\.p" = "_pval", "\\.ci_" = "_", "\\.coef" = "")
    colnames(res) <- str_replace_all(colnames(res), lookup)

    # Compute adjusted p-values and add them to dataframe
    cols <- endsWith(colnames(res), "pval")
    pcols <- lapply(res[ , cols], p.adjust, method = p.adj.method)
    pcols <- do.call(cbind, pcols)
    colnames(pcols) <- gsub("pval", "padj", colnames(pcols))
    res <- cbind(res, pcols)

    # Store model for each mediator into the model_metadata attribute
    if( add.metadata ){
        attr(res, "model_metadata") <- models
    }
    # Order output dataframe by ACME p-values
    if( sort ){
        res <- res[order(res[["acme_padj"]]), ]
    }
    return(res)
}
