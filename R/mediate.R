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
#' returned. (Default: \code{FALSE})
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
#' "mediation". Its columns include:
#' 
#' \describe{
#'   \item{Mediator}{the mediator variable}
#'   \item{ACME_estimate}{the Average Causal Mediation Effect (ACME) estimate}
#'   \item{ADE_estimate}{the Average Direct Effect (ADE) estimate}
#'   \item{ACME_pval}{the adjusted p-value for the ACME estimate}
#'   \item{ADE_pval}{the adjusted p-value for the ADE estimate}
#' }
#'
#' @name getMediation
#'
#' @examples
#' \dontrun{
#' # Import libraries
#' library(mia)
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
#' med_df <- getMediation(tse,
#'                        outcome = "bmi_group",
#'                        treatment = "nationality",
#'                        mediator = "diversity",
#'                        covariates = c("sex", "age"),
#'                        treat.value = "Scandinavia",
#'                        control.value = "CentralEurope",
#'                        boot = TRUE, sims = 100,
#'                        add.metadata = TRUE)
#'  
#'  # Visualise model statistics for 1st mediator
#'  plot(attr(med_df, "metadata")[[1]])
#' 
#' # Apply clr transformation to counts assay
#' tse <- transformAssay(tse,
#'                       method = "clr",
#'                       pseudocount = 1)
#'
#' # Analyse mediated effect of nationality on BMI via clr-transformed features
#' # 100 permutations were done to speed up execution, but ~1000 are recommended
#' tse <- addMediation(tse, name = "assay_mediation",
#'                     outcome = "bmi_group",
#'                     treatment = "nationality",
#'                     assay.type = "clr",
#'                     covariates = c("sex", "age"),
#'                     treat.value = "Scandinavia",
#'                     control.value = "CentralEurope",
#'                     boot = TRUE, sims = 100,
#'                     p.adj.method = "fdr")
#' 
#' # Show results for first 5 mediators
#' head(metadata(tse)$assay_mediation, 5)
#'  
#' # Perform ordination
#' tse <- runMDS(tse, name = "MDS",
#'               method = "euclidean",
#'               assay.type = "clr",
#'               ncomponents = 3)
#' 
#' # Analyse mediated effect of nationality on BMI via NMDS components
#' # 100 permutations were done to speed up execution, but ~1000 are recommended
#' tse <- addMediation(tse, name = "reddim_mediation",
#'                     outcome = "bmi_group",
#'                     treatment = "nationality",
#'                     dimred = "MDS",
#'                     covariates = c("sex", "age"),
#'                     treat.value = "Scandinavia",
#'                     control.value = "CentralEurope",
#'                     boot = TRUE, sims = 100,
#'                     p.adj.method = "fdr")
#' 
#' # Show results for first 5 mediators
#' head(metadata(tse)$reddim_mediation, 5)
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
            add.metadata = FALSE, verbose = TRUE, ...) {
            
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
            add.metadata = FALSE, verbose = TRUE, ...) {
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
        x <- .check.mediate.args(
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
        
        # Create template list of results
        results <- list(
            Mediator = c(), ACME_estimate = c(), ADE_estimate = c(),
            ACME_pval = c(), ADE_pval = c(), Model = list()
        )
        
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
            # Update list of results
            results <- .update.results(results, med_out, mediator)
        }
        # Combine results into dataframe
        med_df <- .make.output(results, p.adj.method, add.metadata)
        return(med_df)
    }
)

# Check that arguments can be passed to mediate and remove unused samples
#' @importFrom stats na.omit
.check.mediate.args <- function(
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
#' @importFrom mediation mediate
#' @importFrom stats lm formula glm
.run_mediate <- function(x, outcome, treatment, mediator = NULL, mat = NULL,
                        family = gaussian(), covariates = NULL, ...) {
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
    med_out <- mediate(
        fit_m, fit_dv,
        treat = "Treatment", mediator = "Mediator",
        covariates = covariates, ...
    )
    return(med_out)
}


# Update list of results
.update.results <- function(results, med_out, mediator) {
    # Update model variables
    results[["Mediator"]] <- c(results[["Mediator"]], mediator)
    # Update stats of ACME (average causal mediation effect)
    results[["ACME_estimate"]] <- c(results[["ACME_estimate"]], med_out$d.avg)
    results[["ACME_pval"]] <- c(results[["ACME_pval"]], med_out$d.avg.p)
    # Update stats of ADE (average direct effect)
    results[["ADE_estimate"]] <- c(results[["ADE_estimate"]], med_out$z.avg)
    results[["ADE_pval"]] <- c(results[["ADE_pval"]], med_out$z.avg.p)
    # Add current model to metadata
    results[["Model"]][[length(results[["Model"]]) + 1]] <- med_out
    return(results)
}

# Combine results into output dataframe
.make.output <- function(results, p.adj.method, add.metadata) {
    
    # Create dataframe with model variables, effect sizes and p-values
    med_df <- do.call(data.frame, results[seq_len(length(results) - 1)])
    
    # Compute adjusted p-values and add them to dataframe
    med_df[["ACME_pval"]] <- p.adjust(
        med_df[["ACME_pval"]],
        method = p.adj.method
    )
    med_df[["ADE_pval"]] <- p.adjust(
        med_df[["ADE_pval"]],
        method = p.adj.method
    )
    
    if( add.metadata ){
        # models for every mediator are saved into the metadata attribute
        attr(med_df, "metadata") <- results[["Model"]]
    }
    
    # Order output dataframe by ACME p-values
    med_df <- med_df[order(med_df[["ACME_pval"]]), ]
    
    return(med_df)
}
