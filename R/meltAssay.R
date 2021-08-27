#' Converting a \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object into a long data.frame
#'
#' \code{meltAssay} Converts a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} object into a
#' long data.frame which can be used for \code{tidyverse}-tools.
#'
#' @details
#' If the \code{colData} contains a column \dQuote{SampleID} or the
#' \code{rowData} contains a column \dQuote{FeatureID}, they will be renamed to
#' \dQuote{SampleID_col} and \dQuote{FeatureID_row}, if row names or column
#' names are set.
#'
#' @param x A numeric matrix or a
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#'   
#' @param abund_values a \code{character} value to select an
#'   \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assayNames}}
#'
#' @param add_col_data \code{NULL}, \code{TRUE} or a \code{character} vector to
#'   select information from the \code{colData} to add to the molten assay data.
#'   If \code{add_col_data = NULL} no data will be added, if
#'   \code{add_col_data = TRUE} all data will be added and if
#'   \code{add_col_data} is a \code{character} vector, it will be used to subset
#'   to given column names in \code{colData}. (default:
#'   \code{add_col_data = NULL})
#'
#' @param add_row_data \code{NULL}, \code{TRUE} or a \code{character} vector to
#'   select information from the \code{rowData} to add to the molten assay data.
#'   If \code{add_row_data = NULL} no data will be added, if
#'   \code{add_row_data = TRUE} all data will be added and if
#'   \code{add_row_data} is a \code{character} vector, it will be used to subset
#'   to given column names in \code{rowData}. (default:
#'   \code{add_row_data = NULL})
#'
#' @param feature_name a \code{character} scalar to use as the output's name
#'   for the feature identifier. (default: \code{feature_name = "FeatureID"})
#'
#' @param sample_name a \code{character} scalar to use as the output's name
#'   for the sample identifier. (default: \code{sample_name = "SampleID"})
#'
#' @param ... optional arguments:
#' \itemize{
#'   \item{check_names}{ A boolean value passed to data.frame function's check.name
#'   argument. Determines if sample names are checked that they are syntactically 
#'   valid variable names and are not duplicated. If they are not, sample names 
#'   are modified. (default: \code{check_names = TRUE})}
#' }
#'
#' @return A \code{tibble} with the molten data. The assay values are given in a
#' column named like the selected assay \code{abund_values}. In addition, a
#' column \dQuote{FeatureID} will contain the rownames, if set, and analogously
#' a column \dQuote{SampleID} with the colnames, if set
#'
#' @name meltAssay

#' @author
#' Sudarshan A. Shetty
#'
#' @examples
#' data(GlobalPatterns)
#' molten_se <- meltAssay(GlobalPatterns,
#'                        abund_values = "counts",
#'                        add_row_data = TRUE,
#'                        add_col_data = TRUE
#'                        )
#' molten_se
NULL

#' @rdname meltAssay
#' @export
setGeneric("meltAssay",
           signature = "x",
           function(x,
                    abund_values = "counts",
                    add_row_data = NULL,
                    add_col_data = NULL,
                    feature_name = "FeatureID",
                    sample_name = "SampleID",
                    ...)
               standardGeneric("meltAssay")
)

.norm_add_row_data <- function(add_row_data, x, feature_name){
    if(is.null(add_row_data)){
        return(NULL)
    }
    if(anyNA(add_row_data)){
        stop("'add_row_data' contains NA.", call. = FALSE)
    }
    cn <- colnames(rowData(x))
    if(is.logical(add_row_data) && length(add_row_data) == 1L && add_row_data){
        add_row_data <- cn
    } else if (isFALSE(all(add_row_data %in% cn))) {
        stop("Please provide valid column names with 'add_row_data' matching ",
             "those in 'rowData(x)'", call. = FALSE)
    }
    if(!is.null(rownames(x)) && feature_name %in% add_row_data){
        warning("'x' contains a column '",feature_name,"' in its ",
                "rowData(), which will ",
                "be renamed to '",feature_name,"_row'", call. = FALSE)
    }
    add_row_data
}

.norm_add_col_data <- function(add_col_data, x, sample_name){
    if(is.null(add_col_data)){
        return(NULL)
    }
    if(anyNA(add_col_data)){
        stop("'add_col_data' contains NA.", call. = FALSE)
    }
    cn <- colnames(colData(x))
    if(is.logical(add_col_data) && length(add_col_data) == 1L && add_col_data){
        add_col_data <- cn
    } else if (isFALSE(all(add_col_data %in% cn))) {
        stop("Please provide valid column names with 'add_col_data' matching ",
             "those in 'colData(x)'", call. = FALSE)
    }
    if(!is.null(colnames(x)) && sample_name %in% add_col_data){
        warning("'x' contains a column '",sample_name,"' in its ",
                "colData(), which will ",
                "be renamed to '",sample_name,"_col'", call. = FALSE)
    }
    add_col_data
}

.col_switch_name <- function(name){
    paste0(name,"_col")
}

.row_switch_name <- function(name){
    paste0(name,"_row")
}

#' @importFrom dplyr mutate select
.format_molten_assay <- function(molten_assay, x,
                                 feature_name,
                                 sample_name){
    if(is.null(rownames(x)) &&
       .row_switch_name(feature_name) %in% colnames(molten_assay) &&
       !anyNA(molten_assay[,.row_switch_name(feature_name)]) &&
       !anyDuplicated(rowData(x)[,feature_name])){
        molten_assay <- molten_assay %>%
            select(!sym(feature_name)) %>%
            dplyr::rename(!!sym(feature_name) := !!sym(.row_switch_name(feature_name)))
    }
    if(is.null(colnames(x)) &&
       .col_switch_name(sample_name) %in% colnames(molten_assay) &&
       !anyNA(molten_assay[,.col_switch_name(sample_name)]) &&
       !anyDuplicated(colData(x)[,sample_name])){
        molten_assay %>%
            select(!sym(sample_name)) %>%
            dplyr::rename(!!sym(sample_name) := !!sym(.col_switch_name(sample_name)))
    }
    molten_assay %>%
        mutate(!!sym(feature_name) := factor(!!sym(feature_name)),
               !!sym(sample_name) := factor(!!sym(sample_name)))
}


#' @rdname meltAssay
#'
#' @export
setMethod("meltAssay", signature = c(x = "SummarizedExperiment"),
    function(x,
             abund_values = "counts", 
             add_row_data = NULL,
             add_col_data = NULL,
             feature_name = "FeatureID",
             sample_name = "SampleID",
             ...) {
        # input check
        .check_assay_present(abund_values, x)
        if(!.is_a_string(feature_name)){
            stop("'feature_name' must be a single non-empty character value.",
                 call. = FALSE)
        }
        if(!.is_a_string(sample_name)){
            stop("'sample_name' must be a single non-empty character value.",
                 call. = FALSE)
        }
        # check selected colnames
        add_row_data <- .norm_add_row_data(add_row_data, x, feature_name)
        add_col_data <- .norm_add_col_data(add_col_data, x, sample_name)
        molten_assay <- .melt_assay(x, abund_values, feature_name, sample_name, ...)
        if(!is.null(add_row_data)){
            molten_assay <-
                .add_row_data_to_molten_assay(molten_assay, x, add_row_data,
                                              feature_name)
        }
        if(!is.null(add_col_data)){
            molten_assay <-
                .add_col_data_to_molten_assay(molten_assay, x, add_col_data,
                                              sample_name, ...)
        }
        .format_molten_assay(molten_assay, x, feature_name, sample_name)
    }
)

# Melts assay in SummarizedExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
.melt_assay <- function(x, abund_values, feature_name, sample_name, check_names = FALSE) {
    assay(x, abund_values) %>%
        data.frame(check.names = check_names) %>%
        rownames_to_column(feature_name) %>%
        # SampleID is unique sample id
        pivot_longer(!sym(feature_name),
                     values_to = abund_values,
                     names_to = sample_name)
}

# Combines molten assay with rowData i.e. taxonomy table.
#' @importFrom SummarizedExperiment rowData
#' @importFrom rlang sym
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr rename left_join
.add_row_data_to_molten_assay <- function(molten_assay, x, add_row_data,
                                          feature_name) {
    rd <- SummarizedExperiment::rowData(x)[,add_row_data,drop=FALSE] %>%
        data.frame()
    if(feature_name %in% add_row_data){
        rd <- rd %>%
            dplyr::rename(!!sym(.row_switch_name(feature_name)) := !!sym(feature_name))
    }
    rd <- rd %>%
        rownames_to_column(feature_name)
    molten_assay %>%
        left_join(rd, by = feature_name)
}

# Combines molten assay and rowData i.e. taxonomy table with
#' @importFrom SummarizedExperiment colData
#' @importFrom rlang sym
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr rename left_join
.add_col_data_to_molten_assay <- function(molten_assay, x, add_col_data,
                                          sample_name, check_names = FALSE) {
    cd <- SummarizedExperiment::colData(x)[,add_col_data,drop=FALSE] %>%
        data.frame()
    # This makes sure that sample names match
    if(check_names == TRUE){
        rownames(cd) <- make.names(rownames(cd))
    }
    if(sample_name %in% add_col_data){
        cd <- cd %>%
            dplyr::rename(!!sym(.col_switch_name(sample_name)) := !!sym(sample_name))
    }
    cd <- cd %>%
        rownames_to_column(sample_name)
    molten_assay %>%
        left_join(cd, by = sample_name)
}
