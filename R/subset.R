.get_subset_args <- function(x, subset = NULL, select = NULL, ...){
    rows <- subset
    columns <- select
    if(is.null(rows)){
        rows <- rep(TRUE,nrow(x))
    }
    if(is.null(columns)){
        columns <- rep(TRUE,ncol(x))
    }
    return(list(rows = rows, columns = columns))
}
