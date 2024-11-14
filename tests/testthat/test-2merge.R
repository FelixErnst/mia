context("merge")
test_that("merge", {
    # .check_f
    expect_error(mia:::.norm_f(),
                 'argument "f" is missing')
    expect_error(mia:::.norm_f(6),
                 'argument "f" is missing')
    expect_error(mia:::.norm_f(6,5),
                 "'f' must be a factor or character vector")
    f <- factor(c(rep("a",3),rep("b",3)))
    expect_true(is.factor(mia:::.norm_f(6,f)))
    # .check_archetype
    expect_error(mia:::.norm_archetype(),
                 'argument "archetype" is missing')
    expect_error(mia:::.norm_archetype(f),
                 'argument "archetype" is missing')
    expect_error(mia:::.norm_archetype(f),
                 'argument "archetype" is missing')
    expect_equal(mia:::.norm_archetype(f, 1),c(1,1))
    expect_equal(mia:::.norm_archetype(f, c(1,2)),c(1,2))
    expect_error(mia:::.norm_archetype(f, c(1,2,3)),
                 "length of 'archetype' must have the same length as levels")
    expect_error(mia:::.norm_archetype(f, c(5)),
                 "'archetype' out of bounds for some levels of 'f'")
    # .norm_archetype
    expect_error(mia:::.norm_archetype(),
                 'argument "archetype" is missing')
    expect_error(mia:::.norm_archetype(f),
                 'argument "archetype" is missing')
    actual <- mia:::.norm_archetype(f, c(1,2))
    expect_equal(actual, c(1,2))
    actual <- mia:::.norm_archetype(f, c(1))
    expect_equal(actual, c(1,1))

    # .get_element_pos
    expect_error(mia:::.get_element_pos(),
                 'argument "archetype" is missing')
    expect_error(mia:::.get_element_pos(f),
                 'argument "archetype" is missing')

    actual <- mia:::.get_element_pos(f, archetype = mia:::.norm_archetype(f, 1))
    expect_equal(actual,c(a = 1, b = 4))
    actual <- mia:::.get_element_pos(f, archetype = mia:::.norm_archetype(f, 2))
    expect_equal(actual,c(a = 2, b = 5))
    actual <- mia:::.get_element_pos(f, archetype = c(2,1))
    expect_equal(actual,c(a = 2, b = 4))

    # .merge_rows_or_cols
    mat <- matrix(1:60, nrow = 6)
    gr <- GRanges("chr1",rep("1-6",6))
    df <- DataFrame(n = c(1:6))
    mcols(gr) <- df
    grl <- splitAsList(gr,1:6)
    expect_error(mia:::.merge_rows_or_cols(by = 1L),
                 'argument "f" is missing')
    x <- SummarizedExperiment(assays = list(mat = mat))
    xr <- SummarizedExperiment(assays = list(mat = mat),
                               rowRanges = gr)
    xrl <- SummarizedExperiment(assays = list(mat = mat),
                                rowRanges = unname(grl))
    expect_error(mia:::.merge_rows_or_cols(x, by = 1L),
                 'argument "f" is missing')
    FUN_check_x <- function(x,archetype=1){
        actual <- agglomerateByVariable(x, by = "rows", f, archetype)
        expect_s4_class(actual,class(x))
        expect_equal(dim(actual),c(2,10))
    }
    lapply(list(x,xr,xrl),FUN_check_x)
    lapply(list(x,xr,xrl),FUN_check_x,archetype=2)
    #
    f <- factor(c(rep("a",3),rep("b",3)))
    mat <- matrix(1:60, nrow = 6)
    gr <- GRanges("chr1",rep("1-6",6))
    df <- DataFrame(n = c(1:6))
    mcols(gr) <- df
    grl <- splitAsList(gr,1:6)
    xtse <- TreeSummarizedExperiment(assays = list(mat = mat),
                                     rowRanges = unname(grl))
    FUN_check_x <- function(x,archetype=1){
        actual <- agglomerateByVariable(x, by = "rows", f, archetype, 
            update.tree = FALSE)
        expect_s4_class(actual,class(x))
        expect_equal(dim(actual),c(2,10))
    }
    lapply(list(xtse),FUN_check_x)
    lapply(list(xtse),FUN_check_x,archetype=2)
    
    # Check that average works as expected. average parameter controls whether
    # to calculate mean or sum. Check that mean is correctly calculated when
    # there are NAs
    #
    # Calculate average and sum for each row group
    summary_FUN_rows <- function(x, col.var){
        # Loop through groups and calculate statistics
        groups <- unique(rowData(x)[[col.var]]) |> sort()
        res <- lapply(groups, function(group) {
            mat_sub <- assay(x[rowData(x)[[col.var]] == group, ])
            list(
                sum = colSums(mat_sub, na.rm = FALSE),
                sum_na = colSums(mat_sub, na.rm = TRUE),
                mean = colMeans(mat_sub, na.rm = FALSE),
                mean_na = colMeans(mat_sub, na.rm = TRUE)
            )
        })
        # Combine results for each statistic across groups
        res <- lapply(c("sum", "sum_na", "mean", "mean_na"), function(stat) {
            do.call(rbind, lapply(res, `[[`, stat))
        })
        names(res) <- c("sum", "sum_na", "mean", "mean_na")
        return(res)
    }
    # Generate data
    tse <- mockSCE()
    rowData(tse)[["group"]] <- sample(LETTERS, nrow(tse), replace = TRUE)
    colData(tse)[["group"]] <- sample(LETTERS, ncol(tse), replace=TRUE)
    # Create a data with NAs
    n_value <- nrow(tse)*ncol(tse)
    assay(tse)[c(1, 5, 3, 6)] <- NA
    # Test with NAs
    res_sum <- agglomerateByVariable(tse, by = 1, group = "group", average = FALSE, na.rm = FALSE)
    res_sum_na <- agglomerateByVariable(tse, by = 1, group = "group", average = FALSE, na.rm = TRUE)
    res_mean <- agglomerateByVariable(tse, by = 1, group = "group", average = TRUE, na.rm = FALSE)
    res_mean_na <- agglomerateByVariable(tse, by = 1, group = "group", average = TRUE, na.rm = TRUE)
    ref <- summary_FUN_rows(tse, "group")
    #
    expect_equal(assay(res_sum), ref[["sum"]], check.attributes = FALSE)
    expect_equal(assay(res_sum_na), ref[["sum_na"]], check.attributes = FALSE)
    expect_equal(assay(res_mean), ref[["mean"]], check.attributes = FALSE)
    expect_equal(assay(res_mean_na), ref[["mean_na"]], check.attributes = FALSE)
    # Calculate average and sum for each column group
    summary_FUN_cols <- function(x, col.var){
        # Loop through groups and calculate statistics
        groups <- unique(colData(x)[[col.var]]) |> sort()
        res <- lapply(groups, function(group){
            mat_sub <- assay(x[, colData(x)[[col.var]] == group ])
            list(
                sum = rowSums(mat_sub, na.rm = FALSE),
                sum_na = rowSums(mat_sub, na.rm = TRUE),
                mean = rowMeans(mat_sub, na.rm = FALSE),
                mean_na = rowMeans(mat_sub, na.rm = TRUE)
            )
        })
        # Combine results for each statistic across groups
        res <- lapply(c("sum", "sum_na", "mean", "mean_na"), function(stat){
            do.call(cbind, lapply(res, `[[`, stat))
        })
        names(res) <- c("sum", "sum_na", "mean", "mean_na")
        return(res)
    }
    # Test with NAs
    res_sum <- agglomerateByVariable(tse, by = 2, group = "group", average = FALSE, na.rm = FALSE)
    res_sum_na <- agglomerateByVariable(tse, by = 2, group = "group", average = FALSE, na.rm = TRUE)
    res_mean <- agglomerateByVariable(tse, by = 2, group = "group", average = TRUE, na.rm = FALSE)
    res_mean_na <- agglomerateByVariable(tse, by = 2, group = "group", average = TRUE, na.rm = TRUE)
    ref <- summary_FUN_cols(tse, "group")
    #
    expect_equal(assay(res_sum), ref[["sum"]], check.attributes = FALSE)
    expect_equal(assay(res_sum_na), ref[["sum_na"]], check.attributes = FALSE)
    expect_equal(assay(res_mean), ref[["mean"]], check.attributes = FALSE)
    expect_equal(assay(res_mean_na), ref[["mean_na"]], check.attributes = FALSE)
    
    # Check that agglomerateByRank and agglomerateByVariable work correctly
    # with na.rm
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns
    col_idx <- sample(seq_len(ncol(tse)), 1)
    row_idx <- sample(seq_len(nrow(tse)), 1)
    tse_mod <- tse
    assay(tse_mod)[row_idx, col_idx] <- NA
    group <- colData(tse)[col_idx, "SampleType"]
    # na.rm = FALSE
    tse_sub <- agglomerateByVariable(tse_mod, by = "cols", group = "SampleType", na.rm = FALSE)
    test_mat <- tse_sub |> assay()
    # na.rm = TRUE
    tse_sub <- agglomerateByVariable(tse_mod, by = "cols", group = "SampleType", na.rm = TRUE)
    ref_mat <- tse_sub |> assay()
    group_idx <- which(colnames(tse_sub) == group)
    expect_true( is.na(test_mat[row_idx, group_idx]) )
    expect_true( !is.na(ref_mat[row_idx, group_idx]) )
    expect_equal(test_mat[-row_idx, -group_idx], ref_mat[-row_idx, -group_idx])
    #
    # na.rm = FALSE
    group <- rowData(tse)[row_idx, "Kingdom"]
    tse_sub <- agglomerateByVariable(tse_mod, by = "rows", group = "Kingdom", na.rm = FALSE)
    test_mat <- tse_sub |> assay()
    # na.rm = TRUE
    tse_sub <- agglomerateByVariable(tse_mod, by = "rows", group = "Kingdom", na.rm = TRUE)
    ref_mat <- tse_sub |> assay()
    group_idx <- which(rownames(tse_sub) == group)
    expect_true( is.na(test_mat[group_idx, col_idx]) )
    expect_true( !is.na(ref_mat[group_idx, col_idx]) )
    expect_equal(test_mat[-group_idx, -col_idx], ref_mat[-group_idx, -col_idx])
    #
    # na.rm = FALSE
    tse_sub <- agglomerateByRank(tse_mod, rank = "Kingdom", na.rm = FALSE)
    test_mat2 <- tse_sub |> assay()
    # na.rm = TRUE
    tse_sub <- agglomerateByRank(tse_mod, rank = "Kingdom", na.rm = TRUE)
    ref_mat2 <- tse_sub |> assay()
    expect_equal(test_mat, test_mat2)
    expect_equal(ref_mat, ref_mat2)
    
    # Check multiple rowTrees
    data(esophagus, package="mia")
    # Add arbitrary groups
    rowData(esophagus)$group <- c(rep(c("A", "B", "C"), each = nrow(esophagus)/3),
                                  rep("A", nrow(esophagus)-round(nrow(esophagus)/3)*3) )
    rowData(esophagus)$group2 <- c(rep(c("A", "B", "C"), each = nrow(esophagus)/3),
                                   rep("A", nrow(esophagus)-round(nrow(esophagus)/3)*3) )
    rowData(GlobalPatterns)$group <- c(rep(c("C", "D", "E"), each = nrow(GlobalPatterns)/3),
                                       rep("C", nrow(GlobalPatterns)-round(nrow(GlobalPatterns)/3)*3) )
    # Merge
    tse <- mergeSEs(esophagus, GlobalPatterns, assay.type="counts")
    # Reorder data since mergeSEs does not order the data based on original order
    # (trees are pruned differently --> first instance represent specific branch)
    tse <- tse[c(rownames(esophagus), rownames(GlobalPatterns)), ]
    # Only esophagus has these groups --> the merge should contain only esophagus
    merged  <- agglomerateByVariable(
        tse, by = "rows", group = rowData(tse)$group2, update.tree=TRUE)
    merged2 <- agglomerateByVariable(
        tse, by = "rows", group = rowData(tse)$group2, update.tree = FALSE)
    merged3 <- agglomerateByVariable(
        esophagus, by = "rows", group = rowData(esophagus)$group2, update.tree = TRUE)
    merged4 <- .merge_features(tse, merge.by = rowData(tse)$group2, update.tree = TRUE)
    merged5 <- agglomerateByVariable(
        tse, by = "rows", group = rowData(tse)$group2, update.tree = TRUE)
    expect_equal( rowLinks(merged)$whichTree,
                  rowLinks(merged2)$whichTree )
    expect_false( all(rowLinks(merged) == rowLinks(merged2)) )
    expect_equal(rowTree(tse), rowTree(merged2))
    expect_equal(merged4, merged5)
    expect_equal(
        agglomerateByVariable(tse, by = "rows", group = rowData(tse)$group2),
        agglomerateByVariable(tse, by = "rows", group = rowData(tse)$group2))

    # Both datasets have group variable
    merged <- agglomerateByVariable(
        tse, by = "rows", group = rowData(tse)$group, update.tree = TRUE)
    merged2 <- agglomerateByVariable(
        tse, by = "rows", group = rowData(tse)$group, update.tree = FALSE)
    expect_equal( rowLinks(merged)$whichTree,
                  rowLinks(merged2)$whichTree )
    expect_false( all(rowLinks(merged) == rowLinks(merged2)) )
    expect_true( rowTree(merged, "phylo")$Nnode < rowTree(merged2, "phylo")$Nnode )
})
