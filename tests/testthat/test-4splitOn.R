context("splitOn")
test_that("splitOn", {
    data(GlobalPatterns, package="mia")
    x <- GlobalPatterns

    ################################## splitOn #################################
    # Test that throughs an error
    expect_error(splitOn(x, "test"))
    expect_error(splitOn(x[1:10, 1:10], x$SampleType))
    set.seed(103)
    rowData(x)$group <- sample(1:10, nrow(x), replace = TRUE)
    set.seed(374)
    colData(x)$group <- sample(1:10, ncol(x), replace = TRUE)
    expect_error(splitOn(x, "group"))
    expect_error(splitOn(x, "SampleType", by = 1))
    expect_error(splitOn(x, "Phylum", by = 2))
    expect_error(splitOn(x, x$SampleType, by = 1))
    expect_error(splitOn(x, rowData(x)$Phylum, by = 2))
    expect_error(splitOn(x))
    expect_error(splitOn(assay(x), x$SampleType))
    expect_error(splitOn(x, "SampleType", use.names = 1))
    expect_error(splitOn(x, "SampleType", use.names = "TRUE"))
    expect_error(splitOn(x, "SampleType", update.tree = 1))
    expect_error(splitOn(x, "SampleType", update.tree = "TRUE"))

    # Test that names of elemetns are correct
    list <- splitOn(x, "SampleType")
    expect_equal(names(list), as.character(unique(x$SampleType)) )
    list <- splitOn(x, "SampleType", use.names = FALSE)
    expect_equal(names(list), NULL )

    # Test that col-wie split is done correctly
    list <- splitOn(x, "group", by = 1)
    expect_equal( colnames(list[[1]]), colnames(x) )
    expect_true( length(list) == 10 )

    # Test that row-wise split is done correctly
    list <- splitOn(x, "group", by = "cols")
    expect_equal( rownames(list[[1]]), rownames(x) )
    expect_true( length(list) == 10 )

    # Test that number of tips of updated rowTree equals number of rows for
    # each tse in the list returned
    list <- splitOn(x, "SampleType", update.tree = TRUE)
    for (k in length(list)){
        expect_equal( length(rowTree(list[[k]], "phylo")$tip.label),
                      nrow(list[[k]]) )
    }

    ################################# unsplitOn ################################
    # Test that error occurs
    expect_error( unsplitOn(x) )
    mod_list <- list
    mod_list[[1]] <- mod_list[[1]][1:2, 1:2]
    expect_error( unsplitOn(mod_list) )
    expect_error(unsplitOn(list, update.tree = 1))
    expect_error(unsplitOn(list, update.tree = "TRUE"))

    # Test that works
    x_sub <- x[1:100, 1:10]
    # Split
    list <- splitOn(x_sub, "group", by = "features")
    # Unsplit
    unsplitted <- unsplitOn(list)
    # Order the data
    unsplitted <- unsplitted[ rownames(x_sub), colnames(x_sub) ]
    # Convert delayed matrix to normal
    assay(unsplitted) <- as.matrix( assay(unsplitted) )
    expect_equal(assay(x_sub), assay(unsplitted) )

    x_sub <- x[1:100, 1:10]
    # Split
    list <- splitOn(x_sub, "group", by = 2)
    # Unsplit
    unsplitted <- unsplitOn(list)
    # Order the data
    unsplitted <- unsplitted[ rownames(x_sub), colnames(x_sub) ]
    # Convert delayed matrix to normal
    assay(unsplitted) <- as.matrix( assay(unsplitted) )
    expect_equal(assay(x_sub), assay(unsplitted) )

    list <- splitOn(x, "SampleType")
    unsplitted <- unsplitOn(list)
    # Order the data
    unsplitted <- unsplitted[ rownames(x), colnames(x) ]
    # Convert delayed matrix to normal
    assay(unsplitted) <- as.matrix( assay(unsplitted) )
    expect_equal(assay(x), assay(unsplitted) )

    # Split
    list <- splitOn(x, "Phylum")
    # Unsplit
    unsplitted <- unsplitOn(list)
    # Order the data
    unsplitted <- unsplitted[ rownames(x), colnames(x) ]
    # Convert delayed matrix to normal
    assay(unsplitted) <- as.matrix( assay(unsplitted) )
    expect_equal(assay(x), assay(unsplitted) )
})
