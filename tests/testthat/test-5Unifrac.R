context("Unifrac beta diversity")
test_that("Unifrac beta diversity", {
    data(esophagus, package="mia")
    tse <- esophagus
    tse <- transformAssay(tse, assay.type="counts", method="relabundance")

    expect_error(
        getDissimilarity(tse, method = "unifrac",assay.type = "test",
                          tree.name = "phylo", weighted = FALSE)
    )
    expect_error(
        getDissimilarity(tse, method = "unifrac", assay.type = 2,
                          tree.name = "phylo", weighted = FALSE)
    )
    expect_error(
        getDissimilarity(tse, method = "unifrac", assay.type = TRUE,
                          tree.name = "phylo", weighted = FALSE)
    )
    expect_error(
      getDissimilarity(tse, method = "unifrac", assay.type = "counts",
                       tree.name = "test", weighted = FALSE)
    )
    expect_error(
      getDissimilarity(tse, method = "unifrac", assay.type = "counts",
                       tree.name = 1, weighted = FALSE)
    )
    expect_error(
      getDissimilarity(tse, method = "unifrac", assay.type = "counts",
                       tree.name = TRUE, weighted = FALSE)
    )
    expect_error(
      getDissimilarity(tse, method = "unifrac", assay.type = "counts",
                       tree.name = "phylo", weighted = 1)
    )

    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns
    # Calculate unweighted unifrac
    unifrac_mia <- as.matrix(getDissimilarity(tse, method = "unifrac",
                                              weighted = FALSE))
    unifrac_rbiom <- as.matrix(rbiom::unifrac(assay(tse), weighted = FALSE,
                                              rowTree(tse)))
    expect_equal(unifrac_mia, unifrac_rbiom)
    # Calculate weighted unifrac. Allow tolerance since weighted unifrac
    # calculation in rbiom has some stochasticity. That is most likely due
    # multithreading and complex structure of tree (loops).
    unifrac_mia <- as.matrix(getDissimilarity(tse, method = "unifrac",
                                              weighted = TRUE))
    unifrac_rbiom <- as.matrix(rbiom::unifrac(assay(tse), weighted = TRUE,
                                              rowTree(tse)))
    expect_equal(unifrac_mia, unifrac_rbiom, tolerance = 5e-3)

    # Test that the function works correctly when there are multiple trees.
    # The function should subset the data based on tree.
    tse <- GlobalPatterns
    trees <- list(phylo = rowTree(tse), tree2 = rowTree(tse))
    links <- rowLinks(tse)
    links[ 1:500 , "whichTree"] <- "tree2"
    tse@rowTree <- trees
    tse@rowLinks <- links
    tse_ref <- tse
    tse_ref <- tse_ref[ rowLinks(tse_ref)[["whichTree"]] == "tree2", ]
    # Calculate unweighted unifrac
    expect_warning(
    unifrac_mia <- getDissimilarity(tse, method = "unifrac",
                                    weighted = FALSE, tree.name = "tree2")
    )
    unifrac_mia <- as.matrix(unifrac_mia)
    unifrac_rbiom <- as.matrix(rbiom::unifrac(assay(tse_ref), weighted = FALSE,
                                              rowTree(tse_ref)))
    expect_equal(unifrac_mia, unifrac_rbiom)
    # Calculate weighted unifrac
    expect_warning(
    unifrac_mia <- getDissimilarity(tse, method = "unifrac",
                                    weighted = TRUE, tree.name = "tree2")
    )
    unifrac_mia <- as.matrix(unifrac_mia)
    unifrac_rbiom <- as.matrix(rbiom::unifrac(assay(tse_ref), weighted = TRUE,
                                              rowTree(tse_ref)))
    expect_equal(unifrac_mia, unifrac_rbiom, tolerance = 1e-3)

    # Test the function with agglomerated data. .get_unifrac renames
    # rownames based on tips and links to them. Then it also prunes the tree so
    # that rows are in tips.
    tse <- GlobalPatterns
    tse <- agglomerateByRank(tse, rank = "Species")
    tse_ref <- tse
    rownames(tse_ref) <- rowLinks(tse_ref)[["nodeLab"]]
    rowTree(tse_ref) <- .prune_tree(rowTree(tse_ref), rowLinks(tse_ref)[["nodeLab"]])
    # Calculate unweighted unifrac
    unifrac_mia <- as.matrix(getDissimilarity(tse, method = "unifrac",
                                              weighted = FALSE))
    unifrac_rbiom <- as.matrix(rbiom::unifrac(assay(tse_ref), weighted = FALSE,
                                              rowTree(tse_ref)))
    # Calculate weighted unifrac. No tolerance needed since the tree has
    # simpler structure after pruning.
    unifrac_mia <- as.matrix(getDissimilarity(tse, method = "unifrac",
                                              weighted = TRUE))
    unifrac_rbiom <- as.matrix(rbiom::unifrac(assay(tse_ref), weighted = TRUE,
                                              rowTree(tse_ref)))
    expect_equal(unifrac_mia, unifrac_rbiom)
})
