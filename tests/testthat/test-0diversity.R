context("diversity estimates")
test_that("diversity estimates", {

    skip_if_not(requireNamespace("vegan", quietly = TRUE))
    data(esophagus, package="mia")

    tse <- esophagus
    tse <- transformAssay(tse, method="relabundance")
    indices <- c("coverage", "fisher", "gini_simpson", "faith",
                 "inverse_simpson", "log_modulo_skewness",
                 "shannon")
    tse_idx <- addAlpha(tse, index = indices, threshold = 0.473)

    # Checks that the type of output is the same as the type of input.
    expect_true(typeof(tse_idx) == typeof(tse))

    # Check that every index is calculated by checking the column names from
    # colData.
    # Check that the order of indices is right / the same as the order
    # in the input vector.
    expect_named(colData(tse_idx), indices)

    lambda <- unname(colSums(assays(tse_idx)$relabundance^2))
    ginisimpson <- 1 - lambda
    invsimpson <- 1 / lambda

    expect_equal(lambda, unname(.simpson_lambda(assays(tse_idx)$relabundance)))
    expect_equal(ginisimpson, unname(colData(tse_idx)$gini_simpson))
    expect_equal(ginisimpson, unname(.calc_gini_simpson(assays(tse_idx)$relabundance)))

    expect_equal(invsimpson, unname(colData(tse_idx)$inverse_simpson))
    expect_equal(invsimpson, unname(.calc_inverse_simpson(assays(tse_idx)$relabundance)))

    cd <- colData(tse_idx)
    expect_equal(unname(round(cd$shannon, 5)), c(2.24937, 2.76239, 2.03249))
    expect_equal(unname(round(cd$gini_simpson, 6)),
                 c(0.831372, 0.903345, 0.665749))
    expect_equal(unname(round(cd$inverse_simpson, 5)),
                 c(5.93021, 10.34606, 2.99177))
    expect_equal(unname(round(cd$coverage, 0)), c(2,3,1))
    expect_equal(unname(round(cd$fisher, 4)), c(8.8037, 10.0989, 13.2783))
    expect_equal(unname(round(cd$log_modulo_skewness, 6)), c(2.013610, 1.827198, 2.013695))

    # Tests that 'quantile' and 'nclasses' are working
    expect_equal(unname(round(colData(
        addAlpha(
            tse, index="log_modulo_skewness", quantile=0.855,nclasses=32)
        )$log_modulo_skewness, 6)), c(1.814770, 1.756495, 1.842704))

    # Tests that .calc_skewness returns right value
    mat <- assay(tse, "counts")
    nclasses <- 61
    quantile <- 0.35

    quantile_point <- quantile(max(mat), quantile)
    cutpoints <- c(seq(0, quantile_point, length=nclasses), Inf)

    freq_table <- table(cut(mat, cutpoints), col(mat))
    test1 <- mia:::.calc_skewness(freq_table)

    test2 <- mia:::.calc_skewness(apply(mat, 2, function(x) {
        table(cut(x, cutpoints))
        }))

    expect_equal(unname(test1), unname(test2))
    expect_equal(unname(round(test1, 6)), c(7.256706, 6.098354, 7.278894))

    # Tests faith index with esophagus data
    for( i in c(1:(length(colnames(tse_idx)))) ){
        # Gets those taxa that are present/absent in the sample
        present <- rownames(tse)[assays(tse)$counts[,i] > 0]
        absent <- rownames(tse)[assays(tse)$counts[,i] == 0]

        # Absent taxa are dropped from the tree
        sub_tree <- ape::drop.tip(rowTree(tse_idx), absent)
        # Faith is now calculated based on the sub tree
        faith <- sum(sub_tree$edge.length)

        expect_equal(cd$faith[i], faith)
    }

    ########## Check that .estimate_faith works correctly ##########
    ########## with different SE object types ##########

    # Creates SE from TSE by dropping, e.g., rowTree
    se <- as(tse, "SummarizedExperiment")

    # Add rownames because they are not included when SE is created from TSE
    rownames(se) <- rownames(tse)

    # Calculates "faith" TSE
    tse_only <- addAlpha(tse, index = "faith")

    # tse_only should be TSE object
    expect_true(class(tse_only)== "TreeSummarizedExperiment")
    # tse_only should include "faith"
    expect_equal(colnames(colData(tse_only)), c(colnames(colData(tse)), "faith"))

    # Calculates "faith" TSE + TREE
    tse_tree <- addAlpha(tse, index = "faith", tree = rowTree(tse))

    # tse_tree should be TSE object
    expect_true(class(tse_tree)== "TreeSummarizedExperiment")
    # tse_tree should include "faith"
    expect_equal(
        colnames(colData(tse_tree)), c(colnames(colData(tse)), "faith"))


    # Calculates "faith" SE + TREE
    se_tree <- addAlpha(se, index = "faith", tree = rowTree(tse))

    # se_tree should be SE object
    expect_true(class(se_tree)== "SummarizedExperiment")
    # se_tree should include "faith"
    expect_equal(colnames(colData(se_tree)), c(colnames(colData(se)), "faith"))

    # Expect error
    expect_error(addAlpha(tse, index = "faith", tree.name = "test"))
    expect_error(addAlpha(tse, index = c("shannon", "faith"), tree.name = "test"))

    data(GlobalPatterns, package="mia")
    data(esophagus, package="mia")
    tse <- mergeSEs(GlobalPatterns, esophagus,  join = "full", assay.type = "counts")
    expect_error(addAlpha(tse, index = c("shannon", "faith"), tree.name = "phylo.1", assay.type="counts"))
    expect_error(addAlpha(tse, index = c("faith"), tree.name = "test"))
    expect_error(addAlpha(tse, index = c("shannon", "faith"), tree.name = TRUE))
    expect_error(addAlpha(tse, index = c("shannon", "faith"), tree.name = 1))
    expect_error(addAlpha(tse, index = c("shannon", "faith"), tree.name = c("phylo", "phylo.1")))
    expect_no_error(addAlpha(tse, index = c("shannon", "faith"), tree.name = "phylo"))

    # Test Faith with picante packages results (version 1.8.2)
    picante_res <- c(
        250.5354, 262.2629, 208.4578, 117.8762, 119.8247, 135.7673, 159.3715,
        123.3516, 143.7972, 111.7095, 156.4513, 147.9323, 247.2830, 253.2101,
        245.1008, 127.2336, 167.7246, 155.5872, 142.3473, 197.6823, 197.2321,
        124.6510, 121.2056, 179.9377, 140.8096, 126.5695)
    tse <- GlobalPatterns
    res <- .estimate_faith(tse, assay(tse), rowTree(tse))
    expect_equal(res, picante_res, tolerance=1e-5)
    # Check only tips paramater
    expect_error(.estimate_faith(tse, assay(tse), rowTree(tse), only.tips = 1))
    expect_error(.estimate_faith(
        tse, assay(tse), rowTree(tse), only.tips = "TRUE"))
    expect_error(.estimate_faith(
        tse, assay(tse), rowTree(tse), only.tips = c(TRUE, FALSE)))
})
