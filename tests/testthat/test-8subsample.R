context("rarefyAssay")
test_that("rarefyAssay", {
    seed = 1938
    set.seed(seed)
    data(GlobalPatterns, package="mia")
    
    tse.subsampled <- rarefyAssay(
        GlobalPatterns, 
        sample = 60000,
        name = "subsampled",
        replace = TRUE)
    # check class 
    expect_s4_class(tse.subsampled, "TreeSummarizedExperiment")
    expect_equal(nrow(tse.subsampled), 12403)
    expect_equal(ncol(tse.subsampled), 25)
    # check number of features removed is correct
    expnFeaturesRemoved <- 6813
    obsnFeaturesRemoved <- nrow(GlobalPatterns) - nrow(tse.subsampled)
    expect_equal(obsnFeaturesRemoved, expnFeaturesRemoved)
    
    # check if same Features removed
    obsFeaturesRemoved <- rownames(GlobalPatterns)[!rownames(GlobalPatterns) %in% rownames(tse.subsampled)]
    
    expFeaturesRemoved <- c("951","244423","586076","246140","143239",
                            "244960","144887","141782","215972","31759")
    
    expect_equal(obsFeaturesRemoved[1:10], expFeaturesRemoved)
    
    # check which sample is removed
    expSampleRemoved <- "TRRsed1"
    obsSampleRemoved <- colnames(GlobalPatterns)[!colnames(GlobalPatterns) %in% colnames(tse.subsampled)]
    expect_equal(obsSampleRemoved, expSampleRemoved)
    
    # check if all samples subsampled to even depth
    expColSums <- rep(60000, 25)
    expect_equal(unname(colSums2(assay(tse.subsampled, "subsampled"))), expColSums)
    
    # When replace = FALSE
    sample <- 60000
    seed = 1938
    set.seed(seed)
    tse.subsampled.rp <- rarefyAssay(
        GlobalPatterns, 
        sample = sample, 
        name = "subsampled",
        replace = FALSE)
    
    # check number of features removed is correct
    set.seed(seed)
    mat <- assay(GlobalPatterns, "counts")
    # When replace = FALSE, the function utilize vegan::rrarefy
    suppressWarnings(
    sub_mat <- t( vegan::rrarefy(t(mat), sample = sample) )
    )
    # The function removes those samples whose library size is lower than
    # 'sample'
    sub_mat <- sub_mat[, colSums(mat) > sample]
    # Moreover, it removes those features that are not present in any of the
    # samples
    sub_mat <- sub_mat[rowSums(sub_mat) > 0, ]
    expect_equal(assay(tse.subsampled.rp, "subsampled"), sub_mat, check.attributes = FALSE)
    
    # check if all samples subsampled to even depth
    expColSumsRp <- rep(sample, 25)
    expect_equal(unname(colSums2(assay(tse.subsampled.rp, "subsampled"))), expColSumsRp)
    
    # check if same Features removed
    expect_true( all(rownames(tse.subsampled.rp) == rownames(sub_mat)) )
})
