context("prevalence")

test_that("getPrevalence", {

    data(GlobalPatterns, package="mia")
    expect_error(getPrevalence(GlobalPatterns, detection="test"),
                 "'detection' must be a single numeric value or coercible to one")
    expect_error(getPrevalence(GlobalPatterns, include.lowest="test"),
                 "'include.lowest' must be TRUE or FALSE")
    expect_error(getPrevalence(GlobalPatterns, sort="test"),
                 "'sort' must be TRUE or FALSE")
    expect_error(getPrevalence(GlobalPatterns, as.relative="test"),
                 "'as.relative' must be TRUE or FALSE")
    expect_error(getPrevalence(GlobalPatterns, assay.type="test"),
                 "'assay.type' must be a valid name")
    # Output should be always a frequency between 0 to 1
    pr <- getPrevalence(GlobalPatterns, detection=0.1/100, as.relative=TRUE)
    expect_true(min(pr) >= 0 && max(pr) <= 1)
    pr <- getPrevalence(GlobalPatterns, detection=0.1/100, as.relative=FALSE)
    expect_true(min(pr) >= 0 && max(pr) <= 1)

    # Same prevalences should be returned for as.relative T/F in certain cases.
    pr1 <- getPrevalence(GlobalPatterns, detection=1, include.lowest=TRUE, as.relative=FALSE)
    pr2 <- getPrevalence(GlobalPatterns, detection=0/100, include.lowest=FALSE, as.relative=TRUE)
    expect_true(all(pr1 == pr2))

    # Same prevalences should be returned for as.relative T/F in certain cases.
    pr1 <- getPrevalence(GlobalPatterns, detection=1, include.lowest=TRUE, as.relative=FALSE)
    pr2 <- getPrevalence(GlobalPatterns, detection=0, include.lowest=FALSE, as.relative=FALSE)
    expect_true(all(pr1 == pr2))

    # Different ways to use relative abundance should yield the same output
    pr2 <- getPrevalence(GlobalPatterns, as.relative=TRUE, assay.type = "counts")
    GlobalPatterns <- transformAssay(GlobalPatterns, method="relabundance")
    pr1 <- getPrevalence(GlobalPatterns, as.relative=FALSE, assay.type = "relabundance")
    expect_true(all(pr1 == pr2))

    # Sorting should put the top values first
    pr <- getPrevalence(GlobalPatterns, sort=TRUE, detection = 0.1/100)
    expect_equal(as.vector(which.max(pr)), 1)
    pr <- names(head(getPrevalence(GlobalPatterns, sort=TRUE,  include.lowest = TRUE), 5L))
    actual <- getTop(GlobalPatterns,
                         method="prevalence",
                         top=5,
                         assay.type="counts")
    expect_equal(pr, actual)
    # Test alias
    alias<- getTop(GlobalPatterns,
                       method="prevalence",
                       top=5,
                       assay.type="counts")
    expect_equal(alias, actual)
    # Check that works also when rownames is NULL
    gp_null <- GlobalPatterns
    rownames(gp_null) <- NULL

    pr1 <- unname(getPrevalence(GlobalPatterns, detection=0.004, as.relative=TRUE))
    # If there are no rownames, getPrevalence returns names "featurex" where
    # x denotes the index.
    pr2 <- unname(getPrevalence(gp_null, detection=0.004, as.relative=TRUE))
    expect_equal(pr1, pr2)

    pr1 <- getPrevalence(GlobalPatterns, detection=0.004, as.relative=TRUE, rank = "Family")
    pr2 <- getPrevalence(gp_null, detection=0.004, as.relative=TRUE, rank = "Family")
    expect_equal(pr1, pr2)

    # Check that na.rm works correctly
    tse <- GlobalPatterns
    # Get reference value
    ref <- getPrevalence(tse, assay.type = "counts")
    # Add NA values to matrix
    remove <- c(1, 3, 10)
    assay(tse, "counts")[remove, ] <- NA
    # There should be 3 NA values if na.rm = FALSE. Otherwise there should be 0
    expect_warning(
        res <- getPrevalence(tse, assay.type = "counts", na.rm = FALSE) )
    expect_true( sum(is.na(res)) == 3)
    expect_warning(
        res <- getPrevalence(tse, assay.type = "counts", na.rm = TRUE) )
    expect_true( sum(is.na(res)) == 0)
    # Expect that other than features with NA values are the same as in reference
    expect_warning(
        res <- getPrevalence(tse, assay.type = "counts", na.rm = TRUE))
    res <- res[ !names(res) %in% remove ]
    ref <- ref[ !names(ref) %in% remove ]
    expect_equal( res[ names(ref) ], res[ names(ref) ] )

    # Now test that the number of samples where feature was detected is correct
    tse <- GlobalPatterns
    ref <- getPrevalence(tse, assay.type = "counts")
    # Add NA values to specific feature that has non-zero value
    feature <- rownames(tse)[[7]]
    assay(tse, "counts")[feature, 1] <- NA
    expect_warning(
        res <- getPrevalence(tse, assay.type = "counts", na.rm = TRUE))
    # Get the feature values and check that they have correct number of samples
    res <- res[ feature ]
    ref <- ref[ feature ]
    expect_true(res*ncol(tse) == 2)
    expect_true(ref*ncol(tse) == 3)

    #
    tse <- GlobalPatterns
    rank <- "Genus"
    # Add NA values to matrix
    remove <- c(15, 200)
    assay(tse, "counts")[remove, ] <- NA
    # Check that agglomeration works
    tse_agg <- agglomerateByRank(tse, ignore.taxonomy = FALSE, empty.rm = FALSE, rank = rank)
    expect_warning(ref <- getPrevalence(tse_agg, na.rm = FALSE))
    expect_warning(res <- getPrevalence(tse, rank = "Genus", empty.rm = FALSE))
    expect_true( all(res == ref, na.rm = TRUE) )
    #
    tse_agg <- agglomerateByRank(
        tse, ignore.taxonomy = FALSE, empty.rm = TRUE, rank = rank)
    ref <- getPrevalence(tse_agg, na.rm = TRUE)
    res <- getPrevalence(
        tse, na.rm = TRUE, rank = "Genus", empty.rm = TRUE)
    expect_true( all(res == ref, na.rm = TRUE) )

    # Test that results of addPrevalence equals to getPrevalence
    # sort should be disabled
    tse <- GlobalPatterns
    tse2 <- addPrevalence(tse, sort = TRUE)
    prev <- getPrevalence(tse, sort = FALSE)
    expect_equal(rowData(tse2)[["prevalence"]], prev)
    tse2 <- addPrevalence(tse, rank = "Family")
    prev <- getPrevalence(tse, rank = "Family")
    expect_equal(rowData(tse2)[["prevalence"]], prev)
    expect_error(addPrevalence(tse, name = 1))
    expect_error(addPrevalence(tse, name = NULL))
    expect_error(addPrevalence(tse, name = c("asd", "test")))
})


test_that("getPrevalent", {

    data(GlobalPatterns, package="mia")
    expect_error(getPrevalent(GlobalPatterns, prevalence="test"),
                 "'prevalence' must be a single numeric value or coercible to one")
    # Results compatible with getPrevalence
    pr1 <- getPrevalent(GlobalPatterns, detection=0.1/100, as.relative=TRUE, sort=TRUE)
    pr2 <- names(getPrevalence(GlobalPatterns, rank = "Kingdom", detection=0.1/100, as.relative=TRUE, sort=TRUE))
    expect_true(all(pr1 == pr2))

    # Same sorting for toptaxa obtained in different ways
    pr1 <- getPrevalent(GlobalPatterns, detection=0.1/100, as.relative=TRUE, sort=TRUE)
    pr2 <- names(getPrevalence(GlobalPatterns, rank = "Kingdom", detection=0.1/100, as.relative=TRUE, sort=TRUE))
    expect_true(all(pr1 == pr2))

    # Retrieved taxa are the same for counts and relative abundances
    pr1 <- getPrevalent(GlobalPatterns, prevalence=0.1/100, as.relative=TRUE)
    pr2 <- getPrevalent(GlobalPatterns, prevalence=0.1/100, as.relative=FALSE)
    expect_true(all(pr1 == pr2))

    # Prevalence and detection threshold at 0 has the same impact on counts and relative abundances
    pr1 <- getPrevalent(GlobalPatterns, detection=0, prevalence=0, as.relative=TRUE)
    pr2 <- getPrevalent(GlobalPatterns, detection=0, prevalence=0, as.relative=FALSE)
    expect_true(all(pr1 == pr2))

    # Check that works also when rownames is NULL
    gp_null <- GlobalPatterns
    rownames(gp_null) <- NULL

    pr1 <- getPrevalent(GlobalPatterns, detection=0.0045, prevalence = 0.25, as.relative=TRUE)
    pr1 <- which(rownames(GlobalPatterns) %in% pr1)
    pr2 <- getPrevalent(gp_null, detection=0.0045, prevalence = 0.25, as.relative=TRUE)
    expect_equal(pr1, pr2)

    # Test alias
    alias <- getPrevalent(gp_null, detection=0.0045, prevalence = 0.25, as.relative=TRUE)
    expect_equal(pr1, alias)

    pr1 <- getPrevalent(GlobalPatterns, detection=0.004, prevalence = 0.1,
                            as.relative=TRUE, rank = "Family")
    pr2 <- getPrevalent(gp_null, detection=0.004, prevalence = 0.1,
                            as.relative=TRUE, rank = "Family")
    expect_equal(pr1, pr2)

})

test_that("getRare", {

    data(GlobalPatterns, package="mia")
    expect_error(getRare(GlobalPatterns, prevalence="test"),
                 "'prevalence' must be a single numeric value or coercible to one")

    ############# Test that output type is correct #############
    expect_type(taxa <- getRare(GlobalPatterns,
                                    detection = 130,
                                    prevalence = 90/100), "character")

    ##### Test that getPrevalent and getRare has all the taxa ####

    # Gets rownames for all the taxa
    all_taxa <- rownames(GlobalPatterns)

    # Gets prevalent taxa
    prevalent_taxa <- getPrevalent(GlobalPatterns,
                                       detection = 0,
                                       prevalence = 90/100,
                                       include.lowest = FALSE)
    # Gets rare taxa
    rare_taxa <- getRare(GlobalPatterns,
                             detection = 0,
                             prevalence = 90/100,
                             include.lowest = FALSE)

    # Concatenates prevalent and rare taxa
    prevalent_and_rare_taxa <- c(prevalent_taxa, rare_taxa)

    # If all the elements are in another vector, vector of TRUEs
    # Opposite --> negative of FALSEs
    # If every element is FALSE, any is FALSE --> opposite --> TRUE
    expect_true( !any( !(all_taxa %in% prevalent_and_rare_taxa) ) &&
                     !any( !( prevalent_and_rare_taxa %in% all_taxa) ) )

    ##### Test that getPrevalent and getRare has all the taxa, ####
    ##### but now with detection limit and with rank #####

    # Check that it works with all the ranks
    ranks <- taxonomyRanks(GlobalPatterns)

    for( rank in ranks ){

        # Agglomerates data by rank
        se <- agglomerateByRank(GlobalPatterns, rank = rank)

        # Gets rownames for all the taxa
        all_taxa <- rownames(se)

        # Gets prevalent taxa
        prevalent_taxa <- getPrevalent(GlobalPatterns,
                                           prevalence = 0.05,
                                           detection = 0.1,
                                           rank = rank,
                                           include.lowest = TRUE, as.relative = TRUE)
        # Gets rare taxa
        rare_taxa <- getRare(GlobalPatterns,
                                 prevalence = 0.05,
                                 detection = 0.1,
                                 rank = rank,
                                 include.lowest = TRUE, as.relative = TRUE)

        # Concatenates prevalent and rare taxa
        prevalent_and_rare_taxa <- c(prevalent_taxa, rare_taxa)

        # If all the elements are in another vector, vector of TRUEs
        # Opposite --> negative of FALSEs
        # If every element is FALSE, any is FALSE --> opposite --> TRUE
        expect_true( !any( !(all_taxa %in% prevalent_and_rare_taxa) ) &&
                         !any( !( prevalent_and_rare_taxa %in% all_taxa) ) )

        # Expect that there are no duplicates
        expect_true(!anyDuplicated(prevalent_taxa))
        expect_true(!anyDuplicated(rare_taxa))

    }

    # Check that works also when rownames is NULL
    gp_null <- GlobalPatterns
    rownames(gp_null) <- NULL

    pr1 <- getRare(GlobalPatterns, detection=0.0045, prevalence = 0.25, as.relative=TRUE)
    pr1 <- which(rownames(GlobalPatterns) %in% pr1)
    pr2 <- getRare(gp_null, detection=0.0045, prevalence = 0.25, as.relative=TRUE)
    expect_equal(pr1, pr2)

    # Test lias
    alias <- getRare(gp_null, detection=0.0045, prevalence = 0.25, as.relative=TRUE)
    expect_equal(pr1, alias)

    pr1 <- getRare(GlobalPatterns, detection=0.004, prevalence = 0.1,
                       as.relative=TRUE, rank = "Family")
    pr2 <- getRare(gp_null, detection=0.004, prevalence = 0.1,
                       as.relative=TRUE, rank = "Family")
    expect_equal(pr1, pr2)

})

test_that("subsetByPrevalent", {
    data(GlobalPatterns, package="mia")
    expect_error(subsetByPrevalent(GlobalPatterns, prevalence="test"),
                 "'prevalence' must be a single numeric value or coercible to one")
    # Expect TSE object
    expect_equal(class(subsetByPrevalent(GlobalPatterns)), class(GlobalPatterns))

    # Results compatible with getPrevalent
    pr1 <- rownames(subsetByPrevalent(
        GlobalPatterns, rank = "Class", detection=0.1/100,
        as.relative=TRUE, sort=TRUE))
    pr2 <- getPrevalent(GlobalPatterns, rank = "Class", detection=0.1/100,
                            as.relative=TRUE, sort=TRUE)
    expect_true(all(pr1 == pr2))

    # Retrieved taxa are the same for counts and relative abundances
    pr1 <- assay(subsetByPrevalent(GlobalPatterns, prevalence=0.1/100, as.relative=TRUE), "counts")
    pr2 <- assay(subsetByPrevalent(GlobalPatterns, prevalence=0.1/100, as.relative=FALSE), "counts")
    expect_true(all(pr1 == pr2))

    # Prevalence and detection threshold at 0 has the same impact on counts and relative abundances
    pr1 <- rownames(subsetByPrevalent(GlobalPatterns, detection=0, prevalence=0, as.relative=TRUE))
    pr2 <- rownames(subsetByPrevalent(GlobalPatterns, detection=0, prevalence=0, as.relative=FALSE))
    expect_true(all(pr1 == pr2))

    # Check that works also when rownames is NULL
    gp_null <- GlobalPatterns
    rownames(gp_null) <- NULL

    pr1 <- subsetByPrevalent(GlobalPatterns, detection=12, prevalence = 0.33)
    pr1 <- unname(assay(pr1, "counts"))
    pr2 <- subsetByPrevalent(gp_null, detection=12, prevalence = 0.33)
    pr2 <- unname(assay(pr2, "counts"))
    expect_equal(pr1, pr2)

    pr1 <- subsetByPrevalent(GlobalPatterns, detection=5, prevalence = 0.33, rank = "Phylum")
    pr1 <- unname(assay(pr1, "counts"))
    pr2 <- subsetByPrevalent(gp_null, detection=5, prevalence = 0.33, rank = "Phylum")
    pr2 <- unname(assay(pr2, "counts"))
    expect_equal(pr1, pr2)

    # Test alias
    alias <- subsetByPrevalent(gp_null, detection=5, prevalence = 0.33, rank = "Phylum")
    alias <- unname(assay(alias, "counts"))
    expect_equal(alias, pr2)

    # Check that tree subsetting works
    expect_error(subsetByPrevalent(GlobalPatterns, update.tree = 1))
    expect_error(subsetByPrevalent(GlobalPatterns, update.tree = NULL))
    expect_error(subsetByPrevalent(GlobalPatterns, update.tree = c(TRUE, FALSE)))
    tse_sub <- subsetByPrevalent(GlobalPatterns, prevalence = 0.4, rank = "Genus", update.tree = TRUE)
    expect_equal(length(rowTree(tse_sub)$tip.label), nrow(tse_sub))

})

test_that("subsetByRare", {
    data(GlobalPatterns, package="mia")
    expect_error(subsetByRare(GlobalPatterns, prevalence="test"),
                 "'prevalence' must be a single numeric value or coercible to one")
    # Expect TSE object
    expect_equal(class(subsetByRare(GlobalPatterns)), class(GlobalPatterns))

    # Results compatible with getRare
    pr1 <- rownames(subsetByRare(
        GlobalPatterns, rank = "Phylum", detection=0.1/100,
        as.relative=TRUE, sort=TRUE))
    pr2 <- getRare(GlobalPatterns, rank = "Phylum", detection=0.1/100,
                       as.relative=TRUE, sort=TRUE)
    expect_true(all(pr1 == pr2))

    # Retrieved taxa are the same for counts and relative abundances
    pr1 <- assay(subsetByRare(GlobalPatterns, prevalence=0.1/100, as.relative=TRUE), "counts")
    pr2 <- assay(subsetByRare(GlobalPatterns, prevalence=0.1/100, as.relative=FALSE), "counts")
    expect_true(all(pr1 == pr2))

    # Prevalence and detection threshold at 0 has the same impact on counts and relative abundances
    pr1 <- rownames(subsetByRare(GlobalPatterns, detection=0, prevalence=0, as.relative=TRUE))
    pr2 <- rownames(subsetByRare(GlobalPatterns, detection=0, prevalence=0, as.relative=FALSE))
    expect_true(all(pr1 == pr2))

    # subsetByRare + subsetByPrevalent should include all the taxa in OTU level
    d <- runif(1, 0.0001, 0.1)
    p <- runif(1, 0.0001, 0.5)
    rare <- rownames(subsetByRare(GlobalPatterns, detection=d, prevalence=p,
                                      as.relative=TRUE))

    prevalent <- rownames(subsetByPrevalent(GlobalPatterns, detection=d, prevalence=p,
                                                as.relative=TRUE))

    all_taxa <- c(rare, prevalent)

    expect_true( all(all_taxa %in% rownames(GlobalPatterns)) )

    # Check that works also when rownames is NULL
    gp_null <- GlobalPatterns
    rownames(gp_null) <- NULL

    pr1 <- subsetByRare(GlobalPatterns, detection=12, prevalence = 0.33)
    pr1 <- unname(assay(pr1, "counts"))
    pr2 <- subsetByRare(gp_null, detection=12, prevalence = 0.33)
    pr2 <- unname(assay(pr2, "counts"))
    expect_equal(pr1, pr2)

    pr1 <- subsetByRare(GlobalPatterns, detection=5, prevalence = 0.33, rank = "Phylum")
    pr1 <- unname(assay(pr1, "counts"))
    pr2 <- subsetByRare(gp_null, detection=5, prevalence = 0.33, rank = "Phylum")
    pr2 <- unname(assay(pr2, "counts"))
    expect_equal(pr1, pr2)

    # Test alias
    alias <- subsetByRare(gp_null, detection=5, prevalence = 0.33, rank = "Phylum")
    alias <- unname(assay(alias, "counts"))
    expect_equal(alias, pr2)

    # Check that tree subsetting works
    expect_error(subsetByRare(GlobalPatterns, update.tree = 1))
    expect_error(subsetByRare(GlobalPatterns, update.tree = NULL))
    expect_error(subsetByRare(GlobalPatterns, update.tree = c(TRUE, FALSE)))
    tse_sub <- subsetByRare(GlobalPatterns, rank = "Class", update.tree = TRUE)
    expect_equal(length(rowTree(tse_sub)$tip.label), nrow(tse_sub))
})

test_that("agglomerateByPrevalence", {

    data(GlobalPatterns, package="mia")
    expect_error(agglomerateByPrevalence(GlobalPatterns, other.name=TRUE),
                 "'other.name' must be a single character value")
    actual <- agglomerateByPrevalence(GlobalPatterns, rank = "Kingdom")
    expect_s4_class(actual,class(GlobalPatterns))
    expect_equal(dim(actual),c(2,26))

    actual <- agglomerateByPrevalence(
        GlobalPatterns,
        rank = "Phylum",
        detection = 1/100,
        prevalence = 50/100,
        as.relative = TRUE,
        other.name = "test",
        update.tree = FALSE)
    expect_s4_class(actual,class(GlobalPatterns))
    expect_equal(dim(actual),c(6,26))
    expect_equal(rowData(actual)[6,"Phylum"],"test")
    expect_false(length(rowTree(actual)$tip.label) == length(rownames(actual)))

    actual <- agglomerateByPrevalence(GlobalPatterns,
                                      rank = NULL,
                                      detection = 0.0001,
                                      prevalence = 50/100,
                                      as.relative = TRUE,
                                      other.name = "test",
                                      update.tree = TRUE)
    expect_equal(agglomerateByPrevalence(GlobalPatterns,
                                           rank = NULL,
                                           detection = 0.0001,
                                           prevalence = 50/100,
                                           as.relative = TRUE,
                                           other.name = "test"),
                 agglomerateByPrevalence(GlobalPatterns,
                                           rank = NULL,
                                           detection = 0.0001,
                                           prevalence = 50/100,
                                           as.relative = TRUE,
                                           other.name = "test"))
    expect_s4_class(actual,class(GlobalPatterns))
    expect_equal(dim(actual),c(6,26))
    expect_true(all(is.na(rowData(actual)[6,])))
    expect_equal(length(rowTree(actual)$tip.label), length(rownames(actual)))
    actual <- agglomerateByPrevalence(GlobalPatterns,
                                      rank = "Phylum",
                                      detection = 1/100,
                                      prevalence = 50/100,
                                      as_relative = TRUE,
                                      other_label = "test",
                                      update.tree = TRUE)
    expect_equal(length(rowTree(actual)$tip.label), length(rownames(actual)))

    # Load data from miaTime package
    skip_if_not(require("miaTime", quietly = TRUE))
    data(SilvermanAGutData)
    se <- SilvermanAGutData

    # checking reference consensus sequence generation
    actual <- agglomerateByPrevalence(se,"Genus", update.refseq = FALSE)
    # There should be only one exact match for each sequence
    seqs_test <- as.character( referenceSeq(actual) )
    seqs_ref <- as.character( referenceSeq(se) )
    expect_true(all(vapply(
    seqs_test, function(seq) sum(seqs_ref %in% seq) == 1,
    FUN.VALUE = logical(1) )) )

    # Merging creates consensus sequences.
    th <- runif(1, 0, 1)
    actual <- agglomerateByPrevalence(
      se, "Genus", update.refseq = TRUE, threshold = th)
    seqs_test <- referenceSeq(actual)
    # Get single taxon as reference. Merge those sequences and test that it
    # equals to one that is output of agglomerateByPrevalence
    seqs_ref <- referenceSeq(se)
    feature <- sample(na.omit(rowData(actual)[["Genus"]]), 1)
    seqs_ref <- seqs_ref[ rowData(se)[["Genus"]] %in% feature ]
    seqs_ref <- .merge_refseq(
      seqs_ref, factor(rep(feature, length(seqs_ref))), rownames(seqs_ref),
      threshold = th)
    seqs_test <- seqs_test[ names(seqs_test) %in% feature ]
    expect_equal(seqs_test, seqs_ref)

    # checking reference consensus sequence generation using 'Genus:Alistipes'
    actual <- agglomerateByPrevalence(se,"Genus", update.refseq = FALSE)
    expect_equal(as.character(referenceSeq(actual)[["Alistipes"]]),
                 paste0("TCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTTTGATAA",
                        "GTTAGAGGTGAAATCCCGGGGCTTAACTCCGGAACTGCCTCTAATACTGTTAG",
                        "ACTAGAGAGTAGTTGCGGTAGGCGGAATGTATGGTGTAGCGGTGAAATGCTTA",
                        "GAGATCATACAGAACACCGATTGCGAAGGCAGCTTACCAAACTATATCTGACG",
                        "TTGAGGCACGAAAGCGTGGGG"))
    actual <- agglomerateByPrevalence(se,"Genus", update.refseq = TRUE)
    expect_equal(as.character(referenceSeq(actual)[["Alistipes"]]),
                 paste0("SCRAGCGTTRTCCGGAWTTAYTGGGYKTAAAGSGMGCGYAGGYGGHBDNKYAA",
                        "GTCWGWWGTGAAAKYYYGSGGCTCAACCSYRRRMBKSCWKTKGAAACTGBVHK",
                        "RCTWGAKTKYVKDWGAGGWRRGYGGAATKCSWVGTGTAGCGGTGAAATGCKTA",
                        "GAKATBWSGARGAACWCCRRTKGCGAAGGCRRCTYWCTRGWCKGWVAMTGACG",
                        "CTGAKGCKCGAAAGYGTGGGK"))
})
