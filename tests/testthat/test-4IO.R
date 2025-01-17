
test_that("Importing biom files yield SummarizedExperiment objects", {
    skip_if_not(require("biomformat", quietly = TRUE))
    rich_dense_file  = system.file("extdata", "rich_dense_otu_table.biom",
                                    package = "biomformat")
    me <- importBIOM(rich_dense_file)
    expect_s4_class(me, "SummarizedExperiment")
    # load from object
    x1 <- biomformat::read_biom(rich_dense_file)
    me2 <- convertFromBIOM(x1)
    expect_s4_class(me2, "SummarizedExperiment")
    expect_equal(dim(me), dim(me2))
    expect_equal(rowData(me), rowData(me2))

    biom_object <- biomformat::read_biom(
        system.file("extdata", "Aggregated_humanization2.biom",
                    package="mia")
    )
    tse <- convertFromBIOM(biom_object,
                                prefix.rm = FALSE,
                                rank.from.prefix = FALSE,
                                artifact.rm = TRUE,
                                pattern = "\"")
    # Testing no prefixes removed
    expect_true(rowData(tse) %>%
                    apply(2,grepl,pattern="sk__|([dkpcofgs]+)__") %>%
                    all())
    # Testing no taxonomy ranks parsed
    expect_false(
        sapply(tolower(colnames(rowData(tse))),
                function(x) x %in% TAXONOMY_RANKS) %>% all())
    # Testing the artifact.rm, since the original artifact in the biom file
    # is '\"'
    expect_false(apply(rowData(tse), 2, grepl, pattern="^\"") %>% all())

    # Testing prefixes removed
    tse <- convertFromBIOM(biom_object,
                                prefix.rm=TRUE,
                                rank.from.prefix=TRUE,
                                artifact.rm = TRUE,
                                pattern = "\"")
    expect_false(rowData(tse) %>%
                    apply(2,grepl,pattern="sk__|([dkpcofgs]+)__") %>%
                    all())

    # Testing parsing taxonomy ranks from prefixes
    tse <- convertFromBIOM(biom_object,
                                prefix.rm=FALSE,
                                rank.from.prefix=TRUE,
                                artifact.rm = TRUE,
                                pattern = "\"")
    expect_true(
        sapply(tolower(colnames(rowData(tse))),
                function(x) x %in% TAXONOMY_RANKS) %>% all())

    # Testing the artifact.rm, the original artifact in the biom file
    # is '\"', as a test we rather try remove a non existing pattern.
    tse <- convertFromBIOM(biom_object,
                                prefix.rm=FALSE,
                                rank.from.prefix=FALSE,
                                artifact.rm = TRUE,
                                pattern = "\\*|\\?")
    # with wrong pattern artifact not cleaned
    expect_true(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())
    # Testing the artifact.rm, with the value 'auto' to automatically
    # detect the artifact and remove it (in our case the artifact is '\"').
    tse <- convertFromBIOM(biom_object,
                                prefix.rm=FALSE,
                                rank.from.prefix=FALSE,
                                artifact.rm = TRUE)
    # Checking if 'auto' has detected and cleaned the artifact
    expect_false(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())
    # Testing the artifact.rm, with the value NULL to not detect or clean
    # anything.
    tse <- convertFromBIOM(biom_object,
                                prefix.rm=FALSE,
                                rank.from.prefix=FALSE,
                                artifact.rm = FALSE)
    # Checking if the '\"' artifact still exists.
    expect_true(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())

    # General final test
    tse <- convertFromBIOM(biom_object,
                                prefix.rm=TRUE,
                                rank.from.prefix=TRUE,
                                artifact.rm = TRUE)
    # check if '\"' cleaned
    expect_false(apply(rowData(tse), 2, grepl, pattern="\"") %>% any())
    # check if taxa prefixes removed
    expect_false(rowData(tse) %>%
                    apply(2,grepl,pattern="sk__|([dkpcofgs]+)__") %>%
                    all())
    # Check if rank names were parsed correctly
    expect_true(
        sapply(tolower(colnames(rowData(tse))),
                function(x) x %in% TAXONOMY_RANKS) %>% all())

    # General final test with another biom file
    biom_object <- biomformat::read_biom(
        system.file("extdata", "rich_dense_otu_table.biom",
                    package = "biomformat")
    )
    tse <- convertFromBIOM(biom_object,
                                prefix.rm=TRUE,
                                rank.from.prefix=TRUE,
                                artifact.rm = TRUE)
    # check if taxa prefixes removed
    expect_false(rowData(tse) %>%
                    apply(2,grepl,pattern="sk__|([dkpcofgs]+)__") %>%
                    all())
    # Check if rank names were parsed correctly
    expect_true(
        sapply(tolower(colnames(rowData(tse))),
                function(x) x %in% TAXONOMY_RANKS) %>% all())

    # Check that convertToBIOM works
    # Get errors if input is incorrect
    expect_error( convertToBIOM() )
    expect_error( convertToBIOM(c(1, 2, 3)) )
    expect_error( convertToBIOM(tse, assay.type = "test") )
    # Check that the output is correct
    biom <- convertToBIOM(tse)
    # Test type
    expect_equal(biomformat::matrix_element_type(biom), "int")
    # Assay
    test <- as.matrix( biomformat::biom_data(biom) )
    ref <- assay(tse)
    expect_equal(test, ref)
    # rowData
    test <- as.data.frame( biomformat::observation_metadata(biom) )
    ref <- as.data.frame( rowData(tse) )
    expect_equal(test, ref)
    # colData
    test <- as.data.frame( biomformat::sample_metadata(biom) )
    ref <- as.data.frame( colData(tse) )
    expect_equal(test, ref)
    # Test with empty rowData and colData and relative abundance (float) assay
    rowData(tse) <- NULL
    colData(tse) <- NULL
    tse <- transformAssay(tse, method = "relabundance")
    biom <- convertToBIOM(tse, assay.type = "relabundance")
    # Test type
    expect_equal(biomformat::matrix_element_type(biom), "float")
    # Assay
    test <- as.matrix( biomformat::biom_data(biom) )
    ref <- assay(tse, "relabundance")
    expect_equivalent(test, ref)
    # rowData has one empty column (only NA values)
    test <- as.data.frame( biomformat::observation_metadata(biom) )
    expect_true( all(is.na(test)) && colnames(test) == "V1" )
    # colData has one empty column (only NA values)
    test <- as.data.frame( biomformat::sample_metadata(biom) )
    expect_true( all(is.na(test)) && colnames(test) == "V1" )

})

test_that("Importing phyloseq objects yield TreeSummarizedExperiment objects", {
    skip_if_not_installed("phyloseq")
    data(GlobalPatterns, package="phyloseq")
    me <- convertFromPhyloseq(GlobalPatterns)
    expect_s4_class(me, "TreeSummarizedExperiment")
    expect_equal(dim(me),c(19216,26))
    data(enterotype, package="phyloseq")
    me <- convertFromPhyloseq(enterotype)
    expect_s4_class(me, "TreeSummarizedExperiment")
    expect_equal(dim(me),c(553,280))
    data(esophagus, package="phyloseq")
    me <- convertFromPhyloseq(esophagus)
    expect_s4_class(me, "TreeSummarizedExperiment")
    expect_equal(dim(me),c(58,3))
    esophagus2 <- esophagus
    phyloseq::otu_table(esophagus2) <- t(phyloseq::otu_table(esophagus))
    me2 <- convertFromPhyloseq(esophagus2)
    expect_equal(me, me2)
})

test_that("Importing dada2 objects yield TreeSummarizedExperiment objects", {
    skip_if_not_installed("dada2")
    fnF <- system.file("extdata", "sam1F.fastq.gz", package="dada2")
    fnR = system.file("extdata", "sam1R.fastq.gz", package="dada2")
    dadaF <- dada2::dada(fnF, selfConsist=TRUE)
    dadaR <- dada2::dada(fnR, selfConsist=TRUE)

    me <- convertFromDADA2(dadaF, fnF, dadaR, fnR)
    expect_s4_class(me, "TreeSummarizedExperiment")
})

test_that("Importing Mothur files yield SummarizedExperiment objects", {

    counts <- system.file("extdata", "mothur_example.shared", package = "mia")
    taxa <- system.file("extdata", "mothur_example.cons.taxonomy", package = "mia")
    taxa2 <- system.file("extdata", "mothur_example.taxonomy", package = "mia")
    meta <- system.file("extdata", "mothur_example.design", package = "mia")
    se <- importMothur(counts)
    se2 <- importMothur(counts)
    expect_s4_class(se, "SummarizedExperiment")
    expect_s4_class(se2, "SummarizedExperiment")
    se <- importMothur(counts, taxa)
    se2 <- importMothur(counts, taxa2)
    expect_s4_class(se, "SummarizedExperiment")
    expect_s4_class(se2, "SummarizedExperiment")
    se <- importMothur(counts, col.file = meta)
    se2 <- importMothur(counts, col.file = meta)
    expect_s4_class(se, "SummarizedExperiment")
    expect_s4_class(se2, "SummarizedExperiment")
    se <- importMothur(assay.file = counts, row.file = taxa, col.file = meta)
    se2 <- importMothur(assay.file = counts, row.file = taxa2, col.file = meta)
    expect_s4_class(se, "SummarizedExperiment")
    expect_s4_class(se2, "SummarizedExperiment")

    # Checks dimensions, rownames, and colnames of assay
    expect_equal(nrow(assays(se)$counts), 100)
    expect_equal(rownames(assays(se)$counts)[1:10],
                        c("Otu001", "Otu002",  "Otu003",  "Otu004",  "Otu005",
                            "Otu006", "Otu007", "Otu008", "Otu009", "Otu010"))
    expect_equal(ncol(assays(se)$counts), 100)
    expect_equal(colnames(assays(se)$counts)[1:10],
                        c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5",
                            "Sample6", "Sample7", "Sample8", "Sample9", "Sample10"))
    expect_equal(nrow(assays(se2)$counts), 100)
    expect_equal(rownames(assays(se2)$counts)[1:10],
                c("Otu001", "Otu002",  "Otu003",  "Otu004",  "Otu005",
                    "Otu006", "Otu007", "Otu008", "Otu009", "Otu010"))
    expect_equal(ncol(assays(se2)$counts), 100)
    expect_equal(colnames(assays(se)$counts)[1:10],
                        c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5",
                            "Sample6", "Sample7", "Sample8", "Sample9", "Sample10"))

    # Checks that rowData has right dimensions, rownames, and colnames
    expect_equal(nrow(rowData(se)), 100)
    expect_equal(rownames(rowData(se))[1:10],
                        c("Otu001", "Otu002",  "Otu003",  "Otu004",  "Otu005",
                            "Otu006", "Otu007", "Otu008", "Otu009", "Otu010"))
    expect_equal(colnames(rowData(se)),
                c("OTU", "Size", "Kingdom", "Phylum", "Order", "Class", "Family", "Genus"))
    expect_equal(nrow(rowData(se2)), 100)
    expect_equal(rownames(rowData(se2))[1:10],
                        c("Otu001", "Otu002",  "Otu003",  "Otu004",  "Otu005",
                            "Otu006", "Otu007", "Otu008", "Otu009", "Otu010"))
    expect_equal(colnames(rowData(se2)),
                c("OTU", "Kingdom", "Phylum", "Order", "Class", "Family", "Genus"))

    expect_equal(nrow(colData(se)), 100)
    expect_equal(rownames(colData(se))[1:10],
                        c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5",
                            "Sample6", "Sample7", "Sample8", "Sample9", "Sample10"))

    # Checks colData's dimensions and names of columns and rows
    expect_equal(colnames(colData(se)),
                        c("group", "sex", "age", "drug", "label", "numOtus", "Group"))
    expect_equal(nrow(colData(se2)), 100)
    expect_equal(rownames(colData(se2))[1:10],
                        c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5",
                            "Sample6", "Sample7", "Sample8", "Sample9", "Sample10"))
    expect_equal(colnames(colData(se2)),
                        c("group", "sex", "age", "drug", "label", "numOtus", "Group"))
})

assay.file <- system.file("extdata", "table.qza", package = "mia")
row.file <- system.file("extdata", "taxonomy.qza", package = "mia")
col.file <- system.file("extdata", "sample-metadata.tsv", package = "mia")
refseq.file <- system.file("extdata", "refseq.qza", package = "mia")

test_that("make TSE worked properly while no sample or taxa data", {
    ## no sample data or taxa data
    expect_silent(tse <- importQIIME2(assay.file))
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_equal(dim(tse), c(770,34))
})

test_that("reference sequences of TSE", {
    skip_if_not(require("biomformat", quietly = TRUE))
    # 1. fasta file of refseq
    tse <- importQIIME2(
        assay.file,
        refseq.file = refseq.file
    )
    tse2 <-  importQIIME2(
        assay.file,
        refseq.file = refseq.file,
        featureNamesAsRefseq = FALSE
    )
    expect_identical(tse@referenceSeq, importQZA(refseq.file)[rownames(tse), ])
    expect_identical(tse2@referenceSeq, importQZA(refseq.file)[rownames(tse), ])

    # 2. row.names of feature table as refseq
    # 2.1 element of row.names of feature table is not DNA sequence
    tse <- importQIIME2(
        assay.file,
        featureNamesAsRefseq = TRUE
    )
    expect_null(tse@referenceSeq)

    # 2.2 element of row.names of feature table is a DNA sequence
    #
    # codes used for create sample data (donot run)
    if (FALSE) {
        .require_package("biomformat")
        feature_tab <- importQZA(assay.file)
        n_feature <- nrow(feature_tab)
        random_seq <- sapply(
            rep(20, n_feature),
            function(size) {
                paste(
                    sample(Biostrings::DNA_BASES, size, replace = TRUE),
                    collapse=""
                )
            }
        )
        row.names(feature_tab) <- random_seq
        obj <- biomformat::make_biom(feature_tab)
        biomformat::write_biom(obj, "data/table-rownamesInSeq.biom")
        # create qza file using QIIME2
        # qiime tools import \
        #     --input-path data/table-rownamesInSeq.biom
        #     --type 'FeatureTable[Frequency]'
        #     --input-format BIOMV100Format
        #     --output-path inst/extdata/table-rownamesInSeq.qza
        unlink("data/table-rownamesInSeq.biom")
    }
    featureTableFile2 <- system.file(
        "extdata",
        "table-rownamesInSeq.qza",
        package = "mia"
    )

    # featureNamesAsRefseq is TRUE, refseq.file is NULL, set row.names of
    # feature table as reference sequences
    tse <- importQIIME2(
        featureTableFile2,
        featureNamesAsRefseq = TRUE
    )
    feature_tab <- importQZA(featureTableFile2)
    names_seq <- Biostrings::DNAStringSet(row.names(feature_tab))
    names(names_seq) <- names_seq
    expect_identical(tse@referenceSeq, names_seq)

    # refseq.file is not NULL, featureNamesAsRefseq is TRUE,
    # set the sequences from refseq.file as reference sequences.
    tse <- importQIIME2(
        assay.file,
        featureNamesAsRefseq = TRUE,
        refseq.file = refseq.file
    )
    expect_identical(tse@referenceSeq, importQZA(refseq.file)[rownames(tse)])

    # 3. refseq.file = NULL, featureNamesAsRefseq = FALSE
    tse <- importQIIME2(
        assay.file,
        refseq.file = NULL,
        featureNamesAsRefseq = FALSE
    )
    expect_null(tse@referenceSeq)
})


test_that("`.parse_taxonomy` work with any combination of taxonomic ranks", {
    test_taxa <- matrix(
        c("a", "k__Bacteria; c__Bacteroidia; s__", 0.88,
            "b", "k__Bacteria; c__Clostridia; s__", 0.9),
        ncol = 3,
        byrow = TRUE,
        dimnames = list(c("a", "b"), c("Feature.ID", "Taxon", "Confidence"))
    )
    expect_silent(mia:::.parse_taxonomy(test_taxa))

    # at certain rank (e.g. species): some taxa can not be determined which
    # species it assigned to (NA)
    test_taxa <- matrix(
        c("a", "k__Bacteria; c__Bacteroidia; s__", 0.88,
            "b", "k__Bacteria; c__Clostridia;", 0.9),
        ncol = 3,
        byrow = TRUE,
        dimnames = list(c("a", "b"), c("Feature.ID", "Taxon", "Confidence"))
    )
    expect_true(is.na(mia:::.parse_taxonomy(test_taxa)[2,"species"]))

    # if the expexted order is not present it will return a correct result
    test_taxa <- matrix(
        c("a", "k__Bacteria; s__test; c__Bacteroidia", 0.88,
            "b", "k__Bacteria; c__Clostridia;", 0.9),
        ncol = 3,
        byrow = TRUE,
        dimnames = list(c("a", "b"), c("Feature.ID", "Taxon", "Confidence"))
    )
    expect_equal(mia:::.parse_taxonomy(test_taxa)[,"species"],c("s__test",NA))
    expect_equal(mia:::.parse_taxonomy(test_taxa, prefix.rm = TRUE)[,"species"],
                c("test",NA))
    # Expect that output have only taxonomy levels that have information
    expect_equal(ncol(mia:::.parse_taxonomy(test_taxa)), 3L)
})

test_that("`.read_q2sample_meta` remove  the row contained `#q2:types`", {
    expect_false(any(as(.read_q2sample_meta(col.file), "matrix") == "#q2:types"))
})

test_that('get file extension', {
    expect_identical(.get_ext("a.b.c"), "c")
    expect_identical(.get_ext("abc/a.b.c"), "c")
})

test_that('read qza file', {
    expect_error(importQZA("abc"), "does not exist")
    expect_error(importQZA(col.file), "must be in `qza` format")
})

test_that("Confidence of taxa is numberic", {
    skip_if_not(require("biomformat", quietly = TRUE))
    tse <- importQIIME2(
        assay.file,
        row.file = row.file
    )
    expect_true(is.numeric(S4Vectors::mcols(tse)$Confidence))
})

test_that("dimnames of feature table is identicle with meta data", {
    skip_if_not(require("biomformat", quietly = TRUE))
    feature_tab <- importQZA(assay.file)

    sample_meta <- .read_q2sample_meta(col.file)
    taxa_meta <- importQZA(row.file)
    data_list <- .set_feature_tab_dimnames(
        feature_tab,
        sample_meta,
        taxa_meta
   )
    expect_identical(colnames(data_list[[1]]), colnames(feature_tab))
    expect_identical(rownames(data_list[[1]]), rownames(feature_tab))
    expect_identical(rownames(data_list[[1]]), rownames(data_list[[2]]))
    expect_identical(colnames(data_list[[1]]), rownames(data_list[[3]]))
    expect_identical(colnames(data_list[[1]]), colnames(data_list[[1]]))

    # sample_meta or feature meta is NULL
    sample_meta2 <- S4Vectors::make_zero_col_DFrame(ncol(feature_tab))
    rownames(sample_meta2) <- colnames(feature_tab)
    taxa_meta2 <- S4Vectors::make_zero_col_DFrame(nrow(feature_tab))
    rownames(taxa_meta2) <- rownames(feature_tab)
    expect_silent(.set_feature_tab_dimnames(feature_tab, sample_meta2, taxa_meta))

    # sample meta or feature meta without any information, only contains sample/feature
    # ID in its rownames
    feature_tab3 <- S4Vectors::DataFrame(
        sample1 = 1:3,
        sample2 = 4:6,
        sample3 = 7:9,
        row.names = paste0("feature", 1:3)
    )
    sample_meta3 <- S4Vectors::DataFrame(row.names = paste0("sample", 3:1))
    feature_meta3 <- S4Vectors::DataFrame(row.names = paste0("feature", c(2, 3, 1)))
    data_list <- .set_feature_tab_dimnames(
        feature_tab3,
        sample_meta3,
        feature_meta3
    )
    expect_identical(rownames(data_list[[1]]), rownames(feature_tab3))
    expect_identical(colnames(data_list[[1]]), colnames(feature_tab3))
    expect_identical(rownames(data_list[[2]]), rownames(feature_tab3))
    expect_identical(rownames(data_list[[3]]), colnames(feature_tab3))
})


test_that("convertToPhyloseq", {

    skip_if_not_installed("phyloseq")

    # TSE object
    data(GlobalPatterns, package="mia")
    tse <- GlobalPatterns

    phy <- convertToPhyloseq(GlobalPatterns)

    # Test that assay is in otu_table
    expect_equal(as.data.frame(phyloseq::otu_table(phy)@.Data), as.data.frame(assays(tse)$counts))

    # Test that rowData is in tax_table
    expect_equal(as.data.frame(phyloseq::tax_table(phy)@.Data), as.data.frame(rowData(tse)))

    # Test that colData is in sample_table
    expect_equal(phyloseq::sample_data(phy),
                phyloseq::sample_data(data.frame(colData(tse))))

    # Test that rowTree is in phy_tree
    expect_identical(phyloseq::phy_tree(phy), rowTree(tse))

    # Test that referenceSeq is in refseq. Expect error, because there should not be
    # reference sequences.
    expect_error(phyloseq::refseq(phy))

    # Test with agglomeration that that pruning is done internally
    test1 <- agglomerateByRank(tse, rank = "Phylum", update.tree = FALSE)
    test2 <- agglomerateByRank(tse, rank = "Phylum", update.tree = TRUE)
    test1_phy <- expect_warning(convertToPhyloseq(test1))
    test2_phy <- convertToPhyloseq(test2)

    expect_equal(length(phyloseq::phy_tree(test1_phy)$node),
                length(ape::keep.tip(rowTree(test1), rowLinks(test1)$nodeLab)$node))
    expect_equal(phyloseq::phy_tree(test1_phy)$tip.label, rownames(test2))
    # The tip labels do not match because of renaming
    expect_identical(phyloseq::phy_tree(test2_phy)$edge, rowTree(test2)$edge)

    # Check that everything works also with agglomerated data
    for (level in colnames(rowData(tse)) ){
        temp <- agglomerateByRank(tse, rank = level)
        expect_no_warning(convertToPhyloseq(temp))
    }

    tse2 <- tse
    # Concerts data frame to factors
    rowData(tse2) <- DataFrame(lapply(rowData(tse2), as.factor))
    phy <- convertToPhyloseq(tse)
    phy2 <- convertToPhyloseq(tse2)
    expect_equal(phyloseq::tax_table(phy2), phyloseq::tax_table(phy))

    # TSE object
    data(esophagus, package="mia")
    tse <- esophagus

    phy <- convertToPhyloseq(tse, assay.type="counts")

    # Test that assay is in otu_table
    expect_equal(as.data.frame(phyloseq::otu_table(phy)@.Data), as.data.frame(assays(tse)$counts))

    # Test that rowTree is in phy_tree
    expect_identical(phyloseq::phy_tree(phy), rowTree(tse))

    # Test that merging objects lead to correct phyloseq
    tse <- mergeSEs(GlobalPatterns, esophagus, assay.type="counts", missing.values = 0)
    pseq <- convertToPhyloseq(tse, assay.type="counts")

    # Include rownames from both trees
    tse_compare <- tse[ c(rownames(GlobalPatterns), rownames(esophagus)), ]
    pseq_compare <- convertToPhyloseq(tse_compare, assay.type="counts")

    expect_equal(phyloseq::otu_table(pseq), phyloseq::otu_table(pseq_compare))
})

test_that("Import HUMAnN file", {
    file_path <- system.file("extdata", "humann_output.tsv", package = "mia")
    tse <- importHUMAnN(file_path)
    #
    expect_true( length(taxonomyRanks(tse)) == 2 )
    expect_true( ncol(rowData(tse)) == 4 )
    expect_true( ncol(tse) == 3 )
    expect_true( nrow(tse) == 12 )
    expect_true( all(!is.na(assay(tse))) )
    #
    expect_error(importHUMAnN("test"))
    expect_error(importHUMAnN(file_path, remove.suffix="test"))
    expect_error(importHUMAnN(file_path, remove.suffix=1))
    expect_error(importHUMAnN(file_path, remove.suffix=c(FALSE, TRUE)))
    #
    tse2 <- importHUMAnN(file_path, remove.suffix=TRUE)
    # There is no suffix, should be equal
    expect_equal(tse, tse2)
    #
    cd <- colData(tse)
    tse2 <- importHUMAnN(file_path, colData = cd)
    # colData should not change the result
    expect_equal(tse, tse2)
})
