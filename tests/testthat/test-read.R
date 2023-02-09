# This tests the readH5AD function (and by implication, SCE2AnnData).
library(SummarizedExperiment)
file <- system.file("extdata", "krumsiek11.h5ad", package = "zellkonverter")

test_that("Reading H5AD works", {
    sce <- readH5AD(file)
    expect_s4_class(sce, "SingleCellExperiment")

    expect_identical(assayNames(sce), "X")
    expect_identical(colnames(colData(sce)), "cell_type")
})

test_that("Reading H5AD works with version 0.7.6", {
    sce <- readH5AD(file, version = "0.7.6")
    expect_s4_class(sce, "SingleCellExperiment")

    expect_identical(assayNames(sce), "X")
    expect_identical(colnames(colData(sce)), "cell_type")
})

test_that("Reading H5AD works with verbose=TRUE", {
    sce <- readH5AD(file, verbose = TRUE)
    expect_s4_class(sce, "SingleCellExperiment")

    expect_identical(assayNames(sce), "X")
    expect_identical(colnames(colData(sce)), "cell_type")
})

test_that("Reading H5AD works with HDF5Arrays", {
    sce <- readH5AD(file, use_hdf5 = TRUE)
    expect_s4_class(sce, "SingleCellExperiment")
    expect_s4_class(DelayedArray::seed(assay(sce)), "HDF5ArraySeed")

    ref <- readH5AD(file)
    expect_identical(as.matrix(assay(ref)), as.matrix(assay(sce)))

    # Properly sleeps to wait for the process to shut down.
    expect_s4_class(
        sce <- readH5AD(file, use_hdf5 = TRUE),
        "SingleCellExperiment"
    )
})

test_that("Reading H5AD works with a mixture of sparse and HDF5Arrays", {
    sce <- readH5AD(file)
    assay(sce, "more") <- as(assay(sce, "X"), "CsparseMatrix")

    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp)

    backed <- readH5AD(temp, use_hdf5 = TRUE)
    expect_s4_class(DelayedArray::seed(assay(backed)), "HDF5ArraySeed")
    expect_s4_class(assay(backed, "more"), "CsparseMatrix")
})

test_that("readH5AD works in a separate process", {
    oldshare <- basilisk::getBasiliskShared()
    basilisk::setBasiliskShared(FALSE)
    oldfork <- basilisk::getBasiliskFork()
    basilisk::setBasiliskFork(FALSE)

    sce <- readH5AD(file)
    expect_s4_class(sce, "SingleCellExperiment")

    basilisk::setBasiliskShared(oldshare)
    basilisk::setBasiliskFork(oldfork)
})

test_that("Reading H5AD works with native reader", {
    sce <- readH5AD(file, reader = "R")
    expect_s4_class(sce, "SingleCellExperiment")

    expect_identical(assayNames(sce), "X")
    expect_identical(colnames(colData(sce)), "cell_type")
})

test_that("Reading v0.8 H5AD works with native reader", {
    sce <- readH5AD(file)

    # Test float with NA
    colData(sce)$dummy_num <- 42.42
    colData(sce)$dummy_num[1] <- NA

    # Test nullable ints
    # HACK: ints work, but only by accident. Should follow the v0.8
    # spec and read in the NA mask properly
    colData(sce)$dummy_int <- as.integer(42)
    colData(sce)$dummy_int[1] <- NA

    ## Test simple strings
    ## NOTE: writeH5AD seems to write strings as factors so we can't
    ## test this yet
    #colData(sce)$dummy_str <- "hello"
    #colData(sce)$dummy_str[1] <- "world"

    ## Test nullable booleans
    ## FIXME: native R is reading in the bools as chars
    #colData(sce)$dummy_bool <- TRUE
    #colData(sce)$dummy_bool[1] <- FALSE
    #colData(sce)$dummy_bool[2] <- NA

    temp <- tempfile(fileext = ".h5ad")
    writeH5AD(sce, temp, version = "0.8")

    sce2 <- readH5AD(temp, reader='R')

    # NOTE: colData columns not in same order. Is this a problem?
    # Just re-sort them for now.
    expect_true(all(colnames(colData(sce2)) %in% colnames(colData(sce))))
    colData(sce2) <- colData(sce2)[,colnames(colData(sce))]

    # sce already contains a factor column (cell_type) so we didn't
    # add a dummy one. But double check it's there as expected
    expect_identical(colnames(colData(sce2))[1], "cell_type")
    expect_identical(class(colData(sce2)$cell_type), "factor")

    expect_identical(rownames(sce), rownames(sce2))
    expect_identical(colnames(sce), colnames(sce2))

    expect_identical(assays(sce), assays(sce2))
    expect_identical(metadata(sce), metadata(sce2))
    expect_identical(colData(sce), colData(sce2))
    expect_identical(rowData(sce), rowData(sce2))
})

test_that("Skipping slot conversion works", {
    sce <- readH5AD(file, layers = FALSE, uns = FALSE, var = FALSE, obs = FALSE,
                    varm = FALSE, obsm = FALSE, varp = FALSE, obsp = FALSE)

    expect_identical(assayNames(sce), "X")
    expect_identical(metadata(sce), list())
    expect_equal(ncol(rowData(sce)), 0)
    expect_equal(ncol(colData(sce)), 0)
    expect_equal(length(reducedDims(sce)), 0)
    expect_equal(length(rowPairs(sce)), 0)
    expect_equal(length(colPairs(sce)), 0)
})

test_that("Selective slot conversion works", {
    sce <- readH5AD(file, uns = "iroot")

    expect_identical(names(metadata(sce)), "iroot")
})

test_that("Selective DF conversion works", {
    sce <- readH5AD(file, obs = "cell_type")

    expect_identical(names(colData(sce)), "cell_type")
})

test_that("Conversion of raw works", {
    skip_if_offline()

    cache <- BiocFileCache::BiocFileCache(ask = FALSE)
    example_file <- BiocFileCache::bfcrpath(
        cache, "https://ndownloader.figshare.com/files/30462915"
    )

    sce <- readH5AD(example_file, raw = TRUE)

    names <- list(
        assays = c("X"),
        colData = c("n_genes", "n_genes_by_counts", "total_counts",
                    "total_counts_mt", "pct_counts_mt", "leiden"),
        rowData = c("gene_ids", "n_cells", "mt", "n_cells_by_counts",
                    "mean_counts", "pct_dropout_by_counts", "total_counts",
                    "highly_variable", "means", "dispersions",
                    "dispersions_norm", "mean", "std"),
        metadata = c("hvg", "leiden", "neighbors", "pca", "rank_genes_groups",
                     "umap"),
        redDim = c("X_pca", "X_umap"),
        varm = c("PCs"),
        colPairs = c("connectivities", "distances"),
        raw_rowData = c("gene_ids", "n_cells", "mt", "n_cells_by_counts",
                        "mean_counts", "pct_dropout_by_counts", "total_counts",
                        "highly_variable", "means", "dispersions",
                        "dispersions_norm")
    )

    missing <- list()

    validateH5ADSCE(sce, names, missing)
})
