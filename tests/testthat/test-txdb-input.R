
test_that("makeTx2geneFromTxdb returns expected columns and content", {
   skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
   skip_if_not_installed("org.Mm.eg.db")
   skip_if_not_installed("genejam")

   txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
   result <- makeTx2geneFromTxdb(txdb, ann_lib = "org.Mm.eg.db")

   expect_s3_class(result, "data.frame")
   expect_true(all(c("transcript_id", "gene_id", "gene_name") %in% colnames(result)))
   expect_gt(nrow(result), 0)
   # Gria1 should be present as a gene_name
   expect_true("Gria1" %in% result$gene_name)
   # No NA values in any column
   expect_false(any(is.na(result$transcript_id)))
   expect_false(any(is.na(result$gene_id)))
   expect_false(any(is.na(result$gene_name)))
})

test_that("splicejamDataFromTxDb returns expected environment structure", {
   skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
   skip_if_not_installed("org.Mm.eg.db")
   skip_if_not_installed("genejam")

   txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
   env <- splicejamDataFromTxDb(
      txdb         = txdb,
      ann_lib      = "org.Mm.eg.db",
      detectedGenes = "Gria1")

   expect_true(is.environment(env))
   # Required objects
   expect_true(exists("flatExonsByGene", envir = env))
   expect_true(exists("flatExonsByTx",   envir = env))
   expect_true(exists("tx2geneDF",       envir = env))
   expect_true(exists("detectedTx",      envir = env))
   expect_true(exists("detectedGenes",   envir = env))
   # Gria1 should be the only gene (and present)
   expect_identical(env$detectedGenes, "Gria1")
   expect_true("Gria1" %in% names(env$flatExonsByGene))
   # tx2geneDF should only contain Gria1 rows
   expect_true(all(env$tx2geneDF$gene_name == "Gria1"))
   # detectedTx should match tx2geneDF
   expect_setequal(env$detectedTx, env$tx2geneDF$transcript_id)
})

test_that("splicejamDataFromTxDb handles detectedTx subsetting", {
   skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
   skip_if_not_installed("org.Mm.eg.db")
   skip_if_not_installed("genejam")

   txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
   # Build full tx2geneDF first to get valid transcript IDs for Gria1
   full_df  <- makeTx2geneFromTxdb(txdb, ann_lib = "org.Mm.eg.db")
   gria1_tx <- full_df$transcript_id[full_df$gene_name == "Gria1"]
   # Use only the first two transcripts
   subset_tx <- head(gria1_tx, 2)

   env <- splicejamDataFromTxDb(
      txdb       = txdb,
      ann_lib    = "org.Mm.eg.db",
      detectedTx = subset_tx)

   expect_setequal(env$detectedTx, subset_tx)
   expect_true("Gria1" %in% env$detectedGenes)
})

test_that("splicejamDataFromTxDb derives color_sub from filesDF", {
   skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
   skip_if_not_installed("org.Mm.eg.db")
   skip_if_not_installed("genejam")
   skip_if_not_installed("colorjam")

   txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene

   # Reuse filesDF from sjenvtest (same mm10 genome)
   data(sjenvtest, package = "splicejam")
   fDF <- sjenvtest$filesDF

   env <- splicejamDataFromTxDb(
      txdb          = txdb,
      ann_lib       = "org.Mm.eg.db",
      detectedGenes = "Gria1",
      filesDF       = fDF)

   expect_true(exists("filesDF",   envir = env))
   expect_true(exists("color_sub", envir = env))
   # color_sub names should cover all sample_ids in filesDF
   expect_true(all(unique(fDF$sample_id) %in% names(env$color_sub)))
})

test_that("splicejamDataFromTxDb silently drops unknown detectedTx/detectedGenes", {
   skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
   skip_if_not_installed("org.Mm.eg.db")
   skip_if_not_installed("genejam")

   txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
   env <- splicejamDataFromTxDb(
      txdb          = txdb,
      ann_lib       = "org.Mm.eg.db",
      detectedGenes = c("Gria1", "NOT_A_REAL_GENE_XYZ"))

   # Unknown gene should be dropped; Gria1 should still be present
   expect_false("NOT_A_REAL_GENE_XYZ" %in% env$detectedGenes)
   expect_true("Gria1" %in% env$detectedGenes)
})
