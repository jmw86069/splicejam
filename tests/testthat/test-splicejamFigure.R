# splicejamFigure
test_that("splicejamFigure-farrisdata-Gria1", {
   #
   testthat::skip_if_not_installed("farrisdata")
   testthat::skip_if_not_installed("vdiffr")

   # use test data
   data(sjenvtest)
   detectedGenes <- sjenvtest$detectedGenes;
   detectedTx <- sjenvtest$detectedTx;

   # Assert names(flatExonsByGene)
   testthat::expect_contains(
      detectedGenes,
      names(sjenvtest$flatExonsByGene))
   # Assert names(flatExonsByTx)
   testthat::expect_contains(
      detectedTx,
      names(sjenvtest$flatExonsByTx))

   Gria1_default <- function() {
      splicejamFigure(sjenv=sjenvtest,
         use_memoise=TRUE,
         gene="Gria1")
   }
   vdiffr::expect_doppelganger("splicejamFigure-Gria1-default",
      Gria1_default)

})
