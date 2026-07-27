# splicejamFigure
test_that("splicejamFigure-farrisdata-Gria1", {
   #
   testthat::skip_if_not_installed("farrisdata")
   testthat::skip_if_not_installed("vdiffr")

   # use test data
   data(sjenvtest)
   detectedGenes <- c("Gria1", "Ntrk3");
   detectedTx <- c(
      'ENSMUST00000036315.15',
      'ENSMUST00000039431.13',
      'ENSMUST00000039438.8',
      'ENSMUST00000094179.10',
      'ENSMUST00000151885.2',
      'ENSMUST00000193002.5',
      'ENSMUST00000195262.5',
      'ENSMUST00000205354.1',
      'ENSMUST00000206268.1',
      'ENSMUST00000206949.1')

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
         use_memoise=FALSE,
         gene="Gria1")
   }
   vdiffr::expect_doppelganger("splicejamFigure-Gria1-default",
      Gria1_default)

})
