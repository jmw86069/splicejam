## code to prepare `sjenvtest`
##
detectedGenes <- c("Gria1", "Ntrk3")
detectedTx <- subset(
   SummarizedExperiment::rowData(farrisdata::farrisTxSE),
   TxDetectedByTPM %in% TRUE &
      geneSymbol %in% detectedGenes)$transcript_id;
sjenvtest <- new.env();
sjenvtest <- sashimiDataConstants(envir=sjenvtest,
   detectedGenes=detectedGenes,
   detectedTx=detectedTx,
   empty_uses_farrisdata=TRUE,
   verbose=FALSE)


sjenvtest$cdsByTx <- sjenvtest$cdsByTx[names(sjenvtest$cdsByTx) %in% detectedTx]
sjenvtest$exonsByTx <- sjenvtest$exonsByTx[names(sjenvtest$exonsByTx) %in% detectedTx]
sjenvtest$tx2geneDF <- subset(sjenvtest$tx2geneDF, transcript_id %in% detectedTx)
rownames(sjenvtest$tx2geneDF) <- sjenvtest$tx2geneDF$transcript_id;
sjenvtest$color_sub <- farrisdata::colorSub;
sjenvtest$filesDF <- subset(farrisdata::farris_sashimi_files_df,
   grepl("CA[12]", sample_id))

usethis::use_data(sjenvtest, overwrite=TRUE)

######################################
## Minimalist testing
sjenvtest1 <- new.env()
keep_minimal <- c(
   "tx2geneDF",
   "filesDF",
   "flatExonsByGene",
   "flatExonsByTx"
   # "detectedGenes", # optional?
   # "detectedTx",    # optional?
   # "color_sub"      # optional
)
for (i in keep_minimal) {
   assign(x=i,
      value=get(i, envir=sjenvtest),
      envir=sjenvtest1)
}
jamba::sdim(sjenvtest1)

sjenvtest2 <- sashimiDataConstants(envir=sjenvtest1)
jamba::sdim(sjenvtest2)
system.time(sjfigtest2 <- splicejamFigure(sjenv=sjenvtest2, gene="Gria1", use_memoise=TRUE))
# 30 sec

system.time(sjfigtest2b <- splicejamFigure(sjenv=sjenvtest2, gene="Gria1", use_memoise=TRUE))
# 7 sec
