
# 30apr2026 new comments

* `flatExonsByGene`,`flatExonsByTx` could be done per-gene as needed,
instead of adding to upfront startup cost for R-shiny.
* Ideal world: SashimiData S4 object

   * derive from gtf, or txdb, or EnsDb, or whatever
   * API interface for other components - it handles cache, conversions, etc.
   * obscures whether underlying data is EnsDb, txdb/gtf, etc.
   * Then `prepareSashimi()` only needs this object as input.


# splicejam compatibility with TxDb, EnsDb

* TxDb can "intuit" the gene_name, flatExons, etc.

   * knownGene uses ENTREZID as 'GENEID' (I did not know that.)
   * Default could use gene_id also as gene_name if no other option exists.

* EnsDb can create tx2geneDF, flatExons, etc.

   * `EnsDb` object is enough to support gene/tx/exon needs of splicejam.

* GTF file can support splicejam

   * GTF -> `TxDb`
   * tx2geneDF
   * txdb -> exonsByTx, cdsByTx -> flatExonsByGene,flatExonsByTx


# splicejam refactoring ideas


`Splicejam` S4 object?

   * contains "whatever it needs" to create sashimi plots
   * Gene-exon data
   
      * flatExonsByGene
      * flatExonsByTx
      * tx2geneDF
      * detectedGenes
      * detectedTx

* `SplicejamCoverage` S4 object

   * Purpose is to describe coverage and junction data
   * filesDF - is it enough as-is?
   * Future goal to support BAM alignments to produce coverage
   * BAM -> GenomicAlignments -> covGR `GRanges` objects with coverages
   I had prototyped it long ago.


* `SplicejamGenes` S4 object?

   * Purpose is to handle the gene-transcript-exon data in one place
   
   * Required data for sashimi plots:
   
      * flatExonsByGene - derive from txdb, or user
      * flatExonsByTx - derive from txdb, or user
      * tx2geneDF - derive from gtf, or user
      * detectedGenes (optional) - from user
      * detectedTx (optional) - from user

   * Paths to derive:
   
      * tx2geneDF - gtf, `data.frame`, `valules(flatExonsByGene)`?
      * flatExonsByGene <- txdb, `GRangesList`, gtf
      * flatExonsByTx <- txdb, `GRangesList`, gtf
      * txdb, tx2geneDF - sufficient to create: flatExonsByGene, flatExonsByTx
      * gtf - sufficient to create: tx2geneDF, txdb, flatExonsByGene, flatExonsByTx

   * Gene-exon data
   
      * flatExonsByGene
      * flatExonsByTx
      * tx2geneDF
      * detectedGenes
      * detectedTx

   * Create using txdb?

