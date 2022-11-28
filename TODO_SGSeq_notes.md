
* SGSeq: https://bioconductor.org/packages/3.16/bioc/vignettes/SGSeq/inst/doc/SGSeq.html

* `plotCoverage()` shows coverage by exon, and arcs per splice
* `plotSpliceGraph()` plots schematic splice graph, flattened per gene.
   
* TxFeatures class: `importTranscripts()` from GFF/GTF file, or `GRangesList`

   * code converts `TxDb` using `convertToTxFeatures()`
   * this function converts exons to `GRanges`
   * it creates junction regions between exons, also `GRanges`
   * exons have one or more `txName`, (CharacterList),
   and one or more `geneName` (CharacterList, maybe multiple are allowed?)
   * exons have a type:
   
      * J (splice junction)
      * I (internal exon)
      * F (first/5'-terminal exon)
      * L (last/5'-terminal exon)
      * U (unspliced transcript).

* SGFeatures class:

   * `TxFeatures` can be converted with `convertToSGFeatures()`
   * model of disjoint exons spliced into transcripts
   * exons have a type:
   
      * J (splice junction) - this junction could be used in splicejam
      * E (disjoint exon bin) - this exon could be used in splicejam
      * D (splice donor site) - this point feature is not used in splicejam
      * A (splice acceptor site) - this point feature is not used in splicejam
   
   * `featureID` unique identifier
   * `txName` CharacterList of transcripts; `geneName` CharacterList of genes

* Input typically by BAM file, STAR or GSNAP that custom tag `XS` for spliced reads

   * `si` (sample information) is a `data.frame` with these columns:
   
      * sample_name: `character` unique sample label
      * file_bam: `character` BAM alignment file
      * paired_end: `logical` indicating whether reads were paired
      * read_length: `integer` length of each read
      * frag_length: `integer` fragment length
      * lib_size: `integer` total aligned fragments

   * `getBamInfo()` can be used to create this `data.frame`:
   
      * `path <- system.file("extdata", package = "SGSeq")`
      * `si$file_bam <- file.path(path, "bams", si$file_bam)`
