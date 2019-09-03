#' @import jamba
NULL

#' Make tx2gene data.frame from a GTF file
#'
#' Make tx2gene data.frame from a GTF file
#'
#' Create a transcript-to-gene data.frame from a GTF file, which is required
#' by a number of transcriptome analysis methods such as those in
#' the DEXseq package, and the limma package functions such as
#' \code{diffSplice()}.
#'
#' This function also only uses \code{data.table::fread()} and does not
#' import the full GTF file using something like Bioconductor
#' \code{GenomicFeatures}, simply because the data.table method is markedly
#' faster when importing only the transcript-to-gene relationship. Also, this
#' method allows the import of more annotations than are supported by the
#' typical Bioconductor \code{rtracklayer::import()} for GTF data.
#'
#' This function is intended to help keep all transcript data consistent by
#' using the same GTF file that is also used by other analysis tools, whether
#' those tools be based in R or more likely, outside R in a terminal
#' environment.
#'
#' For example, the GTF file could be used:
#'
#' \itemize{
#'    \item{to run STAR sequence alignment
#'    then `Rsubread::featureCounts()` to generate a matrix of read
#'    counts per gene, transcript, or exon; or}
#'    \item{to generate a transcript
#'    FASTA sequence file then run a kmer quantitation tool such as
#'    Salmon or Kallisto, then using `tximport::tximport()` to import
#'    results into R for downstream processing.}
#' }
#'
#' @param GTF character file name sent to `data.table::fread()`. When the
#'    file ends with ".gz", the `R.utils` package is recommended, otherwise
#'    the fallback option is to make a system call to `gzcat`
#'    to gunzip the file during the import step. Note this process fails
#'    when `gzcat` is not available in the path of the user environment.
#'    In general, the `R.utils` package is the best solution.
#' @param geneAttrNames character vector of recognized attribute names
#'    as they appear in column 9 of the GTF file, for gene rows.
#' @param txAttrNames character vector of recognized attribute names
#'    as they appear in column 9 of the GTF file, for transcript rows.
#' @param geneFeatureType character value to match column 3 of the GTF
#'    file, used to define gene rows, by default "gene".
#' @param txFeatureType character value to match column 3 of the GTF
#'    file, used to define gene rows, by default "transcript". In some
#'    GTF files, "mRNA" is used, so either is accepted by default.
#' @param nrows integer number of rows to read from the GTF file, by default
#'    -1 means all rows are imported. This parameter is useful to check the
#'    results of a large GTF file using only a subset portion of the file.
#' @param verbose logical whether to print verbose output during processing.
#'
#' @return
#' `data.frame` with colnames indicated by the values in
#' `geneAttrNames` and `txAttrNames`.
#'
#' @family jam RNA-seq functions
#'
#' @export
makeTx2geneFromGtf <- function
(GTF,
 geneAttrNames=c("gene_id", "gene_name", "gene_type"),
 txAttrNames=c("transcript_id", "transcript_type"),
 geneFeatureType="gene",
 txFeatureType=c("transcript","mRNA"),
 nrows=-1L,
 verbose=FALSE,
 ...)
{
   ## Purpose is to use a GTF file to create a transcript to gene
   ## annotation table, tx2geneDF
   ##
   ## geneAttrNames a set of recognized attribute names in column 9 of
   ## the GTF for gene rows.
   ##
   ## txAttrNames a set of recognized attribute names in column 9 of
   ## the GTF for transcript rows.
   ##
   ## geneFeatureType is a value in column 3 of the GTF, used for
   ## gene rows which have gene annotations, by default "gene".
   ##
   ## txFeatureType is a value in column 3 of the GTF, used for
   ## transcript rows which have transcript annotations and no
   ## gene annotations, but will have a link to the parent gene feature.
   ## Usually only "transcript" or "mRNA" are used, so either is acceptable
   ## by default.
   ##
   if (suppressPackageStartupMessages(!require(data.table))) {
      stop(paste0("makeTx2geneFromGtf() requires the data.table package ",
         "to load the GTF file rapidly."));
   }
   if (suppressPackageStartupMessages(!require(jamba))) {
      stop(paste0("makeTx2geneFromGtf() requires the jamba package."));
   }
   if (verbose) {
      printDebug("makeTx2geneFromGtf() :",
         "reading GTF file:",
         GTF);
   }
   if (igrepHas("[.]gz$", GTF) &&
      !suppressPackageStartupMessages(require(R.utils))) {
      if (verbose) {
         jamba::printDebug("makeTx2geneFromGtf() :",
            "using commandline call to gzcat to decompress gtf.gz. ",
            "Install the 'R.utils' package to avoid this step.");
      }
      gtfDF <- data.table::fread(cmd=paste0("gzcat ", GTF),
         sep="\t",
         autostart=20,
         nrows=nrows,
         data.table=FALSE);
   } else {
      if (verbose) {
         jamba::printDebug("makeTx2geneFromGtf() :",
            "using native data.table::fread().");
      }
      gtfDF <- data.table::fread(GTF,
         sep="\t",
         autostart=20,
         nrows=nrows,
         data.table=FALSE);
   }
   ## Subset to clear some memory
   gtfDF <- subset(gtfDF,
      gtfDF[[3]] %in% c(geneFeatureType,
         txFeatureType));

   ## Determine which rows are gene and transcript
   geneRows <- (gtfDF[[3]] %in% geneFeatureType);
   txRows <- (gtfDF[[3]] %in% txFeatureType);


   ## gene attributes
   geneM <- do.call(cbind, lapply(nameVector(geneAttrNames),
      function(attrName){
         if (verbose) {
            jamba::printDebug("makeTx2geneFromGtf() :",
               "gene attributes:",
               attrName);
         }
         attrGrep <- paste0('^.*', attrName, ' ["]([^"]+)["].*$');
         if (jamba::igrepHas(attrGrep, gtfDF[geneRows,][[9]])) {
            attrValues <- gsub(attrGrep,
               "\\1",
               gtfDF[geneRows,,drop=FALSE][[9]]);
         } else {
            if (verbose) {
               jamba::printDebug("makeTx2geneFromGtf(): ",
                  "Note: No gene attributes found for:",
                  attrName);
            }
            attrValues <- NULL;
         }
      }));

   ## transcript attributes
   txM <- do.call(cbind, lapply(nameVector(c(txAttrNames,geneAttrNames)),
      function(attrName){
         if (verbose) {
            jamba::printDebug("makeTx2geneFromGtf(): ",
               "tx attributes:",
               attrName);
         }
         attrGrep <- paste0('^.*', attrName, ' ["]([^"]+)["].*$');
         if (jamba::igrepHas(attrGrep, gtfDF[txRows,][[9]])) {
            attrValues <- gsub(attrGrep,
               "\\1",
               gtfDF[txRows,,drop=FALSE][[9]]);
         } else {
            if (verbose) {
               jamba::printDebug("makeTx2geneFromGtf(): ",
                  "Note: No tx attributes found for:",
                  attrName);
            }
            attrValues <- NULL;
         }
      }));
   if (verbose) {
      jamba::printDebug("makeTx2geneFromGtf() :",
         "Merging gene and transcript annotations.");
   }
   return(merge(geneM,
      txM,
      all.x=TRUE,
      all.y=TRUE));
}

#' tx2ale: detect alternative last exons (ALE) from transcript data
#'
#' detect alternative last exons (ALE) from transcript data
#'
#' This function is intended to encapsulate logic used to define
#' alternative last exons (ALE) based upon gene-transcript models.
#' Specifically, it uses either the 5' end of all 3'UTR regions
#' \code{aleMethod="first"}, or the full 3'UTR range
#' \code{aleMethod="range"}, optionally using only detected transcripts
#' if given by \code{detectedTx}. Any transcripts overlapping the criteria
#' above are merged together by gene, to define a unique set of
#' ALEs for each gene.
#'
#' Note that when using \code{aleMethod="first"},
#' 3'UTRs will be combined if the 5'end of the
#' 3'UTR is identical, which means the length of the 3'UTR is not considered
#' in terms of defining an ALE. The goal is to determine whether the ALE is
#' or is not maintained during transcript processing.
#'
#' Note that when using \code{aleMethod="range"}, 3'UTRs are converted first
#' to a range, which converts multi-exon 3'UTRs to the full range covered by
#' the 3'UTR exons, and not specific exon ranges contained. It is mainly
#' intended to focus on the specific scenario where alternative 3'UTRs for
#' a give gene do not overlap in any way. In reality, several sources of
#' GTF gene models showed unusual transcript isoforms that contained premature
#' stop codons, thus annotating the remaining exons all as "3'UTR" and
#' therefore causing all subsequent 3'UTR regions to be combined into one
#' large spanning range. For this reason, we recommend using
#' \code{aleMethod="first"} by default.
#'
#' We noted substantially improved results when supplying detected transcripts
#' via the parameter \code{detectedTx}, which greatly reduces the full set of
#' possible ALE regions to use only those regions relevant and observed in the
#' given RNA-seq dataset. In other words, there are a large number of possible
#' ALE regions which have no supporting measured data in a particular RNA-seq
#' experiment. It is possible that generating a superset of ALE regions
#' may not be practically problematic in downstream analysis steps, however
#' the impact may be seen predominantly during fase discovery rate (FDR)
#' adjustment for multiple testing, especially if over half the annotated
#' ALE regions have zero reported measurements. If they have no supporting
#' measurements, they are effectively not being tested. As GTF files become
#' increasingly annotated to cover every possible scenario, this issue is
#' likely to become more impactful over time.
#'
#' Input data can be in one of three forms:
#'    * \code{gtf} a GTF file such as those GTF files provided by Gencode
#'    * \code{txdb} a TxDb R object as produced by the Bioconductor
#'    package \code{GenomicFeatures}, suitable for use in calling
#'    \code{GenomicFeatures::threeUTRsByTranscript()}
#'    * \code{threeUtrGRL} a \code{GenomicRanges::GRangesList} object
#'    which is a list of 3'UTR exons by transcript. This input also requires
#'    \code{tx2geneDF} which is used to relate transcripts to genes.
#'
#' Two important outputs are \code{tx2ale} which converts transcripts to
#' ALE names, and if a matrix of transcript expression is supplied with
#' \code{iMatrix}, then \code{iMatrixByALE} is a matrix of expression
#' aggregated at the ALE level, suitable for differential testing in
#' \code{limma::diffSplice()} or \code{DEXSeq}.
#'
#' This function depends upon custom functions: \code{annotateGRLfromGRL()},
#' \code{annotateGRfromGR()}, and \code{assignGRLexonNames()}.

#' @param gtf character path to a GTF file, used to import gene models
#'    using Bioconductor \code{GenomicFeatures::makeTxDbFromGFF()}
#'    to create a transcriptDb object. If supplied, and if \code{tx2geneDF}
#'    is NULL, then \code{tx2geneDF} will be created using
#'    \code{makeTx2geneFromGtf}.
#' @param txdb a \code{TxDb} object as defined by Bioconductor
#'    package \code{GenomicFeatures}. Note that when reloading an R session,
#'    these objects usually become defunct, since they refer internally to
#'    a sqlite database stored in a temporary file. Therefore, when saving
#'    and reloading an R session, it is recommended to use the relevant
#'    functions \code{AnnotationDbi::saveDb} and
#'    \code{AnnotationDbi::loadDb}.
#' @param threeUtrGRL GRangesList object of 3'UTR exons by transcript, whose
#'    names are available in the \code{tx2geneDF} data.frame with colname
#'    \code{transcript_id}.
#' @param detectedTx optional character vector of transcripts detected,
#'    whose values match the transcript names in \code{gtf}, \code{txdb},
#'    \code{threeUtrGRL}, and/or \code{tx2geneDF} as relevant to the
#'    supplied input data.
#' @param tx2geneDF optional data.frame with colnames \code{gene_id} and
#'    \code{transcript_id} which are used to relate transcripts to genes.
#'    If either \code{gtf} or \code{txdb} are supplied, this data.frame
#'    can be determined from that data.
#' @param iMatrix optional numeric matrix of transcript rows, and sample
#'    columns, containing expression abundance data. When supplied, a
#'    matrix will be created \code{iMatrixByALE} where the transcript
#'    abundance values have been aggregated per ALE as defined.
#' @param aleMethod character value describing the method used to define
#'    ALE regions.
#' @param verbose logical whether to print verbose output during processing.
#'
#' @return aleGRL GRangesList object containing the ALE ranges after
#'    processing. Compare to transcript models to see if 3'UTR regions have
#'    been properly handled.
#' @return multiALEgenes vector of gene_name values which contain multiple
#'    ALEs after processing.
#' @return tx2ale vector containing ALE_name values, named by transcript_id,
#'    used to aggregate per-transcript data into per-ALE data.
#' @return iMatrixByALE optional data matrix created only if iMatrix was
#'    supplied as input. The data matrix rows are ALE_name values after
#'    aggregating transcripts into ALEs. Note that only multi-ALE genes
#'    are included, all single-ALE genes are excluded.
#'
#' @family jam ALE-specific RNA-seq functions
#'
#' @export
tx2ale <- function
(gtf=NULL,
 txdb=NULL,
 threeUtrGRL=NULL,
 detectedTx=NULL,
 tx2geneDF=NULL,
 iMatrix=NULL,
 aleMethod=c("first", "range"),
 verbose=FALSE,
 ...)
{
   ## Purpose is to encapsulate the logic of using Gencode GTF gene transcript
   ## models to infer alternative last exons (ALEs).
   ##
   ## Basic logic flow:
   ## - determine 3'UTR for each transcript
   ## - convert 3'UTR from multiple exons to a spanning range
   ## - melt 3'UTR ranges by gene, the separate 3'UTR ranges become ALEs
   ## - assign transcripts to ALEs, many-to-one relationship
   ## - require multiple ALEs per gene
   ## - aggregate transcript count matrix by ALE
   ## - test ALE matrix for differential ALEs using limma::diffSplice()
   ##
   ## Input required:
   ## - either gtf (Gencode GTF file), or txdb, or threeUtrGRL. The gtf or
   ##   txdb are used to create a GRangesList of 3'UTR elements, with
   ##   GenomicFeatures::threeUTRsByTranscript(). However, the GRangesList
   ##   can be supplied directly, bypassing the use of gtf or txdb.
   ## - detectedTx, vector of transcript_id values matching transcript names
   ##   from the gtf or txdb data.
   ##   This data is very helpful in removing
   ##   noisy transcript entries which are not expressed, but whose annotations
   ##   sometimes cause problems when trying to recognize ALEs. Specifically,
   ##   there are transcripts with alternate translation stop codon early in
   ##   the transcript, which cases the entire remaining set of exons to be
   ##   marked as 3'UTR. Some transcripts thus cover 80% or more of the exons
   ##   with 3'UTR, which during the step that melts 3'UTR ranges results
   ##   in one large 3'UTR covering almost the entire gene. Thus, if there
   ##   were separate and real 3'UTR ranges downstream, they will have been
   ##   melted into one large range, rendering the gene invisible to this ALE
   ##   analysis.
   ## - tx2geneDF is optionally a data.frame with gene_name and transcript_id
   ##   columns, but will be created from either the txdb or the gtf file
   ##   if not supplied.
   ## - iMatrix optional matrix containing log2 expression values per
   ##   transcript, whose rownames match transcript_id values from the gtf
   ##   or txdb data. If supplied, it will be used to aggregate the expression
   ##   values by ALE, creating iMatrixByALE. If not supplied, nothing will
   ##   be produced.
   ##
   ## Output:
   ## - aleGRL GRangesList object containing the ALE ranges after
   ##   processing. Compare to transcript models to see if 3'UTR regions have
   ##   been properly handled.
   ## - multiALEgenes vector of gene_name values which contain multiple ALEs
   ##   after processing.
   ## - tx2ale vector containing ALE_name values, named by transcript_id,
   ##   used to aggregate per-transcript data into per-ALE data.
   ## - iMatrixAle data matrix created only if iMatrix was supplied as input.
   ##   The data matrix rows are ALE_name values after aggregating transcripts
   ##   into ALEs. Note that only multi-ALE genes are included.
   ##
   ## aleMethod "range" will convert 3'UTR to the full range, and will
   ##    combine any other 3'UTRs from the same gene which overlap any part
   ##    of that range.
   ## aleMethod "first" will use only the first exon from a 3'UTR, which has
   ##    the effect of only combining 3'UTRs which begin at or near the same
   ##    point in the transcript. It therefore does not combine transcripts
   ##    that overlap some downstream 3'UTR exons, which can happen sometimes
   ##    with transcripts that contain alternate annotated stop codons far
   ##    upstream, turning all downstream CDS exons into 3'UTR, and thus
   ##    combining all transcripts into one huge ALE.
   ##
   ## Other R custom function dependencies:
   ## - annotateGRLfromGRL(), annotateGRfromGR()
   ## - assignGRLexonNames()
   if (suppressPackageStartupMessages(!require(jamba))) {
      stop("gencode2ale() requires the jamba package, installable with: devtools::install_github('jmw86069/jamba')");
   }
   if (suppressPackageStartupMessages(!require(GenomicRanges))) {
      stop("gencode2ale() requires the GenomicRanges package.");
   }
   if (suppressPackageStartupMessages(!require(GenomicFeatures))) {
      stop("gencode2ale() requires the GenomicFeatures package, GenomicFeatures::threeUTRsByTranscript()");
   }
   retVals <- list();
   aleMethod <- match.arg(aleMethod);

   ####################################################
   ## Determine appropriate input, ultimately create 3'UTR GRL
   if (length(threeUtrGRL) > 0 &&
         igrepHas("GRangesList", class(threeUtrGRL))) {
      ## use threeUtrGRL as-is
      if (verbose) {
         printDebug("gencode2ale(): ",
            "Using threeUtrGRL as-is.");
      }
   } else if (length(txdb) > 0) {
      ## use txdb
      if (verbose) {
         printDebug("gencode2ale(): ",
            "Creating threeUtrGRL from txdb.");
      }
      threeUtrGRL <- GenomicFeatures::threeUTRsByTranscript(txdb,
         use.names=TRUE);
      retVals$threeUtrGRL <- threeUtrGRL;
   } else if (length(gtf) > 0) {
      ## use GTF file
      if (verbose) {
         printDebug("gencode2ale(): ",
            "Creating txdb from gtf.");
      }
      #{startTimer();
      txdb <- makeTxDbFromGFF(gtf);
      retVals$txdb <- txdb;
      #stopTimer();}
      ## Extract 3'UTR
      if (verbose) {
         printDebug("gencode2ale(): ",
            "Creating txdb from gtf.");
      }
      threeUtrGRL <- GenomicFeatures::threeUTRsByTranscript(txdb,
         use.names=TRUE);
      retVals$threeUtrGRL <- threeUtrGRL;
   }

   ####################################################
   ## create tx2geneDF
   if (length(tx2geneDF) == 0) {
      if (length(gtf) == 0) {
         stop("gencode2ale() requires either a GTF file or tx2geneDF data.frame");
      }
      if (verbose) {
         printDebug("gencode2ale(): ",
            "Creating tx2geneDF from gtf file.");
      }
      tx2geneDF <- makeTx2geneFromGtf(GTF=gtf);
      rownames(tx2geneDF) <- tx2geneDF[,"transcript_id"];
      retVals$tx2geneDF <- tx2geneDF;
   }

   ####################################################
   ## create detectedTx
   if (length(detectedTx) == 0) {
      detectedTx <- names(threeUtrGRL);
      if (verbose) {
         printDebug("gencode2ale(): ",
            "No detectedTx supplied, keeping all ",
            formatInt(length(detectedTx)),
            " transcripts.");
      }
   } else {
      detectedTx <- intersect(detectedTx,
         names(threeUtrGRL));
      if (verbose) {
         printDebug("gencode2ale(): ",
            "Keeping ",
            formatInt(length(detectedTx)),
            " transcripts found in detectedTx.");
      }
   }

   ####################################################
   ## Apply filter to threeUtrGRL transcripts
   threeUtrGRLdet <- threeUtrGRL[names(threeUtrGRL) %in% detectedTx];

   ####################################################
   ## Convert 3'UTR exons to ranges
   if (aleMethod %in% "range") {
      if (verbose) {
         printDebug("gencode2ale(): ",
            "Using ",
            "full 3'UTR range",
            " to determine ALEs.");
      }
      threeUtrGRLdetRange <- range(threeUtrGRLdet);
   } else if (aleMethod %in% "first") {
      if (verbose) {
         printDebug("gencode2ale(): ",
            "Using ",
            "first 3'UTR stranded exon",
            " to determine ALEs.");
      }
      threeUtrGRLdetRange <- getFirstStrandedFromGRL(threeUtrGRLdet,
         method="direct",
         verbose=verbose);
   }
   ## add transcript_id annotation
   GenomicRanges::values(threeUtrGRLdetRange@unlistData)[,"transcript_id"] <- rep(
      names(threeUtrGRLdetRange),
      S4Vectors::elementNROWS(threeUtrGRLdetRange));
   ## add gene_name annotation
   GenomicRanges::values(threeUtrGRLdetRange@unlistData)[,"gene_name"] <- tx2geneDF[
      match(GenomicRanges::values(threeUtrGRLdetRange@unlistData)[,"transcript_id"],
         tx2geneDF[,"transcript_id"]),"gene_name"];

   ####################################################
   ## Re-aggregate by gene
   if (verbose) {
      printDebug("gencode2ale(): ",
         "Splitting ranges by gene.");
   }
   threeUtrGRLdetGeneGRL <- GenomicRanges::split(
      threeUtrGRLdetRange@unlistData,
      f=GenomicRanges::values(threeUtrGRLdetRange@unlistData)[,"gene_name"]);
   ## Reduce (melt) 3'UTR ranges, combining overlapping ranges per gene
   if (verbose) {
      printDebug("gencode2ale(): ",
         "Reducing ranges.");
   }
   threeUtrGRLdetGeneGRLred <- GenomicRanges::reduce(threeUtrGRLdetGeneGRL);

   ####################################################
   ## Annotate transcripts to the reduced 3'UTR ranges
   if (verbose) {
      printDebug("gencode2ale(): ",
         "Annotating reduced ranges with transcript_id per gene.");
   }
   threeUtrGRLdetGeneGRLred2 <- annotateGRLfromGRL(
      threeUtrGRLdetGeneGRLred,
      threeUtrGRLdetGeneGRL,
      verbose=verbose);

   ####################################################
   ## Assign ALE numbers in stranded order for each gene
   if (verbose) {
      printDebug("gencode2ale(): ",
         "Assigning stranded numbers to the ranges.");
      printDebug("head(threeUtrGRLdetGeneGRLred2):");
      print(head(threeUtrGRLdetGeneGRLred2));
   }
   threeUtrGRLdetGeneGRLred2 <- assignGRLexonNames(
      exonNameColname="ALE_name",
      geneSymbolColname="gene_name",
      threeUtrGRLdetGeneGRLred2,
      suffix="_ale",
      verbose=verbose);
   if (verbose) {
      printDebug("gencode2ale(): ",
         "Assigned stranded numbers to the ranges.");
   }
   names(threeUtrGRLdetGeneGRLred2@unlistData) <-
      GenomicRanges::values(threeUtrGRLdetGeneGRLred2@unlistData)[,"ALE_name"];
   GenomicRanges::values(threeUtrGRLdetGeneGRLred2@unlistData)[,"score"] <- 0.5;
   retVals$aleGRL <- threeUtrGRLdetGeneGRLred2;

   ####################################################
   ## Subset ALE containing 2 or more ALEs per gene
   GencodeALEmin2 <- threeUtrGRLdetGeneGRLred2[
      S4Vectors::elementNROWS(threeUtrGRLdetGeneGRLred2) > 1];
   if (verbose) {
      printDebug("gencode2ale(): ",
         "Filtered ",
         formatInt(length(threeUtrGRLdetGeneGRLred2)),
         " genes to ",
         formatInt(length(GencodeALEmin2)),
         " genes having multiple ranges.");
   }
   ## vector of genes containing multiple ALEs
   GencodeALEmin2genes <- jamba::mixedSort(unique(GenomicRanges::values(GencodeALEmin2@unlistData)[,"gene_name"]));
   retVals$multiALEgenes <- GencodeALEmin2genes;

   ####################################################
   ## Create a tx-to-ALE xref using multi-ALE genes
   if (verbose) {
      printDebug("gencode2ale(): ",
         "Creating transcript-to-ALE xref for multi-range genes.");
   }
   ale2txL <- strsplit(jamba::nameVector(
      GenomicRanges::values(subset(threeUtrGRLdetGeneGRLred2@unlistData,
         gene_name %in% GencodeALEmin2genes))[,c("transcript_id","ALE_name")]),
      ",");
   tx2ale <- jamba::nameVector(list2df(ale2txL)[,c("item","value")]);
   retVals$tx2ale <- tx2ale;

   ####################################################
   ## Aggregate expression values by ALE
   if (length(iMatrix) > 0) {
      ## Note: data is exponentiated before taking the sum,
      ## then is log2-transformed
      if (verbose) {
         printDebug("gencode2ale(): ",
            "Aggregating transcript abundances by ranges.");
      }
      if (all(names(tx2ale) %in% rownames(iMatrix))) {
         iMatrixAle <- log2(1+
               shrinkMatrix(
                  2^(iMatrix[names(tx2ale),,drop=FALSE])-1,
                  groupBy=tx2ale,
                  shrinkFunc=function(x){sum(x, na.rm=TRUE)},
                  returnClass="matrix"));
      } else {
         printDebug("gencode2ale(): ",
            "Warning: ",
            "names(tx2ale) were not present in all rownames(iMatrix).",
            fgText=c("orange","red"));
         iMatrixUse <- (2^(iMatrix[match(names(tx2ale), rownames(iMatrix)),,drop=FALSE])-1);
         rownames(iMatrixUse) <- names(tx2ale);
         iMatrixAle <- log2(1+shrinkMatrix(iMatrixUse,
            groupBy=tx2ale,
            shrinkFunc=function(x){sum(x, na.rm=TRUE)},
            returnClass="matrix"));
         #iMatrixAle <- NULL;
      }
   } else {
      iMatrixAle <- NULL;
   }
   retVals$iMatrixAle <- iMatrixAle;
   return(retVals);
}

#' Define detected transcripts
#'
#' Define detected transcripts
#'
#' This function aims to combine evidence from RNA-seq sequence
#' read counts (or pseudocounts from a kmer tool such as Salmon or
#' Kallisto), along with alternative TPM quantitation, to determine
#' the observed "detected" transcript space for a given experiment.
#'
#' Each input data matrix is assumed to be appropriately log-transformed,
#' typically using `log2(1+x)`. If any value is `>= 50` then the data
#' matrix will be log2-transformed using `log2(1+x)`.
#'
#' The criteria must be met in at least one sample group, but all
#' criteria must be met in the same sample group for an isoform to
#' be considered "detected".
#'
#' In our experience the use of TPM values appears more robust and
#' is conceptually the best approach for comparing the relative
#' quantity of one transcript isoform to another. Our reasoning is
#' that TPM is intended to be roughly a molar quantity of transcript
#' molecules, independent of the transcript length, and the potential
#' for overlapping regions between isoforms. We also recommend the use
#' of a kmer quantitation method, such as Salmon or Kallisto, which
#' estimates isoform abundances not by specific read counts, but by
#' quantifying kmers unique to particular isoforms for a given gene.
#'
#' In all cases, the thresholds for detection can be modified, however
#' from our experiences thus far the default values perform reasonably
#' well at identifying expressed isoforms, while filtering out isoforms
#' that we considered to be spuriously expressed.
#'
#' There are three default requirements for a transcript to be considered
#' "detected".
#'
#' 1. An isoform must be expressed at least 10% of the max isoform for
#' a given gene, using TPM values.
#' 2. An isoform must have at least log2(32) pseudocounts to be considered
#' detected, based upon our view of Salmon pseudocount data using MA-plots.
#' 3. An isoform must have at least log2(2) TPM units to be considered
#' detected, based upon our view of Salmon TPM values using MA-plots.
#'
#' Each experiment is likely to be different in terms of total
#' sequenced reads, quality of read alignment or quantitation to the
#' transcriptome, etc. We suggest observing MA-plots for the counts and
#' TPM values, for the point at which the signal substantially increases
#' from baseline zero. We also plotted the TPM versus count per sample,
#' noted the point at which the two signals began to correlate. These
#' observations along with careful review of numerous gene model
#' transcript isoforms supported our selection of these criteria.
#'
#' Lastly, the requirement for 10 percent of max isoform expression
#' was motivated by observing highly expressed genes, which sometimes had
#' alternative isoforms with extremely low abundance compared to the
#' most abundant isoform, but which was notably higher than the minimum
#' for detection. For example Gapdh expression above 100,000 pseudocounts,
#' may have an isoform with 120 pseudocounts. When we reviewed the sequence
#' coverage, we could find no compelling evidence to support the minor
#' isoform, and theorized that the pseudocounts arose from the stochastic
#' nature of rebalancing relative expression among isoforms.
#'
#' Note the argument `zeroAsNA=TRUE`, which by default treats any
#' expression value of zero (or less than zero) as `NA`,
#' thus removing them from group
#' mean calculations. When `iMatrixTxGrp` and `iMatrixTxTPMGrp`
#' are not supplied, this option is helpful in calculating a more
#' appropriate group mean expression value, notably when a
#' value of zero represents absence of data. Any group mean that is
#' `NA` as a result is converted to zero for the purpose of applying
#' filters.
#'
#' @return
#' List with the following elements:
#' \describe{
#'    \item{txExprGrpTx}{Numeric matrix representing the expression
#'       counts per transcript, grouped by `"gene_name"`.}
#'    \item{txPctMaxGrpAll}{Numeric matrix representing the
#'       percent expression of each transcript isoform per gene, as
#'       compared to the highest expression of isoforms for that gene.}
#'    \item{txExprGrpAll}{Numeric matrix of sample group counts,
#'       exponentiated and rounded to integer values.}
#'    \item{txTPMExprGrpAll}{Numeric matrix of sample group TPM values,
#'       exponentiated and rounded to integer values.}
#'    \item{txFilterM}{Numeric matrix indicating whether each isoform met
#'       the criteria to be considered detected. The criteria must be
#'       met in the same group for an isoform to be considered detected.}
#'    \item{detectedTx}{Character vector of transcripts, as defined by
#'       the `rownames(iMatrixTx)`.}
#' }
#'
#' @param iMatrixTx numeric matrix of read counts (or pseudocounts) with
#'    transcript rows and sample columns. This data is assumed to be
#'    log2-transformed, and if any value is higher than 50, it will
#'    be log2-transformed with `log2(x+1)`.
#' @param iMatrixTxGrp numeric matrix of read counts averaged by sample
#'    group. If this matrix is not provided, it will be calculated
#'    from `iMatrixTx`
#'    using `jamba::rowGroupMeans()` and the `groups` parameter.
#'    This data is assumed to be
#'    log2-transformed, and if any value is higher than 50, it will
#'    be log2-transformed with `log2(x+1)`.
#' @param iMatrixTxTPM numeric matrix of TPM values, with sample columns
#'    and transcript rows. Note that if this parameter is not supplied,
#'    the counts in `iMatrixTx` will be used to determine the percent
#'    max isoform expression.
#' @param iMatrixTxTPMGrp numeric matrix of TPM values averaged by sample
#'    group. If this matrix is not provided, it will be calculated
#'    from `iMatrixTxTPM`
#'    using `jamba::rowGroupMeans()` and the `groups` parameter.
#' @param groups vector of group labels, either as character vector
#'    or factor. It should be named by `colnames(iMatrixTx)`.
#' @param tx2geneDF data.frame with colnames including
#'    `c("transcript_id","gene_name")`, where the values in the
#'    `"transcript_id"` column must match the `rownames(iMatrixTx)`.
#' @param cutoffTxPctMax numeric value scaled from 0 to 100
#'    indicating the percentage of
#'    the maximum isoform expression per gene, for an alternate isoform
#'    to be considered for detection.
#' @param cutoffTxExpr numeric value indicating the minimum group mean
#'    counts in `iMatrixTxGrp` for a transcript to be considered for
#'    detection.
#' @param cutoffTxTPMExpr numeric value indicating the minimum group mean
#'    TPM in `iMatrixTxTPMGrp` for a transcript to be considered for
#'    detection.
#' @param txColname,geneColname the `colnames(tx2geneDF)` representing
#'    the `rownames(iMatrixTx)` matched by `tx2geneDF[,txColname]`,
#'    and the associated genes given by `tx2geneDF[,geneColname]`.
#'    Note that `detectedTx` must also contain values in
#'    `rownames(iMatrixTx)` and `tx2geneDF[,txColname]`.
#' @param zeroAsNA logical indicating whether values of zero
#'    (or less than zero) should
#'    be treated as `NA` values, thus removing them from mean
#'    calculations. This argument is only relevant when `iMatrixTxGrp`
#'    and `iMatrixTxTPMGrp` are not supplied.
#'    Argument `zeroAsNA=TRUE` is recommended when using kmer
#'    quantitation tools such as Salmon or Kallisto, which can
#'    sometimes allocate all expression to one or another transcript
#'    isoform when two isoforms are nearly identical. Also use
#'    `TRUE` when a value of zero represents the absense of data.
#' @param useMedian logical indicating whether to use group median
#'    values instead of group mean values.
#' @param verbose logical indicating whether to print verbose output.
#'
#' @family jam RNA-seq functions
#'
#' @export
defineDetectedTx <- function
(iMatrixTx=NULL,
 iMatrixTxGrp=NULL,
 iMatrixTxTPM=NULL,
 iMatrixTxTPMGrp=NULL,
 groups=NULL,
 tx2geneDF=NULL,
 cutoffTxPctMax=10,
 cutoffTxExpr=5,
 cutoffTxTPMExpr=2,
 txColname="transcript_id",
 geneColname="gene_name",
 zeroAsNA=TRUE,
 useMedian=FALSE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to define detected transcripts given some thresholds
   ## for minimum observed counts, and fraction relative expression
   ## of each isoform per gene.
   ##
   ## cutoffTxExpr is a threshold of counts, after exponentiating iMatrixTx
   ## via (2^iMatrixTx - 1).
   ##
   ## iMatrixTxTPM is a data matrix containing TPM values. When supplied,
   ## it is used during the step to filter isoforms by percent max expression.
   ## We noticed that when isoforms were substantially shorter length than
   ## the isoform with max counts, it had much lower TPM value proportionally,
   ## as expected. As a result, the suggested method involved using counts
   ## to filter for minimum counts above noise, and TPM for percent
   ## max isoform expression.
   ##
   ## TODO:
   ## - Detect whether data needs log2 transformation, and apply it as
   ##   necessary.

   retVals <- list();

   ## group mean values for counts
   if (length(iMatrixTxGrp) == 0) {
      if (length(iMatrixTx) == 0) {
         stop("defineDetectedTx() requires either iMatrixTx or iMatrixTxGrp.");
      }
      if (length(groups) == 0) {
         stop("defineDetectedTx() requires groups, named by colnames(iMatrixTx).");
      }
      ## Optionally apply log2(1+x) transformation
      if (max(iMatrixTx, na.rm=TRUE) >= 50) {
         if (verbose) {
            printDebug("defineDetectedTx(): ",
               "Applying log2(1+x) transform to:",
               "iMatrixTx");
         }
         iMatrixTx <- jamba::noiseFloor(log2(1+iMatrixTx),
            minimum=0,
            adjustNA=TRUE);
      }
      ## Optionally convert zero to NA
      if (length(zeroAsNA) > 0 && zeroAsNA && any(iMatrixTx <= 0)) {
         if (verbose) {
            printDebug("defineDetectedTx(): ",
               "Converting Tx zero to NA prior to group mean calculations");
         }
         iMatrixTx[iMatrixTx <= 0] <- NA;
      }
      ## Calculate group mean values
      if (verbose) {
         printDebug("defineDetectedTx(): ",
            "Calculating iMatrixTxGrp.");
      }
      iMatrixTxGrp <- jamba::rowGroupMeans(iMatrixTx,
         useMedian=useMedian,
         groups=groups);
   }
   ## Process group mean counts
   ## Convert any remaining NA to zero
   if (any(is.na(iMatrixTxGrp))) {
      iMatrixTxGrp[is.na(iMatrixTxGrp)] <- 0;
   }
   ## Optionally apply log2(1+x) transformation
   if (max(iMatrixTxGrp, na.rm=TRUE) >= 50) {
      if (verbose) {
         printDebug("defineDetectedTx(): ",
            "Applying log2(1+x) transform to:",
            "iMatrixTxGrp");
      }
      iMatrixTxGrp <- jamba::noiseFloor(log2(1+iMatrixTxGrp),
         minimum=0,
         adjustNA=TRUE);
   }

   ## group mean values for TPM
   if (length(iMatrixTxTPMGrp) == 0) {
      if (length(iMatrixTxTPM) > 0) {
         if (length(groups) == 0) {
            stop("defineDetectedTx() requires groups, named by colnames(iMatrixTxTPM).");
         }
         ## Optionally apply log2(1+x) transformation
         if (max(iMatrixTxTPM, na.rm=TRUE) >= 50) {
            if (verbose) {
               printDebug("defineDetectedTx(): ",
                  "Applying log2(1+x) transform to:",
                  "iMatrixTxTPM");
            }
            iMatrixTxTPM <- jamba::noiseFloor(log2(1+iMatrixTxTPM),
               minimum=0,
               adjustNA=TRUE);
         }
         ## Optionally convert zero to NA
         if (length(zeroAsNA) > 0 && zeroAsNA && any(iMatrixTxTPM <= 0)) {
            if (verbose) {
               printDebug("defineDetectedTx(): ",
                  "Converting TxTPM zero to NA prior to group mean calculations");
            }
            iMatrixTxTPM[iMatrixTxTPM <= 0] <- NA;
         }
         ## Calculate group mean values
         if (verbose) {
            printDebug("defineDetectedTx(): ",
               "Calculating iMatrixTxTPMGrp.");
         }
         iMatrixTxTPMGrp <- rowGroupMeans(iMatrixTxTPM,
            useMedian=useMedian,
            groups=groups);
      }
   }
   ## Process group mean TPM
   ## Convert any remaining NA to zero
   if (any(is.na(iMatrixTxTPMGrp))) {
      iMatrixTxTPMGrp[is.na(iMatrixTxTPMGrp)] <- 0;
   }
   ## Optionally apply log2(1+x) transformation
   if (max(iMatrixTxTPMGrp, na.rm=TRUE) >= 50) {
      if (verbose) {
         printDebug("defineDetectedTx(): ",
            "Applying log2(1+x) transform to:",
            "iMatrixTxTPMGrp");
      }
      iMatrixTxTPMGrp <- jamba::noiseFloor(log2(1+iMatrixTxTPMGrp),
         minimum=0,
         adjustNA=TRUE);
   }

   ######################################################################
   ## matrix associating gene to transcript_id, mainly useful since it
   ## returns transcripts in the same order as each matrix below, which
   ## otherwise has rownames based upon gene and not transcript.
   iRows <- match(rownames(iMatrixTxGrp), tx2geneDF[,txColname]);
   if (verbose) {
      printDebug("defineDetectedTx(): ",
         "shrinkMatrix Tx names.");
   }
   txExprGrpTx <- shrinkMatrix(rownames(iMatrixTxGrp),
      groupBy=tx2geneDF[iRows,geneColname],
      shrinkFunc=c,
      returnClass="matrix",
      verbose=verbose);
   retVals$txExprGrpTx <- txExprGrpTx;


   ######################################################################
   ## Percent max isoform expression per isoform per gene
   if (length(iMatrixTxTPMGrp) > 0) {
      if (verbose) {
         printDebug("defineDetectedTx(): ",
            "Calculating percent max isoform expression by TPM.");
      }
      txPctMaxGrpAll <- shrinkMatrix(2^iMatrixTxTPMGrp-1,
         groupBy=tx2geneDF[iRows,geneColname],
         shrinkFunc=function(i){
            round(i/max(c(1, max(i)))*100)
         },
         returnClass="matrix");
      attr(txPctMaxGrpAll, "TxMeasurement") <- "TPM";
      rownames(txPctMaxGrpAll) <- txExprGrpTx[,1];
      #retVals$iMatrixTxTPMGrp <- iMatrixTxTPMGrp;
   } else {
      if (verbose) {
         printDebug("defineDetectedTx(): ",
            "Calculating percent max isoform expression by counts.");
      }
      txPctMaxGrpAll <- shrinkMatrix((2^iMatrixTxGrp-1),
         groupBy=tx2geneDF[iRows,geneColname],
         shrinkFunc=function(i){
            round(i/max(c(1, max(i, na.rm=TRUE)), na.rm=TRUE)*100);
         },
         returnClass="matrix");
      attr(txPctMaxGrpAll, "TxMeasurement") <- "counts";
      rownames(txPctMaxGrpAll) <- txExprGrpTx[,1];
      #retVals$iMatrixTxGrp <- iMatrixTxGrp;
   }
   retVals$txPctMaxGrpAll <- txPctMaxGrpAll;


   ######################################################################
   ## Exponentiated expression counts per isoform per gene
   if (verbose) {
      printDebug("defineDetectedTx(): ",
         "Creating exponentiated expression counts, rounded to integers.");
   }
   txExprGrpAll <- shrinkMatrix(2^iMatrixTxGrp-1,
      groupBy=tx2geneDF[iRows,geneColname],
      shrinkFunc=function(i){
         round(i)
      },
      returnClass="matrix");
   rownames(txExprGrpAll) <- txExprGrpTx[,1];
   retVals$txExprGrpAll <- txExprGrpAll;


   ######################################################################
   ## Exponentiated expression TPM per isoform per gene
   if (length(iMatrixTxTPMGrp) > 0) {
      if (verbose) {
         printDebug("defineDetectedTx(): ",
            "Creating exponentiated expression TPM.");
      }
      txTPMExprGrpAll <- shrinkMatrix(2^iMatrixTxTPMGrp-1,
         groupBy=tx2geneDF[iRows,geneColname],
         shrinkFunc=function(i){
            round(i)
         },
         returnClass="matrix");
      rownames(txTPMExprGrpAll) <- txExprGrpTx[,1];
      retVals$txTPMExprGrpAll <- txTPMExprGrpAll;
   }


   ######################################################################
   ## Apply the filters for detected, expressed Tx
   ## Tx must be at least 10% of the max TX abundance
   ## Tx must have expression at least 2^5 (32 counts)
   ## These conditions must occur in the same sample group,
   ## but if it occurs in any sample group the Tx is "detected"
   if (length(iMatrixTxTPMGrp) > 0) {
      if (verbose) {
         printDebug("defineDetectedTx(): ",
            "Applying filters for percentMaxIsoformTPM, Expr counts, Expr TPM.");
      }
      txFilterM <- (
         txPctMaxGrpAll >= cutoffTxPctMax &
            txExprGrpAll >= cutoffTxExpr &
            txTPMExprGrpAll >= cutoffTxTPMExpr
      )*1;
   } else {
      if (verbose) {
         printDebug("defineDetectedTx(): ",
            "Applying filters for percentMaxIsoform, Expr counts.");
      }
      txFilterM <- (
         txPctMaxGrpAll >= cutoffTxPctMax &
            txExprGrpAll >= cutoffTxExpr
      )*1;
   }
   retVals$txFilterM <- txFilterM;
   whichTx <- which(rowMaxs(txFilterM, na.rm=TRUE) > 0);
   detectedTx <- txExprGrpTx[whichTx,1];
   retVals$detectedTx <- detectedTx;
   return(retVals);
}

#' Shrink numeric matrix by groups of rows
#'
#' Shrink numeric matrix by groups of rows
#'
#' This function is mainly a wrapper around the very efficient
#' functions in the `data.table` package, with the ability to
#' provide a custom function to shrink row values.
#'
#' @param x numeric matrix
#' @param groupBy vector of group labels, whose length equals `nrow(x)`.
#'    These values will become rownames in the output data.
#' @param shrinkFunc function that takes vector input and returns
#'    vector output. The vector class can be checked, in order to
#'    call a function on numeric or character data separately, as
#'    needed.
#' @param returnClass character string indicating the return data type,
#'    `"data.frame"` returns a `data.frame` whose first column contains
#'    entries from `groups`; `"matrix"` returns a numeric matrix whose
#'    rownames are entries from `groups`.
#' @param verbose logical indicating whether to print verbose output.
#'
#' @import data.table
#'
#' @family jam matrix functions
#'
#' @export
shrinkMatrix <- function
(x,
 groupBy,
 shrinkFunc=function(x){.Internal(mean(x))},
 returnClass=c("data.frame", "matrix"),
 verbose=FALSE,
 ...)
{
   ## Purpose is to wrapper the data.table() fast function to apply numerical
   ## operations to grouped rows of a matrix; using matrix because this method
   ## is written to use only numeric/integer matrices and not text values.
   ##
   ## However, in principle this method will work fine as long as all columns
   ## are to be treated the same way, thus class can be anything valid for the
   ## function.
   ##
   ## Note: using .Internal(mean(x)) is 5x faster than using mean(x)
   ##
   ## Timings which show data.table and sqldf are fastest at apply() functions
   ## on groups:
   ## http://zvfak.blogspot.com/2011/03/applying-functions-on-groups-sqldf-plyr.html
   ##
   ## Another alternative is package sqldf, example syntax:
   ## n <- 100000;
   ## grp1 <- sample(1:750, n, replace=TRUE);
   ## grp2 <- sample(1:750, n, replace=TRUE);
   ## d <- data.frame(x=rnorm(n), y=rnorm(n), grp1=grp1, grp2=grp2, n, replace=TRUE);
   ## rsqldf <- system.time(sqldf("select grp1, grp2, avg(x), avg(y) from dgroup by grp1, grp2"));
   ##
   ## DT <- data.table(d);
   ## rdataT <- system.time(DT[,list(.Internal(mean(x)), .Internal(mean(y))), by=list(grp1,grp2)]);
   if (!suppressPackageStartupMessages(require(data.table))) {
      stop("This method requires the data.table package.");
   }
   returnClass <- match.arg(returnClass);

   ## Create DT object
   if (verbose) {
      t1 <- Sys.time();
   }
   DT <- data.table(
      data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         x,
         groupBy=groupBy),
      key="groupBy");
   if (verbose) {
      t2 <- Sys.time();
   }

   ## Operate on the DT object
   byDT <- DT[,lapply(.SD, shrinkFunc),
      by="groupBy"];
   if (verbose) {
      t3 <- Sys.time();
   }

   if (verbose) {
      jamba::printDebug("shrinkMatrix(): ",
         "Duration for data.table DT creation: ",
         format(t2-t1));
      jamba::printDebug("shrinkMatrix(): ",
         "Duration for data.table shrinkMatrix: ",
         format(t3-t2));
   }
   retData <- as(byDT, "data.frame");
   if (returnClass %in% "matrix") {
      retData <- matrix(ncol=ncol(retData)-1,
         data=(as.matrix(retData[,-1,drop=FALSE])),
         dimnames=list(retData[,"groupBy"], colnames(retData)[-1]));
   }
   return(retData);
}

#' Make codon usage data.frame
#'
#' Make codon usage data.frame
#'
#' This function imports a text file containing codon usage, and
#' creates a data.frame with columns containing the `codon`,
#' `fraction`, and `count`. The colname `fraction` represents
#' the fraction observed counts per codon
#' compared to the codon with the highest counts, for each
#' amino acid. For example, an amino acid with only one codon
#' would be `c(1.0)`. An amino acid with two codons, with observed
#' counts `c(100, 200)`, would have fraction values `c(0.5, 1.0)`.
#'
#' The format is derived from published files, one example is
#' included in the package `extdata/Mouse_codon_usage.txt`,
#' which can be accessed using the command:
#'
#' `system.file("extdata", "Mouse_codon_usage.txt", package="farrisdata")`
#'
#' @param file path to the codon usage file
#'
#' @return data.frame containing codon usage data, with colnames
#' `codon`, `frequency`, and `count`.
#'
#' @examples
#' codonFile <- system.file("extdata", "Mouse_codon_usage.txt", package="splicejam");
#' codonDF <- codonUsage2df(codonFile);
#'
#' @family jam codon usage functions
#'
#' @export
codonUsage2df <- function
(file,
 ...)
{
   ## Purpose is to import a text codon usage file and return
   ## a data.frame.
   codonTxt <- read.table(file,
      comment.char="#",
      header=FALSE,
      sep="\t",
      stringsAsFactors=FALSE)[,1];

   ## Split into vector per codon
   codonTxt <- gsub("[ \t]+([ACGTUN]{3})[ \t]+",
      "!\\1 ",
      codonTxt);
   codonV <- unlist(strsplit(codonTxt, "[!]+"));
   codonDF <- data.frame(stringsAsFactors=FALSE,
      jamba::rbindList(strsplit(codonV, "[() ]+")));
   ncolDF <- seq_len(min(c(ncol(codonDF), 3)));
   colnames(codonDF)[ncolDF] <- c("codon", "frequency", "count")[ncolDF];

   ## Repair numeric columns, ensure it doesn't convert everything to NA
   for (i in tail(colnames(codonDF), -1)) {
      if (!all(is.na(as.numeric(codonDF[,i])))) {
         codonDF[,i] <- as.numeric(codonDF[,i]);
      }
   }
   ## Create rownames substituting 't' for 'u'
   rownames(codonDF) <- gsub("u", "t", tolower(codonDF[,1]));
   codonDF;
}

#' Convert DNA to 3-base codons
#'
#' Convert DNA to 3-base codons
#'
#' This function takes character string input, either as one large
#' string or vector of strings such as those read from a FASTA file,
#' and returns a vector of 3-base codons. The final codon must contain
#' all three character values otherwise it is removed.
#'
#' @return character vector where each element contains three characters,
#' e.g. `all(nchar(dna2codon(x)) == 3)`.
#'
#' @param x character vector assumed to contain DNA or RNA nucleotides
#' @param ... additional arguments are ignored.
#'
#' @examples
#' dna <- c("atgggattatag");
#' dna2codon(dna);
#'
#' # example adding one extra base, which does not make a full codon
#' dnaext1 <- c("atgggattataga");
#' dna2codon(dnaext1);
#'
#' @family jam codon usage functions
#'
#' @export
dna2codon <- function
(x,
 ...)
{
   ## Purpose is to convert a text string containing DNA sequence
   ## into 3-base codons.
   ##
   ## Split into a vector
   y <- unlist(strsplit(x, ""));
   h <- floor(length(y) / 3);
   ## The purpose of head() is to make sure we do not use the last
   ## codon unless it contains all three nucleotides. Otherwise
   ## matrix() will recycle text to fill the last row, creating a
   ## false codon.
   head(
      jamba::pasteByRow(
         matrix(ncol=3,
            byrow=TRUE,
            unlist(strsplit(x, ""))),
            #seqinr::s2c(x)),
      sep=""),
      h);
}

#' Calculate Codon Adaptation Index
#'
#' Calculate Codon Adaptation Index
#'
#' This function is equivalent to `seqinr::cai()` except that it
#' runs in vectorized form and is markedly more efficient.
#'
#' @param x character vector containing DNA or RNA sequence, intended
#'    to include amino acid codons in frame beginning with the first
#'    character. It is sent to `dna2codon()` which creates a character
#'    vector whose elements all contain three characters. The input
#'    `x` represents sequence data for one protein-coding transcript
#'    sequence.
#' @param codonDF data.frame containing amino acid codons per row,
#'    with colname `"codon"` containing three-character codon sequences,
#'    and colname `"multifreq"` containing the relative frequency
#'    versus max of codon use per amino acid, as produced by
#'    `codonUsage2df()`.
#' @param ... additional arguments are ignored.
#'
#' @return numeric vector with length `length(x)`, named with
#'    `names(x)` containing codon adaptation index (cai) values
#'    equivalent to those produced by `seqinr::cai()`.
#'
#' @family jam codon usage functions
#'
#' @export
jamCai <- function
(x,
 codonDF,
 ...)
{
   ## Purpose is to run the equivalent of cai() to test speed
   ##
   ## Todo:
   ## * validation checks on codonDF data.frame columns,
   ## potentially allow setting the colname to use.
   ##
   sapply(x, function(i){
      jamGeomean(rmNA(codonDF[dna2codon(i),"multiFreq"]));
   });
}

#' Modified geometric mean for positive and negative values
#'
#' Modified geometric mean for positive and negative values
#'
#' This function calculates a geometric mean using a formula which
#' tolerates positive and negative values. It also tolerates zeros
#' without resulting in zero. The intent is to provide
#' a mean summary value which closely models the classical
#' geometric mean, while retaining information present in
#' vectors that contain either zeros or negative values.
#'
#' The classical geometric mean is defined as
#' the exponentiated mean of log-transformed values. Said another
#' way, it is the `n`th root of the product of `n` numeric values.
#' This formula is analogous to geometric distance. The formula
#' does not allow negative values, however, and if any value is
#' zero the result is also zero.
#'
#' There are several proposed alternatives to address negative numbers,
#' or zeros. This function makes the following adjustments:
#'
#' * Add `1` to absolute value of the input, so the numeric sign
#' is not flipped during log transformation: `i <- log2(1+ abs(x))`
#' * Multiply that vector by the `sign(x)` to retain the original
#' positive and negative directionality: `j <- i * sign(x)`
#' * Take the mean: `k <- mean(j)`
#' * Exponentiate the absolute value: `m <- 2^(abs(k))`
#' * Multiply by the sign:
#' `n <- m * sign(k)`
#' * Subtract `1`: `o <- n - 1;`
#'
#' The properties of the calculation:
#'
#' * Symmetry around zero, for example `jamGeomean(x) = -jamGeomean(-x)`
#' * Results are slightly different than classical geometric mean values,
#' as a result of adding `1` prior to log transformation. The difference
#' is larger with increasing `range(x)` and is most noticeable when one
#' input value in `x` is zero.
#'
#' @return numeric value representing the modified geometric mean
#' of input values.
#'
#' @param x numeric vector which may contain positive and negative values.
#' @param na.rm logical indicating whether to ignore `NA` values.
#' @param ... additional parameters are ignored.
#'
#' @examples
#' x <- c(1, 10, 40);
#' jamGeomean(x);
#' # compare to classical geometric mean
#' geomean(x);
#'
#' # Positive and negative values should offset
#' x <- c(-20, 20);
#' jamGeomean(x);
#'
#' x <- c(-20,10,40);
#' jamGeomean(x);
#'
#' @family jam numeric functions
#'
#' @export
jamGeomean <- function
(x,
 na.rm=TRUE,
 ...)
{
   ## Purpose is to calculate geometric mean while allowing for
   ## positive and negative values
   x2 <- mean(log2(1+abs(x))*sign(x));
   sign(x2)*(2^abs(x2)-1);
}

#' Classical geometric mean
#'
#' Classical geometric mean
#'
#' This function calculates the classical geometric mean.
#'
#' The classical geometric mean is defined as
#' the exponentiated mean of log-transformed values. Said another
#' way, it is the `n`th root of the product of `n` numeric values.
#' This formula is analogous to geometric distance. The formula
#' does not allow negative values, however, and if any value is
#' zero the result is also zero.
#'
#' @return numeric value representing the geometric mean of input values
#'
#' @param x numeric vector containing only positive values
#' @param na.rm logical indicating whether to ignore `NA` values. Note that
#'    `NA` values are removed prior to log-transformation, to avoid negative
#'    numbers being dropped completely. To drop negative numbers, do so
#'    prior to calling `geomean()`.
#' @param offset numeric value added to input `x` prior to log
#'    transformation, intended only when values between 0 and 1 should be
#'    retained. Note that the offset makes the result slightly different
#'    than classical geometric mean.
#' @param ... additional parameters are ignored.
#'
#' @examples
#' x <- c(2, 4);
#' geomean(x);
#'
#' x <- c(-2, 2, 4);
#' geomean(x);
#'
#' x <- c(0, 4000, 200000);
#' geomean(x);
#'
#' @family jam numeric functions
#'
#' @export
geomean <- function
(x,
 na.rm=TRUE,
 naValue=NA,
 offset=0,
 ...)
{
   ## Purpose is to calculate the classical geometric mean
   if (na.rm && any(is.na(x))) {
      x <- rmNA(x,
         naValue=naValue);
   }
   2^mean(log2(x + offset)) - offset;
}

#' Summarize detected transcript results
#'
#' Summarize detected transcript results
#'
#' This function provides a simple summary of the results of
#' `defineDetectedTx()`, typically for a given gene of interest.
#'
#' @return list of `data.frame` objects, each containing one
#'    summary table of data used to support whether each transcript
#'    were called "detected."
#'
#' @param detectedTxL list output from `defineDetectedTx()`
#' @param Gene optional character vector of one or more genes
#'    of interest, used to find transcript_id entries in the
#'    `detectedTxL` input data. If not supplied, `Tx` is expected.
#' @param Tx optional character vector used to subset summary data,
#'    usually intended to keep rows in a specific order for a
#'    given set of transcripts.
#' @param ... additional arguments are ignored.
#'
#' @family jam RNA-seq functions
#'
#' @export
detectedTxInfo <- function
(detectedTxL,
 Gene=NULL,
 Tx=NULL,
 ...)
{
   ## Purpose is to summarize detectedTx supporting data for a gene
   ## where detectedTxL is the output of defineDetectedTx()
   ##
   ## detectedTxInfoByGene(detectedTxTPML, "Actb")
   ## detectedTxInfo(detectedTxTPML, Tx=i1)
   if (length(Gene) == 0 && length(Tx) == 0) {
      stop("detectedTxInfoByGene() requires Gene or TX to be supplied.");
   }

   if (length(Gene) == 0) {
      iTx <- detectedTxL$txExprGrpTx[detectedTxL$txExprGrpTx[,1] %in% Tx,1];
      iGene <- rownames(detectedTxL$txExprGrpTx)[detectedTxL$txExprGrpTx[,1] %in% Tx];
   } else {
      iTx <- detectedTxL$txExprGrpTx[rownames(detectedTxL$txExprGrpTx) %in% Gene,1];
      iGene <- rownames(detectedTxL$txExprGrpTx)[rownames(detectedTxL$txExprGrpTx) %in% Gene];
   }

   txDatNames <- setdiff(vigrep("^tx", names(detectedTxL)), "txExprGrpTx");
   lapply(jamba::nameVector(txDatNames), function(i){
      i1 <- match(iTx, rownames(detectedTxL[[i]]));
      iM <- detectedTxL[[i]][i1,,drop=FALSE];
      colnames(iM) <- gsub("^x[.]", "", colnames(iM));
      iMetric <- gsub("^tx|GrpAll$|M$", "", i);
      iDF <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         transcript_id=iTx,
         gene_name=iGene,
         metric=i,
         iM);
      iDF;
   });
}

#' Convert factor to a factor label
#'
#' Convert factor to a factor label
#'
#' This function is intended to take a vector of factor levels
#' and convert to labels using summary statistics. For example
#' by default it will add the number of items for each factor
#' level.
#'
#' This function is intended to help create useful ordered
#' factor labels that can be used in a ggplot2 visualization.
#'
#' @return factor vector with the same length as input `x` but
#'    where the levels also include summary information, such
#'    as the count of each factor level.
#'
#' @param x factor vector
#' @param valuesL optional list of numeric values, where each vector
#'    has the same length as `x`. If it is not a factor, it is
#'    converted to factor, using `jamba::mixedSort()` to order
#'    the factor levels.
#' @param types character value indicating the summary to use for
#'    the input factor `x`.
#' @param aggFun summary function used for each vector of numeric values
#'    in `valuesL` when supplied.
#' @param digits,big.mark arguments passed to `base::format()` to
#'    create a suitable text label.
#' @param itemName character vector indicating the name to associate to
#'    counts.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' x <- factor(rep(letters[1:5], c(2,4,3,2,1)));
#' levels(factor2label(x));
#' factor2label(x);
#'
#' @family jam plot functions
#'
#' @export
factor2label <- function
(x,
 valuesL=NULL,
 types=c("count","none"),
 aggFun=mean,
 digits=3,
 big.mark=",",
 itemName="items",
 ...)
{
   ## Purpose is to convert a factor vector into a factor using labels
   ## to summarize the factor. For example to count the number of
   ## entries for each factor level.
   ## x is a factor vector
   ##
   ## valuesL is either a named list or data.frame
   types <- match.arg(types);
   if (!igrepHas("factor", class(x))) {
      x <- factor(x,
         levels=mixedSort(unique(x)));
   }
   xNames <- levels(x);
   if (igrepHas("count", types)) {
      xVals1 <- paste(format(table(x),
         digits=digits,
         big.mark=big.mark,
         trim=TRUE),
         itemName);
   } else {
      xVals1 <- NULL;
   }
   if (length(valuesL) > 0) {
      xValsL <- lapply(nameVectorN(valuesL), function(i){
         iL <- split(valuesL[[i]], x);
         sapply(iL, function(j){
            paste(format(aggFun(rmNA(j)),
               digits=digits,
               big.mark=big.mark,
               trim=TRUE),
               i);
         });
      });
   } else {
      xValsL <- NULL;
   }
   x1L <- rmNULL(c(list(xVals1), xValsL));
   x1 <- do.call(paste, c(x1L, sep="; "));
   x2 <- paste0(xNames, " (", x1, ")");
   names(x2) <- xNames;
   x2f <- factor(x2, levels=x2);
   x2f[x];
}

#' Get first stranded GRanges feature per GRangesList
#'
#' Get first stranded GRanges feature per GRangesList
#'
#' This function returns the first feature per GRangesList,
#' relative to the strand of the GRangesList. It assumes each
#' GRangesList element has only one strand and one seqname,
#' and will stop otherwise.
#'
#' @return GRangesList containing only the first stranded
#'    GRanges feature per input GRangesList element. When
#'    `method="flank"` the output contains no metadata values,
#'    but this method is deprecated.
#'
#' @family jam ALE-specific RNA-seq functions
#' @family jam GRanges functions
#'
#' @param grl GRangesList
#' @param method character value in `c("direct", "endoapply", "flank")`
#'    representing which method to use to define the first feature.
#'    The `"direct"` method uses `IRanges::heads()` or
#'    `IRanges::tails()` depending upon the strandedness;
#'    the `"endoapply"` method uses `S4Vectors::endoapply()` to
#'    iterate each GRangesList element, then returns the result
#'    from `head()` or `tail()`. The `"flank"` method is
#'    intended to be equivalent but uses `IRanges::heads()` or
#'    `IRanges::tails()` then `GenomicRanges::flank()`, which
#'    is vectorized. The `"flank"` method is deprecated and
#'    `"direct"` is recommended as the preferred vectorized
#'    replacement.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' gr12 <- GenomicRanges::GRanges(seqnames=rep("chr1", 9),
#'    ranges=IRanges::IRanges(
#'       start=c(100, 200, 400, 500, 300, 100, 200, 400, 600),
#'       width=c(50,50,50, 50,50,50, 50,50,50)
#'    ),
#'    strand=rep(c("+", "-", "+"), c(3,3,3)),
#'    gene_name=rep(c("GeneA", "GeneB", "GeneC"), each=3)
#' )
#' grl <- GenomicRanges::split(gr12, GenomicRanges::values(gr12)$gene_name)
#' getFirstStrandedFromGRL(grl)
#'
#' @export
getFirstStrandedFromGRL <- function
(grl,
 method=c("direct", "endoapply", "flank"),
 verbose=FALSE,
 ...)
{
   ## Purpose is to take a GRangesList and return the first element
   ## in stranded order.
   ## It assumes each GRangesList element contains only one strand, such
   ## as with exons of a transcript.
   method <- match.arg(method);
   if (suppressPackageStartupMessages(!require(GenomicRanges))) {
      stop("getFirstStrandedFromGRL() requires the GenomicRanges Bioconductor package.");
   }
   ## First, ensure the GRangesList is sorted by strand, since we use that
   ## assumption in downstream steps
   if (verbose) {
      jamba::printDebug("getFirstStrandedFromGRL(): ",
         "sorting grl");
   }
   grl <- sortGRL(grl,
      verbose=verbose);

   # First check that each GRL has only one strand and seqname
   if (any(lengths(unique(strand(grl))) > 1)) {
      stop("Input GRangesList must contain only one strand per element.");
   }
   if (any(lengths(unique(GenomicRanges::seqnames(grl))) > 1)) {
      stop("Input GRangesList must contain only one seqname per element.");
   }

   if (method %in% "direct") {
      if (verbose) {
         jamba::printDebug("getFirstStrandedFromGRL(): ",
            "performing direct logic");
      }
      grl2 <- IRanges::heads(grl, 1);
      is_minus <- as.vector(unlist(strand(range(grl)))) %in% "-";
      if (any(is_minus)) {
         grl2[is_minus] <- IRanges::tails(grl[is_minus], 1);
      }
      return(grl2);
   } else if (method %in% "endoapply") {
      if (verbose) {
         jamba::printDebug("getFirstStrandedFromGRL(): ",
            "performing endoapply() logic");
      }
      grl2 <- S4Vectors::endoapply(grl, function(iGR){
         if ("-" %in% as.vector(strand(iGR[1]))) {
            tail(iGR, 1);
         } else {
            head(iGR, 1);
         }
      });
      return(grl2);
   } else {
      if (verbose) {
         jamba::printDebug("getFirstStrandedFromGRL(): ",
            "performing flank() logic");
      }
      grlWidths1 <- width(IRanges::heads(grl, 1));
      grlWidths2 <- width(IRanges::tails(grl, 1));
      if (verbose) {
         jamba::printDebug("getFirstStrandedFromGRL(): ",
            "applying range() function, class(grl):",
            class(grl));
      }
      grlRanges <- GenomicRanges::range(grl);
      if (verbose) {
         jamba::printDebug("getFirstStrandedFromGRL(): ",
            "Validating one strand per range(grl).");
      }
      if (length(grlRanges) != length(grlRanges@unlistData)) {
         stop("getFirstStrandedFromGRL() requires each GRanges element to have only one strand and one seqname.");
      }
      if (verbose) {
         jamba::printDebug("getFirstStrandedFromGRL(): ",
            "applying ifelse() logic");
      }
      grlFlankWidth <- ifelse(as.vector(unlist(strand(grlRanges))) %in% "+",
         unlist(grlWidths1),
         unlist(grlWidths2));
      if (verbose) {
         jamba::printDebug("getFirstStrandedFromGRL(): ",
            "applying flank() function");
      }
      grl3 <- GenomicRanges::flank(grlRanges,
         start=TRUE,
         width=-grlFlankWidth);
      if (verbose) {
         jamba::printDebug("getFirstStrandedFromGRL(): ",
            "flank() function complete.");
      }
      return(grl3);
   }
}

#' Sort GRangesList elements by chromosome and position
#'
#' Sort GRangesList elements by chromosome and position
#'
#' This function sorts GRangesList elements by chromosome and
#' position, using the default sort routine for GRanges objects.
#' It accomplishes the task by sorting the GRanges elements, then
#' splitting that GRanges into a GRangesList based upon the
#' original `names(GrangesList)`.
#'
#' @return GRangesList sorted by chromosome and position.
#'
#' @family jam GRanges functions
#'
#' @param GRL GRangesList to be sorted.
#' @param splitColname intermediate colname used to split values
#'    from GRanges to GRangesList. It only needs to be defined
#'    if for some reason the default "splitColname" is already
#'    in `colnames(values(GRL))`.
#' @param keep_order logical indicating whether to maintain the
#'    input order of GRanges; `FALSE` will return GRangesList
#'    ordered by the first sorted GRanges occurrence.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' gr12 <- GenomicRanges::GRanges(
#'    seqnames=rep(c("chr1", "chr2", "chr1"), c(3,3,3)),
#'    ranges=IRanges::IRanges(
#'       start=c(100, 200, 400, 500, 300, 100, 200, 400, 600),
#'       width=c(50,50,50, 50,50,50, 50,50,50)
#'    ),
#'    strand=rep(c("+", "-", "+"), c(3,3,3)),
#'    gene_name=rep(c("GeneA", "GeneB", "GeneC"), each=3)
#' )
#' grl <- GenomicRanges::split(gr12, GenomicRanges::values(gr12)$gene_name);
#' grl;
#' GenomicRanges::sort(grl);
#' sortGRL(grl);
#' sortGRL(grl, keep_order=FALSE);
#'
#' @export
sortGRL <- function
(GRL,
 splitColname="splitColname",
 keep_order=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to sort GRangesList by chromosome, using default sort(...)
   ## on the GRanges entries, then split back into GRangesList in consistent
   ## order
   if (!suppressPackageStartupMessages(require(GenomicRanges))) {
      stop("The GenomicRanges package is required for sortGRL()");
   }
   ## Ensure input GRL contains names
   if (is.null(names(GRL))) {
      names(GRL) <- jamba::makeNames(rep("GRL", length(GRL)));
   }

   GenomicRanges::values(GRL@unlistData)[,splitColname] <- factor(rep(names(GRL), IRanges::elementNROWS(GRL)),
      levels=names(GRL));
   GR1 <- GenomicRanges::sort(GRL@unlistData);
   if (!keep_order) {
      GenomicRanges::values(GR1)[[splitColname]] <- factor(GenomicRanges::values(GR1)[[splitColname]],
         levels=as.character(unique(GenomicRanges::values(GR1)[[splitColname]])));
   }
   if (verbose) {
      jamba::printDebug("sortGRL(): ",
         "splitting GRanges into list");
   }
   splitColnum <- match(splitColname, colnames(GenomicRanges::values(GR1)));
   GRL <- GenomicRanges::split(GR1[,-splitColnum],
      f=GenomicRanges::values(GR1)[[splitColname]]);

   return(GRL);
}

#' Annotate GRanges using another GRanges object
#'
#' Annotate GRanges using another GRanges object
#'
#' This function adds annotations to features in the given GRanges
#' object, from overlapping features in a second GRanges object.
#' It is especially useful after performing a manipulation that
#' results in a GRanges object without any `values` annotation,
#' for example `GenomicRanges::reduce()` or `GenomicRanges::intersect()`.
#'
#' In theory this function is relatively simple, it applies annotations
#' to overlapping entries. In practice, it gets complicated when multiple
#' annotations overlap the same GRange entry. In that case, numerical
#' values by default return the `base::mean()` while string values call
#' `jamba::cPasteUnique()` by default, and return the unique, sorted,
#' comma-delimited set of string values. For example, overlapping several
#' exons from the same gene, the resulting annotation might include
#' just one gene symbol.
#'
#' The numeric values can be aggregated using another function, controlled
#' by `numShrinkFunction` named using the colname. For example:
#' `numShrinkFunction="min"` would return the `base::min()` value for all
#' numeric columns. But `numShrinkFunction=list(default=mean, score=max)`
#' would return the `base::mean()` for all numeric columns except the
#' `"score"` column which would report the `base::max()`.
#'
#' Numeric values currently use the `data.table` package workflow, which
#' provides extremely efficient calculations for numeric values. However,
#' vectorized functions have the potential to be notably faster, as is
#' true with `jamba::cPasteUnique()` especially when applying
#' `jamba::mixedSort()` to sort alphanumeric values. In future, the
#' implementation may change to accomodate vectorized numeric functions,
#' instead of using `data.table`.
#'
#' TODO: This function is a specific case where another function
#' `shrinkDataFrame()` may be a more general purpose use. That function
#' is yet to be ported to the Jam package suite.
#'
#' @return GRanges object with colnames added to `values`, with length
#' and order equal to the input `GR1` GRanges object.
#'
#' @family jam GRanges functions
#'
#' @param GR1 GRanges object, the reference object to which annotations
#'    are added.
#' @param GR2 GRanges object used for annotations
#' @param grOL overlaps object, optionally used when the
#'    `GenomicRanges::findOverlaps()` function has already been run, mostly
#'    intended to save re-processing the overlaps for large objects.
#' @param numAsStrings logical indicating whether numerical values should
#'    be treated as strings for the purpose of added annotations. When
#'    `TRUE`, numerical values will be converted to character strings and
#'    retained as if they were labels. This argument treats all numeric
#'    columns as string, to treat only a subset use the `addStringCols`
#'    argument.
#' @param stringShrinkFunc function used to shrink string annotations
#'    by overlapping feature. This function should take a list as input,
#'    and `sep` as an argument indicating the delimiter if applicable,
#'    and return a character vector of the same length as the list. By
#'    default, `jamba::cPasteUnique()` is used. Control over the sort
#'    and uniqueness should be applied with this function.
#'    Note: `stringShrinkFunc` can be a single function, or can be a
#'    named list of function calls, whose names match the colnames
#'    of `values`. If a named list is provided, the first entry is used
#'    for any any names in `values` which are not in `names(stringShrinkFunc)`.
#' @param numShrinkFunc function used to shrink numerical values
#'    by overlapping feature. For numeric values, the `data.table` package
#'    is used to apply the function, which enables custom functions not
#'    otherwise available via the `GenomicRanges` methods. Therefore, the
#'    function does not take a list as input, instead takes a numeric
#'    vector as input.
#'    Note: `numShrinkFunc` can be a single function, or can be a
#'    named list of function calls, whose names match the colnames
#'    of `values`. If a named list is provided, the first entry is used
#'    for any any names in `values` which are not in `names(numShrinkFunc)`.
#' @param addStringCols character vector, optional, of numeric colnames that
#'    should be treated as string values. This option is a subset of the
#'    `numAsStrings`.
#' @param type,ignore.strand,select arguments sent to
#'    `GenomicRanges::findOverlaps()`. Note these options are only used
#'    when `grOL` is not supplied.
#' @param sep character value indicating the delimiter to use between
#'    string values.
#' @param verbose logical indicating whether to print verbose output.
#' @param DEBUG logical indicating whether to print debugging output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' gr12 <- GenomicRanges::GRanges(
#'    seqnames=rep(c("chr1", "chr2", "chr1"), c(3,3,3)),
#'    ranges=IRanges::IRanges(
#'       start=c(100, 200, 400, 500, 300, 100, 200, 400, 600),
#'       width=c(100,150,50, 50,50,100, 50,200,50)
#'    ),
#'    strand=rep(c("+", "-", "+"), c(3,3,3)),
#'    gene_name=rep(c("GeneA", "GeneB", "GeneC"), each=3)
#' )
#' gr1 <- gr12[,0];
#'
#' # Say for example you have a GRanges object with no annotation
#' gr1;
#'
#' # And had another GRanges object with annotations
#' gr12;
#'
#' # To add annotations
#' annotateGRfromGR(gr1, gr12);
#'
#' # Notice certain features overlap multiple annotations,
#' # which may be acceptable.
#'
#' # If you want to keep annotations distinct,
#' # use annotateGRLfromGRL()
#' grl1 <- GenomicRanges::split(gr12[,0],
#'    GenomicRanges::values(gr12)$gene_name);
#' grl2 <- GenomicRanges::split(gr12,
#'    GenomicRanges::values(gr12)$gene_name);
#'
#' # The first object is a GRangesList with no annotations
#' grl1;
#'
#' # The second object is a GRangesList with annotation,
#' # assumed to be in the same order
#' grl2;
#'
#' annotateGRLfromGRL(grl1, grl2);
#'
#' @export
annotateGRfromGR <- function
(GR1,
 GR2,
 grOL=NULL,
 numAsStrings=FALSE,
 stringShrinkFunc=function(...){jamba::cPasteUnique(..., doSort=TRUE)},
 numShrinkFunc=sum,
 addStringCols=NULL,
 type=c("any", "start", "end", "within", "equal"),
 ignore.strand=FALSE,
 select="all",
 sep=",",
 verbose=FALSE,
 DEBUG=FALSE,
 ...)
{
   ## Purpose is to try a new faster method of annotating GR1 with values(GR2)
   ## than annotateGRfromGR()
   ##
   ## Note: this function assumes data is stored as character vectors
   ## TODO: handle different types, e.g. numeric, CharacterList (tricky)
   ##
   ## TODO: for entries overlapping only one other feature, skip the shrinkMatrix()
   ## step and perform a straight paste(), which should be notably faster.
   ##
   ## useMixedSort=TRUE will use mixedSort() and mixedOrder() to sort
   ## string values, so the end result will be, for example:
   ## chr1, chr2, chr3, chr10, chr11
   ## and not:
   ## chr1, chr10, chr11, chr2, chr3
   ##
   ## grOL is an optional findOverlaps object, intended to allow
   ## external logic to be applied

   ## First make sure the GRanges objects both have names
   if (is.null(names(GR1)) || any(table(names(GR1)) > 1)) {
      names(GR1) <- paste0("GR1_", padInteger(seq_along(GR1)));
   }
   if (is.null(names(GR2)) || any(names(GR2) %in% c(NA, "")) || any(table(names(GR2)) > 1)) {
      names(GR2) <- paste0("GR2_", padInteger(seq_along(GR2)));
   }

   ## Added type argument, to allow specific types of overlaps
   type <- match.arg(type);

   ## Run findOverlaps()
   if (verbose) {
      printDebug("annotateGRfromGR(): ",
         "Running findOverlaps(GR1, GR2)");
   }
   if (is.null(grOL)) {
      grOL <- GenomicRanges::findOverlaps(query=GR1,
         subject=GR2,
         type=type,
         select=select,
         ignore.strand=ignore.strand);
   }
   if (verbose) {
      printDebug("annotateGRfromGR(): ",
         "Completed findOverlaps(GR1, GR2)");
      printDebug("annotateGRfromGR(): ",
         "length(grOL):",
         length(grOL));
   }
   if (length(grOL) == 0) {
      return(GR1);
   }

   grOLm <- as.matrix(grOL);

   ## Test for query entries with only one overlap in the subject
   grOLtable <- table(S4Vectors::from(grOL));
   grOLm1 <- grOLm[grOLm[,1] %in% which(grOLtable == 1),,drop=FALSE];
   grOLq1 <- grOLm1[,"queryHits"];
   grOLs1 <- grOLm1[,"subjectHits"];

   grOLmUse <- grOLm[!grOLm[,1] %in% which(grOLtable == 1),,drop=FALSE];
   grOLq <- grOLmUse[,"queryHits"];
   grOLs <- grOLmUse[,"subjectHits"];
   if (verbose) {
      printDebug("annotateGRfromGR(): ",
         "   grOLq1:", head(grOLq1, 10));
      printDebug("annotateGRfromGR(): ",
         "   grOLs1:", head(grOLs1, 10));
      printDebug("annotateGRfromGR(): ",
         "   grOLq:", head(grOLq, 10));
      printDebug("annotateGRfromGR(): ",
         "   grOLs:", head(grOLs, 10));
   }

   ## Shrink the values
   if (verbose) {
      printDebug("annotateGRfromGR(): ",
         "Running shrinkMatrix on ",
         formatInt(sum(grOLtable > 1)),
         " entries out of ",
         formatInt(length(grOLtable)));
   }

   ## Define the column types by inspecting data, and
   ## applying function arguments
   colClasses <- sapply(colnames(GenomicRanges::values(GR2)), function(iCol){
      class(GenomicRanges::values(GR2)[,iCol])
   });
   if (numAsStrings) {
      stringCols <- names(colClasses)[sapply(colClasses, function(iClass){
         igrepHas("integer|numeric|character|factor|ordered", iClass);
      })];
      numCols <- character(0);
      if (!is.list(stringShrinkFunc)) {
         stringShrinkFunc <- nameVector(
            rep(list(stringShrinkFunc),
               length(stringCols)),
            stringCols);
      } else if (!all(stringCols %in% names(stringShrinkFunc))) {
         iNewCols <- setdiff(stringCols,
            names(stringShrinkFunc));
         stringShrinkFunc <- c(stringShrinkFunc,
            nameVector(
               rep(stringShrinkFunc,
                  length.out=length(iNewCols)),
               iNewCols));
      }
   } else {
      numCols <- names(colClasses)[sapply(colClasses, function(iClass){
         igrepHas("integer|numeric|float", iClass);
      })];
      stringCols <- names(colClasses)[sapply(colClasses, function(iClass){
         igrepHas("character|factor|ordered", iClass);
      })];
      if (length(addStringCols) > 0 && any(addStringCols %in% numCols)) {
         moveNumCols <- intersect(addStringCols, numCols);
         numCols <- setdiff(numCols, moveNumCols);
         stringCols <- c(stringCols, moveNumCols);
         if (verbose) {
            printDebug("annotateGRfromGR(): ",
               "Moving ",
               cPaste(moveNumCols),
               " from numCols to stringCols.");
         }
      }
      if (length(numCols) > 0 && !is.list(numShrinkFunc)) {
         numShrinkFunc <- nameVector(
            rep(list(numShrinkFunc),
               length(numCols)),
            numCols);
      } else if (length(numCols) > 0 && !all(numCols %in% names(numShrinkFunc))) {
         iNewCols <- setdiff(numCols,
            names(numShrinkFunc));
         numShrinkFunc <- c(numShrinkFunc,
            nameVector(
               rep(numShrinkFunc,
                  length.out=length(iNewCols)),
               iNewCols));
      }
      if (length(stringCols) > 0 && !is.list(stringShrinkFunc)) {
         stringShrinkFunc <- nameVector(
            rep(list(stringShrinkFunc),
               length(stringCols)),
            stringCols);
      } else if (length(stringCols) > 0 && !all(stringCols %in% names(stringShrinkFunc))) {
         iNewCols <- setdiff(stringCols,
            names(stringShrinkFunc));
         stringShrinkFunc <- c(stringShrinkFunc,
            nameVector(
               rep(stringShrinkFunc,
                  length.out=length(iNewCols)),
               iNewCols));
      }
   }
   if (verbose) {
      printDebug("annotateGRfromGR(): ",
         "numCols:   ",
         numCols);
      printDebug("annotateGRfromGR(): ",
         "stringCols:",
         stringCols);
   }

   #####################################
   ## Shrink numeric columns
   if (length(numCols) > 0) {
      if (verbose) {
         printDebug("annotateGRfromGR(): ",
            "Shrinking num columns.");
      }
      if (nrow(grOLm1) > 0) {
         if (verbose) {
            printDebug("annotateGRfromGR(): ",
               "   grOLm1>0");
         }
         numShrunk1 <- lapply(nameVector(numCols), function(iCol){
            if (verbose) {
               printDebug("annotateGRfromGR(): ",
                  "   ",
                  iCol);
            }
            GenomicRanges::values(GR2[grOLs1])[,iCol];
         });
      }
      numShrunk <- lapply(nameVector(numCols), function(iCol){
         if (verbose) {
            printDebug("annotateGRfromGR(): ",
               "   ",
               iCol);
            printDebug("annotateGRfromGR(): ",
               "   groupBy:",
               head(names(GR1)[grOLq]));
            printDebug("annotateGRfromGR(): ",
               "   grOLq:",
               head(grOLq));
         }
         grOLi <- shrinkMatrix(as.data.frame(GenomicRanges::values(GR2)[grOLs,iCol,drop=FALSE]),
            groupBy=seq_along(GR1)[grOLq],
            shrinkFunc=numShrinkFunc[[iCol]],
            returnClass="matrix");
      });
   }

   #####################################
   ## Shrink string columns
   if (length(stringCols) > 0) {
      if (verbose) {
         printDebug("annotateGRfromGR(): ",
            "Shrinking string columns:",
            stringCols);
         printDebug("annotateGRfromGR(): ",
            "nrow(grOLm1):",
            nrow(grOLm1));
         printDebug("annotateGRfromGR(): ",
            "nrow(grOLmUse):",
            nrow(grOLmUse));
      }
      if (nrow(grOLm1) > 0) {
         ## Note: Decided to convert to character vector to avoid
         ## factors being converted to numeric values, then to string,
         ## and generally to be consistent with the type produced by
         ## the string shrink function below
         stringShrunk1 <- lapply(nameVector(stringCols), function(iCol){
            if (verbose) {
               printDebug("annotateGRfromGR(): ",
                  "   ",
                  iCol);
            }
            as.data.frame(GenomicRanges::values(GR2)[grOLs1,iCol,drop=FALSE]);
         });
      }
      if (nrow(grOLmUse) > 0) {
         stringShrunk <- lapply(nameVector(stringCols), function(iCol){
            if (verbose) {
               printDebug("annotateGRfromGR(): ",
                  "   ",
                  iCol);
            }
            if (igrepHas("list", class(GenomicRanges::values(GR2)[grOLs,iCol]))) {
               ## list column
               ## What was all this effort doing?
               #iVals <- values(GR2)[grOLs,iCol];
               #names(iVals) <- names(GR2)[grOLs];
               #iValsX <- unlist(iVals);
               #iValsXnames1 <- rep(grOLq,
               #   S4Vectors::elementNROWS(iVals));
               #iValsSplit <- split(iValsX, iValsXnames1);
               #iValsSplit <- iVals;
               iValsSplit <- GenomicRanges::values(GR2)[grOLs,iCol];

               iXnonNA <- stringShrinkFunc[[iCol]](iValsSplit,
                  sep=sep);
               iX[names(iXnonNA)] <- iXnonNA;
               iX;
            } else {
               ## Non-list column
               iValsX <- GenomicRanges::values(GR2)[grOLs,iCol];
               iValsXnames1 <- grOLq;

               ## Split the values
               nonNA <- which(!is.na(iValsX));
               stringLnonNA <- split(iValsX[nonNA], names(GR1)[iValsXnames1[nonNA]]);
               #stringLnonNA <- split(iValsX[nonNA], iValsXnames1[nonNA]);
               iXnonNA <- stringShrinkFunc[[iCol]](stringLnonNA,
                  sep=sep);
               iX <- nameVector(rep(NA, length(unique(grOLq))), names(GR1)[unique(grOLq)]);
               iX[names(iXnonNA)] <- iXnonNA;
               iX;
            }
         });
      }
      if (DEBUG) {
         return(stringShrunk);
      }
   }
   ## Note: we keep track of grOLqAll to represent the single- and multi-entry rows
   ## from the original GRanges object
   grOLqAll <- c(unique(grOLq), grOLq1);
   grOLsAll <- c(grOLs, grOLs1);


   #####################################
   ## Put it all back together
   stringShrunkDF <- NULL;
   numShrunkDF <- NULL;
   if (length(numCols) > 0) {
      if (verbose) {
         printDebug("annotateGRfromGR(): ",
            "length(numShrunk):",
            length(numShrunk));
         if (length(numShrunk) > 0) {
            printDebug("annotateGRfromGR(): ",
               "class(numShrunk):",
               class(numShrunk));
            print(head(numShrunk, 3));
         }
      }
      if (nrow(grOLmUse) > 0) {
         numShrunkDF <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            do.call(cbind, numShrunk));
         if (verbose) {
            printDebug("annotateGRfromGR(): ",
               "   nrow(numShrunkDF):",
               formatInt(nrow(numShrunkDF)),
               " (before)");
         }
      }
      if (nrow(grOLm1) > 0) {
         ## Create data.frame using the original entries
         numShrunkDF1 <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            do.call(cbind, numShrunk1));
         ## Append the shrunken multi-entry and the single-entry data.frames
         if (nrow(grOLmUse) > 0 && !is.null(numShrunkDF)) {
            numShrunkDF <- rbind(numShrunkDF, numShrunkDF1);
         } else {
            numShrunkDF <- numShrunkDF1;
         }
      }
   }
   if (length(stringCols) > 0) {
      if (nrow(grOLmUse) > 0 && !is.null(stringShrunk)) {
         stringShrunkDF <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            do.call(cbind, stringShrunk));
      }
      if (nrow(grOLm1) > 0) {
         if (verbose) {
            printDebug("annotateGRfromGR(): ",
               "Appending multi- and single-entry ",
               "string",
               " data.frames");
         }
         ## Create data.frame using the original entries
         stringShrunkDF1 <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            do.call(cbind, stringShrunk1));
         ## Append the shrunken multi-entry and the single-entry data.frames
         if (nrow(grOLmUse) > 0 && !is.null(stringShrunkDF)) {
            stringShrunkDF <- rbind(stringShrunkDF,
               stringShrunkDF1);
         } else {
            stringShrunkDF <- stringShrunkDF1;
         }
      }
      if (length(numCols) > 0 && !is.null(stringShrunkDF)) {
         grOL1 <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            stringShrunkDF,
            numShrunkDF);
      } else {
         grOL1 <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            stringShrunkDF);
      }
   } else if (length(numCols) > 0) {
      grOL1 <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         numShrunkDF);
   } else {
      grOL1 <- data.frame(row.names=unique(names(GR1)[grOLqAll]),
         olNames=unique(names(GR1)[grOLqAll]))[,-1,drop=FALSE];
   }

   ## Make duplicate colnames uniquely named
   ## TODO: review make.unique() for making column renaming robust,
   ## currently does not handle renaming duplicated versioned names.
   if (ncol(GenomicRanges::values(GR1)) > 0 &&
         any(colnames(grOL1) %in% colnames(GenomicRanges::values(GR1)))) {
      newColnames <- make.unique(c(colnames(GenomicRanges::values(GR1)), colnames(grOL1)),
         sep="_v");
      newColnames2 <- tail(newColnames, ncol(grOL1));
      colnames(grOL1) <- newColnames2;
   }

   ## Now append each column, meanwhile fix some issues with ,-delimiters
   for (iCol in colnames(grOL1)) {
      if (igrepHas("integer|numeric", class(grOL1[,iCol]))) {
         GenomicRanges::values(GR1)[,iCol] <- numeric(0);
      } else {
         grOL1[,iCol] <- gsub(", ", ",", grOL1[,iCol]);
         GenomicRanges::values(GR1)[,iCol] <- "";
      }
      if (verbose) {
         printDebug("annotateGRfromGR(): ",
            "head(grOL1):");
         print(head(grOL1, 5));
      }
      GenomicRanges::values(GR1)[grOLqAll,iCol] <- grOL1[,iCol];
      blankRows <- seq_along(GR1)[-grOLqAll];
      if (length(blankRows) > 0) {
         GenomicRanges::values(GR1[blankRows])[,iCol] <- NA;
      }
   }
   GR1;
}

#' Annotate GRangesList from GRangesList objects
#'
#' Annotate GRangesList from GRangesList objects
#'
#' This function extends `annotateGRfromGR()` for the special case
#' of GRangesList objects. It requires both GRangesList objects have
#' identical length, and assumes both are in equivalent order.
#' It then restricts all overlapping annotations to those where the
#' query and subject are the same original GRangesList index.
#'
#' This function is particularly useful following an operation on
#' a GRangesList object that otherwise removes all annotations in
#' `values(GRL1)`, for example `GenomicRanges::reduce()` or
#' `GenomicRanges::flank()`. This function can be used to re-annotate
#' the resulting features using the original GRangesList object.
#'
#' Note that annotations are added at the level of individual GRanges
#' entries, equivalent to `values(GRL1@unlistData)`. This function does
#' not currently apply annotations at the GRangesList level, thus
#' it does not use `values(GRL2)` if they exist.
#'
#' @return GRangesList object with the same length and lengths as
#'    the input `GRL1`, with annotation columns added from `GRL2`.
#'
#' @family jam GRanges functions
#'
#' @param GRL1 GRangesList query
#' @param GRL2 GRangesList subject, used to add annotations to `GRL1`
#' @param annoName1 character value indicating either the colname of
#'    `values(GRL1)` to use as the name, or if `"name"` then it uses
#'    `names(GRL1)`.
#' @param annoName2 character value indicating either the colname of
#'    `values(GRL2)` to use as the name, or if `"name"` then it uses
#'    `names(GRL2)`.
#' @param grlOL overlap result (optional) from `GenomicRanges::findOverlaps()`
#'    for these same GRangesList objects, used to save time by not re-running
#'    `GenomicRanges::findOverlaps()` again.
#' @param add_grl_name logical indicating whether to add the names of each
#'    GRangesList object to the output object, useful for tracking the
#'    annotations to the source data.
#' @param returnType character value indicating whether to return
#'    GRangesList "GRL" or GRange "GR" object.
#' @param splitColname character value used internally to indicate how
#'    to split the resulting GRanges annotated data back to GRangesList.
#'    Almost always, this value should be identical to `annoName1` which
#'    will split the resulting GRanges back into the identical input
#'    GRangesList.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to `annotateGRfromGR()`.
#'    To customize the aggregation functions, supply `numShrinkFunc` or
#'    `stringShrinkFunc` as described in `annotateGRfromGR()`.
#'
#' @examples
#' gr12 <- GenomicRanges::GRanges(
#'    seqnames=rep(c("chr1", "chr2", "chr1"), c(3,3,3)),
#'    ranges=IRanges::IRanges(
#'       start=c(100, 200, 400, 500, 300, 100, 200, 400, 600),
#'       width=c(100,150,50, 50,50,100, 50,200,50)
#'    ),
#'    strand=rep(c("+", "-", "+"), c(3,3,3)),
#'    gene_name=rep(c("GeneA", "GeneB", "GeneC"), each=3)
#' )
#'
#' # Now split into GRangesList
#' grl1 <- GenomicRanges::split(gr12[,0],
#'    GenomicRanges::values(gr12)$gene_name);
#' grl2 <- GenomicRanges::split(gr12,
#'    GenomicRanges::values(gr12)$gene_name);
#'
#' # The first object is a GRangesList with no annotations
#' grl1;
#'
#' # The second object is a GRangesList with annotation,
#' # assumed to be in the same order
#' grl2;
#'
#' annotateGRLfromGRL(grl1, grl2);
#'
#' @export
annotateGRLfromGRL <- function
(GRL1,
 GRL2,
 annoName1="name",
 annoName2="name",
 grlOL=NULL,
 add_grl_name=FALSE,
 returnType=c("GRL", "GR"),
 splitColname=annoName1,
 verbose=FALSE,
 ...)
{
   ## Purpose is to run annotateGRfromGR() except allow for GRangesList
   ## objects as input.  It accomplishes the task by running
   ## findOverlapsGRL() which ensures GRangesList overlaps are only
   ## allowed for the same annotated entries, defined in annoName1,
   ## and annoName2.
   returnType <- match.arg(returnType);
   ## Assign names to GRL1 and GRL2 as needed
   if (length(names(GRL1@unlistData)) == 0) {
      names(GRL1@unlistData) <- jamba::makeNames(rep("grl1",
         length(GRL1@unlistData)));
   }
   if (length(names(GRL2@unlistData)) == 0) {
      names(GRL2@unlistData) <- jamba::makeNames(rep("grl2",
         length(GRL2@unlistData)));
   }

   ## Validate annoName1
   annoName1 <- head(annoName1, 1);
   if ("name" %in% annoName1 || "name" %in% splitColname) {
      GenomicRanges::values(GRL1@unlistData)[,"grl_name1"] <- rep(names(GRL1),
         S4Vectors::elementNROWS(GRL1));
   }
   if ("name" %in% annoName1) {
      annoName1 <- "grl_name1";
   }
   if ("name" %in% splitColname) {
      splitColname[splitColname %in% "name"] <- "grl_name1";
   }
   if (!annoName1 %in% colnames(GenomicRanges::values(GRL1@unlistData))) {
      stop(paste0("Supplied annoName1:'",
         annoName1,
         "' was not found in colnames(values(GRL1@unlistData))."));
   }
   ## Validate annoName2
   annoName2 <- head(annoName2, 1);
   if ("name" %in% annoName2) {
      GenomicRanges::values(GRL2@unlistData)[,"grl_name2"] <- rep(names(GRL2),
         S4Vectors::elementNROWS(GRL2));
      annoName2 <- "grl_name2";
   }
   if (!annoName2 %in% colnames(GenomicRanges::values(GRL2@unlistData))) {
      stop(paste0("Supplied annoName2:'",
         annoName2,
         "' was not found in colnames(values(GRL2@unlistData))."));
   }
   annoNames2 <- setdiff(colnames(GenomicRanges::values(GRL2@unlistData)), annoName2);
   if (verbose) {
      printDebug("annotateGRLfromGRL(): ",
         "annoNames2:",
         annoNames2);
   }

   ## Find overlaps in GRL fashion
   if (is.null(grlOL)) {
      grlOL <- findOverlapsGRL(GRL1,
         GRL2,
         annoName1=annoName1,
         annoName2=annoName2);
   }

   GR12 <- annotateGRfromGR(GRL1@unlistData,
      GRL2@unlistData[,annoNames2],
      grOL=grlOL,
      verbose=verbose,
      ...);
   if (returnType %in% "GR") {
      if (!add_grl_name) {
         keepCols <- setdiff(colnames(GenomicRanges::values(GR12)),
            c("grl_name1", "grl_name2"));
         GenomicRanges::values(GR12) <- GenomicRanges::values(GR12)[,keepCols];
      }
      return(GR12);
   } else {
      GRL12 <- GenomicRanges::split(GR12,
         GenomicRanges::values(GR12)[,splitColname]);
      if (!add_grl_name) {
         keepCols <- setdiff(colnames(GenomicRanges::values(GRL12@unlistData)),
            c("grl_name1", "grl_name2"));
         GenomicRanges::values(GRL12@unlistData) <- GenomicRanges::values(GRL12@unlistData)[,keepCols];
      }
      return(GRL12);
   }
}

#' Find overlaps between two GRangesList objects
#'
#' Find overlaps between two GRangesList objects
#'
#' This function implements `GenomicRanges::findOverlaps()` for the special
#' case of two GRangesList objects, restricting results to those including
#' the same GRangesList index in the subject and query.
#'
#' @family jam GRanges functions
#'
#' @return Hits object, or the natural output from
#'    `GenomicRanges::findOverlaps()` dependent upon the `...` arguments,
#'    subsetted for only entries with matching values defined
#'    by `annoName1` and `annoName2`.
#'
#' @family jam GRanges functions
#'
#' @param GRL1 GRangesList query
#' @param GRL2 GRangesList subject
#' @param annoName1 character value indicating either the colname of
#'    `values(GRL1)` to use as the name, or if `"name"` then it uses
#'    `names(GRL1)`.
#' @param annoName2 character value indicating either the colname of
#'    `values(GRL2)` to use as the name, or if `"name"` then it uses
#'    `names(GRL2)`.
#' @param check_names logical indicating whether the values defined
#'    by `annoName1` should match the values defined by `annoName2`.
#'    Note that when `check_names=FALSE` the overlaps returned will
#'    depend upon GRL1 and GRL2 being in identical order. Otherwise,
#'    the values in `annoName1` and `annoName2` are used, which means
#'    GRL1 and GRL2 do not have to be in identical order.
#' @param ... additional arguments are passed to
#'    `GenomicRanges::findOverlaps()`, useful for customizing the overlap
#'    criteria.
#'
#' @export
findOverlapsGRL <- function
(GRL1,
 GRL2,
 annoName1="name",
 annoName2="name",
 check_names=TRUE,
 ...)
{
   ## Purpose is to run findOverlaps() on GRangesList objects,
   ## where overlaps are required to share the same annotation,
   ## in the annoName column.
   ##
   ## Specifically, it helps run findOverlaps() on GRangesList objects
   ## separated by gene, then only returning entries which match the same
   ## gene.
   if (annoName1 %in% "name") {
      GRLnames1 <- rep(names(GRL1), elementNROWS(GRL1));
   } else {
      GRLnames1 <- GenomicRanges::values(GRL1@unlistData)[,annoName1];
   }
   if (annoName2 %in% "name") {
      GRLnames2 <- rep(names(GRL2), elementNROWS(GRL2));
   } else {
      GRLnames2 <- GenomicRanges::values(GRL2@unlistData)[,annoName2];
   }
   if (check_names && !any(GRLnames1 %in% GRLnames2)) {
      warning("No names are shared between GRL1 and GRL2. Please try again.");
      return(GRL1);
   }

   ## Perform overlap between of all GRanges
   grOL <- GenomicRanges::findOverlaps(GRL1@unlistData,
      GRL2@unlistData,
      ...);
   grOLdf <- data.frame(as.data.frame(grOL));
   grOLdf[,"queryName"] <- GRLnames1[grOLdf[,"queryHits"]];
   grOLdf[,"subjectName"] <- GRLnames2[grOLdf[,"subjectHits"]];
   ## Subset the overlap for matching GRL entries
   if (check_names) {
      ## Match by name
      grOLdfUse <- which(grOLdf[,"queryName"] == grOLdf[,"subjectName"]);
   } else {
      ## Match by the index, assuming GRL1 and GRL2 are in identical order
      grOLdfUse <- which(grOLdf[,"queryHits"] == grOLdf[,"subjectHits"]);
   }
   return(grOL[grOLdfUse]);
}

#' Assign exon names to GRangesList
#'
#' Assign exon names to GRangesList
#'
#' This function takes a GRangesList object with an annotated gene symbol
#' column, and defines exon numbers for each distinct (non-adjacent)
#' range per gene. When multiple ranges overlap, one exon number is
#' applied to them all, and disjoint ranges are denoted using a letter
#' suffix.
#'
#' For example the exon labels below:
#'
#' `|======|......|======|=======|======|......|======|=======|`
#' `.exon1.........exon2a.exon2b..exon2c........exon3a..exon3b.`
#'
#' The full name for each feature will become:
#' * Gene_exon1
#' * Gene_exon2a
#' * Gene_exon2b
#' * Gene_exon2c
#' * Gene_exon3a
#' * Gene_exon3b
#'
#' The reasoning, is to preserve the gene symbol for readability, but
#' to number exons to indicate the numbered contiguous exon, with suffix
#' to indicate the sub-section of each exon.
#'
#' @return GRangesList object with additional columns indicating the
#' exon name.
#'
#' @family jam GRanges functions
#' @family jam RNA-seq functions
#'
#' @param GRL GRangesList input object. Ideally, the input GRanges are
#'    disjoint, meaning no two exons overlap, but instead are represented
#'    by non-overlapping regions that are sometimes adjacent.
#' @param geneSymbolColname character string indicating with colname
#'    of `values(GRL)` contains the unique gene symbol. If this value
#'    does not already exist, it is created and populated using
#'    `names(GRL)`.
#' @param exonNameColname character string to be used for the resulting
#'    exon name. By default `"Exon"` is appended to the end of the
#'    `geneSymbolColname`.
#' @param suffix character value indicating the suffix to add to the
#'    `geneSymbolColname` to indicate individual exon numbers.
#' @param renameOnes logical indicating whether to name exon sub-sections
#'    when the exon is not subdivided. For example, when `renameOnes=FALSE`,
#'    an exon with one section would be named, `"exon1"`, otherwise
#'    when `renameOnes=TRUE` an exon with one section would be named,
#'    `"exon1a"`. This distinction can be helpful to recognize exons that
#'    contain no subsections.
#' @param filterTwoStrand logical indicating whether to filter out genes
#'    occurring on two different strands.
#' @param checkDisjoin character value indicating how to handle non-disjoint
#'    input GRL ranges. When `checkDisjoin="stop"` then any non-disjoint
#'    GRanges in an element of `GRL` will cause the function to fail.
#' @param assignGRLnames logical indicating whether names for the
#'    resulting GRangesList should use the exon names.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
assignGRLexonNames <- function
(GRL,
 geneSymbolColname="geneSymbol",
 exonNameColname=paste0(geneSymbolColname, "Exon"),
 suffix="_exon",
 renameOnes=FALSE,
 filterTwoStrand=TRUE,
 checkDisjoin=c("warn","none","stop"),
 assignGRLnames=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to assign exon names using numbers
   ## to represent contiguous segments, and lowercase
   ## letters to represent subsections of each exon.
   ##
   ## filterTwoStrand=TRUE will remove entries which have two strands
   ## for the same GRL entry
   ##
   ## This function is a light wrapper for renumberGRanges()
   ##
   ## checkDisjoin="stop" will check to make sure exons in each set
   ## of GRanges are disjoint, otherwise the numbering can be
   ## problematic.
   ## The most common symptom is negative strand embedded exons receive
   ## sub-numbers in opposite order.
   ##
   checkDisjoin <- match.arg(checkDisjoin);

   ## First verify that incoming data is valid per assumptions
   ## that exons for a transcript would all appear only on one strand
   if (verbose) {
      printDebug("assignGRLexonNames(): ",
         "class(GRL):",
         class(GRL));
   }
   GRLstrandL <- unique(GenomicRanges::strand(GRL));
   if (filterTwoStrand && any(S4Vectors::elementNROWS(GRLstrandL) > 1)) {
      if (verbose) {
         printDebug("assignGRLexonNames(): ",
            "removing some multi-stranded exon entries.");
      }
      iRemove <- which(S4Vectors::elementNROWS(GRLstrandL) > 1);
      GRL <- GRL[-iRemove];
   } else {
      if (verbose) {
         printDebug("assignGRLexonNames(): ",
            "No multi-stranded exon entries.");
      }
   }

   ## check disjoint GRanges
   if (checkDisjoin %in% c("warn","stop")) {
      if (verbose) {
         printDebug("assignGRLexonNames(): ",
            "Checking disjoint ranges.");
      }
      GRLdis <- GenomicRanges::disjoin(GRL);
      if (!all(S4Vectors::elementNROWS(GRLdis) == S4Vectors::elementNROWS(GRL))) {
         if (checkDisjoin %in% "stop") {
            stop("assignGRLexonNames() detected overlapping GRanges, stopping.");
         } else {
            printDebug("assignGRLexonNames(): ",
               "detected overlapping GRanges, continuing.",
               fgText=c("red","orange"));
         }
      }
   }

   ## Reduce entries
   if (verbose) {
      printDebug("assignGRLexonNames(): ",
         "Reducing ranges.");
   }
   GRLred <- GenomicRanges::reduce(GRL);

   ## Add geneSymbolColname if it does not already exist
   if (!geneSymbolColname %in% colnames(GenomicRanges::values(GRLred@unlistData))) {
      GenomicRanges::values(GRLred@unlistData)[,geneSymbolColname] <- rep(names(GRLred),
         S4Vectors::elementNROWS(GRLred));
   }
   if (verbose) {
      jamba::printDebug("assignGRLexonNames(): ",
         "head(GRLred):");
      print(head(GRLred));
   }
   if (verbose) {
      jamba::printDebug("assignGRLexonNames(): ",
         "geneSymbolColname:",
         geneSymbolColname);
   }
   if (verbose) {
      jamba::printDebug("assignGRLexonNames(): ",
         "geneSymbolColname values:",
         head(GenomicRanges::values(GRLred@unlistData)[,geneSymbolColname], 10));
   }
   GRLredStrand <- unlist(unique(strand(GRLred)));
   GRLredStrandP <- which(GRLredStrand %in% "+");
   GRLredStrandN <- which(GRLredStrand %in% "-");

   ## Stranded exon numbering
   GenomicRanges::values(GRLred@unlistData)[,exonNameColname] <- "";
   if (verbose) {
      printDebug("assignGRLexonNames(): ",
         "head(GRLred):");
      print(head(GRLred));
      printDebug("assignGRLexonNames(): ",
         "head(GRLredStrandP):");
      print(head(GRLredStrandP));
      printDebug("assignGRLexonNames(): ",
         "head(GRLredStrandN):");
      print(head(GRLredStrandN));
   }
   if (length(GRLredStrandP) > 0) {
      GenomicRanges::values(GRLred[GRLredStrandP]@unlistData)[,exonNameColname] <- jamba::makeNames(
         GenomicRanges::values(GRLred[GRLredStrandP]@unlistData)[,geneSymbolColname],
         suffix=suffix,
         renameOnes=TRUE);
   }
   if (length(GRLredStrandN) > 0) {
      GenomicRanges::values(GRLred[GRLredStrandN]@unlistData)[,exonNameColname] <- rev(jamba::makeNames(
         GenomicRanges::values(GRLred[rev(GRLredStrandN)]@unlistData)[,geneSymbolColname],
         suffix=suffix,
         renameOnes=TRUE));
   }
   if (verbose) {
      printDebug("assignGRLexonNames(): ",
         "exonNameColname values:",
         head(GenomicRanges::values(GRLred@unlistData)[,exonNameColname], 10));
   }

   ## Add lowercase letter suffix
   GRLcolnames <- unvigrep(paste0(exonNameColname, "(_v[0-9]|)$"),
      colnames(GenomicRanges::values(GRL@unlistData)));
   if (verbose) {
      printDebug("assignGRLexonNames(): ",
         "head(GRL[,GRLcolnames]):");
      print(head(GRL[,GRLcolnames]));
      printDebug("assignGRLexonNames(): ",
         "head(GRLred[,exonNameColname]):");
      print(head(GRLred[,exonNameColname]));
   }
   GRLnew <- annotateGRLfromGRL(GRL1=GRL[,GRLcolnames],
      GRL2=GRLred[,exonNameColname],
      verbose=verbose);
   if (verbose) {
      printDebug("assignGRLexonNames(): ",
         "Completed annotateGRLfromGRL().");
   }
   GRLnewStrand <- unlist(unique(strand(GRLnew)));
   GRLnewStrandP <- which(GRLnewStrand %in% "+");
   GRLnewStrandN <- which(GRLnewStrand %in% "-");
   GRLnewStrandNn <- names(GRLnew[GRLnewStrandN]@unlistData);
   subFeatureNumberStyle <- "letters";
   subFeatureSuffix <- "";
   exonNameColname1 <- paste0(exonNameColname, "1");

   GenomicRanges::values(GRLnew@unlistData)[,exonNameColname] <-
      GenomicRanges::values(GRLnew@unlistData)[,exonNameColname];
   GenomicRanges::values(GRLnew[GRLnewStrandP]@unlistData)[,exonNameColname] <- (
      jamba::makeNames(
         GenomicRanges::values(GRLnew[GRLnewStrandP]@unlistData)[,exonNameColname],
         numberStyle=subFeatureNumberStyle,
         suffix=subFeatureSuffix,
         renameOnes=renameOnes));
   GenomicRanges::values(GRLnew@unlistData[rev(GRLnewStrandNn),])[,exonNameColname] <- (
      jamba::makeNames(
         GenomicRanges::values(GRLnew@unlistData[rev(GRLnewStrandNn),])[,exonNameColname],
         numberStyle=subFeatureNumberStyle,
         suffix=subFeatureSuffix,
         renameOnes=renameOnes));

   if (assignGRLnames) {
      names(GRLnew@unlistData) <- jamba::makeNames(
         GenomicRanges::values(GRLnew@unlistData)[,exonNameColname]);
   }

   return(GRLnew);

}

#' Prepare ALE data for violin plots
#'
#' Prepare ALE data for violin plots
#'
#' This function takes output from `tx2ale()`, applies some filtering
#' to the output data, then returns a tall data.frame sufficient
#' for viewing as a violin plot using `ggplot2::geom_violin()`.
#'
#' @param iMatrixALE numeric matrix of expression data containing
#'    ALE rows, as output from `tx2ale()`. Each rowname is expected
#'    to have a suffix "_ale1" where the number represents the stranded
#'    order of ALE elements in a given gene. Therefore "_ale1" is the
#'    shortest form, closest to the 5-prime end of the transcript, and
#'    "_ale2" is the next ALE element downstream, and so on.
#' @param groups vector of group labels, named by `colnames(iMatrixALE)`.
#' @param facet_groups vector of group labels, named by `colnames(iMatrixALE)`,
#'    as a possible alternative to using the `groups`, for example
#'    for higher level grouping.
#' @param facet_name character string used to label the `facet_groups`.
#' @param maxGroupMeanALE numeric value indicating the threshold for
#'    including an ALE in the output data, where the max group mean
#'    (the highest group mean) is at least this value.
#' @param removeAboveAleNum character string of the ALE colname to remove,
#'    restricting data to include only ALE values below this number. For
#'    example "ale3" would remove "ale3" and all higher ALE numbers,
#'    thereby restricting data to "ale1" and "ale2".
#' @param maxGroupMeanFloor numeric threshold used as a noise floor.
#' @param returnAll logical indicating whether to return intermediate
#'    data formats in the output results list.
#' @param geneLists list containing vectors of genes, or a data.frame
#'    with two columns "gene_name" and "geneList", used
#'    to assign genes to one or more lists in the resulting violin plot.
#' @param lineAlpha numeric value of alpha transparency, scaled between
#'    0 and 1, used to draw lines from "_ale1" to "_ale2" on the
#'    violin plot.
#' @param subsetFunc optional function that takes the tall format data used in the
#'    primary violin plot, and returns data in the same format after
#'    applying logic specific for filtering this data. Intended to
#'    restrict display of genes to `groups` or `facet_groups` that are
#'    relevant to each `geneLists` entry.
#' @param make_ggplots logical indicating whether to create plot objects
#'    using ggplot2.
#' @param verbose logical indicating whether to print verbose output
#' @param ... additional arguments are ignored.
#'
#' @family jam ALE-specific RNA-seq functions
#'
#' @export
ale2violin <- function
(iMatrixAle=NULL,
 iMatrixAleGrp=NULL,
 groups,
 facet_groups=groups,
 facet_name="Group",
 maxGroupMeanALE=2,
 removeAboveAleNum="ale3",
 maxGroupMeanFloor=0,
 returnAll=FALSE,
 geneLists,
 lineAlpha=0.1,
 subsetFunc=NULL,
 make_ggplots=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to take an ALE expression matrix, and create
   ## a tall version suitable for plotting in ggplot2 to create
   ## a violin plot, connected by gene_Region
   if (suppressPackageStartupMessages(!require(ggplot2))) {
      stop("ale2violin() requires the ggplot2 package.");
   }
   retVals <- list();

   if (length(iMatrixAleGrp) == 0) {
      if (!all(names(groups) %in% colnames(iMatrixAle))) {
         stop("ale2violin() requires all(names(groups) %in% colnames(iMatrixAle))");
      }
      iMatrixAle <- iMatrixAle[,names(groups),drop=FALSE];
      ## Calculate group mean expression values
      iDiffAle <- data.frame(check.names=FALSE,
         rowGroupMeans(iMatrixAle,
            useMedian=FALSE,
            groups=groups));
   } else {
      if (!all(names(groups) %in% colnames(iMatrixAleGrp))) {
         stop("ale2violin() requires all(names(groups) %in% colnames(iMatrixAleGrp))");
      }
      iDiffAle <- iMatrixAleGrp[,names(groups),drop=FALSE];
   }

   ## Get groupMeans for each ALE for each sample group
   iDiffAle[,"aleNum"] <- gsub("^.+_(ale[0-9]+)$", "\\1",
      rownames(iDiffAle));
   iDiffAle[,"gene_name"] <- gsub("_ale.+", "",
      rownames(iDiffAle));

   ## Create tall version of the data
   iDiffAleTall <- data.table::melt(iDiffAle,
         variable.name="Group",
         id.vars=c("aleNum", "gene_name"),
         value.name="intensity");
   ## Convert back to wide format, using ale number per column
   iDiffAleWide <- data.table::dcast(iDiffAleTall,
      formula=gene_name + Group ~ aleNum,
      value.var="intensity")
   aleCols <- vigrep("^ale[0-9]+$", colnames(iDiffAleWide));
   iDiffAleWideM <- as.matrix(iDiffAleWide[,aleCols,drop=FALSE]);
   rownames(iDiffAleWideM) <- pasteByRow(iDiffAleWide[,c("gene_name","Group")]);

   if (returnAll) {
      retVals$iDiffAle <- iDiffAle;
      retVals$iDiffAleTall <- iDiffAleTall;
      retVals$iDiffAleWide <- iDiffAleWide;
   }

   ## Filter rows where:
   ## - the max ale value is below threshold
   ## - two UTRs only, removing genes with only 1 or with 3+
   if (verbose) {
      printDebug("ale2violin(): ",
         "filtering ALE maxGroupMean>=", maxGroupMeanALE,
         " and removeAboveAleNum:", removeAboveAleNum);
   }

   ## Filter for:
   ## - rowMax at or above threshold
   ## - non-empty ale2
   ## - empty ale3
   iDiffAleWideMGMkeep <- (
      rowMaxs(iDiffAleWideM, na.rm=TRUE) > maxGroupMeanALE &
      !is.na(iDiffAleWideM[,"ale2"]) &
      (!removeAboveAleNum %in% colnames(iDiffAleWideM) |
         length(removeAboveAleNum) == 0 |
         is.na(iDiffAleWideM[,"ale3"]))
      );

   iDiffAleWideUse <- iDiffAleWide[iDiffAleWideMGMkeep,,drop=FALSE];
   iDiffAleWideUse[,aleCols] <- noiseFloor(iDiffAleWideUse[,aleCols,drop=FALSE],
      minimum=maxGroupMeanFloor);
   #minimum=maxGroupMeanALE);

   if (verbose) {
      printDebug("ale2violin(): ",
         "keeping ", formatInt(sum(iDiffAleWideMGMkeep)),
         " out of ", formatInt(length(iDiffAleWideMGMkeep)),
         " rows.");
   }
   if (returnAll) {
      retVals$iDiffAleWideUse <- iDiffAleWideUse;
   }

   ## Create a data.frame from the geneLists
   if (jamba::igrepHas("data.frame", class(geneLists))) {
      if (!all(c("gene_name", "geneList") %in% colnames(geneLists))) {
         stop("ale2violin() requires geneList to have colnames gene_name and geneList.");
      }
      geneListsDF <- geneLists[,c("gene_name","geneList"),drop=FALSE];
   } else {
      geneListsDF <- renameColumn(list2df(geneLists),
         from=c("item","value"),
         to=c("geneList","gene_name"));
      geneListsDF$geneList <- factor(geneListsDF$geneList,
         levels=names(geneLists));
   }

   ## Merge the geneList name for each gene, making one geneList label per row
   ## Also subset iDiffAleWideUse to include genes in geneList
   iDiffAleWide2 <- merge(all.x=TRUE,
      all.y=TRUE,
      x=subset(iDiffAleWideUse, gene_name %in% geneListsDF$gene_name),
      subset(geneListsDF, gene_name %in% iDiffAleWideUse$gene_name));
   if (returnAll) {
      retVals$iDiffAleWide2 <- iDiffAleWide2;
   }

   ## Also filter to make sure gene is present in both "CB" and "DE"
   ## for each geneList
   if (1 == 2) {
      iDiffAleWide2$Region <- gsub("_.+", "", iDiffAleWide2$sampleGroup);
      iDiffAleWide2$CellType <- gsub("^.+_", "", iDiffAleWide2$sampleGroup);
      iGeneListRegionL <- lapply(split(iDiffAleWide2[,c("gene_name","CellType")],
         pasteByRow(iDiffAleWide2[,c("geneList","Region")])), function(i){
            names(tcount(i$gene_name, minCount=2))
         });
      iGeneListRegionDF <- list2df(iGeneListRegionL);
      iGeneListRegionDF$geneList <- gsub("_.+", "", iGeneListRegionDF$item);
      iGeneListRegionDFuniq <- unique(iGeneListRegionDF[,c("name","geneList")]);
      iGeneListRegionDFuniqL <- split(iGeneListRegionDFuniq$name, iGeneListRegionDFuniq$geneList);
   }

   ## Get ale1 groupMean value to use in centering
   #centerAle1Group <- sapply(split(iDiffAleWide2, iDiffAleWide2$gene_name), function(iDF){
   #   mean(unlist(iDF[,"ale1"]), na.rm=TRUE);
   #});


   ## Make a version which centers by ale1
   iDiffAleWide2ctr <- iDiffAleWide2;
   iDiffAleWide2ctr[,aleCols] <- (iDiffAleWide2[,aleCols] - iDiffAleWide2[,"ale1"]);
   ## Remove rows with NA values

   ## Convert wide to tall
   if (returnAll) {
      retVals$iDiffAleWide2ctr <- iDiffAleWide2ctr;
   }
   iDiffAleTall2ctr <- data.table::melt(
      iDiffAleWide2ctr[,c("gene_name","Group",aleCols,"geneList")],
      variable.name="ALE",
      id.vars=c("gene_name", "Group", "geneList"),
      value.name="intensity"
   );
   ## Remove NA rows
   iDiffAleTall2ctr <- subset(iDiffAleTall2ctr, !is.na(intensity));

   ## Labels to include the number of rows and average value
   #iDFsub <- subset(iDiffAleTall2ctr, variable %in% c("ale1","ale2"));
   if (1 == 2) {
      iDFsub <- subset(iDiffAleTall2ctr,
         !ALE %in% "ale1");
      listGroups <- c("Region","geneList","ALE");
      listSizes <- sdim(split(iDFsub,
         pasteByRowOrdered(iDFsub[,listGroups,drop=FALSE],
         #pasteByRow(iDFsub[,c("Group","geneList","ALE")],
            sep="!")));
      listSizesDF <- data.frame(check.names=FALSE,
         listSizes[,c("rows"),drop=FALSE],
         rbindList(strsplit(rownames(listSizes), "!"),
            newColnames=listGroups));
      listSizesDF$geneList <- factor(listSizesDF$geneList,
         levels=levels(iDiffAleTall2ctr$geneList));
   }

   ## Add column indicating the facet_group
   facet_groupsV <- cPasteUnique(split(facet_groups[names(groups)], groups));
   iDiffAleTall2ctr[[facet_name]] <- facet_groupsV[iDiffAleTall2ctr$Group];

   ## Optionally subset data at this step, to ensure geneLists entries
   ## are only represented in relevant groups or facet_groups.
   if (length(subsetFunc) > 0 && is.function(subsetFunc)) {
      iDiffAleTall2ctr <- subsetFunc(iDiffAleTall2ctr);
   }

   ## Violin plot where each gene line is drawn for each sample group
   iDiffAleTall2ctr$Line <- pasteByRowOrdered(iDiffAleTall2ctr[,c("gene_name","Group")]);
   retVals$iDiffAleTall2ctr <- iDiffAleTall2ctr;
   if (make_ggplots) {
      if (verbose) {
         printDebug("ale2violin(): ",
            "Preparing violin including all groups.");
      }
      ggAle1 <- ggplot(iDiffAleTall2ctr,
            aes(y=intensity, x=ALE, fill=geneList)) +
         geom_violin(draw_quantiles=c(0.5), scale="width") +
         #geom_violin(draw_quantiles=c(0.5), scale="count") +
         geom_line(aes(group=Line), alpha=lineAlpha) +
         facet_grid(as.formula(paste0(facet_name, "~geneList"))) +
         ylab("log2 difference from ale1") +
         #ylim(c(-10,10)) +
         theme_jam() + scale_fill_manual(values=colorSubGeneLists);
      retVals$ggAle1 <- ggAle1;
   }

   ## Take the average value per facet_group
   iDiffAleTall2ctrFacet1 <- splicejam::shrinkMatrix(iDiffAleTall2ctr[,"intensity"],
      groupBy=pasteByRow(iDiffAleTall2ctr[,c("gene_name",facet_name,"ALE","geneList")],
         sep="!"),
      shrinkFunc=mean);
   iDF2 <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      rbindList(strsplit(iDiffAleTall2ctrFacet1$groupBy, "!"),
      newColnames=c("gene_name",facet_name,"ALE","geneList")));
   iDF2$geneList <- factor(iDF2$geneList,
      levels=levels(iDiffAleTall2ctr$geneList));
   iDiffAleTall2ctrFacet <- data.frame(iDF2,
      intensity=iDiffAleTall2ctrFacet1$x);
   retVals$iDiffAleTall2ctrFacet <- iDiffAleTall2ctrFacet;

   ## Violin plot where each gene line is drawn using mean per CB and DE
   if (make_ggplots) {
      if (verbose) {
         printDebug("ale2violin(): ",
            "Preparing violin using mean region values.");
      }
      ggAle2 <- ggplot(iDiffAleTall2ctrFacet,
         aes(y=intensity, x=ALE, fill=geneList)) +
         geom_violin(draw_quantiles=c(0.5), scale="width") +
         geom_line(aes(group=gene_name), alpha=lineAlpha) +
         facet_grid(as.formula(paste0(facet_name, "~geneList"))) +
         ylab("log2 mean difference from ale1") +
         xlab("Alternative last 3'UTR") +
         theme_jam() + scale_fill_manual(values=colorSubGeneLists);
      retVals$ggAle2 <- ggAle2;
   }
   return(retVals);
}

#' Perform differential isoform analysis using diffSplice
#'
#' Perform differential isoform analysis using diffSplice
#'
#' This function is intended to be a convenient method
#' to call `limma::diffSplice()` and return the output of
#' `limma::topSplice()` in helpful formats for downstream
#' use.
#'
#' The basic input required:
#'
#' * `iMatrixTx` a numeric matrix of transcript expression
#' * `detectedTx` (optional) subset of `rownames(iMatrixTx)`, sometimes
#'    determined by `defineDetectedTx()`.
#' * `tx2geneDF` data.frame with "transcript_id" and "gene_name" columns.
#'    This data.frame is often the same one used with `tximport::tximport()`
#'    when importing transcriptome data to the gene level.
#' * `iDesign`,`iContrasts` design and contrast matrix, as created by
#'    `groups2contrasts()`.
#'
#' These steps are run in order:
#'
#' * `limma:voom()` (Optional.) This step is enabled with
#'    the argument `useVoom=TRUE` and should only be used when `iMatrixTx`
#'    contains count or pseudocount data. When using TPM or FPKM values,
#'    set `useVoom=FALSE`.
#' * `limma::lmFit()`
#' * `limma::contrasts.fit()`
#' * `limma::diffSplice()`
#' * `limma::topSplice()` This function is called on each contrast
#'       in order to return a list of data.frames.
#'
#' Each contrast is tested for differential transcript expression.
#' No other contrasts are tested.
#'
#' The output is a list with an element `"statsDFs"` that is itself
#' list of data.frames for each contrast. By default when
#' `collapseByGene=TRUE` each row is collapsed to gene level,
#' using the best statistical hit per gene as an exemplar.
#'
#' The rows  of `iMatrixTx` are expected to contain expression values
#' per transcript isoform, but may contain alternative measurements
#' such as: junction counts per gene; exon expression per gene.
#'
#' When `useVoom=TRUE`, the `iMatrixTx` data is exponentiated
#' prior to running `limma::voom()`, for the purpose of
#' calculating a weights matrix.
#' The original `iMatrixTx` data is used in `limma::lmFit()` as-is,
#' alongside the voom weights. The voom-normalized data is not
#' used. Therefore, the input data is assumed to be normalized.
#'
#' Statistical results can be summarized at the gene level, after
#' applying thresholds for statistical hits in the form of
#' required adjusted P-value and/or fold change. It may be helpful
#' to review results per gene, alongside the specific transcript
#' isoforms which are called statistical hits.
#'
#' @return list containing:
#' * `detectedTx`, `detectedTxUse` the detectedTx
#' values representing genes with multiple transcripts;
#' * `fit` the initial model fit;
#' * `fit2` the contrast model fit;
#' * `splice` the output from `limma::diffSplice()`;
#' * `statsDFs` list of data.frame output from `limma::topSplice()`
#' either at the transcript level or the gene level.
#'
#' Note that when `collapseByGene=TRUE` the results will return the first
#' transcript entry per gene that has the best P-value. Often one gene
#' will have two transcripts with identical P-value, and the order
#' that the transcripts appear is inconsistent. Therefore, the direction
#' of fold change is not meaningful by itself, except with respect to
#' the specific isoform returned. See `limma::topSplice()` argument
#' `test` for more information about transcript- and gene-level
#' summaries.
#'
#' @family jam RNA-seq functions
#'
#' @param iMatrixTx numeric matrix of expression, with transcripts as
#'    rows and samples as columns. The data is assumed to be log2-transformed
#'    using the format `log2(1 + x)`. This data should be normalized using
#'    appropriate methods, outside the scope of this function.
#' @param detectedTx character vector containing all or a subset of
#'    `rownames(iMatrix)` used for statistical testing.
#' @param tx2geneDF data.frame with colnames `c(txColname, geneColname)`,
#'    where all entries of `rownames(iMatrix)` are represented
#'    in `tx2geneDF[,txColname]`.
#' @param txColname,geneColname the `colnames(tx2geneDF)` representing
#'    the `rownames(iMatrixTx)` matched by `tx2geneDF[,txColname]`,
#'    and the associated genes given by `tx2geneDF[,geneColname]`.
#'    Note that `detectedTx` must also contain values in `rownames(iMatrixTx)`
#'    and `tx2geneDF[,txColname]`.
#' @param iDesign numeric matrix representing the design matrix for
#'    the experiment design. For example, `limma::model.matrix(~0+group)`
#'    will represent each group. Typically, `rownames(iDesign)` should
#'    be defined to match the `colnames(iMatrix)` even if it requires
#'    extra processing. The `colnames(iDesign)` should represent
#'    group names used in `rownames(iContrasts)`.
#' @param iContrasts numeric matrix representing the contrasts used
#'    in statistical comparisons. This matrix can be generated by
#'    running `limma::makeContrasts()` using a format similar to the
#'    following: `limma::makeContrasts(contrasts="group1-group2", levels=iDesign)`,
#'    see "Examples" for more info.
#' @param cutoffFDR numeric value indicating a statistical threshold
#'    on the FDR (adjusted P-value). Values should be between 0 and 1, where
#'    `cutoffFDR=1` would impose no threshold on the adjusted P-value.
#' @param cutoffFold numeric value indicating the minimum normal space
#'    fold change allowed for statistical hits. For example `cutoffFold=2`
#'    would require a 2-fold change, equivalent to log2 fold change >= 1.
#' @param collapseByGene logical indicating whether results should be
#'    summarized at the gene level after filtering statistical hits.
#' @param spliceTest character value described in `limma::topSplice()`
#'    which defines the statistical test to return. The default `"t"`
#'    returns the t-test result for each isoform, which is mainly beneficial
#'    because it also includes fold change that can be filtered. The
#'    `"F"` returns F-test per gene, and `"simes"` returns the per-gene
#'    t-test P-value after Simes adjustment per gene.
#' @param sep character value used as a delimiter in output data.frame
#'    colnames, such that each stats is followed by the contrast name,
#'    separated by this delimiter.
#' @param useVoom logical indicating whether to apply the `limma::voom()`
#'    adjustment prior to running `limma::diffSplice()`. This value
#'    should be `TRUE` when analyzing count or pseudocount data.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @references
#' Law, CW, Chen, Y, Shi, W, and Smyth, GK (2014).
#' Voom: precision weights unlock linear model analysis
#' tools for RNA-seq read counts.
#' *Genome Biology* **15**, R29.
#'
#' @examples
#' # example for defining iDesign
#' # first define a vector of sample groups
#' iGroups <- jamba::nameVector(paste(rep(c("WT", "KO"), each=6),
#'    rep(c("Control", "Treated"), each=3),
#'    sep="_"));
#' iGroups <- factor(iGroups, levels=unique(iGroups));
#' iGroups;
#' # next define sample identifiers
#' iSamples <- names(iGroups);
#'
#' # given a vector of groups, make iDesign
#' iDesign <- stats::model.matrix(~0+iGroups);
#'
#' # It is good practice to rename colnames(iDesign) and rownames(iDesign),
#' # that is, for the love of all that is good, use colnames and rownames
#' # that help confirm that these matrices are consistent.
#' colnames(iDesign) <- levels(iGroups);
#' rownames(iDesign) <- names(iGroups);
#' iDesign;
#'
#' # define contrasts
#' # the example below includes a two-way contrast, which is a test
#' # of the pairwise fold changes
#' iContrasts <- limma::makeContrasts(contrasts=c(
#'    "WT_Treated-WT_Control", "KO_Treated-KO_Control",
#'    "(KO_Treated-KO_Control)-(WT_Treated-WT_Control)"),
#'    levels=iDesign);
#' iContrasts;
#'
#' # for validation, verify these constraints:
#' # - all(rownames(iDesign) == iSamples)
#' # - colnames(iDesign) == the actual group names
#' # - all(colnames(iDesign) == rownames(iContrasts))
#' # you can see which samples are included in each test with crossproduct:
#' iDesign %*% iContrasts;
#'
#' ## Another efficient way to define iDesign and iContrasts:
#' iDC <- splicejam::groups2contrasts(iGroups, returnDesign=TRUE);
#' iDesign <- iDC$iDesign;
#' iContrasts <- iDC$iContrasts;
#'
#' @export
runDiffSplice <- function
(iMatrixTx,
 detectedTx=rownames(iMatrixTx),
 tx2geneDF,
 txColname="transcript_id",
 geneColname="gene_name",
 iDesign,
 iContrasts,
 cutoffFDR=0.05,
 cutoffFold=1.5,
 collapseByGene=TRUE,
 spliceTest=c("t", "simes", "F"),
 sep=" ",
 useVoom=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to provide a wrapper around the steps required to run
   ## limma::diffSplice() and the associated results export steps
   retVals <- list();

   ## Validate the input data
   if (!all(colnames(iMatrixTx) %in% rownames(iDesign))) {
      stop("runDiffSplice() requires all(colnames(iMatrixTx) %in% rownames(iDesign))");
   }
   if (!all(colnames(iDesign) %in% rownames(iContrasts))) {
      stop("runDiffSplice() requires all(colnames(iDesign) %in% rownames(iContrasts))");
   }
   spliceTest <- match.arg(spliceTest);

   ########################################################
   ## starting with detected transcripts
   ## determine which genes have multiple transcripts
   if (length(detectedTx) == 0) {
      if (verbose) {
         printDebug("runDiffSplice(): ",
            "Using all rownames(iMatrixTx) as ",
            "detectedTx");
      }
      detectedTx <- rownames(iMatrixTx);
   }
   ## Match detectedTx to the tx2geneDF data.frame,
   ## then drop NA unmatched values
   iTxMatch <- jamba::rmNA(match(detectedTx,
      tx2geneDF[,txColname]));
   ## Count the total unique genes represented
   iTxGeneCt <- length(unique(tx2geneDF[iTxMatch,geneColname]));
   ## Look for genes represented more than once
   iGeneCt2 <- names(jamba::tcount(tx2geneDF[iTxMatch,geneColname],
      minCount=2));

   detectedTxUse <- subset(tx2geneDF[iTxMatch,,drop=FALSE],
      tx2geneDF[iTxMatch,geneColname] %in% iGeneCt2)[[txColname]];
   iTxMatchUse <- match(detectedTxUse,
      tx2geneDF[[txColname]]);
   retVals$detectedTx <- detectedTx;
   retVals$detectedTxUse <- detectedTxUse;

   if (verbose) {
      printDebug("runDiffSplice(): ",
         "Started with ",
         formatInt(length(detectedTx)),
         " entries representing ",
         formatInt(iTxGeneCt),
         " genes.");
      printDebug("runDiffSplice(): ",
         "  Ended with ",
         formatInt(length(detectedTxUse)),
         " entries representing ",
         formatInt(length(iGeneCt2)),
         " multi-entry genes.");
   }


   ########################################################
   ## Define ExpressionSet
   iMatrixTxES <- ExpressionSet(
      assayData=2^iMatrixTx[detectedTxUse,,drop=FALSE]-1,
      featureData=new("AnnotatedDataFrame",
         data=data.frame(
            check.names=FALSE,
            stringsAsFactors=FALSE,
            probes=detectedTxUse,
            tx2geneDF[iTxMatchUse,geneColname,drop=FALSE],
            row.names=detectedTxUse)));

   ########################################################
   ## Ensure iDesign is consistently ordered with iMatrixTx
   iDesign <- iDesign[colnames(iMatrixTx),,drop=FALSE];
   ## Ensure that iContrasts is consistently ordered with iDesign
   iContrasts <- iContrasts[colnames(iDesign),,drop=FALSE];


   ########################################################
   ## Optionally run voom prior to the initial model fit
   if (useVoom) {
      if (verbose) {
         printDebug("runDiffSplice(): ",
            "Running voom().");
      }
      v <- voom(iMatrixTxES,
         design=iDesign,
         plot=FALSE);
      if (verbose) {
         printDebug("runDiffSplice(): ",
            "Running lmFit().");
      }
      #fit <- lmFit(v,
      #   design=iDesign);
      exprs(iMatrixTxES) <- log2(1+exprs(iMatrixTxES));
      fit <- lmFit(iMatrixTxES,
         weights=v$weights,
         design=iDesign);
   } else {
      exprs(iMatrixTxES) <- log2(1+exprs(iMatrixTxES));
      if (verbose) {
         printDebug("runDiffSplice(): ",
            "Running lmFit().");
      }
      fit <- lmFit(iMatrixTxES,
         design=iDesign);
   }
   ## Add transcript back to the output data
   fit$genes[,txColname] <- fit$genes[,"probes"];
   retVals$fit <- fit;

   ########################################################
   ## Fit contrasts
   if (verbose) {
      printDebug("runDiffSplice(): ",
         "Running contrasts.fit().");
   }
   fit2 <- contrasts.fit(fit,
      iContrasts);
   retVals$fit2 <- fit2;

   ########################################################
   ## Run diffSplice()
   if (verbose) {
      printDebug("runDiffSplice(): ",
         "Running diffSplice().");
   }
   splice <- diffSplice(fit2,
      geneid=geneColname,
      exonid=txColname);
   retVals$splice <- splice;

   ## Clean up each contrast into a data.frame
   if (verbose) {
      printDebug("runDiffSplice(): ",
         "Preparing statsDFs.");
   }
   iCoefs <- colnames(coefficients(splice));
   statsDFs <- lapply(nameVector(iCoefs), function(iCoef){
      iCoefGroups <- unlist(strsplit(gsub("^[(]|[)]$", "", iCoef), "[-()]+"));
      iGroupMeans <- coefficients(fit)[,iCoefGroups,drop=FALSE];
      iGroupMeansDF <- data.frame(
         setNames(data.frame(rownames(iGroupMeans)), txColname),
         renameColumn(iGroupMeans,
            from=colnames(iGroupMeans),
            to=paste("groupMean",
               colnames(iGroupMeans),
               sep=sep)),
         check.names=FALSE,
         stringsAsFactors=FALSE);
      spliceDF <- limma::topSplice(splice,
         coef=iCoef,
         test=spliceTest,
         #test="t",#"simes","F"
         n=Inf);
      if (verbose) {
         printDebug("runDiffSplice(): ",
            "head(spliceDF):");
         print(head(spliceDF));
      }
      if (collapseByGene) {
         pvColname <- head(vigrep("^P.Value", colnames(spliceDF)), 1);
         bestPbyGene <- shrinkMatrix(spliceDF[,pvColname],
            groupBy=spliceDF[,geneColname],
            shrinkFunc=min,
            returnClass="matrix");
         countTxByGeneM <- shrinkMatrix(data.frame(numTx=spliceDF[,"probes"]),
            groupBy=spliceDF[,geneColname],
            shrinkFunc=length,
            returnClass="matrix");

         spliceDF[,"best P.Value"] <- bestPbyGene[spliceDF[,geneColname],1];
         spliceDFsub <- subset(spliceDF,
            `P.Value` == `best P.Value`);
         iMatchUniqGene <- match(unique(spliceDFsub[,geneColname]),
            spliceDFsub[,geneColname]);
         keepCols <- setdiff(colnames(spliceDFsub), "probes");
         spliceDFuse <- spliceDFsub[iMatchUniqGene,keepCols,drop=FALSE];
         spliceDFuse$numTx <- countTxByGeneM[match(spliceDFuse[,geneColname],
            rownames(countTxByGeneM)),"numTx"];
         iMatch2 <- match(spliceDFuse[,txColname],
            rownames(iGroupMeansDF));
         spliceDFuse[,colnames(iGroupMeansDF)] <- iGroupMeansDF[iMatch2,,drop=FALSE];
      } else {
         spliceDFuse <- spliceDF;
      }
      ## Add a hit flag
      spliceDFuse[,"hit"] <- (spliceDFuse[,"FDR"] <= cutoffFDR &
         (spliceTest %in% c("F","simes") |
            abs(spliceDFuse[,"logFC"]) >= log2(cutoffFold) )
         )*1;

      ## Rename columns to include the contrast
      renameCols <- c("logFC", "t", "P.Value", "FDR",
         "best P.Value", txColname, "hit");
      spliceDFuse <- jamba::renameColumn(spliceDFuse,
         from=renameCols,
         to=paste(renameCols, iCoef));
      jamba::mixedSortDF(spliceDFuse,
         byCols=provigrep(c("FDR","P.Value"),
            colnames(spliceDFuse)));
   });
   retVals$statsDFs <- statsDFs;

   return(retVals);
}

#' Get gaps in GRanges
#'
#' Get gaps in GRanges
#'
#' This function returns the gaps between GRanges regions, calculated
#' for each chromosome (using `GenomicRanges::seqnames(gr)`), and when `strandSpecific=TRUE`
#' it determines gaps in stranded fashion.
#'
#' @family jam GRanges functions
#'
#' @param gr GRanges object
#' @param strandSpecific logical indicating whether to convert strand
#'    to `"*"` prior to determining gaps between features.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to `getGRLgaps()`.
#'
#' @export
getGRgaps <- function
(gr,
 strandSpecific=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to wrapper the gaps() function from GenomicRanges
   ## except to return only the gaps between features on the same
   ## chromosome and strand
   ##
   ## keepValues=TRUE will keep values(GR) colnames, but will fill with NA
   ## so the resulting object can be appended to the original GR.
   ##
   #GRDF <- as.data.frame(GR);
   ##
   ## This function is essentially a wrapper around getGRLgaps()
   #if (!strandSpecific) {
   #   strand(gr) <- "*";
   #}
   #grl <- GRangesList(split(gr,
   #   pasteByRowOrdered(as.data.frame(gr)[,c("seqnames", "strand")])));
   gapsGRL <- getGRLgaps(grl=GenomicRanges::GRangesList(list(`gr`=gr)),
      strandSpecific=strandSpecific,
      verbose=verbose,
      ...);
   return(gapsGRL@unlistData);
}

#' Get gaps in GRangesList objects
#'
#' Get gaps in GRangesList objects
#'
#' This function returns gaps between GRanges regions in a GRangesList
#' object. When `strandSpecific=TRUE` is determines gaps per strand,
#' otherwise strands are converted to `"*"`. It will also determine
#' gaps within chromosome for each GRanges entry in GRangesList.
#'
#' @family jam GRanges functions
#'
#' @return GRangesList object with gaps for each chromosome and strand
#'    present in each GRanges entry. It does not return gap sequence
#'    at the edges of GRanges regions to the chromosome ends.
#'
#' @param grl GRangesList object.
#' @param strandSpecific logical indicating whether to determine gaps
#'    within strand.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
getGRLgaps <- function
(grl,
 strandSpecific=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to wrapper the gaps() function from GenomicRanges
   ## except to return only the gaps between features on the same
   ## chromosome and strand
   ##
   ## Check for one strand per grl
   ## grl <- grlGria1;strand(grl@unlistData)[3:4] <- "-";
   ## seqnames(grl@unlistData)[3:4] <- "chr7";
   if (!strandSpecific) {
      strand(grl) <- "*";
   }
   isMultiStrand <- any(lengths(unique(GenomicRanges::strand(grl))) > 1);
   isMultiSeqname <- any(lengths(unique(GenomicRanges::seqnames(grl))) > 1);

   if (length(names(grl)) == 0) {
      names(grl) <- jamba::makeNames(rep("grl", length(grl)));
   }
   grlNames <- names(grl);

   ## pre-process GRangesList to split each GRanges by seqnames_strand
   ## TODO: check whether any GRanges entry has multiple seqnames or
   ## multiple strands -- if not then skip this step.
   ## If so, then split each GRanges by seqnames_strand
   if (isMultiStrand || isMultiSeqname) {
      GenomicRanges::values(grl@unlistData)[,"grl_name"] <- rep(names(grl),
         S4Vectors::elementNROWS(grl));
      grl <- GenomicRanges::GRangesList(GenomicRanges::split(grl@unlistData,
         pasteByRowOrdered(sep=":!:",
            as.data.frame(grl@unlistData)[,c("grl_name","seqnames","strand")])
         )
      );
   }

   ## Use gaps() method directly
   if (verbose) {
      printDebug("getGRLgaps(): ",
         "Began gaps logic.");
   }
   IRL <- as(grl, "IRangesList");
   gapsIRL <- IRanges::gaps(IRL);
   irlName <- rep(names(gapsIRL),
      S4Vectors::elementNROWS(gapsIRL));
   grlName <- gsub(":!:.*$", "", irlName);
   grlName <- factor(grlName,
      levels=unique(c(grlNames, grlName)));
   grNew <- GRanges(seqnames=rep(as.character(GenomicRanges::seqnames(range(grl))),
      S4Vectors::elementNROWS(gapsIRL)),
      range=IRanges::IRanges(start=IRanges::start(gapsIRL@unlistData),
         end=IRanges::end(gapsIRL@unlistData)),
      strand=rep(GenomicRanges::strand(range(grl)@unlistData),
         S4Vectors::elementNROWS(gapsIRL)),
      irl_name=irlName,
      grl_name=grlName
   );
   ## Split by the original grl name, which will combine different
   ## seqnames and strands if needed
   grlNew <- GenomicRanges::split(grNew[,0],
      GenomicRanges::values(grNew)[,"grl_name"]);
   keepColnames <- setdiff(colnames(GenomicRanges::values(grlNew@unlistData)), "grl_name");
   GenomicRanges::values(grlNew@unlistData) <- GenomicRanges::values(grlNew@unlistData)[,keepColnames];
   return(grlNew);
}

#' Flatten exons by gene or transcript
#'
#' Flatten exons by gene or transcript
#'
#' This function takes as input:
#'
#'  * `exonsByTx` as a `GRangesList` object
#' of transcript exons named by the `transcript_id`,
#' * `tx2geneDF` a `data.frame` with transcript-gene cross-reference,
#' * `detectedTx` an optional character vector of `transcript_id`
#' values, used to subset the overall transcripts
#' * `cdsByTx` an optional `GRangesList` object, similar to `exonsByTx`
#' except that it only contains the CDS portion of exons
#'
#' This function groups exons together by gene, producing a flattened,
#' disjoint (non-overlapping) set of exons, where exons are subdivided
#' when there are multiple boundaries.
#'
#' Finally, it labels each exon using a defined naming scheme:
#'
#' * Each contiguous exon is numbered in order, starting at `1` for the
#' first stranded exon for the gene, for example `exon1`, `exon2`,
#' `exon3`.
#' * When an exon is sub-divided, each section is labeled with
#' an alphabetic suffix to indicate the order within that exon,
#' for example `exon1a`, `exon1b`, `exon1c`.
#'
#' A text schematic is shown below:
#'
#' `|=======|======|......|=======|.....|=======|=======|=======|`
#'
#' `|.exon1a|exon1b|......|.exon2.|.....|.exon3a|.exon3b|.exon3c|`
#'
#' Where
#'
#' * `|====|` represents an exon,
#' * `|====|====|` represents one contiguous exon with two
#' sub-divided parts, and
#' * `.....` represents an intron.
#'
#' It is recommended but not required to supply `detectedTx`,
#' since it can greatly reduce the total number of transcripts.
#' This step has two benefits: It can greatly simplify the
#' resulting exon model based upon observed data, and it
#' has the by-product of removing potentially erroneous
#' transcripts which often has no supporting observed data.
#'
#' @family jam RNA-seq functions
#' @family GRanges functions
#'
#' @return GRangesList named by gene when `by="gene"` or transcript when
#'    `by="tx"`, containing non-overlapping GRanges with exon names
#'    as described above.
#'
#' @param exonsByTx GRangesList named by transcript, containing one or
#'    more GRanges representing exons. This data is often produced
#'    from `TxDb` data using `GenomicFeatures::exonsBy(...,by="tx")`.
#' @param tx2geneDF data.frame containing at least two columns with
#'    transcript and gene annotation, whose colnames are defined by
#'    arguments `txColname` and `geneColname` respectively. When
#'    using a GTF file, `makeTx2geneFromGtf()` can be used to
#'    create a `tx2geneDF` in `data.frame` format.
#' @param by character string to group exons, `"gene"` groups multiple
#'    transcripts per gene, and `"tx"` groups exons per transcript.
#'    Note that in both cases, it combines `exonsByTx` and `cdsByTx`
#'    when `cdsByTx` is also supplied.
#' @param detectedTx character vector of detected transcripts, used to
#'    subset the overall transcripts prior to producing a flattened gene
#'    exon model.
#' @param genes optional character vector, representing a subset of
#'    genes for which flattened exons will be prepared. This argument
#'    is useful when focusing on only one or a subset of genes.
#' @param txColname character string indicating a column from
#'    `colnames(tx2geneDF)` used to identify transcripts.
#' @param geneColname character string indicating a column from
#'    `colnames(tx2geneDF)` used to identify gene name, or gene symbol.
#' @param cdsByTx `GRangesList` named by transcript, containing `GRanges`
#'    exons that only include CDS regions. This data is often produced
#'    from `TxDb` data using `GenomicFeatures::cdsBy(...,by="tx")`.
#' @param cdsByGene `GRangesList` named by gene, containing `GRanges`
#'    exons that only include CDS regions. This data is often produced
#'    from `TxDb` data using `GenomicFeatures::cdsBy(...,by="gene")`.
#'    Note this input is only used when `by="gene"`.
#' @param verbose logical indicating whether to print verbose output.
#'
#' @export
flattenExonsBy <- function
(exonsByTx,
 tx2geneDF,
 by=c("gene", "tx"),
 detectedTx=NULL,
 genes=NULL,
 txColname="transcript_id",
 geneColname="gene_name",
 cdsByTx=NULL,
 cdsByGene=NULL,
 verbose=FALSE)
{
   ##
   if (!suppressPackageStartupMessages(require(GenomicRanges))) {
      stop("The GenomicRanges package is required.");
   }
   if (!igrepHas("data.frame|tibble|tbl|dataframe", class(tx2geneDF))) {
      stop("tx2geneDF must be a form of data.frame or related class.");
   }
   if (!all(c(txColname, geneColname) %in% colnames(tx2geneDF))) {
      stop("colnames(tx2geneDF) must contain txColname and geneColname.");
   }
   by <- match.arg(by);
   if (length(detectedTx) > 0) {
      iTxs <- intersect(
         c(names(exonsByTx),
            names(cdsByTx)),
         detectedTx);
   } else {
      iTxs <- unique(c(names(exonsByTx),
         names(cdsByTx)));
   }
   ## Optionally subset tx2geneDF by genes
   if (length(genes) > 0) {
      tx2geneDF <- subset(tx2geneDF,
         tx2geneDF[[geneColname]] %in% genes);
      if (nrow(tx2geneDF) == 0) {
         stop("tx2geneDF[[geneColname]] contains no values matching the supplied genes.");
      }
   }
   ## Validate iTxs in tx2geneDF
   iTxs <- intersect(iTxs,
      tx2geneDF[[txColname]]);
   tx2geneDF <- subset(tx2geneDF,
      tx2geneDF[[txColname]] %in% iTxs);
   if (length(iTxs) == 0) {
      stop("There are no Tx entries shared by: names(exonsByTx), tx2geneDF[,txColname], detectedTx.");
   }
   if (length(exonsByTx) > 0) {
      exonsByTx <- exonsByTx[names(exonsByTx) %in% tx2geneDF[[txColname]]];
   }
   if (length(cdsByTx) > 0) {
      cdsByTx <- cdsByTx[names(cdsByTx) %in% tx2geneDF[[txColname]]];
   }
   if (verbose) {
      printDebug("flattenExonsBy(): ",
         "Flattening ", length(iTxs), " transcripts from ",
         length(unique(tx2geneDF[[geneColname]])),
         " unique genes.");
   }

   ## Subset exonsByTx and add gene annotations
   iTxExonsGRL <- exonsByTx[iTxs];
   iTxMatch <- match(names(iTxExonsGRL),
      tx2geneDF[[txColname]]);
   GenomicRanges::values(iTxExonsGRL@unlistData)[,geneColname] <- rep(
      as.character(tx2geneDF[iTxMatch,geneColname]),
      S4Vectors::elementNROWS(iTxExonsGRL));

   ## split exons by gene
   if (verbose) {
      printDebug("flattenExonsBy(): ",
         "Splitting tx exons by gene.");
   }
   if ("gene" %in% by) {
      exonsByGene <- GenomicRanges::GRangesList(
         GenomicRanges::split(
            iTxExonsGRL@unlistData,
            GenomicRanges::values(iTxExonsGRL@unlistData)[[geneColname]])
         );
   } else {
      GenomicRanges::values(iTxExonsGRL@unlistData)[,txColname] <- rep(
         names(iTxExonsGRL),
         S4Vectors::elementNROWS(iTxExonsGRL));
      exonsByGene <- iTxExonsGRL[,c(txColname,geneColname)];
   }


   ## Disjoin exons within each gene GRL
   if ("gene" %in% by) {
      if (verbose) {
         printDebug("flattenExonsBy(): ",
            "Preparing disjoint gene exons.");
      }
      iGeneExonsDisGRL <- GenomicRanges::disjoin(exonsByGene);
   } else {
      iGeneExonsDisGRL <- exonsByGene;
   }
   if (verbose) {
      printDebug("flattenExonsBy(): ",
         "Completed disjoint gene exons.");
   }
   ## Add gene annotation to each entry
   if ("gene" %in% by) {
      GenomicRanges::values(iGeneExonsDisGRL@unlistData)[,geneColname] <- rep(
         names(iGeneExonsDisGRL),
         elementNROWS(iGeneExonsDisGRL));
   #} else {
      #GenomicRanges::values(iGeneExonsDisGRL@unlistData)[,txColname] <- rep(
      #   names(iGeneExonsDisGRL),
      #   elementNROWS(iGeneExonsDisGRL));
      #txMatch <- match(names(iGeneExonsDisGRL),
      #   tx2geneDF[[txColname]]);
      #GenomicRanges::values(iGeneExonsDisGRL)[,geneColname] <- tx2geneDF[txMatch, geneColname]
   }

   ## Optionally subdivide by CDS boundary if supplied
   if (length(cdsByTx) > 0) {
      if (verbose) {
         printDebug("flattenExonsBy(): ",
            "Creating cdsByGene from cdsByTx.");
      }
      cdsByTx <- cdsByTx[names(cdsByTx) %in% iTxs];
      if (length(cdsByTx) > 0) {
         if (!geneColname %in% colnames(GenomicRanges::values(cdsByTx))) {
            txMatch <- match(names(cdsByTx), tx2geneDF[[txColname]]);
            GenomicRanges::values(cdsByTx@unlistData)[,geneColname] <- rep(
               tx2geneDF[txMatch,geneColname],
               elementNROWS(cdsByTx));
         }
         if ("gene" %in% by) {
            cdsByGene <- GenomicRanges::reduce(GenomicRanges::GRangesList(
               GenomicRanges::split(cdsByTx@unlistData,
                  GenomicRanges::values(cdsByTx@unlistData)[[geneColname]])));
         } else {
            cdsByGene <- cdsByTx;
            GenomicRanges::values(cdsByGene@unlistData)[,txColname] <- rep(
               names(cdsByGene),
               elementNROWS(cdsByGene)
            );
         }
         if (verbose) {
            printDebug("flattenExonsBy(): ",
               "length(cdsByGene):",
               length(cdsByGene));
         }
      }
   }
   if (length(cdsByGene) > 0 && any(names(iGeneExonsDisGRL) %in% names(cdsByGene))) {
      if (verbose) {
         printDebug("flattenExonsBy(): ",
            "Adding CDS exon boundary information.");
      }
      cdsByGene <- cdsByGene[names(cdsByGene) %in% names(iGeneExonsDisGRL)];
      ## Use subset of exonsByGene that have cds exons
      exonsByGeneSub <- iGeneExonsDisGRL[names(cdsByGene)];
      exonsByGeneSubCds <- GenomicRanges::intersect(exonsByGeneSub, cdsByGene);
      GenomicRanges::values(exonsByGeneSubCds@unlistData)[,"subclass"] <- "cds";
      if ("gene" %in% by) {
         GenomicRanges::values(exonsByGeneSubCds@unlistData)[,geneColname] <- rep(
            names(exonsByGeneSubCds),
            S4Vectors::elementNROWS(exonsByGeneSubCds));
         exonsByGeneCds <- sort(GenomicRanges::disjoin(GenomicRanges::GRangesList(
            GenomicRanges::split(
            c(exonsByGeneSub@unlistData[,geneColname],
               exonsByGeneSubCds@unlistData[,geneColname]),
            c(GenomicRanges::values(exonsByGeneSub@unlistData)[,geneColname],
               GenomicRanges::values(exonsByGeneSubCds@unlistData)[,geneColname])))));
         GenomicRanges::values(exonsByGeneCds@unlistData)[,geneColname] <- rep(
            names(exonsByGeneCds),
            S4Vectors::elementNROWS(exonsByGeneCds));
      } else {
         GenomicRanges::values(exonsByGeneSubCds@unlistData)[,txColname] <- rep(
            names(exonsByGeneSubCds),
            S4Vectors::elementNROWS(exonsByGeneSubCds));
         #txMatch <- match(names(exonsByGeneSubCds),
         #   tx2geneDF[[txColname]]);
         #values(exonsByGeneSubCds@unlistData)[,geneColname] <- rep(
         #   tx2geneDF[txMatch, geneColname],
         #   elementNROWS(exonsByGeneSubCds)
         #);
         if (verbose) {
            printDebug("flattenExonsBy(): ",
               "split()");
         }
         exonsByGeneCds <- GenomicRanges::split(
               c(exonsByGeneSub@unlistData,
                  exonsByGeneSubCds@unlistData),
               c(GenomicRanges::values(exonsByGeneSub@unlistData)[,txColname],
                  GenomicRanges::values(exonsByGeneSubCds@unlistData)[,txColname]));
         if (verbose) {
            printDebug("flattenExonsBy(): ",
               "disjoin()");
         }
         exonsByGeneCds <- GenomicRanges::disjoin(exonsByGeneCds);
         if (verbose) {
            printDebug("flattenExonsBy(): ",
               "Sorting.");
         }
         exonsByGeneCds <- GenomicRanges::sort(exonsByGeneCds);
         if (verbose) {
            printDebug("flattenExonsBy(): ",
               "Adding txColname,geneColname to disjoint tx exons.");
         }
         GenomicRanges::values(exonsByGeneCds@unlistData)[,txColname] <- rep(
            names(exonsByGeneCds),
            S4Vectors::elementNROWS(exonsByGeneCds));
         txMatch <- match(names(exonsByGeneCds),
            tx2geneDF[[txColname]]);
         GenomicRanges::values(exonsByGeneCds@unlistData)[,geneColname] <- rep(
            tx2geneDF[txMatch, geneColname],
            S4Vectors::elementNROWS(exonsByGeneCds)
         );
      }
      if (verbose) {
         printDebug("flattenExonsBy(): ",
            "Running annotateGRLfromGRL on CDS disjoint exons.");
      }
      exonsByGeneCds <- annotateGRLfromGRL(exonsByGeneCds,
         exonsByGeneSubCds[,"subclass"]);
      naClass <- is.na(GenomicRanges::values(exonsByGeneCds@unlistData)[,"subclass"]);
      GenomicRanges::values(exonsByGeneCds@unlistData)[naClass,"subclass"] <- "noncds";
      GenomicRanges::values(iGeneExonsDisGRL@unlistData)[,"subclass"] <- "noncds";
      iGeneExonsDisGRL[names(exonsByGeneCds)] <- exonsByGeneCds[,c(geneColname, "subclass")];
   }

   ## Assign exon names and numbers
   if (verbose) {
      printDebug("flattenExonsBy(): ",
         "Assigning exon labels to disjoint gene exons.");
   }
   if ("gene" %in% by) {
      iGeneExonsDisGRL <- assignGRLexonNames(iGeneExonsDisGRL,
         geneSymbolColname=geneColname,
         verbose=FALSE);
      GenomicRanges::values(iGeneExonsDisGRL)[,geneColname] <- names(iGeneExonsDisGRL);
   } else {
      iGeneExonsDisGRL <- assignGRLexonNames(iGeneExonsDisGRL,
         geneSymbolColname=txColname,
         verbose=FALSE);
      GenomicRanges::values(iGeneExonsDisGRL)[,txColname] <- names(iGeneExonsDisGRL);
      GenomicRanges::values(iGeneExonsDisGRL@unlistData)[,txColname] <- rep(
         names(iGeneExonsDisGRL),
         S4Vectors::elementNROWS(iGeneExonsDisGRL)
      )
      txMatch <- match(names(iGeneExonsDisGRL),
         tx2geneDF[[txColname]]);
      GenomicRanges::values(iGeneExonsDisGRL)[,geneColname] <- tx2geneDF[txMatch, geneColname];
      GenomicRanges::values(iGeneExonsDisGRL@unlistData)[,geneColname] <- rep(
         GenomicRanges::values(iGeneExonsDisGRL)[,geneColname],
         S4Vectors::elementNROWS(iGeneExonsDisGRL)
      )
   }
   GenomicRanges::values(iGeneExonsDisGRL@unlistData)[,"feature_type"] <- "exon";
   ## TODO: optionally add intron regions between exons of each gene

   return(iGeneExonsDisGRL);
}

#' Add gaps between GRanges regions
#'
#' Add gaps between GRanges regions
#'
#' This function adds gaps between each GRanges region where
#' there is a gap between two GRanges for
#' the same seqnames. When `strandSpecific=TRUE` the gaps are
#' determined per strand.
#'
#' This function is a wrapper around `getGRgaps()`, which is then
#' concatenated to the input `gr` GRanges object using `base::c()`.
#' When the input `gr` has column `S4Vectors::values()` then the
#' gaps GRanges object will have `NA` values used by default. To supply
#' values, use the `newValues` argument, which assigns name-value pairs.
#'
#' @family jam GRanges functions
#'
#' @return GRanges object, sorted when `doSort=TRUE`. When `newValues`
#'    is supplied, the values for gaps GRanges elements will be assigned,
#'    otherwise any column values present in `gr` will be `NA` for
#'    gaps elements. The names of gaps elements are assigned using
#'    `gapname` then are made unique using `jamba::makeNames()`,
#'    unless `gapname is NULL`.
#'
#' @param gr GRanges object
#' @param strandSpecific logical indicating whether the gaps are calculated
#'    per strand, see `getGRgaps()`.
#' @param gapname,suffix character vector supplying the name to assign to new
#'    gap GRanges elements, using `jamba::makeNames()` with `suffix` as
#'    described to define non-duplicated names. If `gapname is NULL` then
#'    no names are assigned to new gap GRanges entries, however when the
#'    input `gr` GRanges object has names, the concatenation of gaps
#'    causes names `""` to be assigned to all gap GRanges elements, which
#'    are duplicated for multiple gaps.
#' @param newValues list of values to add to the resulting gap GRanges,
#'    whose names become `colnames(gr)`, and whose values are used
#'    to populate each column. By default a colname `"feature_type"` is
#'    added, with value `"gap"` added to each row. When `newValues is NULL`
#'    then no values are added to the gaps GRanges.
#' @param doSort logical indicating whether to sort the resulting
#'    GRanges object. When `doSort=FALSE` the gaps are added to the end
#'    of the `gr` input GRanges object.
#' @param ... additional arguments are passed to `getGRgaps()`.
#'
#' @examples
#' gr <- GenomicRanges::GRanges(seqnames=rep(c("chr1","chr2"), c(3,2)),
#'    ranges=IRanges::IRanges(start=c(100, 300, 400, 300, 700),
#'       end=c(199, 450, 500, 600, 800)),
#'    strand=rep(c("+","-"), c(3,2)));
#' gr;
#' getGRLgaps(GenomicRanges::split(gr, GenomicRanges::seqnames(gr)))
#' getGRgaps(gr);
#'
#' @export
addGRgaps <- function
(gr,
 strandSpecific=TRUE,
 gapname="gap",
 suffix="_v",
 newValues=list(feature_type="gap"),
 default_feature_type="exon",
 feature_type_colname="feature_type",
 doSort=TRUE,
 ...)
{
   ## Purpose is to add gaps as GRanges elements between the
   ## GRanges elements given
   if (length(feature_type_colname) > 0 &&
         !feature_type_colname %in% colnames(GenomicRanges::values(gr))) {
      if (length(default_feature_type) == 0) {
         default_feature_type <- c("exon");
      } else {
         default_feature_type <- head(default_feature_type, 1);
      }
      GenomicRanges::values(gr)[[feature_type_colname]] <- default_feature_type;
   }
   grGaps <- getGRgaps(gr,
      strandSpecific=strandSpecific,
      ...);
   if (length(gapname) > 0) {
      gapnames <- jamba::makeNames(
         rep(gapname,
            length.out=length(grGaps)),
         suffix=suffix);
      names(grGaps) <- gapnames;
   }
   if (length(newValues) > 0) {
      for (newValueName in names(newValues)) {
         GenomicRanges::values(grGaps)[[newValueName]] <- rep(newValues[[newValueName]],
            length(grGaps));
      }
   }
   if (doSort) {
      grNew <- sort(c(gr, grGaps));
   } else {
      grNew <- c(gr, grGaps);
   }
   return(grNew);
}

#' Add gaps between GRangesList regions
#'
#' Add gaps between GRangesList regions
#'
#' This function adds gaps between each GRanges region separately
#' for each GRangesList element, where
#' there is a gap between two GRanges for
#' the same seqnames. When `strandSpecific=TRUE` the gaps are
#' determined per strand.
#'
#' This function is a wrapper around `getGRLgaps()`, which is then
#' concatenated to the input `gr` GRanges object using `S4Vectors::pc()`.
#' When the input `grl` GRanges has column `S4Vectors::values()` then the
#' gaps GRanges object will have `NA` values used by default. To supply
#' values, use the `newValues` argument, which assigns name-value pairs.
#'
#' @family jam GRanges functions
#'
#' @return GRangesList object, sorted per GRangesList element
#'    when `doSort=TRUE`. When `newValues`
#'    is supplied, the values for gaps GRanges elements will be assigned,
#'    otherwise any column values present in `gr` will be `NA` for
#'    gaps elements. The names of gaps elements are assigned using
#'    `gapname` then are made unique using `jamba::makeNames()`,
#'    unless `gapname is NULL`.
#'
#' @param grl GRangesList object
#' @param strandSpecific logical indicating whether the gaps are calculated
#'    per strand, see `getGRLgaps()`.
#' @param gapname,suffix character vector supplying the name to assign to new
#'    gap GRanges elements, using `jamba::makeNames()` with `suffix` as
#'    described to define non-duplicated names. If `gapname is NULL` then
#'    no names are assigned to new gap GRanges entries, however when the
#'    input `gr` GRanges object has names, the concatenation of gaps
#'    causes names `""` to be assigned to all gap GRanges elements, which
#'    are duplicated for multiple gaps.
#' @param newValues list of values to add to the resulting gap GRanges,
#'    whose names become `colnames(grl@unlistData)`, and whose values are used
#'    to populate each column. By default a colname `"feature_type"` is
#'    added, with value `"gap"` added to each row. When `newValues` is `NULL`
#'    then no values are added to the gaps GRanges. For every entry in
#'    `names(newValues)`, if the colname does not exist, it is created
#'    and populated with `default_feature_type` as a default value,
#'    prior to adding gaps.
#' @param default_feature_type,feature_type_colname character values,
#'    indicating the type and colname to populate in `values(grl@unlistData)`
#'    prior to adding gaps. When `feature_type_colname` is NULL,
#'    no action is taken. When `feature_type_colname` does not
#'    exist in `colnames(values(grl@unlistData))`, it is created
#'    with values in `default_feature_type`. Also any colname
#'    defined by `newValues` that does not already exist is also
#'    created and populated with `default_feature_type`.
#' @param doSort logical indicating whether to sort the resulting
#'    GRanges objects. When `doSort=FALSE` the gaps are added to the end
#'    of each input `grl` GRanges object. Note that the GrangesList object
#'    is not sorted, only the GRanges objects within the GRangesList
#'    are sorted.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to `getGRLgaps()`.
#'
#' @examples
#' gr <- GenomicRanges::GRanges(seqnames=rep(c("chr1","chr2"), c(3,2)),
#'    ranges=IRanges::IRanges(start=c(100, 300, 400, 300, 700),
#'       end=c(199, 450, 500, 600, 800)),
#'    strand=rep(c("+","-"), c(3,2)),
#'    feature_type=rep("exon", 5));
#' names(gr) <- jamba::makeNames(rep("exon", length(gr)));
#' gr;
#' addGRgaps(gr);
#'
#' grl <- GenomicRanges::split(gr, GenomicRanges::seqnames(gr));
#' grl;
#' addGRLgaps(grl);
#' addGRLgaps(grl, strandSpecific=FALSE);
#'
#' @export
addGRLgaps <- function
(grl,
 strandSpecific=TRUE,
 gapname="gap",
 suffix="_v",
 newValues=list(feature_type="gap"),
 default_feature_type="exon",
 feature_type_colname="feature_type",
 doSort=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to add gaps as GRanges elements between the
   ## GRanges elements given, applied to each element in the
   ## GRangesList.

   ## Populate feature_type_colname with default_feature_type
   ## if the colname does not already exist
   if (length(feature_type_colname) > 0 &&
         !feature_type_colname %in% colnames(GenomicRanges::values(grl@unlistData)) &&
         length(default_feature_type) > 0) {
      GenomicRanges::values(grl@unlistData)[[feature_type_colname]] <- rep(default_feature_type,
         length.out=length(grl@unlistData));
   }
   ## Populate colnames(values(grl@unlistData))
   ## using names(newValues) when they are not already present
   if (length(newValues) > 0 &&
         length(names(newValues)) > 0 &&
         length(default_feature_type)) {
      for (i in names(newValues)) {
         if (!i %in% colnames(GenomicRanges::values(grl@unlistData))) {
            GenomicRanges::values(grl@unlistData)[[i]] <- rep(default_feature_type,
               length.out=length(grl@unlistData));
         }
      }
   }
   if (verbose) {
      printDebug("addGRLgaps(): ",
         "calling getGRLgaps().");
   }
   grlGaps <- getGRLgaps(grl,
      strandSpecific=strandSpecific,
      ...);
   if (length(grlGaps) == 0) {
      return(grl);
   }
   if (length(gapname) > 0) {
      gapnames <- jamba::makeNames(
         rep(gapname,
            length.out=length(grlGaps@unlistData)),
         suffix=suffix);
      names(grlGaps@unlistData) <- gapnames;
   }
   if (length(newValues) > 0) {
      for (newValueName in names(newValues)) {
         GenomicRanges::values(grlGaps@unlistData)[[newValueName]] <- rep(newValues[[newValueName]],
            length(grlGaps@unlistData));
      }
   }
   if (doSort) {
      grlNew <- sort(S4Vectors::pc(grl, grlGaps));
   } else {
      grlNew <- S4Vectors::pc(grl, grlGaps);
   }
   return(grlNew);
}


#' Find closest exon to splice junction ends
#'
#' Find closest exon to splice junction ends
#'
#' This function is used to annotate splice junction GRanges entries
#' based upon the closest compatible stranded exon boundary.
#'
#' This function should usually be called by `spliceGR2junctionDF()`
#' and not called directly.
#'
#' @family jam RNA-seq functions
#' @family jam GRanges functions
#'
#' @param spliceGRgene GRanges representing splice junctions
#' @param exonsGR GRanges representing flattened exons per gene, as
#'    is produced by `flattenExonsBy()`.
#' @param flipNegativeStrand logical indicating whether to flip
#'    the orientation of features on the negative strand.
#' @param sampleColname character value matching the colname that
#'    defines distinct sample identifier, for which junctions will
#'    be kept separate.
#' @param reportActualCoords logical indicating whether to report
#'    genomic coordinates or transcriptome coordinates. (Work in
#'    progress.)
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
closestExonToJunctions <- function
(spliceGRgene,
 exonsGR,
 flipNegativeStrand=TRUE,
 sampleColname="sample_id",
 reportActualCoords=FALSE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to take:
   ## - spliceGRgene: GRanges representing splice junctions where
   ##    start-end represents the first and last base of the junction span;
   ## - exonsGR: GRanges representing exons per gene, flattened so that no
   ##    two exons overlap on the same strand.
   ##    Expected to have names which are useful downstream, typically gene:exonnum
   ##
   ## This function augments distanceToNearest in two ways:
   ## 1. It returns the closest feature along with the distance, and
   ## 2. It returns the stranded distance, so one can tell whether the splice
   ##    start comes before or after the annotated end of an exon, similarly whether
   ##    splice end comes before or after the annotated start of an exon.
   ##
   ## The hope is that the relative position of the splice site can help name
   ## splice junctions when the site lies between two exons, e.g. 5.1, 5.2, 5.3, etc.
   ##
   ## flipNegativeStrand=TRUE will orient the "from" and "to" to follow
   ## the direction of the transcript, instead of being relative to the genome coordinates.
   ##
   ##
   ## First make sure the spliceGRgene supplied already has values in the colnames we
   ## will be propagating, otherwise we may only update the entries per this method which
   ## might be incomplete
   updateColnames <- c("distFrom", "distTo", "nameFrom", "nameTo",
      "genesDiffer", "genesMatch", "tooFarFrom", "tooFarTo", "tooFar");
   if (any(updateColnames %in% colnames(GenomicRanges::values(spliceGRgene)))) {
      GenomicRanges::values(spliceGRgene) <-
         GenomicRanges::values(spliceGRgene)[,setdiff(colnames(spliceGRgene), updateColnames),drop=FALSE];
   }

   ## Distance from splice start to exon end
   ## ignore.strand=TRUE on resize makes the method search the left side of a splice site with the right side of an exon.
   if (verbose) {
      printDebug("closestExonToJunctions(): ",
         "Finding closest exons for splice starts.");
   }
   spliceStartExonEndD1 <- as.data.frame(
      GenomicRanges::distanceToNearest(
         GenomicRanges::resize(spliceGRgene,
            width=1,
            fix="start",
            ignore.strand=TRUE),
         GenomicRanges::resize(exonsGR,
            width=1,
            fix="end",
            ignore.strand=TRUE),
         select="all"));

   ## Calculate stranded distance
   if (verbose) {
      printDebug("closestExonToJunctions(): ",
         "Calculating stranded distance.");
      print(spliceStartExonEndD1);
   }
   spliceStartExonEndDactual1 <- (
      GenomicRanges::start(
         GenomicRanges::resize(spliceGRgene,
            width=1,
            fix="start",
            ignore.strand=TRUE)[spliceStartExonEndD1[,"queryHits"]]) -
         GenomicRanges::start(
            GenomicRanges::resize(exonsGR,
               width=1,
               fix="end",
               ignore.strand=TRUE)[spliceStartExonEndD1[,"subjectHits"]]));

   ## TODO: add the actual end coordinate of the exon matched
   ## end(exonsGR[spliceStartExonEndD1[,"subjectHits"]])
   spliceStartExonEndDactual1end <- GenomicRanges::end(exonsGR[spliceStartExonEndD1[,"subjectHits"]]);
   spliceStartExonEndDactual <- spliceStartExonEndDactual1 - sign(spliceStartExonEndDactual1);

   ## Flip the direction when strand is negative
   spliceGRgeneNeg <- as.vector(strand(spliceGRgene)) %in% "-";
   if (any(spliceGRgeneNeg)) {
      spliceStartExonEndDactual[spliceGRgeneNeg] <- spliceStartExonEndDactual[spliceGRgeneNeg] * -1;
   }

   ## Distance from splice end to exon start
   ## ignore.strand=TRUE on resize makes the method search the right side of a splice site with the left side of an exon.
   if (verbose) {
      printDebug("closestExonToJunctions(): ",
         "Finding closest exons for splice ends.");
   }
   spliceEndExonStartD1 <- as.data.frame(
      GenomicRanges::distanceToNearest(
         GenomicRanges::resize(spliceGRgene,
            width=1,
            fix="end",
            ignore.strand=TRUE),
         GenomicRanges::resize(exonsGR,
            width=1,
            fix="start",
            ignore.strand=TRUE),
         select="all"));

   ## Calculate stranded distance
   if (verbose) {
      printDebug("closestExonToJunctions(): ",
         "Calculating stranded distance.");
   }
   spliceEndExonStartDactual1 <- (
      GenomicRanges::start(
         GenomicRanges::resize(spliceGRgene,
            width=1,
            fix="end",
            ignore.strand=TRUE)[spliceEndExonStartD1[,"queryHits"]]) -
         GenomicRanges::start(
            GenomicRanges::resize(exonsGR,
               width=1,
               fix="start",
               ignore.strand=TRUE)[spliceEndExonStartD1[,"subjectHits"]]));

   ## TODO: add the actual start coordinate of the exon matched
   ## start(exonsGR[spliceEndExonStartD1[,"subjectHits"]])
   spliceEndExonStartDactual1start <- GenomicRanges::start(exonsGR[spliceEndExonStartD1[,"subjectHits"]]);
   spliceEndExonStartDactual <- spliceEndExonStartDactual1 - sign(spliceEndExonStartDactual1);

   ## Flip the direction when strand is negative
   if (any(spliceGRgeneNeg)) {
      spliceEndExonStartDactual[spliceGRgeneNeg] <- spliceEndExonStartDactual[spliceGRgeneNeg] * -1;
   }

   ## Add the stranded distance to the data.frame typically returned by distanceToNearest
   spliceStartExonEndD1[,"strandedDistance"] <- spliceStartExonEndDactual;
   spliceEndExonStartD1[,"strandedDistance"] <- spliceEndExonStartDactual;

   #NAfrom <- (!seq_along(spliceGRgene) %in% spliceStartExonEndD1[,"queryHits"]);
   #NAto <- (!seq_along(spliceGRgene) %in% spliceEndExonStartD1[,"queryHits"]);
   #itherNA <- (NAfrom | NAto);

   ## Update the exon name on the start side (from) and end side (to) of the splice junction
   GenomicRanges::values(spliceGRgene)[spliceStartExonEndD1[,"queryHits"],"nameFrom"] <-
      names(exonsGR[spliceStartExonEndD1[,"subjectHits"]]);
   GenomicRanges::values(spliceGRgene)[spliceEndExonStartD1[,"queryHits"],"nameTo"] <-
      names(exonsGR[spliceEndExonStartD1[,"subjectHits"]]);

   ## Add nameFromTo column
   GenomicRanges::values(spliceGRgene)[,"nameFromTo"] <- pasteByRow(
      GenomicRanges::values(spliceGRgene)[,c("nameFrom", "nameTo")],
      sep=" ",
      na.rm=TRUE);

   ## Update the stranded distance on the start side (from) and end side (to) of the splice junction
   GenomicRanges::values(spliceGRgene)[spliceStartExonEndD1[,"queryHits"],"distFrom"] <-
      spliceStartExonEndD1[,"strandedDistance"];
   GenomicRanges::values(spliceGRgene)[spliceEndExonStartD1[,"queryHits"],"distTo"] <-
      spliceEndExonStartD1[,"strandedDistance"];

   ## TODO: report the actual coordinate boundary being matched
   if (reportActualCoords) {
      if (verbose) {
         printDebug("closestExonToJunctions(): ",
            "Reporting actual coords.");
      }
      GenomicRanges::values(spliceGRgene)[spliceStartExonEndD1[,"queryHits"],"coordFrom"] <-
         spliceStartExonEndDactual1end;
      GenomicRanges::values(spliceGRgene)[spliceEndExonStartD1[,"queryHits"],"coordTo"] <-
         spliceEndExonStartDactual1start;
      ## TODO: flip the from/to for negative strand entries
      if (flipNegativeStrand) {
      }
   }

   ## Optionally flip negative strand "from" and "to" entries
   if (flipNegativeStrand) {
      if (verbose) {
         printDebug("closestExonToJunctions(): ",
            "Flipping negative strand.");
      }
      fromToCols <- paste0(rep(c("dist", "name", "coord", "tooFar"), each=2), c("From", "To"));
      switchCols1 <- intersect(fromToCols, colnames(GenomicRanges::values(spliceGRgene)));
      switchCols2 <- as.vector(matrix(nrow=2, switchCols1)[2:1,])
      if (verbose) {
         printDebug("closestExonToJunctions(): ",
            "switchCols1:",
            switchCols1);
         printDebug("closestExonToJunctions(): ",
            "switchCols2:",
            switchCols2);
      }
      negStrand <- (as.vector(strand(spliceGRgene)) %in% "-");
      if (any(negStrand)) {
         GenomicRanges::values(spliceGRgene)[negStrand,switchCols1] <-
            GenomicRanges::values(spliceGRgene)[negStrand,switchCols2];
      }
   }

   retVal <- list(
      spliceStartExonEndD=spliceStartExonEndD1,
      spliceEndExonStartD=spliceEndExonStartD1,
      spliceGRgene=spliceGRgene);
   return(retVal);
}

#' Splice junction data.frame summary
#'
#' Splice junction data.frame summary
#'
#' This function takes a GRanges object representing multiple splice
#' junction ranges, with associated scores, and returns a data.frame
#' summary of junctions with annotated boundaries using a set
#' of gene exon models. Junctions whose ends are within `spliceBuffer`
#' distance are combined, and the scores are summed.
#'
#' By default, junctions not within `spliceBuffer` of a compatible exon
#' boundary are named by the nearest exon boundary, and the distance
#' upstream or downstream from the boundary.
#'
#' Multiple samples can be processed together, and the results will
#' be aggregated within each sample, using `sampleColname`. The results
#' in that case may be cast to wide format using `nameFromTo` as the
#' row identifier, `score` as the value column, and `sampleColname` as
#' the new column headers.
#'
#' @family jam RNA-seq functions
#' @family jam GRanges functions
#'
#' @param spliceGRgene GRanges object containing splice junctions, where
#'    the `scoreColname` contains numeric scores.
#' @param exonsGR GRanges object containing flattened exons by gene,
#'    as is provided by `flattenExonsBy()`.
#' @param spliceBuffer integer distance allowed from a compatible exon
#'    boundary, for a junction read to be snapped to that boundary.
#' @param useOnlyValidEntries logical indicating whether to remove
#'    junctions that do not align with a compatible exon boundary.
#' @param renameTooFar logical indicating whether junctions are
#'    named by the nearest exon boundary and the distance to that
#'    boundary.
#' @param scoreColname,sampleColname colnames in `values(spliceGRgene)`
#'    to define the score, and `sample_id`.
#' @param flipNegativeStrand logical indicating whether to flip the
#'    orientation of negative strand features when matching exon
#'    boundaries. This argument is passed to `closestExonToJunctions()`.
#' @param returnGRanges logical indicating whether to return GRanges,
#'    or by default, `data.frame`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
spliceGR2junctionDF <- function
(spliceGRgene,
 exonsGR,
 spliceBuffer=3,
 geneExonSep="(:|_exon)",
 useOnlyValidEntries=FALSE,
 renameTooFar=TRUE,
 scoreColname="score",
 sampleColname="sample_id",
 flipNegativeStrand=TRUE,
 returnGRanges=FALSE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to take junctions as GRanges, and exons as GRanges, and
   ## name the junctions by gene, then exonFrom, exonTo, for use in
   ## downstream differential analyses.
   ##
   ## geneExonSep indicates the separater used to delimit the geneSymbol from the exon
   ##    in exonsGR, so the geneSymbol can be compared for the junction start and end.
   ##
   ## useOnlyValidEntries=TRUE means that only entries whose junction start and end
   ## are within spliceBuffer distance of an annotated exon, and where the two exons
   ## are from the same gene.
   ##
   ## flipNegativeStrand=TRUE will flip the exonFrom, exonTo to follow the direction
   ## of the transcript
   ##
   retVals <- list();
   ##
   ## TODO: confirm strandedness is used in edge cases where exons from two genes overlap on opposite
   ## strands; confirm that the junction sites are not ambiguously assigned to these genes without regard
   ## to the strandedness of the splice junction and the exons.
   ##
   ## TODO: snap coordinates to the exon in cases where it is within the spliceBuffer distance
   ##
   ## Quick method to review strandedness -- note that all entries had perfectly matched strandFrom and strandTo
   #   values(spliceGRgene)[!eitherNA,"strandFrom"] <- as.vector(strand(exonsGR[spliceStartExonEnd]));
   #   values(spliceGRgene)[!eitherNA,"strandTo"] <- as.vector(strand(exonsGR[(spliceEndExonStart)]));
   #   table(strandFrom=values(spliceGRgene)[!eitherNA,"strandFrom"], strandTo=values(spliceGRgene)[!eitherNA,"strandTo"]);
   #   table(strandSplice=as.vector(strand(spliceGRgene[!eitherNA])), strandTo=values(spliceGRgene)[!eitherNA,"strandTo"]);
   sampleColname <- intersect(sampleColname, colnames(GenomicRanges::values(spliceGRgene)));

   ## Vectorized logic:
   ## - find closest start/end for each splice start/end
   ## - subset for splice start/end having same gene
   ##
   ## Call a method which encapsulates distanceToNearest(), calculates stranded distance,
   ## and knows to match splice start with exon end, etc.
   #spliceGRgeneVals <- closestExonToJunctions(spliceGRgene=spliceGRgene[,scoreColname], exonsGR=exonsGR,
   #   flipNegativeStrand=flipNegativeStrand, ...);
   spliceGRgene <- closestExonToJunctions(spliceGRgene=spliceGRgene[,c(scoreColname,sampleColname)],
      exonsGR=exonsGR,
      flipNegativeStrand=flipNegativeStrand,
      sampleColname=sampleColname,
      verbose=verbose)$spliceGRgene;
   #...)$spliceGRgene;
   #spliceGRgene <- spliceGRgeneVals$spliceGRgene;

   ## Check when genes match for the two junction sites
   ## Note: we changed from genesDiffer, because in cases where no gene is returned, the two
   ## sides of the junction would be equal, but still not represent the desired outcome.
   genesMatch <- (!is.na(GenomicRanges::values(spliceGRgene)[,"nameFrom"]) &
         !is.na(GenomicRanges::values(spliceGRgene)[,"nameTo"]) &
         (gsub(paste0(geneExonSep, ".*$"), "", GenomicRanges::values(spliceGRgene)[,"nameFrom"]) ==
            gsub(paste0(geneExonSep, ".*$"), "", GenomicRanges::values(spliceGRgene)[,"nameTo"]) ) );
   GenomicRanges::values(spliceGRgene)[,"genesMatch"] <- genesMatch;
   numGenesMatch <- sum(GenomicRanges::values(spliceGRgene)[,"genesMatch"]);
   numGenesDiffer <- (length(spliceGRgene) - numGenesMatch);
   if (verbose && numGenesMatch > 0) {
      printDebug("spliceGR2junctionDF(): ",
         formatInt(numGenesMatch),
         " entries out of ",
         formatInt(length(spliceGRgene)),
         " had the same gene for splice start and end, excepting ",
         formatInt(numGenesDiffer),
         " entries.");
   }

   ## Check when the junction is too far from the nearest exon
   tooFarFrom <- (
      abs(GenomicRanges::values(spliceGRgene)[,"distFrom"]) > spliceBuffer |
      is.na(GenomicRanges::values(spliceGRgene)[,"distFrom"]));
   tooFarTo <- (
      abs(GenomicRanges::values(spliceGRgene)[,"distTo"]) > spliceBuffer |
      is.na(GenomicRanges::values(spliceGRgene)[,"distTo"]));
   tooFar <- (tooFarFrom | tooFarTo);
   GenomicRanges::values(spliceGRgene)[,"tooFarFrom"] <- tooFarFrom;
   GenomicRanges::values(spliceGRgene)[,"tooFarTo"] <- tooFarTo;
   GenomicRanges::values(spliceGRgene)[,"tooFar"] <- tooFar;
   numTooFar <- sum(GenomicRanges::values(spliceGRgene)[,"tooFar"]);
   if (verbose && numTooFar > 0) {
      printDebug("spliceGR2junctionDF(): ",
         formatInt(numTooFar),
         " entries were farther from the nearest exon than spliceBuffer:",
         spliceBuffer);
   }

   ## Rename entries by the distance from the nearest exon.
   ## Note that junctions within spliceBuffer get "snapped" to the exon and combined
   ## with other junctions also within this splice buffer distance.
   ## But entries outside the spliceBuffer of an exon are not "snapped" to other junctions
   ## which might be within spliceBuffer of each other.
   ## The potential negative is that novel splice sites would not have the benefit
   ## of combined counts if junction reads are within spliceBuffer distance of each
   ## other.
   if (renameTooFar && numTooFar > 0) {
      if (verbose) {
         printDebug("spliceGR2junctionDF(): ",
            "Renaming exon sites by distance for entries outside spliceBuffer.");
      }
      if (any(tooFarFrom)) {
         GenomicRanges::values(spliceGRgene)[tooFarFrom,"nameFrom"] <- paste(
            GenomicRanges::values(spliceGRgene)[tooFarFrom,"nameFrom"],
            GenomicRanges::values(spliceGRgene)[tooFarFrom,"distFrom"],
            sep=".");
      }
      if (any(tooFarTo)) {
         GenomicRanges::values(spliceGRgene)[tooFarTo,"nameTo"] <- paste(
            GenomicRanges::values(spliceGRgene)[tooFarTo,"nameTo"],
            GenomicRanges::values(spliceGRgene)[tooFarTo,"distTo"],
            sep=".");
      }
   }

   ## TODO: rename entries which are "tooFar" so the exon number is distinctly different
   ## from entries with are not "tooFar".
   ##
   ## I.e. entries within spliceBuffer for Apoe:003 would be called "Apoe:003" but entries
   ##    farther away than spliceBuffer would be called something like "Apoe:003.1".
   ## This would affect the downstream collapse of splice read counts by exon site, which itself
   ## has the effect of combining read counts which are within spliceBuffer of an annotated
   ## exon site.
   ##
   ## One idea is to run all samples to find the unique global set of start and end sites, then
   ## use this set to define new "exons" which become a new exonsGRextended...
   ## This exonsGRextended would include extra exons, named like this:
   ##    "Apoe:003", "Apoe:003.1", "Apoe:003.2", "Apoe:004", etc.
   ## Then when running this function, use exonsGR=exonsGRextended, and all entries should be
   ## within spliceBuffer distance.
   ##

   ## Optionally return only entries where both ends of the junction snap to an annotated exon,
   ## and where the exons are both from the same gene.
   if (useOnlyValidEntries) {
      if (verbose) {
         printDebug("spliceGR2junctionDF(): ",
            "Only entries whose genes match, and which are within spliceBuffer, will be used for the junction count matrix.");
      }
      toUse <- which(GenomicRanges::values(spliceGRgene)[,"genesMatch"] &
            !GenomicRanges::values(spliceGRgene)[,"tooFar"]);
   } else {
      ## Note that we still only return entries for which there is some nearby exon
      ## TODO: allow for un-named entries to have a temporary name for the purpose of
      ## returning the values.
      if (verbose) {
         printDebug("spliceGR2junctionDF(): ",
            "All entries having any nearest gene will be used, without regard to matching genes or distance from exon.");
      }
      toUse <- which(!is.na(GenomicRanges::values(spliceGRgene)[,"nameFrom"]) &
            !is.na(GenomicRanges::values(spliceGRgene)[,"nameTo"]));
   }
   if (verbose) {
      printDebug("spliceGR2junctionDF(): ",
         "Using ", ifelse(length(toUse) == length(spliceGRgene), "all ", ""),
         formatInt(length(toUse)),
         " entries out of ",
         formatInt(length(spliceGRgene)),
         " to create a data.frame of junction counts by exon.");
   }
   if (verbose && length(toUse) != length(spliceGRgene)) {
      printDebug("spliceGR2junctionDF(): ",
         "Note: ",
         formatInt(length(spliceGRgene) - length(toUse)),
         " entries did not meet the criteria.");
   }

   spliceColnames <- c("seqnames", "start", "end", "nameFrom", "nameTo",
      scoreColname, "strand", sampleColname);
   spliceCountsDF <- as.data.frame(spliceGRgene[toUse])[,spliceColnames,drop=FALSE];
   spliceCountsDF[,"nameFromTo"] <- pasteByRow(spliceCountsDF[,c("nameFrom","nameTo"),drop=FALSE],
      sep=" ",
      na.rm=TRUE);
   spliceCountsDF[,"nameFromToSample"] <- pasteByRow(
      spliceCountsDF[,c("nameFromTo", sampleColname),drop=FALSE],
      sep=":!:");

   if (verbose) {
      printDebug("spliceGR2junctionDF(): ",
         "Shrinking matrix to combine counts spanning the same junctions.");
   }
   spliceCountsDFshrunk <- renameColumn(
      shrinkMatrix(spliceCountsDF[,scoreColname,drop=FALSE],
         groupBy=spliceCountsDF[,"nameFromToSample"],
         shrinkFunc=sum),
      from="groupBy",
      to="nameFromToSample");
   if (length(sampleColname) > 0) {
      spliceCountsDFshrunk[,c("nameFromTo",sampleColname)] <- rbindList(
         strsplit(spliceCountsDFshrunk[,"nameFromToSample"], ":!:"));
   } else {
      spliceCountsDFshrunk[,"nameFromTo"] <- spliceCountsDFshrunk[,"nameFromToSample"];
   }
   spliceCountsDFshrunk[,c("nameFrom", "nameTo")] <- rbindList(
      strsplit(spliceCountsDFshrunk[,"nameFromTo"], " "));
   #spliceCountsDFshrunk <- cbind(spliceCountsDFshrunk,
   #   rbindList(strsplit(spliceCountsDFshrunk[,"groupBy"], " ")));
   #colnames(spliceCountsDFshrunk) <- c("nameFromTo", spliceGRname, "nameFrom", "nameTo");

   ## Now add the start and end coordinates, so the results can be plotted
   if (verbose) {
      printDebug("spliceGR2junctionDF(): ",
         "Adding ref, strand, start, end.");
   }
   matchFromTo <- match(spliceCountsDFshrunk[,"nameFromToSample"],
      spliceCountsDF[,"nameFromToSample"]);
   spliceCountsDFshrunk[,"ref"] <- as.vector(spliceCountsDF[matchFromTo,"seqnames"]);
   spliceCountsDFshrunk[,"start"] <- as.vector(spliceCountsDF[matchFromTo,"start"]);
   spliceCountsDFshrunk[,"end"] <- as.vector(spliceCountsDF[matchFromTo,"end"]);
   spliceCountsDFshrunk[,"strand"] <- as.vector(spliceCountsDF[matchFromTo,"strand"]);
   spliceCountsDFshrunk[,"width"] <- spliceCountsDFshrunk[,"end"] - spliceCountsDFshrunk[,"start"];
   spliceCountsDFshrunk <- mixedSortDF(spliceCountsDFshrunk,
      byCols=c("ref", "start", "end", "strand"));

   ## geneGsub removes the geneExon separator, and everything after it
   geneGsub <- paste0(geneExonSep, ".*$");
   spliceCountsDFshrunk[,"geneSymbolFrom"] <- gsub(geneGsub, "",
      spliceCountsDFshrunk[,"nameFrom"]);
   spliceCountsDFshrunk[,"geneSymbolTo"] <- gsub(geneGsub, "",
      spliceCountsDFshrunk[,"nameTo"]);
   ## For visual ease, add exon numbers to their own columns
   exonGsub <- paste0("^.+", geneExonSep);
   spliceCountsDFshrunk[,"exonFrom"] <- gsub("^[0]+", "",
      gsub(exonGsub, "",
         spliceCountsDFshrunk[,"nameFrom"]));
   ## Remove the padded zeros from the beginning of each exon number
   spliceCountsDFshrunk[,"exonTo"] <- gsub("^[0]+", "",
      gsub(exonGsub, "",
         spliceCountsDFshrunk[,"nameTo"]));

   ## Create junctionID only after we decided how to orient exonFrom and exonTo
   ## for negative strand entries, avoiding re-creating the name for those rows
   spliceCountsDFshrunk[,"junctionID"] <- pasteByRow(
      spliceCountsDFshrunk[,c("exonFrom","exonTo",sampleColname),drop=FALSE],
      sep="_",
      na.rm=TRUE);

   ## Prepare data to return
   if (returnGRanges) {
      retVals$spliceGRgene <- spliceGRgene;
      retVals$spliceCountsDFshrunk <- spliceCountsDFshrunk;
   } else {
      retVals <- spliceCountsDFshrunk;
   }
   return(retVals);
}
