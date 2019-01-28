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
#'    file ends with ".gz", gzcat is used to gunzip the file during the
#'    import step. (TODO: verify that handling of gz files is portable across
#'    architectures, or gracefully exit as needed.)
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
 verbose=TRUE,
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
   if (igrepHas("[.]gz$", GTF)) {
      gtfDF <- fread(paste0("gzcat ", GTF),
         sep="\t",
         autostart=20,
         #col.names,
         nrows=nrows);
   } else {
      gtfDF <- fread(GTF,
         sep="\t",
         autostart=20,
         #col.names,
         nrows=nrows);
   }
   ## Subset to clear some memory
   gtfDF <- subset(gtfDF, gtfDF[[3]] %in% c(geneFeatureType,
      txFeatureType));

   ## Determine which rows are gene and transcript
   geneRows <- (gtfDF[[3]] %in% geneFeatureType);
   txRows <- (gtfDF[[3]] %in% txFeatureType);


   ## gene attributes
   geneM <- do.call(cbind, lapply(nameVector(geneAttrNames),
      function(attrName){
         if (verbose) {
            jamba::printDebug("makeTx2geneFromGtf() :",
               "gene attributes:", attrName);
         }
         attrGrep <- paste0('^.*', attrName, ' ["]([^"]+)["].*$');
         if (jamba::igrepHas(attrGrep, gtfDF[geneRows,][[9]])) {
            attrValues <- gsub(attrGrep,
               "\\1",
               gtfDF[geneRows,,drop=FALSE][[9]]);
         } else {
            jamba::printDebug("Note: No gene attributes found for:", attrName);
            attrValues <- NULL;
         }
      }));

   ## transcript attributes
   txM <- do.call(cbind, lapply(nameVector(c(txAttrNames,geneAttrNames)),
      function(attrName){
         if (verbose) {
            jamba::printDebug("makeTx2geneFromGtf(): ",
               "tx attributes:", attrName);
         }
         attrGrep <- paste0('^.*', attrName, ' ["]([^"]+)["].*$');
         if (jamba::igrepHas(attrGrep, gtfDF[txRows,][[9]])) {
            attrValues <- gsub(attrGrep,
               "\\1",
               gtfDF[txRows,,drop=FALSE][[9]]);
         } else {
            jamba::printDebug("Note: No tx attributes found for:", attrName);
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
#' \code{annotateGRfromGRnew()}, and \code{assignGRLexonNames()}.

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
#' @family jam RNA-seq functions
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
 verbose=TRUE,
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
   ## - annotateGRLfromGRL(), annotateGRfromGRnew()
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
         method="flank",
         verbose=verbose);
   }
   ## add transcript_id annotation
   values(threeUtrGRLdetRange@unlistData)[,"transcript_id"] <- rep(
      names(threeUtrGRLdetRange),
      S4Vectors::lengths(threeUtrGRLdetRange));
   ## add gene_name annotation
   values(threeUtrGRLdetRange@unlistData)[,"gene_name"] <- tx2geneDF[
      match(values(threeUtrGRLdetRange@unlistData)[,"transcript_id"],
         tx2geneDF[,"transcript_id"]),"gene_name"];

   ####################################################
   ## Re-aggregate by gene
   if (verbose) {
      printDebug("gencode2ale(): ",
         "Splitting ranges by gene.");
   }
   threeUtrGRLdetGeneGRL <- GenomicRanges::split(
      threeUtrGRLdetRange@unlistData,
      f=values(threeUtrGRLdetRange@unlistData)[,"gene_name"]);
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
      values(threeUtrGRLdetGeneGRLred2@unlistData)[,"ALE_name"];
   values(threeUtrGRLdetGeneGRLred2@unlistData)[,"score"] <- 0.5;
   retVals$aleGRL <- threeUtrGRLdetGeneGRLred2;

   ####################################################
   ## Subset ALE containing 2 or more ALEs per gene
   GencodeALEmin2 <- threeUtrGRLdetGeneGRLred2[
      S4Vectors::lengths(threeUtrGRLdetGeneGRLred2) > 1];
   if (verbose) {
      printDebug("gencode2ale(): ",
         "Filtered ",
         formatInt(length(threeUtrGRLdetGeneGRLred2)),
         " genes to ",
         formatInt(length(GencodeALEmin2)),
         " genes having multiple ranges.");
   }
   ## vector of genes containing multiple ALEs
   GencodeALEmin2genes <- jamba::mixedSort(unique(values(GencodeALEmin2@unlistData)[,"gene_name"]));
   retVals$multiALEgenes <- GencodeALEmin2genes;

   ####################################################
   ## Create a tx-to-ALE xref using multi-ALE genes
   if (verbose) {
      printDebug("gencode2ale(): ",
         "Creating transcript-to-ALE xref for multi-range genes.");
   }
   ale2txL <- strsplit(nameVector(
      values(subset(threeUtrGRLdetGeneGRLred2@unlistData,
         gene_name %in% GencodeALEmin2genes))[,c("transcript_id","ALE_name")]),
      ",");
   tx2ale <- nameVector(list2df(ale2txL)[,c("item","value")]);
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
         iMatrixAle <- log2(1+shrinkMatrix(2^(iMatrix[names(tx2ale),,drop=FALSE])-1,
            groupBy=tx2ale,
            shrinkFunc=sum,
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
#'    transcript rows and sample columns.
#' @param iMatrixTxGrp numeric matrix of read counts averaged by sample
#'    group. If this matrix is not provided, it will be calculated
#'    from `iMatrixTx`
#'    using `jamba::rowGroupMeans()` and the `groups` parameter.
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
      if (verbose) {
         printDebug("defineDetectedTx(): ",
            "Calculating iMatrixTxGrp.");
      }
      iMatrixTxGrp <- rowGroupMeans(iMatrixTx,
         useMedian=useMedian,
         groups=groups);
   }

   ## group mean values for TPM
   if (length(iMatrixTxTPMGrp) == 0) {
      if (length(iMatrixTxTPM) > 0) {
         if (length(groups) == 0) {
            stop("defineDetectedTx() requires groups, named by colnames(iMatrixTxTPM).");
         }
         if (verbose) {
            printDebug("defineDetectedTx(): ",
               "Calculating iMatrixTxTPMGrp.");
         }
         iMatrixTxTPMGrp <- rowGroupMeans(iMatrixTxTPM,
            useMedian=useMedian,
            groups=groups);
      }
   }

   ######################################################################
   ## matrix associating gene to transcript_id, mainly useful since it
   ## returns transcripts in the same order as each matrix below, which
   ## otherwise has rownames based upon gene and not transcript.
   iRows <- match(rownames(iMatrixTxGrp), tx2geneDF$transcript_id);
   if (verbose) {
      printDebug("defineDetectedTx(): ",
         "shrinkMatrix Tx names.");
   }
   txExprGrpTx <- shrinkMatrix(rownames(iMatrixTxGrp),
      groupBy=tx2geneDF[iRows,"gene_name"],
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
         groupBy=tx2geneDF[iRows,"gene_name"],
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
         groupBy=tx2geneDF[iRows,"gene_name"],
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
      groupBy=tx2geneDF[iRows,"gene_name"],
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
         groupBy=tx2geneDF[iRows,"gene_name"],
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
#' @param returnClass character string indicating the return data type,
#'    `"data.frame"` returns a `data.frame` whose first column contains
#'    entries from `groups`; `"matrix"` returns a numeric matrix whose
#'    rownames are entries from `groups`.
#' @param verbose logical indicating whether to print verbose output.
#'
#' @import data.table
#'
#' @family jam table functions
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
         x=x,
         groupBy=groupBy),
      key="groupBy");
   if (verbose) {
      t2 <- Sys.time();
   }

   ## Operate on the DT object
   byDT <- DT[,lapply(.SD, shrinkFunc), by="groupBy"];
   if (verbose) {
      t3 <- Sys.time();
   }

   if (verbose) {
      printDebug("shrinkMatrix(): ",
         "Duration for data.table DT creation: ",
         format(t2-t1));
      printDebug("shrinkMatrix(): ",
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
#' codonFile <- system.file("extdata", "Mouse_codon_usage.txt", package="farrisdata");
#' codonDF <- codonUsage2df(codonFile);
#'
#' @family jam RNA-seq functions
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
      rbindList(strsplit(codonV, "[() ]+")));
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
#' @family jam RNA-seq functions
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
      pasteByRow(
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
#' @family jam RNA-seq functions
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
   lapply(nameVector(txDatNames), function(i){
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
#' GRangesList element has only one strand, as is true for exons
#' in one transcript.
#'
#' @return GRangesList containing only the first GRanges feature
#'    in stranded order.
#'
#' @family jam GRanges functions
#'
#' @param grl GRangesList
#' @param method character value in `c("endoapply", "flank")`
#'    representing which method to use to define the first feature.
#'    The `"endoapply"` method uses `S4Vectors::endoapply()` to
#'    iterate each GRangesList element. The `"flank"` method is
#'    intended to be equivalent but uses the `GenomicRanges::flank()`
#'    function, which is vectorized.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
getFirstStrandedFromGRL <- function
(grl,
 method=c("endoapply","flank"),
 verbose=TRUE,
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
      printDebug("getFirstStrandedFromGRL(): ",
         "sorting grl");
   }
   grl <- sortGRL(grl,
      verbose=verbose);

   if (method %in% "endoapply") {
      if (verbose) {
         printDebug("getFirstStrandedFromGRL(): ",
            "performing endoapply() logic");
      }
      grl2 <- endoapply(grl, function(iGR){
         if ("-" %in% strand(iGR)[1]) {
            tail(iGR, 1);
         } else {
            head(iGR, 1);
         }
      });
      return(grl2);
   } else {
      if (verbose) {
         printDebug("getFirstStrandedFromGRL(): ",
            "performing flank() logic");
      }
      #grlWidths <- width(grl);
      grlHeads1 <- heads(grl, 1);
      grlWidths1 <- width(grlHeads1);
      grlWidths2 <- width(tails(grl, 1));
      if (verbose) {
         printDebug("getFirstStrandedFromGRL(): ",
            "applying range() function, class(grl):",
            class(grl));
      }
      grlRanges <- range(grl);
      if (verbose) {
         printDebug("getFirstStrandedFromGRL(): ",
            "Validating one strand per range(grl).");
      }
      if (length(grlRanges) != length(grlRanges@unlistData)) {
         stop("getFirstStrandedFromGRL() requires each GRanges element to have only one strand.");
      }
      if (verbose) {
         printDebug("getFirstStrandedFromGRL(): ",
            "applying ifelse() logic");
      }
      grlFlankWidth <- ifelse(as.vector(unlist(strand(grlHeads1))) %in% "+",
         unlist(grlWidths1),
         unlist(grlWidths2));
      if (verbose) {
         printDebug("getFirstStrandedFromGRL(): ",
            "applying flank() function");
      }
      grl3 <- flank(grlRanges,
         start=TRUE,
         width=-grlFlankWidth);
      if (verbose) {
         printDebug("getFirstStrandedFromGRL(): ",
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
#' @param removeSplitColname logical indicating whether to remove
#'    the `splitColname` from the output GRangesList.
#' @param ... additional arguments are ignored.
#'
#' @export
sortGRL <- function
(GRL,
 splitColname="splitColname",
 removeSplitColname=TRUE,
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
      names(GRL) <- makeNames(rep("GRL", length(GRL)));
   }

   values(GRL@unlistData)[,splitColname] <- factor(rep(names(GRL), elementNROWS(GRL)),
      levels=names(GRL));
   GR1 <- sort(GRL@unlistData);
   #values(GR1)[,splitColname] <- factor(values(GR1)[,splitColname], unique(values(GR1)[,splitColname]));
   if (verbose) {
      printDebug("sortGRL(): ",
         "splitting GRanges into list");
   }
   GRL <- GenomicRanges::split(GR1,
      f=values(GR1)[,splitColname]);

   ## Optionally remove splitColname from the result
   if (removeSplitColname) {
      keepNames <- setdiff(names(values(GRL@unlistData)),
         splitColname);
      #GRL@unlistData <- GRL@unlistData[,keepNames];
      values(GRL@unlistData) <- values(GRL@unlistData)[,keepNames,drop=FALSE];
   }
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
#' @export
annotateGRfromGR <- function
(GR1,
 GR2,
 grOL=NULL,
 numAsStrings=FALSE,
 stringShrinkFunc=function(...){cPasteUnique(..., doSort=TRUE)},
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
   grOLtable <- table(from(grOL));
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
   colClasses <- sapply(colnames(values(GR2)), function(iCol){
      class(values(GR2)[,iCol])
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
            values(GR2[grOLs1])[,iCol];
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
         grOLi <- shrinkMatrix(as.data.frame(values(GR2)[grOLs,iCol,drop=FALSE]),
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
            as.data.frame(values(GR2)[grOLs1,iCol,drop=FALSE]);
         });
      }
      if (nrow(grOLmUse) > 0) {
         stringShrunk <- lapply(nameVector(stringCols), function(iCol){
            if (verbose) {
               printDebug("annotateGRfromGR(): ",
                  "   ",
                  iCol);
            }
            if (igrepHas("list", class(values(GR2)[grOLs,iCol]))) {
               ## list column
               iVals <- values(GR2)[grOLs,iCol];
               names(iVals) <- names(GR2)[grOLs];
               iValsX <- unlist(iVals);
               iValsXnames1 <- rep(grOLq,
                  S4Vectors::lengths(iVals));
               if (1 == 2) {
                  iValsDF <- data.frame(stringsAsFactors=FALSE,
                     check.names=FALSE,
                     iValsX=iValsX,
                     iValsXnames1=iValsXnames1);
                  if (useMixedSort) {
                     if (verbose) {
                        printDebug("annotateGRfromGR(): ",
                           "Using mixedSortDF() on iCol:",
                           iCol);
                        printDebug("annotateGRfromGR(): ",
                           "      iValsDF:");
                        print(head(iValsDF));
                     }
                     if (DEBUG) {
                        return(iValsDF);
                     }
                     #iValsDF <- mmixedOrderDF(iValsDF);
                     iValsDF <- mixedSortDF(iValsDF);
                  }
                  iValsDFnonNA <- which(!is.na(iValsDF$iValsX));
                  iX <- nameVector(rep(NA, length(GR1)), names(GR1));
                  iValsSplit <- split(iValsDF[,"iValsX"][iValsDFnonNA],
                     iValsDF[,"iValsXnames1"][iValsDFnonNA]);
               } else {
                  iValsSplit <- split(iValsX, iValsXnames1);
               }
               iXnonNA <- stringShrinkFunc[[iCol]](iValsSplit,
                  sep=sep);
               iX[names(iXnonNA)] <- iXnonNA;
               iX;
            } else {
               ## Non-list column
               iValsX <- values(GR2)[grOLs,iCol];
               iValsXnames1 <- grOLq;

               if (1 == 2) {
                  if (useMixedSort) {
                     if (verbose) {
                        printDebug("      useMixedSort started");
                     }
                     iValsXo <- mixedOrder(iValsX);
                     iValsX <- iValsX[iValsXo];
                     iValsXnames1 <- iValsXnames1[iValsXo];
                     if (verbose) {
                        printDebug("      useMixedSort completed");
                     }
                  }
               }

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
         printDebug("annotateGRfromGR(): ",
            "   nrow(numShrunkDF):",
            formatInt(nrow(numShrunkDF)),
            " (before)");
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
   if (ncol(values(GR1)) > 0 && any(colnames(grOL1) %in% colnames(values(GR1)))) {
      newColnames <- make.unique(c(colnames(values(GR1)), colnames(grOL1)),
         sep="_v");
      newColnames2 <- tail(newColnames, ncol(grOL1));
      colnames(grOL1) <- newColnames2;
   }

   ## Now append each column, meanwhile fix some issues with ,-delimiters
   for (iCol in colnames(grOL1)) {
      if (igrepHas("integer|numeric", class(grOL1[,iCol]))) {
         values(GR1)[,iCol] <- numeric(0);
      } else {
         grOL1[,iCol] <- gsub(", ", ",", grOL1[,iCol]);
         values(GR1)[,iCol] <- "";
      }
      if (verbose) {
         printDebug("annotateGRfromGR(): ",
            "head(grOL1):");
         print(head(grOL1, 5));
      }
      values(GR1)[grOLqAll,iCol] <- grOL1[,iCol];
      blankRows <- seq_along(GR1)[-grOLqAll];
      if (length(blankRows) > 0) {
         values(GR1[blankRows])[,iCol] <- NA;
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
#' @param addGRLnames logical indicating whether to add the names of each
#'    GRangesList object to the output object, useful for tracking the
#'    annotations to the source data.
#' @param returnType character value indicating whether to return
#'    GRangesList "GRL" or GRange "GR" object.
#' @param splitColname character value used internally to indicate how
#'    to split the resulting GRanges annotated data back to GRangesList.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to `annotateGRfromGR()`.
#'    To customize the aggregation functions, supply `numShrinkFunc` or
#'    `stringShrinkFunc` as described in `annotateGRfromGR()`.
#'
#' @export
annotateGRLfromGRL <- function
(GRL1,
 GRL2,
 annoName1="name",
 annoName2="name",
 grlOL=NULL,
 addGRLnames=TRUE,
 returnType=c("GRL", "GR"),
 splitColname=annoName1,
 verbose=verbose,
 ...)
{
   ## Purpose is to run annotateGRfromGR() except allow for GRangesList
   ## objects as input.  It accomplishes the task by running
   ## findOverlapsGRL() which ensures GRangesList overlaps are only
   ## allowed for the same annotated entries, defined in annoName1,
   ## and annoName2.
   returnType <- match.arg(returnType);
   if (is.null(grlOL)) {
      grlOL <- findOverlapsGRL(GRL1,
         GRL2,
         annoName1=annoName1,
         annoName2=annoName2);
   }
   if (addGRLnames) {
      if (annoName1 %in% "name") {
         values(GRL1@unlistData)[,"GRL1name"] <- rep(names(GRL1),
            S4Vectors::lengths(GRL1));
         annoName1 <- "GRL1name";
      }
      if (annoName2 %in% "name") {
         values(GRL2@unlistData)[,"GRL2name"] <- rep(names(GRL2),
            S4Vectors::lengths(GRL2));
         annoName2 <- "GRL2name";
      }
   }
   annoNames2 <- setdiff(colnames(values(GRL2@unlistData)), annoName2);
   GR12 <- annotateGRfromGR(GRL1@unlistData,
      GRL2@unlistData[,annoNames2],
      grOL=grlOL,
      verbose=verbose,
      ...);
   if (returnType %in% "GR") {
      return(GR12);
   } else {
      return(GenomicRanges::split(GR12,
         values(GR12)[,splitColname]));
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
#' @return Hits object, or the natural output from
#'    `GenomicRanges::findOverlaps()` dependent upon the `...` arguments.
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
 ...)
{
   ## Purpose is to run findOverlaps() on GRangesList objects,
   ## where overlaps are required to share the same annotation,
   ## in the annoName column.
   ##
   ## Specifically, it helps run findOverlaps() on GRangesList objects
   ## separated by gene, then only returning entries which match the same
   ## gene.
   grOL <- findOverlaps(GRL1@unlistData,
      GRL2@unlistData,
      ...);
   if (annoName1 %in% "name") {
      GRLnames1 <- rep(names(GRL1), elementNROWS(GRL1));
   } else {
      GRLnames1 <- values(GRL1@unlistData)[,annoName1];
   }
   if (annoName2 %in% "name") {
      GRLnames2 <- rep(names(GRL2), elementNROWS(GRL2));
   } else {
      GRLnames2 <- values(GRL2@unlistData)[,annoName2];
   }
   if (!any(GRLnames1 %in% GRLnames2)) {
      stop("No names are shared between GRL1 and GRL2. Please try again.");
   }
   grOLdf <- data.frame(as.data.frame(grOL));
   grOLdf[,"queryName"] <- GRLnames1[grOLdf[,"queryHits"]];
   grOLdf[,"subjectName"] <- GRLnames2[grOLdf[,"subjectHits"]];
   grOLdfUse <- which(grOLdf[,"queryName"] == grOLdf[,"subjectName"]);
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
#' |======|......|======|=======|======|......|======|=======|
#' .exon1.........exon2a.exon2b..exon2c........exon3a..exon3b.
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
   if (filterTwoStrand && any(S4Vectors::lengths(GRLstrandL) > 1)) {
      if (verbose) {
         printDebug("assignGRLexonNames(): ",
            "removing some multi-stranded exon entries.");
      }
      iRemove <- which(S4Vectors::lengths(GRLstrandL) > 1);
      GRL <- GRL[-iRemove];
   }

   ## check disjoint GRanges
   if (checkDisjoin %in% c("warn","stop")) {
      GRLdis <- disjoin(GRL);
      if (!all(S4Vectors::lengths(GRLdis) == S4Vectors::lengths(GRL))) {
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
   GRLred <- reduce(GRL);

   ## Add geneSymbolColname if it does not already exist
   if (!geneSymbolColname %in% colnames(values(GRLred@unlistData))) {
      values(GRLred@unlistData)[,geneSymbolColname] <- rep(names(GRLred),
         lengths(GRLred));
   }
   if (verbose) {
      printDebug("assignGRLexonNames(): ",
         "head(GRLred):");
      print(head(GRLred));
   }
   if (verbose) {
      printDebug("assignGRLexonNames(): ",
         "geneSymbolColname:",
         geneSymbolColname);
   }
   if (verbose) {
      printDebug("assignGRLexonNames(): ",
         "geneSymbolColname values:",
         head(values(GRLred@unlistData)[,geneSymbolColname], 10));
   }
   GRLredStrand <- unlist(unique(strand(GRLred)));
   GRLredStrandP <- which(GRLredStrand %in% "+");
   GRLredStrandN <- which(GRLredStrand %in% "-");

   ## Stranded exon numbering
   values(GRLred@unlistData)[,exonNameColname] <- "";
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
      values(GRLred[GRLredStrandP]@unlistData)[,exonNameColname] <- makeNames(
         values(GRLred[GRLredStrandP]@unlistData)[,geneSymbolColname],
         suffix=suffix,
         renameOnes=TRUE);
   }
   if (length(GRLredStrandN) > 0) {
      values(GRLred[GRLredStrandN]@unlistData)[,exonNameColname] <- rev(makeNames(
         values(GRLred[rev(GRLredStrandN)]@unlistData)[,geneSymbolColname],
         suffix=suffix,
         renameOnes=TRUE));
   }
   if (verbose) {
      printDebug("assignGRLexonNames(): ",
         "exonNameColname values:",
         head(values(GRLred@unlistData)[,exonNameColname], 10));
   }

   ## Add lowercase letter suffix
   GRLcolnames <- unvigrep(paste0(exonNameColname, "(_v[0-9]|)$"), colnames(values(GRL@unlistData)));
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
   if (1 == 2) {
      values(GRLnew@unlistData)[,exonNameColname1] <- values(GRLnew@unlistData)[,exonNameColname];
      values(GRLnew[GRLnewStrandP]@unlistData)[,exonNameColname1] <- (
         makeNames(values(GRLnew[GRLnewStrandP]@unlistData)[,exonNameColname],
            numberStyle=subFeatureNumberStyle,
            suffix=subFeatureSuffix,
            renameOnes=renameOnes));
      values(GRLnew@unlistData[rev(GRLnewStrandNn),])[,exonNameColname1] <- (
         makeNames(values(GRLnew@unlistData[rev(GRLnewStrandNn),])[,exonNameColname],
            numberStyle=subFeatureNumberStyle,
            suffix=subFeatureSuffix,
            renameOnes=renameOnes));
   } else {
      values(GRLnew@unlistData)[,exonNameColname] <- values(GRLnew@unlistData)[,exonNameColname];
      values(GRLnew[GRLnewStrandP]@unlistData)[,exonNameColname] <- (
         makeNames(values(GRLnew[GRLnewStrandP]@unlistData)[,exonNameColname],
            numberStyle=subFeatureNumberStyle,
            suffix=subFeatureSuffix,
            renameOnes=renameOnes));
      values(GRLnew@unlistData[rev(GRLnewStrandNn),])[,exonNameColname] <- (
         makeNames(values(GRLnew@unlistData[rev(GRLnewStrandNn),])[,exonNameColname],
            numberStyle=subFeatureNumberStyle,
            suffix=subFeatureSuffix,
            renameOnes=renameOnes));
   }

   return(GRLnew);

   ## subsection exon numbering using lowercase letters
   GRnew <- renumberGRanges(GR1=GRL@unlistData,
      groupColname=geneSymbolColname);
   GRLnew <- split(GRnew[,c(colnames(values(GRL@unlistData)),"exon_id")],
      values(GRnew)[,geneSymbolColname]);
   values(GRLnew@unlistData) <- renameColumn(values(GRLnew@unlistData),
      from="exon_id", to=exonNameColname)

   return(GRLnew);
}
