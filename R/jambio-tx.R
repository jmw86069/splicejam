
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
         printDebug("makeTx2geneFromGtf() gene attributes:", attrName);
         attrGrep <- paste0('^.*', attrName, ' ["]([^"]+)["].*$');
         if (igrepHas(attrGrep, gtfDF[geneRows,][[9]])) {
            attrValues <- gsub(attrGrep,
               "\\1",
               gtfDF[geneRows,,drop=FALSE][[9]]);
         } else {
            printDebug("Note: No gene attributes found for:", attrName);
            attrValues <- NULL;
         }
      }));

   ## transcript attributes
   txM <- do.call(cbind, lapply(nameVector(c(txAttrNames,geneAttrNames)),
      function(attrName){
         printDebug("makeTx2geneFromGtf() tx attributes:", attrName);
         attrGrep <- paste0('^.*', attrName, ' ["]([^"]+)["].*$');
         if (igrepHas(attrGrep, gtfDF[txRows,][[9]])) {
            attrValues <- gsub(attrGrep,
               "\\1",
               gtfDF[txRows,,drop=FALSE][[9]]);
         } else {
            printDebug("Note: No tx attributes found for:", attrName);
            attrValues <- NULL;
         }
      }));
   if (verbose) {
      printDebug("makeTx2geneFromGtf() :",
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
         method="flank");
   }
   ## add transcript_id annotation
   values(threeUtrGRLdetRange@unlistData)[,"transcript_id"] <- rep(
      names(threeUtrGRLdetRange),
      lengths(threeUtrGRLdetRange));
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
   threeUtrGRLdetGeneGRL <- split(
      threeUtrGRLdetRange@unlistData,
      values(threeUtrGRLdetRange@unlistData)[,"gene_name"]);
   ## Reduce (melt) 3'UTR ranges, combining overlapping ranges per gene
   if (verbose) {
      printDebug("gencode2ale(): ",
         "Reducing ranges.");
   }
   threeUtrGRLdetGeneGRLred <- reduce(threeUtrGRLdetGeneGRL);

   ####################################################
   ## Annotate transcripts to the reduced 3'UTR ranges
   if (verbose) {
      printDebug("gencode2ale(): ",
         "Annotating reduced ranges with transcript_id per gene.");
   }
   threeUtrGRLdetGeneGRLred2 <- annotateGRLfromGRL(
      threeUtrGRLdetGeneGRLred,
      threeUtrGRLdetGeneGRL);

   ####################################################
   ## Assign ALE numbers in stranded order for each gene
   if (verbose) {
      printDebug("gencode2ale(): ",
         "Assigning stranded numbers to the ranges.");
   }
   threeUtrGRLdetGeneGRLred2 <- assignGRLexonNames(
      exonNameColname="ALE_name",
      threeUtrGRLdetGeneGRLred2, suffix="_ale");
   names(threeUtrGRLdetGeneGRLred2@unlistData) <-
      values(threeUtrGRLdetGeneGRLred2@unlistData)[,"ALE_name"];
   values(threeUtrGRLdetGeneGRLred2@unlistData)[,"score"] <- 0.5;
   retVals$aleGRL <- threeUtrGRLdetGeneGRLred2;

   ####################################################
   ## Subset ALE containing 2 or more ALEs per gene
   GencodeALEmin2 <- threeUtrGRLdetGeneGRLred2[
      lengths(threeUtrGRLdetGeneGRLred2) > 1];
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
   tx2ale <- nameVector(list2df(ale2txL));
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
