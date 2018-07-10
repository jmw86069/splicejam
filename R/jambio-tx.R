
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
#'    * to run STAR sequence alignment
#'    then \code{Rsubread::featureCounts()} to generate a matrix of read
#'    counts per gene, transcript, or exon; or
#'    * to generate a transcript
#'    FASTA sequence file then run a kmer quantitation tool such as
#'    Salmon or Kallisto, then using \code{tximport::tximport()} to import
#'    results into R for downstream processing.
#'
#' @param GTF character file name sent to \code{data.table::fread()}. When the
#'    file ends with ".gz" then gzcat is used to gunzip the file during the
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
#' @return data.frame with colnames indicated by the values in
#' \code{geneAttrNames} and \code{txAttrNames}.
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
         col.names,
         nrows=nrows);
   } else {
      gtfDF <- fread(GTF,
         sep="\t",
         autostart=20,
         col.names,
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


