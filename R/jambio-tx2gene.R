# makeTx2geneFromTxdb()
# makeTx2geneFromGtf()


#' Make tx2gene data.frame from a GTF file
#'
#' Make tx2gene data.frame from a GTF file
#'
#' Create a transcript-to-gene data.frame from a GTF file, which is required
#' by a number of transcriptome analysis methods such as those in
#' the DEXseq package, and the limma package functions such as
#' `diffSplice()`.
#'
#' This function also only uses `data.table::fread()` and does not
#' import the full GTF file using something like Bioconductor
#' `GenomicFeatures`, simply because the data.table method is markedly
#' faster when importing only the transcript-to-gene relationship. Also, this
#' method allows the import of more annotations than are supported by the
#' typical Bioconductor `rtracklayer::import()` for GTF data.
#'
#' This function is intended to help keep all transcript data consistent by
#' using the same GTF file that is also used by other analysis tools, whether
#' those tools be based in R or more likely, outside R in a terminal
#' environment.
#'
#' For example, the GTF file could be used to:
#'
#' * run STAR sequence alignment
#'    then `Rsubread::featureCounts()` to generate a matrix of read
#'    counts per gene, transcript, or exon
#' * generate a transcript
#'    FASTA sequence file then run a kmer quantitation tool such as
#'    Salmon or Kallisto, then using `tximport::tximport()` to import
#'    results into R for downstream processing.
#'
#'
#' @param GTF `character` file name sent to `data.table::fread()`. When the
#'    file ends with ".gz", the `R.utils` package is recommended, otherwise
#'    the fallback option is to make a system call to `gzcat`
#'    to gunzip the file during the import step. Note this process fails
#'    when `gzcat` is not available in the path of the user environment.
#'    In general, the `R.utils` package is the best solution.
#' @param geneAttrNames `character` recognized attribute names
#'    as they appear in column 9 of the GTF file, for gene rows.
#'    The defaults include typical entries in Gencode, plus "range" which
#'    creates one field with format "chromosome:start-end:strand".
#' @param txAttrNames `character` vector of recognized attribute names
#'    as they appear in column 9 of the GTF file, for transcript rows.
#'    The defaults include typical entries in Gencode, plus "range" which
#'    creates one field with format "chromosome:start-end:strand".
#' @param geneFeatureType `character` value to match column 3 of the GTF
#'    file, used to define gene rows, by default "gene".
#' @param txFeatureType `character` value to match column 3 of the GTF
#'    file, used to define gene rows, by default "transcript". In some
#'    GTF files, "mRNA" is used, so either is accepted by default.
#' @param nrows `integer` number of rows to read from the GTF file, by default
#'    -1 means all rows are imported. This parameter is useful to check the
#'    results of a large GTF file using only a subset portion of the file.
#' @param zcat_command `character` name or path to zcat or gzcat executable,
#'    only used when input `GTF` is a file with `".gz"` extension, and when
#'    R package `R.utils` is not available.
#' @param verbose `logical` whether to print verbose output during processing.
#' @param ... additional arguments are ignored.
#'
#' @returns `data.frame` with colnames defined by
#'    `geneAttrNames` and `txAttrNames`.
#'
#' @family GTF functions
#'
#' @import data.table
#'
#' @export
makeTx2geneFromGtf <- function
(GTF,
 geneAttrNames=c("gene_id",
    "gene_name",
    "gene_type",
    "range"),
 txAttrNames=c("transcript_id",
    "transcript_type",
    "range"),
 geneFeatureType="gene",
 txFeatureType=c("transcript",
    "mRNA"),
 nrows=-1L,
 zcat_command="zcat",
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
   gtfDF <- readGtf(GTF=GTF,
      nrows=nrows,
      zcat_command=zcat_command,
      verbose=verbose,
      ...)

   ## Subset to clear some memory
   colnames(gtfDF) <- jamba::makeNames(rep("V", ncol(gtfDF)), suffix="");
   gtfDF <- subset(gtfDF,
      gtfDF[[3]] %in% c(geneFeatureType,
         txFeatureType));
   if (verbose > 1) {
      jamba::printDebug("makeTx2geneFromGtf(): ",
         "nrow for subset gene/tx feature types: ",
         jamba::formatInt(nrow(gtfDF)));
   }

   # Make data unique by chromosome,name,source,annotation
   # data in other columns is not relevant here.
   # This step reduces 90% rows in many files.
   dupe_rows <- duplicated(gtfDF[,c(1,2,3,9),drop=FALSE]);
   gtfDF <- subset(gtfDF, !dupe_rows);

   # Determine which rows are gene and transcript
   # TODO: recognize when geneRows,txRows are identical
   # then process them all in only one step
   geneRows <- (gtfDF[[3]] %in% geneFeatureType);
   txRows <- (gtfDF[[3]] %in% txFeatureType);
   txAttrNames <- unique(c(txAttrNames,
      geneAttrNames));
   if (all(txRows == geneRows)) {
      geneAttrNames <- NULL;
   }

   # 0.0.83.900: special attributes recognized in other columns
   coordAttrNames <- list(
      "chr"=1,
      "start"=4,
      "end"=5,
      "strand"=7,
      "range"=c(1, 4, 5, 7))

   # gene attributes
   names(geneAttrNames) <- ifelse(geneAttrNames %in% names(coordAttrNames),
      paste0("gene_", geneAttrNames),
      geneAttrNames);
   geneM <- NULL;
   if (sum(geneRows) > 0 && length(geneAttrNames) > 0) {
      geneM <- unique(getGtfAttrs(gtfDF,
         useRows=geneRows,
         attrNames=geneAttrNames,
         featureType="gene"))
      # check for duplicate gene_id
      # detect rowname
      if (length(geneM) > 0 && nrow(geneM) > 0) {
         geneid_colname <- head(jamba::provigrep(
            unique(c("^gene[_. ]*id",
               "^gene[_. ]*name",
               names(geneAttrNames))),
            colnames(geneM)), 1);
         dupe_geneid <- duplicated(geneM[[geneid_colname]]);
         if (any(dupe_geneid)) {
            if (TRUE %in% verbose) {
               dupe_geneids <- unique(geneM[[geneid_colname]][dupe_geneid]);
               jamba::printDebug("makeTx2geneFromGtf(): ",
                  "Warning: ",
                  jamba::formatInt(length(dupe_geneids)),
                  " duplicated gene IDs in column ",
                  geneid_colname, ", for example: ",
                  jamba::middle(dupe_geneids, 5));
            }
         }
         rownames(geneM) <- jamba::makeNames(geneM[[geneid_colname]]);
      }
      # head(geneM)
      # geneM2 <- getGtfAttrs(gtfDF, useRows=geneRows, attrNames=geneAttrNames, featureType="gene", matchMethod="multi")
      # head(geneM2)
   }

   ## transcript attributes
   names(txAttrNames) <- ifelse(txAttrNames %in% names(coordAttrNames),
      paste0("tx_", txAttrNames),
      txAttrNames);
   txM <- NULL;
   txid_colname <- NULL;
   if (sum(txRows) > 0 && length(txAttrNames) > 0) {
      txM <- getGtfAttrs(gtfDF, useRows=txRows, attrNames=txAttrNames, featureType="tx")
      # head(txM)
      # detect rowname
      if (length(txM) > 0 && nrow(txM) > 0) {
         txid_colname <- head(jamba::provigrep(
            unique(c("^transcript[_. ]*id",
               "^tx[_. ]*id",
               "^transcript[_. ]*name",
               "^tx[_. ]*name",
               names(txAttrNames))),
            colnames(txM)), 1);
         if (length(txid_colname) == 1) {
            rownames(txM) <- jamba::makeNames(txM[, txid_colname],
               ...);
         }
      }
   }

   # merge gene attributes into the same table
   # which only adds new information if the transcript rows did
   # not already contain this data
   if (length(geneM) > 0) {
      if (length(txM) == 0) {
         txM <- geneM
      } else {
         if (verbose) {
            jamba::printDebug("makeTx2geneFromGtf(): ",
               "Merging gene and transcript annotations.");
         }
         txM <- jamba::mergeAllXY(geneM, txM);
         # detect rowname
         txid_colname <- head(jamba::provigrep(
            unique(c("^transcript[_. ]*id",
               "^tx[_. ]*id",
               "^transcript[_. ]*name",
               "^tx[_. ]*name",
               names(txAttrNames))),
            colnames(txM)), 1);
         if (length(txid_colname) == 1) {
            rownames(txM) <- jamba::makeNames(txM[, txid_colname],
               ...);
         }
      }
   }
   if (length(txM) == 0) {
      stop("No gene or transcript data.")
   }

   return(txM);
}


# makeTx2geneFromTxdb

#' Make tx2gene data.frame from a TxDb object
#'
#' Make tx2gene data.frame from a TxDb object and an org.*.eg.db annotation
#' package.
#'
#' This function converts a Bioconductor `TxDb` annotation package (e.g.
#' `TxDb.Mmusculus.UCSC.mm10.knownGene`) into the three-column `data.frame`
#' expected by splicejam: `transcript_id`, `gene_id`, `gene_name`.
#'
#' Gene IDs in UCSC-style `TxDb` objects are NCBI Entrez IDs (ENTREZID).
#' These are mapped to gene symbols (SYMBOL) using
#' `genejam::freshenGenes()` with the supplied `ann_lib`. When an ENTREZID
#' cannot be resolved to a SYMBOL, the original ENTREZID is retained as the
#' `gene_name` (controlled by `empty_rule="original"` in `freshenGenes()`).
#' 
#' **Exception:** Some `TxDb` packages with suffix 'ensDb' use EnsEMBL
#' gene identifier, in the form 'ENSG00000001' or 'ENSMUSG00000001'.
#' In this case, `genejam::freshenGenes()` will be used to query by
#' the 'ENSEMBL2EG' annotation data, and will fallback to use the
#' EnsEMBL gene_id when not found.
#'
#' @param txdb `TxDb` object, such as those provided by Bioconductor packages
#'    like `TxDb.Mmusculus.UCSC.mm10.knownGene` or
#'    `TxDb.Hsapiens.UCSC.hg38.knownGene`.
#' @param ann_lib `character` string naming an installed Bioconductor organism
#'    annotation package (e.g. `"org.Mm.eg.db"` for mouse or
#'    `"org.Hs.eg.db"` for human). Must match the organism used in `txdb`.
#'    No validation of organism concordance is performed.
#' @param verbose `logical` whether to print verbose output during processing.
#' @param ... additional arguments are ignored.
#'
#' @returns `data.frame` with colnames `"transcript_id"`, `"gene_id"`,
#'    `"gene_name"`. `gene_id` contains the ENTREZID values from the `TxDb`,
#'    and `gene_name` contains the resolved gene SYMBOL.
#'
#' @family GTF functions
#'
#' @seealso `makeTx2geneFromGtf()` for the GTF-based equivalent,
#'    `splicejamDataFromTxDb()` for the higher-level workflow that calls this
#'    function.
#'
#' @export
makeTx2geneFromTxdb <- function
(txdb,
 ann_lib,
 verbose=FALSE,
 ...)
{
   ## Purpose: convert a TxDb + org.*.eg.db package name into the
   ## three-column data.frame (transcript_id, gene_id, gene_name) used
   ## throughout splicejam.

   ## Step 1: get transcripts grouped by gene_id (ENTREZID)
   ## Note: use.names=TRUE is not supported when by="gene", so we use the
   ## default call and extract names manually.
   if (verbose) {
      jamba::printDebug("makeTx2geneFromTxdb(): ",
         "Calling GenomicFeatures::transcriptsBy(by='gene')");
   }
   txByGene <- GenomicFeatures::transcriptsBy(txdb, by="gene")

   ## Step 2: flatten into a two-column data.frame
   tx_unlisted   <- unlist(txByGene)
   gene_ids      <- names(tx_unlisted)          # ENTREZID, repeated per tx
   tx_ids        <- tx_unlisted$tx_name         # transcript ID

   tx2gene_raw <- data.frame(
      transcript_id = tx_ids,
      gene_id       = gene_ids,
      stringsAsFactors = FALSE)
   if (verbose) {
      jamba::printDebug("makeTx2geneFromTxdb(): ",
         "tx2gene_raw nrow: ",
         jamba::formatInt(nrow(tx2gene_raw)));
   }

   ## Step 3: convert ENTREZID → SYMBOL using genejam
   unique_gene_ids <- unique(tx2gene_raw$gene_id)
   if (verbose) {
      jamba::printDebug("makeTx2geneFromTxdb(): ",
         "Calling genejam::freshenGenes() for ",
         jamba::formatInt(length(unique_gene_ids)),
         " unique gene_id values.");
   }
   #
   # Todo: Consider AnnotationDbi::metadata(txdb)
   # name='Type of Gene ID' and value is
   # value == 'Entrez Gene ID' or
   # value == 'Ensembl gene ID'
   #
   if (any(grepl("^ENS", unique_gene_ids))) {
      if (verbose) {
         jamba::printDebug("makeTx2geneFromTxdb(): ",
            "Querying by ENSEMBL");
      }
      # use EnsEMBL gene_id approach
      use_df  <- data.frame(ENSEMBL=unique_gene_ids,
         stringsAsFactors=FALSE)
      gene_df <- genejam::freshenGenes(use_df,
         ann_lib=ann_lib,
         try_list=c("ENSEMBL2EG",
            "ACCNUM2EG",
            "SYMBOL2EG"),
         empty_rule="original")
      # gene_df has columns ENSEMBL + SYMBOL
      symbol_map <- setNames(gene_df$SYMBOL, gene_df$ENSEMBL)
   } else {
      if (verbose) {
         jamba::printDebug("makeTx2geneFromTxdb(): ",
            "Querying by ENTREZID");
      }
      # Use ENTREZID approach
      use_df  <- data.frame(ENTREZID=unique_gene_ids,
         stringsAsFactors=FALSE)
      gene_df <- genejam::freshenGenes(use_df,
         ann_lib=ann_lib,
         empty_rule="original")
      # gene_df has columns ENTREZID + SYMBOL
      symbol_map <- setNames(gene_df$SYMBOL, gene_df$ENTREZID)
   }

   ## Step 4: attach gene_name; fall back to gene_id when SYMBOL is absent
   tx2gene_raw$gene_name <- symbol_map[tx2gene_raw$gene_id]
   na_idx <- is.na(tx2gene_raw$gene_name)
   if (any(na_idx)) {
      if (verbose) {
         jamba::printDebug("makeTx2geneFromTxdb(): ",
            jamba::formatInt(sum(na_idx)),
            " transcripts had no SYMBOL; retaining gene_id as gene_name.");
      }
      tx2gene_raw$gene_name[na_idx] <- tx2gene_raw$gene_id[na_idx]
   }

   ## Return canonical column order
   tx2gene_raw[, c("transcript_id", "gene_id", "gene_name"), drop=FALSE]
}
