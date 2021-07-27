
#' Convert PSL alignment to data.frame
#'
#' Convert PSL alignment to data.frame
#'
#' This function takes PSL alignment format, as produced by BLAT,
#' and converts to data.frame. The driving reason for this function
#' is that current methods to convert PSL to other formats lose
#' some useful information, for example conversion to BED12 format
#' loses the query alignment coordinates, in favor of storing
#' only the reference coordinates.
#'
#' To supply text as input use `base::textConnection()` to wrap
#' a text connection around the input text.
#'
#' @family jam data import functions
#'
#' @param psl file or other connection compatible with `base::readLines()`
#'    of data in PSL alignment format.
#' @param ... additional arguments are ignored.
#'
#' @examples
#'
#' psls <- c(
#'    "2252\t0\t0\t0\t0\t0\t5\t56285\t+\tQuery_Sequence\t2252\t0\t2252\tReference_Sequence\t99095\t1992\t60529\t6\t227,86,71,79,77,1712,\t0,227,313,384,463,540,\t1992,6931,11020,39871,45861,58817,",
#'    "664\t0\t0\t0\t1\t3\t1\t1\t+\tQuery_Sequence\t2252\t1291\t1958\tReference_Sequence\t99095\t58842\t59507\t2\t657,7,\t1291,1951,\t58842,59500,"
#'    );
#' psl2df(textConnection(psls), verbose=TRUE)
#'
#' @export
psl2df <- function
(psl,
 verbose=FALSE,
 ...)
{
   ## Purpose is to convert psl format to data.frame
   #pslLines <- gsub("[' ]+", "",
   #   readLines("Adam19_KO_e5_Tomato_mRNA.psl")[-1*c(1,2,5)]);
   pslHeaderCheck <- readLines(psl, n=1);
   if (jamba::igrepHas("^pslayout", pslHeaderCheck)) {
      if (verbose) {
         printDebug("psl2df(): ",
            "Detected psLayout header line.");
      }
      pslLines <- gsub("[' ]+", "",
         readLines(psl)[-1*c(1,4)]);
      psldf1 <- jamba::rbindList(
         strsplit(
            head(pslLines, 2),
            "\t"));
      pslHeader <- jamba::pasteByRow(t(psldf1), sep="");
      psldf <- read.table(text=tail(pslLines, -2),
         sep="\t",
         stringsAsFactors=FALSE,
         header=FALSE);
   } else {
      if (verbose) {
         printDebug("psl2df(): ",
            "Detected no psLayout header line.");
      }
      pslHeader <- c("match", "mis-match", "rep.match",
         "Ns", "Qgapcount", "Qgapbases",
         "Tgapcount", "Tgapbases", "strand",
         "Qname", "Qsize", "Qstart", "Qend",
         "Tname", "Tsize", "Tstart", "Tend",
         "blockcount", "blockSizes",
         "qStarts", "tStarts");
      pslLines <- gsub("[' ]+", "",
         c(pslHeaderCheck,
            readLines(psl)));
      printDebug("length(pslLines):", length(pslLines));
      if (jamba::igrepHas("^[-]+$", pslLines)) {
         pslN <- max(grep("^[-]+$", pslLines));
         if (verbose) {
            printDebug("psl2df(): ",
               "Detected ----- divider on line:",
               pslN);
         }
         pslLines <- pslLines[-1*seq_len(pslN)];
      }
      psldf <- read.table(text=pslLines,
         sep="\t",
         stringsAsFactors=FALSE,
         header=FALSE);
   }
   colnames(psldf) <- pslHeader;
   return(psldf);
}
