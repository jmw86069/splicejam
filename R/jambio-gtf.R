
#' read gtf/gff3 file
#'
#' read gtf/gff3 file
#'
#' @return `data.frame` representing tab-delimited data stored in the
#'    gtf or gff3 file.
#'
#' @param GTF `character` file name sent to `data.table::fread()`. When the
#'    file ends with ".gz", the `R.utils` package is recommended, otherwise
#'    the fallback option is to make a system call to `gzcat`
#'    to gunzip the file during the import step. Note this process fails
#'    when `gzcat` is not available in the path of the user environment.
#'    In general, the `R.utils` package is the best solution.
#' @param nrows `integer` number of rows to read from the GTF file, by default
#'    -1 means all rows are imported. This parameter is useful to check the
#'    results of a large GTF file using only a subset portion of the file.
#' @param zcat_command `character` name or path to zcat or gzcat executable,
#'    only used when input `GTF` is a file with `".gz"` extension, and when
#'    R package `R.utils` is not available.
#' @param verbose `logical` whether to print verbose output during processing.
#' @param ... additional arguments are ignored.
#'
#' @export
readGtf <- function
(GTF,
 nrows=-1,
 zcat_command="zcat",
 verbose=FALSE,
 ...)
{
   #
   if (verbose) {
      jamba::printDebug("readGtf(): ",
         "reading GTF file:",
         GTF);
   }
   if (jamba::igrepHas("[.]gz$", GTF) &&
         !jamba::check_pkg_installed("R.utils")) {
      if (verbose) {
         jamba::printDebug("readGtf(): ",
            "using commandline call to ",
            zcat_command,
            " to decompress gtf.gz. ",
            "Install the 'R.utils' package to avoid this step.");
      }
      gtfDF <- data.table::fread(cmd=paste(zcat_command, GTF),
         sep="\t",
         header=FALSE,
         nrows=nrows,
         data.table=FALSE);
   } else {
      if (verbose) {
         jamba::printDebug("readGtf(): ",
            "using native data.table::fread().");
      }
      gtfDF <- data.table::fread(GTF,
         sep="\t",
         header=FALSE,
         nrows=nrows,
         data.table=FALSE);
   }
   return(gtfDF);
}


#' describe gtf/gff3 attribute names
#'
#' describe gtf/gff3 attribute names
#'
#' @export
describeGtfAttrNames <- function
(GTF,
 geneFeatureType="gene",
 txFeatureType=c("transcript",
    "mRNA"),
 nrows=10000,
 zcat_command="zcat",
 verbose=FALSE,
 ...)
{
   #
   gtfDF <- readGtf(GTF=GTF,
      nrows=nrows,
      zcat_command=zcat_command,
      verbose=verbose,
      ...)

   source_names <- unique(gtfDF[[2]]);
   feature_names <- unique(gtfDF[[3]]);
   source_attrs <- lapply(split(gtfDF, gtfDF[[3]]), function(idf){
      idf <- head(idf, 10)
      iattrs <- strsplit(idf[[9]], "[ ]*[;][ ]*");
      if (jamba::igrepHas('=', unlist(iattrs))) {
         # gff
         iattrs <- split(
            gsub('"', '',
               sub("[ ]*=[ ]*", "!!", unlist(iattrs))),
            rep(seq_along(iattrs), lengths(iattrs)))
      } else {
         # gtf
         iattrs <- split(
            gsub('"', '',
               sub(" ", "!!", unlist(iattrs))),
            rep(seq_along(iattrs), lengths(iattrs)))
      }
      ilist <- lapply(iattrs, function(iattr){
         jamba::nameVector(jamba::rbindList(head(strsplit(iattr, "!!"), 20))[,2:1,drop=FALSE])
      })
      inames <- unique(unlist(lapply(ilist, names)))
      jdf <- data.frame(check.names=FALSE,
         jamba::rbindList(lapply(ilist, function(ivector){
            jamba::nameVector(ivector[inames], inames)
         })))
      # as.data.frame(lapply(jamba::nameVectorN(jdf), function(icol){
      #    head(c(setdiff(jdf[[icol]], c(NA, "N/A", "NA")), rep(NA, 3)), 3)
      # }))
      jdf
   })
   return(source_attrs[feature_names])
}
