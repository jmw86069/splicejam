
#' read gtf/gff3 file
#'
#' read gtf/gff3 file
#'
#' @family jam gtf functions
#'
#' @returns `data.frame` representing tab-delimited data stored in the
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


#' describe gtf/gff3 attribute names by feature type
#'
#' describe gtf/gff3 attribute names by feature type
#'
#' * Note that when the "name" in a name/value pair is repeated,
#' the first instance retains the name, while subsequent instances
#' are versioned by `jamba::makeNames(x, renameFirst=FALSE)`.
#' For example `"tag"` may appear multiple times, the resulting colnames
#' will become: `c("tag", "tag_v1", "tag_v2")`.
#'
#' @returns `list` named by `c(geneFeatureType, txFeatureType` with
#'    `data.frame` objects which have split the name/value pairs
#'    into columns. Each `data.frame` may have different columns,
#'    using the observed name/value pair data.
#'
#' @family jam gtf functions
#'
#' @param GTF `character` path to GTF or GFF3 file, or `data.frame`
#'    containing GTF or GFF3 data.
#' @param geneFeatureType,txFeatureType `character` vectors with values in
#'    column 3 of the GTF or GFF3 file, used to subset then split the
#'    output data.
#'    * Return all feature types by providing any of these terms:
#'    ".", "any", "all"
#' @param nrows `integer` max number of rows to process. For this purpose,
#'    summarizing the type of data seen for each feature type, a subset
#'    of rows is usually sufficient.
#' @param maxNper `integer` default 10,  number of entries retained
#'    within each feature type. Set to `Inf` to retain all data.
#'    As a brief summary, 10 is sufficient to show typical content.
#' @param maxAttrs `integer` default 50, maximum attributes to retain
#'    for each entry. Only in rare cases are more than 50 attributes
#'    present for one record, and typicaly these are the rare cases
#'    where those attributes were not necessary for annotation purposes.
#' @param zcat_command `character` name or path to the `zcat` command or
#'    equivalent, used only when the R package 'R.utils' is not installed,
#'    and the input GTF has `.gz` file extension.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
describeGtfAttrNames <- function
(GTF,
 geneFeatureType="gene",
 txFeatureType=c("transcript",
    "mRNA"),
 nrows=10000,
 maxNper=10,
 maxAttrs=50,
 zcat_command="zcat",
 verbose=FALSE,
 ...)
{
   #
   if (inherits(GTF, "data.frame")) {
      gtfDF <- GTF;
   } else {
      gtfDF <- readGtf(GTF=GTF,
         nrows=nrows,
         zcat_command=zcat_command,
         verbose=verbose,
         ...)
   }
   if (length(nrows) == 1 && nrow(gtfDF) > nrows) {
      gtfDF <- head(gtfDF, nrows)
   }

   source_names <- unique(gtfDF[[2]]);
   feature_names <- unique(gtfDF[[3]]);

   # subset for matchinf feature types
   subsetTypes <- c(geneFeatureType, txFeatureType);
   if (!any(c(".", "all", "any"), tolower(subsetTypes))) {
      gtfDF <- subset(gtfDF, gtfDF[[3]] %in% c(geneFeatureType,
         txFeatureType))
      gtfDF[[3]] <- factor(gtfDF[[3]],
         levels=unique(subsetTypes));
   }

   # split by feature type
   source_attrs <- lapply(split(gtfDF, gtfDF[[3]]), function(idf){
      idf <- head(idf, maxNper)

      # remove trailing ";" whitespace, newline, linefeed
      idf[[9]] <- gsub("[; \r\n]*$", "", idf[[9]]);
      # split multiple fields using ";" delimiter
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
         jamba::nameVector(renameFirst=FALSE,
            jamba::rbindList(
               head(strsplit(iattr, "!!"), maxAttrs))[, 2:1, drop=FALSE])
      })
      inames <- unique(unlist(lapply(ilist, names)))
      jdf <- data.frame(check.names=FALSE,
         jamba::rbindList(lapply(ilist, function(ivector){
            jamba::nameVector(ivector[inames], inames)
         })))
      jdf
   })
   return(source_attrs);
}

#' Grab GTF or GFF3 attributes into a data.frame
#'
#' Grab GTF or GFF3 attributes into a data.frame
#'
#' The input retrieves data for known attribute names, although when
#' an attribute name is not present, it is silently ignored.
#' However, to see available attribute names, try `describeGtfAttrNames()`.
#'
#' @returns `data.frame` with colnames derived from the GTF or GFF3 data.
#'    When there is no data recognized, it returns `NULL`, for example
#'    when `useRows` is entirely `FALSE` or no `attrNames` are recognized
#'    in the input `gtfDF`.
#'
#' @family jam gtf functions
#'
#' @param gtfDF `data.frame` with GTF or GFF3 formatted data
#' @param useRows `logical`, default NULL, rows in `gtfDF` to use for analysis,
#'    or NULL to use all rows supplied in `gtfDF`.
#' @param attrNames `character` vector of attribute names, matching one
#'    of two types of fields in the `gtfDF` data:
#'    1. Name/value pairs in column 9, typically used for miscellaneous
#'    annotations. The formats are intended to represent accepted formats
#'    according to the GTF and GFF3 specifications, for example:
#'       * GTF: `attrName "value"`
#'       * GTF: `attrName "value";`
#'       * GTF: `attrName value`
#'       * GTF: `attrName value;`
#'       * GFF3: `attrName=value`
#'       * GFF3: `attrName=value;`
#'       * GFF3: `attrName="value"`
#'       * GFF3: `attrName="value";`
#'
#'    2. Recognized column from within the GTF or GFF3 format:
#'       * **chr** - column 1
#'       * **start** - column 4
#'       * **end** - column 5
#'       * **strand** - column 7
#'       * **range** - columns 1,4,5,7 concatenated into one field,
#'       for example: "chr:start-end:strand" or "chr1:10000-11000:+"
#' @param featureType `character`, default NULL, with feature type used
#'    as a prefix for coordinate/range attributes, and otherwise only
#'    used in messaging when `verbose=TRUE`.
#' @param ignore.case `logical` default FALSE, passed to `grep()` and
#'    `gsub()` to enable case-insensitive matching of attribute names.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
getGtfAttrs <- function
(gtfDF,
 useRows=NULL,
 attrNames,
 featureType=NULL,
 ignore.case=FALSE,
 matchMethod=c("once", "multi"),
 verbose=FALSE,
 ...)
{
   #
   if (nrow(gtfDF) == 0) {
      return(NULL)
   }
   if (length(useRows) == 0) {
      useRows <- rep(TRUE, nrow(gtfDF))
   }
   if (all(useRows %in% FALSE)) {
      return(NULL)
   }

   # two closely related approaches
   matchMethod <- match.arg(matchMethod);

   # 0.0.83.900: special attributes recognized in other columns
   coordAttrNames <- list(
      "chr"=1,
      "start"=4,
      "end"=5,
      "strand"=7,
      "range"=c(1, 4, 5, 7))

   # attributes
   usePrefix <- "";
   if (length(featureType) == 1 && nchar(featureType) > 0) {
      usePrefix <- paste0(featureType, "_");
      featureType <- paste0(featureType, " ")
   } else {
      featureType <- ""
   }
   names(attrNames) <- ifelse(attrNames %in% names(coordAttrNames),
      paste0(usePrefix, attrNames),
      attrNames);
   attrM <- NULL;

   # iterate attrNames
   if (sum(useRows) > 0 && length(attrNames) > 0) {
      # Experimental batch pivot processing
      attrDF <- NULL;
      if ("multi" %in% matchMethod) {
         x1 <- gtfDF[useRows, , drop=FALSE][[9]]
         x1l <- strsplit(x1, "[ ]*;[ ]*")
         x1v <- unlist(x1l);
         x1vf <- rep(seq_along(x1l), lengths(x1l))
         x1df <- data.frame(x1vf=x1vf,
            jamba::rbindList(strsplit(sub("[ =]+", "!!", x1v), "!!"),
               newColnames=c("attr", "value")))
         # subset for attrNames
         if (TRUE %in% ignore.case) {
            x1df <- subset(x1df, tolower(attr) %in% tolower(attrNames));
         } else {
            x1df <- subset(x1df, attr %in% attrNames);
         }
         if (nrow(x1df) > 0) {
            x1df$value <- gsub("^\"(.*)\"$", "\\1", x1df$value)
            # handle dupes
            x1nametag <- paste0(x1df[[1]], "!!", x1df[[2]]);
            x1df_dupe <- duplicated(x1nametag);
            if (any(x1df_dupe)) {
               x1df_dupeattr <- unique(gsub("^.+!!", "", x1nametag[x1df_dupe]))
               x1df_whichdupe <- which(x1df$attr %in% x1df_dupeattr)
               x1df_undupe <- gsub("^.+!!", "",
                  jamba::makeNames(x1nametag[x1df_whichdupe],
                     renameFirst=FALSE))
               x1df$attr[x1df_whichdupe] <- x1df_undupe;
            }
            x1df_wide <- tidyr::pivot_wider(x1df,
               names_sort=FALSE,
               # subset(x1df, !attr %in% "tag"),
               id_cols="x1vf",
               names_sep=",", values_fill=NA,
               names_from="attr", values_from="value")
            x1df_wide_df <- data.frame(check.names=FALSE, x1df_wide)
            keep_colnames <- jamba::provigrep(paste0("^", attrNames, "(|_v[0-9]+)$"),
               colnames(x1df_wide_df));
            attrDF <- x1df_wide_df[, keep_colnames, drop=FALSE];
         }
      }
      attrL <- jamba::rmNULL(lapply(attrNames, function(attrName){
         if (verbose) {
            jamba::printDebug("getGtfAttrs(): ",
               c(featureType, "attributes: "),
               attrName,
               sep="")
         }
         # 0.0.83.900: option for coordinate columns
         if (attrName %in% names(coordAttrNames)) {
            useCols <- coordAttrNames[[attrName]];
            if (length(useCols) == 1) {
               attrValues <- gtfDF[useRows, , drop=FALSE][[useCols]];
            } else {
               attrValues <- jamba::pasteByRow(
                  gtfDF[useRows, useCols, drop=FALSE],
                  sep="!");
               attrValues <- sub("!", ":", attrValues)
               attrValues <- sub("!", "-", attrValues)
               attrValues <- sub("!", ":", attrValues)
            }
            return(attrValues);
         }
         if ("multi" %in% matchMethod) {
            return(NULL)
         }
         # 0.0.77.900: alternative grep pattern tolerant to gtf and gff3
         # gtf format example:   attrName "value";
         # gff3 format example:  attrName=value;
         attrGrep <- paste0('(^.+;[ ]*|^[ ]*)',
            attrName,
            '[ =]+["]*([^";]+)[ ]*($|[";].*$)');

         ## only test up to the first 2000 rows
         testNum <- sum(useRows);
         if (testNum > 50000) {
            testNum <- 2000;
         }
         testVals <- head(gtfDF[useRows, , drop=FALSE][[9]], testNum)
         if (jamba::igrepHas(attrGrep, testVals, ignore.case=ignore.case)) {
            attrFound <- grepl(attrGrep,
               gtfDF[useRows, , drop=FALSE][[9]],
               ignore.case=ignore.case)
            attrValues <- rep("", sum(useRows))
            attrValues[attrFound] <- gsub(attrGrep,
               "\\2",
               gtfDF[useRows, , drop=FALSE][[9]][attrFound],
               ignore.case=ignore.case)
            if (TRUE %in% verbose && any(!attrFound)) {
               jamba::printDebug("getGtfAttrs(): ",
                  "Some attrValues (",
                  jamba::formatInt(sum(!attrFound)),
                  ") were empty for attrName:",
                  attrName);
            }
         } else {
            if (verbose) {
               jamba::printDebug("getGtfAttrs(): ",
                  "All attrValues were empty for attrName:",
                  attrName);
            }
            attrValues <- NULL;
         }
         attrValues;
      }));

      if ("multi" %in% matchMethod) {
         if (length(attrDF) == 0) {
            if (length(attrL) == 0) {
               return(NULL)
            }
         }
      }
      if (length(attrDF) == 0 && length(attrL) == 0) {
         return(NULL)
      }
      # assemble into columns
      if (length(attrL) > 0) {
         attrM <- data.frame(check.names=FALSE,
            do.call(cbind, attrL))
         attrM <- unique(attrM);
      }
      if (length(attrDF) > 0) {
         attrM <- data.frame(check.names=FALSE,
            attrDF,
            attrM);
         keep_colnames <- jamba::provigrep(
            paste0("^", names(attrNames), "(|_v[0-9]+)$"),
            colnames(attrM));
         attrM <- attrM[, keep_colnames, drop=FALSE];
      }
   }
   return(attrM)
}

