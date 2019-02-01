
#' Define experimental contrasts from sample groups
#'
#' Define experimental contrasts from sample groups
#'
#' This function is intended to define statistical contrasts
#' that compare one factor at a time. For two-factor designs,
#' it will create two-way contrasts, defined as the contrast
#' of pairwise contrasts.
#'
#' @return list of data matrices: `iDesign` numeric design matrix;
#' `iContrasts` numeric contrast matrix; `contrastNames` data.frame
#' showing the full factor breakdown, with colnames "contrastName"
#' which shows a text contrast suitable for use in `limma::makeContrasts()`.
#' When `returnDesign=FALSE` the output is only the `contrastNames`
#' data.frame.
#'
#' @param iFactors vector of sample groups with one entry per sample,
#'    or data.frame whose colnames are experimental factors, and rows
#'    are samples.
#' @param groupColumns character vector or NULL, to define an optional
#'    subset of colnames when `iFactors` is a data.frame.
#' @param iSamples character vector or NULL, optionally used to subset
#'    the sample identifiers used in subsequent steps. Note that only
#'    groups and contrasts that contain samples will be defined.
#' @param iDesign optional numeric design matrix, an optional method of
#'    defining sample-to-group mapping.
#' @param factorOrder integer vector, optionally used to define the
#'    order of factor contrasts when there are multiple experimental
#'    factors. It can be helpful to force a secondary factor to be
#'    compared before a primary factor especially in two-way contrasts.
#' @param omitGrep character grep pattern used to exclude secondary
#'    factors from contrasts, mainly used internally by this function.
#' @param maxDepth integer value, the maximum number of factor "depth"
#'    to define contrasts, for example `maxDepth=2` will define two-way
#'    contrasts, `maxDepth=1` will only define one-way contrasts.
#' @param currentDepth integer value used internally by `groups2contrasts()`
#'    for iterative operations.
#' @param factorSep,contrastSep character values used as delimiter in
#'    factor and contrast names, respectively.
#' @param renameFirstDepth logical used internally for iterative calls
#'    to `groups2contrasts()`.
#' @param returnDesign logical indicating whether to return the full
#'    set of design (`iDesign`), contrast (`iContrasts`) matrices,
#'    in addition to the `contrastNames` data.frame.
#' @param removePairs list of pairwise vectors of factors which should
#'    not be compared, or NULL to include all comparisons.
#' @param makeUnique logical indicating whether to make output
#'    contrasts unique.
#' @param addContrastNamesDF data.frame or NULL, optionally used to append
#'    to the calculated `contrastNames` data.frame, useful to add custom
#'    contrasts.
#' @param preControlTerms character vector or NULL, optionally used to
#'    help define factor order, for example `preControlTerms=c("WT")` would
#'    help order `"WT"` before `"KO"` when defining control factor levels,
#'    so the resulting contrasts would become `"KO-WT"`. This vector should
#'    contain the factor levels that should be used as the preferred
#'    control term in each contrast, where the earlier terms are preferred.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @family jam RNA-seq functions
#' @family jam design functions
#'
#' @examples
#' # first define a vector of sample groups
#' iGroups <- jamba::nameVector(paste(rep(c("WT", "KO"), each=6),
#'    rep(c("Control", "Treated"), each=3),
#'    sep="_"));
#' iGroups <- factor(iGroups, levels=unique(iGroups));
#' iGroups;
#'
#' iDesignL <- groups2contrasts(iGroups, returnDesign=TRUE);
#' iDesignL$iDesign;
#' iDesignL$iContrasts;
#'
#' # now you can visualize the samples used in each contrast
#' iDesignL$iDesign %*%  iDesignL$iContrasts;
#'
#' # you can adjust the order of factor levels per comparison
#' groups2contrasts(as.character(iGroups))$contrastName
#'
#' groups2contrasts(as.character(iGroups), preControlTerms=c("WT"), factorOrder=2:1)$contrastName
#'
#' @export
groups2contrasts <- function
(iFactors,
 groupColumns=NULL,
 iSamples=NULL,
 iDesign=NULL,
 factorOrder=NULL,
 omitGrep="[-,]",
 maxDepth=2,
 currentDepth=1,
 factorSep="_",
 contrastSep="-",
 renameFirstDepth=TRUE,
 returnDesign=FALSE,
 removePairs=NULL,
 makeUnique=TRUE,
 addContrastNamesDF=NULL,
 preControlTerms=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to take a data.frame, whose rows are groups,
   ## and whose columns are factors with factor levels as column values,
   ## and generate pairwise contrast names where only one factor changes
   ## at a time
   ##
   ## iFactors can be one of the following:
   ##
   ## - data.frame whose columns represent each statistical factor,
   ## whose values are either character, numeric, or factor, the latter
   ## can be ordered in order to provide preference to control groups.
   ##
   ## - vector of character strings representing each group,
   ## where the factors are separated by factorSep, e.g. "WT_Dex", "NT_Veh"
   ##
   ## - iDesign matrix whose colnames represent group names, and rownames
   ## represent samples present in those groups.
   ##
   ## - allNorm list object, with "targets" containing a data.frame of sample
   ## annotations, and groupColumns defines the columns to use for grouping.
   ##
   ## - removePairs is a list of vectors, where each vector is expected to
   ## contain two elements representing two factor levels not to be compared.
   ## For example, an experiment with Control, NTC, Vehicle, Dex, might not
   ## want to compare NTC-Control, Vehicle-Control, Dex-Control,
   ## removePairs <- list(c("NTC","Control"),c("Vehicle","Control"),c("Dex","Control"));
   ##
   ## TODO: enable removePairs to filter out contrasts after they are defined,
   ## for example c("NTC,Control", "d0") would remove the contrast NTC_d0-Control-d0
   ##
   ## makeUnique=TRUE will only return one entry for each set of factors compared,
   ## which will remove cases where factor 2 is tested, then factor 1 tested as an
   ## interaction; if factor 1 and factor 2 are already represented in another
   ## interaction contrast.
   ##
   ## Ultimately a table of experiment design is created, with number of columns
   ## equal to the number of factors. By default the contrasts are applied for
   ## each factor in order of colnames, but factorOrder can be used to specify
   ## a custom order. This change can affect the way two-way contrasts are
   ## defined, by forcing the first/internal contrast to use a particular
   ## factor. In theory the result is simply aesthetic, as the underlying
   ## significance of the two-way comparison will be identical. But if not
   ## for aesthetics, what are we doing?
   ##
   ## TODO: fix issue when one column contains numeric values instead of
   ## character or factor, e.g. when "Time" contains c(15,45).
   ## One solution is convert to factor, then proceed.
   if (!suppressPackageStartupMessages(require(limma))) {
      stop("limma is required for groups2contrasts()).");
   }
   sample2group <- NULL;
   #iDesign <- NULL;

   ## Handle removePairs by expanding to both orientations of contrast
   if (!is.null(removePairs)) {
      removePairsM <- jamba::rbindList(removePairs);
      removePairsFull <- c(
         jamba::pasteByRow(
            removePairsM[,seq_len(ncol(removePairsM)),drop=FALSE], sep=","),
         jamba::pasteByRow(
            removePairsM[,rev(seq_len(ncol(removePairsM))),drop=FALSE], sep=","));
      if (verbose >= 2) {
         jamba::printDebug("groups2contrasts(): ",
            "removePairsFull:");
         print(removePairsFull);
      }
   }

   ## Special case where one data.frame column is sent, which is delimited.
   ## Mainly we treat as a vector, except that we keep the rownames
   ## so we can derive iSamples.
   if (igrepHas("data.frame", class(iFactors)) &&
         ncol(iFactors) == 1) {
      iFactors <- jamba::nameVector(iFactors[,1], rownames(iFactors));
   }

   if (igrepHas("factor|character", class(iFactors))) {
      #####################################################
      ## Vector input
      ##
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "splitting vector into groups");
      }
      if (!is.null(names(iFactors)) && !is.null(iSamples)) {
         iFactors <- iFactors[match(iSamples, names(iFactors))];
      }
      if (igrepHas("factor", class(iFactors))) {
         iFactorsL <- strsplitOrdered(iFactors, factorSep);
         names(iFactorsL) <- names(iFactors);
         iFactorsLevels <- levels(iFactorsL[[1]]);
         iFactors <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            jamba::rbindList(strsplit(as.character(iFactors),
               factorSep)));
         rownames(iFactors) <- names(iFactorsL);
         for (i in seq_len(ncol(iFactors))) {
            iFactors[,i] <- factor(iFactors[,i],
               levels=intersect(iFactorsLevels, iFactors[,i]));
         }
      } else {
         iFactors <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            jamba::rbindList(strsplit(iFactors, factorSep)));
         ## Convert each column to factor for ordering
         for (iCol in seq_len(ncol(iFactors))) {
            if (igrepHas("[a-z]", iFactors[,iCol])) {
               iFactors[,iCol] <- factor(iFactors[,iCol],
                  levels=sortSamples(unique(iFactors[,iCol]),
                     preControlTerms=preControlTerms,
                     ...));
            }
         }
      }
      if (length(groupColumns) > 0) {
         colnames(iFactors) <- jamba::makeNames(rep(groupColumns,
            length.out=ncol(iFactors)),
            renameFirst=FALSE);
      }
      if (length(rownames(iFactors)) == 0) {
         rownames(iFactors) <- jamba::makeNames(
            jamba::pasteByRow(iFactors, sep=factorSep),
            suffix="_rep");
      }
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "iFactors:");
         print(head(iFactors, 40));
      }
      if (returnDesign) {
         ## Assume for now sample rows and group columns
         rowGroups <- jamba::pasteByRowOrdered(iFactors, sep=factorSep);
         sample2group <- split(rownames(iFactors), rowGroups);
         if (length(iDesign) == 0) {
            iDesign <- list2im(sample2group)[rownames(iFactors),levels(rowGroups),drop=FALSE];
         }
      }
   } else if (igrepHas("data.frame", class(iFactors))) {
      #####################################################
      ## data.frame input
      ##
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "using existing data.frame");
      }
      if (length(rownames(iFactors)) > 0) {
         if (length(iSamples) > 0) {
            iFactors <- iFactors[match(iSamples, rownames(iFactors)),,drop=FALSE];
            if (verbose) {
               jamba::printDebug("groups2contrasts(): ",
                  "Specifying iFactors[iSamples,]");
               print(head(iFactors));
            }
         } else {
            iSamples <- rownames(iFactors);
         }
      }
      if (length(groupColumns) > 0 && all(groupColumns %in% colnames(iFactors))) {
         iFactors <- iFactors[,groupColumns,drop=FALSE];
         if (verbose) {
            jamba::printDebug("groups2contrasts(): ",
               "Specifying iFactors[,groupColumns,drop=FALSE]");
            jamba::printDebug("groups2contrasts(): ",
               "groupColumns:",
               groupColumns);
         }
      } else if (length(colnames(iFactors)) == 0) {
         colnames(iFactors) <- jamba::makeNames(rep("groupFactor",
            length.out=ncol(iFactors)),
            suffix="_");
      }
      ## Do quick mixedSortDF() which will order the data.frame,
      ## to influence this sorting, use factors in character columns
      iFactors <- mixedSortDF(iFactors);
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "iFactors:");
         print(head(iFactors));
      }

      ## rowGroups is the unique set of group names, used to keep the original order
      rowGroups <- jamba::pasteByRowOrdered(iFactors,
         sep=factorSep);
      if (length(rownames(iFactors)) == 0) {
         rownames(iFactors) <- jamba::makeNames(rowGroups,
            suffix="_rep");
      }
      ## Assume for now sample rows and group columns
      sample2group <- split(rownames(iFactors), rowGroups);
      if (length(iDesign) == 0) {
         iDesign <- list2im(sample2group)[rownames(iFactors),unique(rowGroups),drop=FALSE];
         if (!is.null(iSamples)) {
            iDesign <- iDesign[match(iSamples, rownames(iDesign)),,drop=FALSE];
         }
      } else {
         if (length(iSamples) > 0) {
            iDesign <- iDesign[match(iSamples, rownames(iDesign)),,drop=FALSE];
         }
      }
   } else if (igrepHas("matrix", class(iFactors)) && all(c(0,1) %in% iFactors)) {
      ##################################
      ## iDesign input
      ##
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "converting iDesign into iFactors data.frame");
      }
      ## Assume for now, iDesign matrix with sample rows and group columns
      sample2group <- split(rownames(iFactors), sapply(seq_len(nrow(iFactors)), function(i){
         colnames(iFactors)[which(iFactors[i,] != 0)];
      }));
      iDesign <- list2im(sample2group)[rownames(iFactors),names(sample2group)];
      iFactorsCols <- colnames(iFactors);
      iFactors <- jamba::rbindList(strsplit(iFactorsCols, factorSep));
      if (!is.null(groupColumns)) {
         colnames(iFactors) <- jamba::makeNames(rep(groupColumns, length.out=ncol(iFactors)),
            renameFirst=FALSE);
      } else {
         colnames(iFactors) <- jamba::makeNames(rep("groupFactor", length.out=ncol(iFactors)),
            suffix="_");
      }
      rownames(iFactors) <- jamba::pasteByRow(iFactors, sep=factorSep);
   }
   if (verbose >= 2) {
      jamba::printDebug("groups2contrasts(): ",
         "iFactors:");
      print(head(iFactors));
      if (!is.null(sample2group)) {
         jamba::printDebug("sample2group:");
         print(head(sample2group));
      }
   }

   ##########################################################
   ## First check to make sure the iFactors values are unique
   ## and if not, use only unique entries
   iContrastGroupsUse <- colnames(iFactors);
   iFactorsV <- jamba::pasteByRow(iFactors, sep=factorSep);
   iKeepRows <- match(unique(iFactorsV), iFactorsV);
   iFactors <- iFactors[iKeepRows,,drop=FALSE];
   if (renameFirstDepth && currentDepth==1) {
      rownames(iFactors) <- jamba::pasteByRow(iFactors, sep=factorSep);
   }

   if (verbose) {
      jamba::printDebug("groups2contrasts(): ",
         "iFactors:");
      print(head(iFactors));
   }


   if (verbose) {
      jamba::printDebug("groups2contrasts(): ",
         "currentDepth:",
         currentDepth);
   }

   ##########################################################
   ## Iterate each factor in order, and create valid contrasts
   ## Note: we allow applying contrasts in a different order than the
   ## columns in iFactor, if !is.null(factorOrder)
   ##
   if (is.null(factorOrder)) {
      factorOrder <- seq_along(colnames(iFactors));
   }
   ##
   iContrastNames <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      jamba::rbindList(lapply(factorOrder, function(iChange){
         if (verbose) {
            jamba::printDebug("groups2contrasts(): ",
               "iChange:",
               iChange);
         }
         iNoChange <- setdiff(seq_len(ncol(iFactors)), iChange);
         ## Optionally omit certain values from consideration,
         ## notably for "," or "-" which already contain changing factors
         iFactorUseRows <- unigrep(omitGrep, iFactors[,iChange]);

         if (length(iNoChange) == 0) {
            iSplit <- rep("", length(iFactorUseRows));
         } else {
            #iSplit <- pasteByRow(iFactors[iFactorUseRows,iNoChange,drop=FALSE],
            #   sep=factorSep);
            iSplit <- jamba::pasteByRowOrdered(iFactors[iFactorUseRows,iNoChange,drop=FALSE],
               sep=factorSep);
         }

         ## Split rows by constant values in non-changing factor columns
         iSplitL <- split(iFactorUseRows, iSplit);
         iSplitL <- iSplitL[lengths(iSplitL) > 1];
         ## Only consider contrasts when there are multiple rows
         if (length(iSplitL) > 0) {
            iDF <- jamba::rbindList(lapply(iSplitL, function(iSplitRows) {
               if (verbose >= 1) {
                  jamba::printDebug("groups2contrasts(): ",
                     "   iSplitRows:",
                     iSplitRows);
               }
               iFactorsSub <- iFactors[iSplitRows,,drop=FALSE];
               iFactorVals <- iFactorsSub[,iChange];
               iMatch <- match(
                  sortSamples(iFactorVals,
                     preControlTerms=preControlTerms),
                  iFactorVals);
               iCombn <- combn(iMatch, 2);
               iGrp1 <- ifelse(grepl("-", rownames(iFactorsSub)[iCombn[2,]]),
                  paste0("(", rownames(iFactorsSub)[iCombn[2,]], ")"),
                  rownames(iFactorsSub)[iCombn[2,]]);
               iGrp2 <- ifelse(grepl("-", rownames(iFactorsSub)[iCombn[1,]]),
                  paste0("(", rownames(iFactorsSub)[iCombn[1,]], ")"),
                  rownames(iFactorsSub)[iCombn[1,]]);
               iContrastName <- paste0(iGrp1, "-", iGrp2);
               iContrastDF <- data.frame(check.names=FALSE,
                  shrinkMatrix(iFactorsSub[intercalate(iCombn[2,], iCombn[1,]),,drop=FALSE],
                     groupBy=rep(iContrastName, each=2),
                     shrinkFunc=function(x){cPasteUnique(x, doSort=FALSE)},
                     returnClass="matrix"),
                  contrastName=iContrastName, row.names=iContrastName);

               ## Create a string representing the combination of factors.
               ## which we will use to prevent re-creating the same contrasts.
               ##
               ## Modified the string to include colname, to ensure that two
               ## factors which may share some levels, will not be confused.
               iContrastDF[,"contrastString"] <- jamba::pasteByRow(
                  iContrastDF[,colnames(iFactorsSub),drop=FALSE],
                  includeNames=TRUE,
                  sep=";",
                  sepName=":");
               iContrastDF;
            }));
            rownames(iDF) <- iDF[,"contrastName"];
            iDF;
         } else {
            NULL;
         }
      })));

   ## Optionally spike in some pre-defined non-standard contrasts
   if (!is.null(addContrastNamesDF)) {
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "Adding custom ",
            "addContrastNamesDF");
      }
      iContrastNames <- rbind(iContrastNames, addContrastNamesDF);
   }

   ## Optionally make each row unique in terms of the factors compared
   if (makeUnique) {
      iDFcomponents <- jamba::pasteByRow(
         iContrastNames[,setdiff(colnames(iContrastNames), "contrastName"),drop=FALSE],
         sep="!");
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "makeUnique=",
            "TRUE");
         jamba::printDebug("iDFcomponents:\n",
            iDFcomponents, sep="\n");
         jamba::printDebug("unique(iDFcomponents):\n",
            unique(iDFcomponents), sep="\n");
      }
      iDFrowKeep <- match(unique(iDFcomponents), iDFcomponents);
      iContrastNames <- iContrastNames[iDFrowKeep,,drop=FALSE];
   }

   if ("contrastName" %in% colnames(iContrastNames)) {
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "tcount(iContrastNames$contrastName):")
         print(jamba::tcount(iContrastNames[,"contrastName"]));
      }
      #rownames(iContrastNames) <- iContrastNames[,"contrastName"];
      rownames(iContrastNames) <- jamba::makeNames(iContrastNames[,"contrastName"]);
   }
   ## Optionally remove contrasts with factor pairs not of interest
   if (!is.null(removePairs)) {
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "Processing any removePairs contrasts.");
      }
      for (iCol in setdiff(colnames(iContrastNames), "contrastName")) {
         if (verbose) {
            jamba::printDebug("groups2contrasts(): ",
               "   Checking for removePairs in column:", iCol);
         }
         if (any(iContrastNames[,iCol] %in% removePairsFull)) {
            if (verbose) {
               jamba::printDebug("      removedPair:\n",
                  head(intersect(iContrastNames[,iCol], removePairsFull), 10),
                  fgText=c("yellow", "purple"), sep="\n");
            }
            iContrastNames <- iContrastNames[which(!iContrastNames[,iCol] %in% removePairsFull),,drop=FALSE];
         }
      }
   }

   if (verbose) {
      jamba::printDebug("groups2contrasts(): ",
         "iContrastNames:");
      print(head(iContrastNames));
   }
   if (length(setdiff(colnames(iContrastNames), "contrastName")) > 1 && currentDepth < maxDepth) {
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "   Defining interactions contrasts.");
         print(head(iContrastNames[,iContrastGroupsUse,drop=FALSE]));
      }
      iContrastNamesInt <- groups2contrasts(iContrastNames[,iContrastGroupsUse,drop=FALSE],
         omitGrep=omitGrep,
         currentDepth=currentDepth+1,
         maxDepth=maxDepth,
         returnDesign=FALSE,
         factorSep=factorSep,
         factorOrder=factorOrder,
         contrastSep=contrastSep,
         renameFirstDepth=renameFirstDepth,
         removePairs=removePairs,
         makeUnique=makeUnique,
         preControlTerms=preControlTerms,
         verbose=verbose,
         ...);
      if (igrepHas("[(]", rownames(iContrastNamesInt[[1]]))) {
         return(iContrastNamesInt);
      }
      if (ncol(iContrastNamesInt) > 1 && any(is.na(iContrastNamesInt[,1]))) {
         iContrastNamesInt <- iContrastNamesInt[!is.na(iContrastNamesInt[,1]),,drop=FALSE];
      }
      if (length(iContrastNamesInt) == 0 || ncol(iContrastNamesInt) > 1) {
         if (verbose) {
            jamba::printDebug("groups2contrasts(): ",
               "begin iContrastNamesInt:");
            print(head(iContrastNamesInt));
            jamba::printDebug("  end iContrastNamesInt:");
         }
         iContrastNames <- jamba::rbindList(list(iContrastNames,
            iContrastNamesInt));
      }
   } else {
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "   Skipping interactions");
         jamba::printDebug("      ncol(iContrastNames):",
            ncol(iContrastNames));
      }
   }
   if ("contrastName" %in% colnames(iContrastNames)) {
      rownames(iContrastNames) <- jamba::makeNames(iContrastNames[,"contrastName"]);
   }
   if (returnDesign && currentDepth == 1) {
      iContrasts <- NULL;
      if (!is.null(iDesign)) {
         iContrasts <- limma::makeContrasts(contrasts=iContrastNames[,"contrastName"],
            levels=iDesign);
      }
      retvals <- list(iContrastNames=iContrastNames,
         iContrasts=iContrasts,
         iDesign=iDesign)
   } else {
      retvals <- iContrastNames;
   }
   return(retvals);
}

#' Sort biological sample labels for experimental design
#'
#' Sort biological sample labels for experimental design
#'
#' This function sorts a vector of sample labels using typical
#' heuristics that order typical control groups terms before
#' test groups. For example, `"Vehicle"` would be returned
#' before `"Treatment"` since `"Vehicle"` is a recognized control
#' term.
#'
#' It also employs `jamba::mixedSort()` for
#' proper alphanumeric sorting, for example so `"Time_5hr"` would
#' be sorted before `"Time_12hr"`.
#'
#' @return character vector ordered such that control terms are
#' preferentially first before non-control terms.
#'
#' @param x character vector or factor
#' @param controlTerms vector of regular expression patterns used to
#'    determine control terms, where the patterns are matched and
#'    returned in order.
#' @param preControlTerms vector or NULL, optional control
#'    terms or regular expressions to use before the `controlTerms`
#'    above. This argument is used as a convenient prefix to the
#'    default terms.
#' @param postControlTerms vector or NULL, optional control
#'    terms or regular expressions to use after the `controlTerms`
#'    above. This argument is used as a convenient suffix to the
#'    default terms.
#' @param keepFactorsAsIs logical indicating whether to maintain
#'    factor level order, if `x` is supplied as a factor. If
#'    `keepFactorsAsIs==TRUE` then only `sort(x)` is returned.
#' @param ... additional arguments are ignored.
#'
#' @family jam string functions
#' @family jam RNA-seq functions
#'
#' @examples
#' # the defaults perform well for clear descriptors
#' sortSamples(c("Trt_12h", "Trt_9h", "Trt_1h", "Trt_9h", "Vehicle"));
#'
#' # custom terms can be added before the usual control terms
#' sortSamples(c("Trt_12h", "Trt_9h", "Trt_1h", "Trt_9h", "Fixated", "Vehicle"),
#'    preControlTerms="fixate");
#'
#' # custom terms can be added after the usual control terms
#' sortSamples(c("Trt_12h", "Trt_9h", "Trt_1h", "Trt_9h", "Fixated", "Vehicle"),
#'    postControlTerms="fixate");
#'
#' @export
sortSamples <- function
(x,
 controlTerms=c("WT", "(^|[-_ ])NT($|[-_ ]|[0-9])",
    "ETOH",
    "control|ctrl|ctl",
    "ntc",
   "Vehicle|veh",
    "none|empty|blank",
    "scramble",
    "ttx",
    "PBS"),
 sortFunc=jamba::mixedSort,
 preControlTerms=NULL,
 postControlTerms=NULL,
 keepFactorsAsIs=TRUE,
 ...)
{
   ## Purpose is to order sample names by typical descriptions
   ## of control groups versus treatment groups
   ##
   ## Test set:
   ## sortSamples(c("Trt_12h", "Trt_9h", "Trt_1h", "Vehicle"))
   ## sortSamples(c("RA_Brg1", "EtOH_WT", "RA_WT", "EtOH_Brg1"))
   ## sortSamples(c("HCTWT_DXR6", "HCTWT_DXR12", "HCTWT_DXR24", "HCTWT_NT24"))
   #order1 <- proigrep(c(controlTerms), x);
   #order2 <- proigrep(c(controlTerms, "."), sortFunc=sortFunc, x);
   ##
   ## keepFactorsAsIs=TRUE will keep factor levels unchanged, and use those levels in the sort
   ## instead of looking for control terms
   if (keepFactorsAsIs && jamba::igrepHas("factor", class(x))) {
      sort(x);
   } else {
      controlTerms <- c(preControlTerms,
         controlTerms,
         postControlTerms);
      xU <- jamba::provigrep(c(controlTerms, "."),
         sortFunc=sortFunc,
         x);
      xOrder <- order(match(x, xU));
      x[xOrder];
   }
}

#' Split the elements of an ordered factor vector
#'
#' Split the elements of an ordered factor vector
#'
#' This function performs `base::strsplit()` while trying to maintain
#' the order of factor levels in the output, based upon the order of
#' factor levels in the input data.
#'
#' @return list of factor vectors, where each factor shares the same
#'    global factor levels based upon the input data.
#'
#' @param x character or factor vector.
#' @param split character split value sent to `base::strsplit()`.
#' @param fixed,perl,useBytes additional arguments sent to `base::split()`.
#' @param sortFunc function used to sort character values when the input
#'    `x` is a character vector. The default `jamba::mixedSort()` applies
#'    alphanumeric sort.
#' @param keepOrder logical indicating whether to keep the order of values
#'    in the input data, for example with character input the values will
#'    be ordered by the first appearance of each term.
#' @param ... additional arguments are ignored.
#'
#' @family jam string functions
#'
#' @examples
#' # first define a vector of sample groups
#' iGroups <- jamba::nameVector(paste(rep(c("WT", "KO"), each=6),
#'    rep(c("Control", "Treated"), each=3),
#'    sep="_"));
#' iGroups <- factor(iGroups, levels=unique(iGroups));
#' iGroups;
#' strsplitOrdered(iGroups, "_");
#'
#' @export
strsplitOrdered <- function
(x,
 split="_",
 fixed=FALSE,
 perl=FALSE,
 useBytes=FALSE,
 sortFunc=jamba::mixedSort,
 keepOrder=TRUE,
 ...)
{
   ## Purpose is to run strsplit() on factors, ordering the new factor
   ## levels consistent with the input
   if (!jamba::igrepHas("factor", class(x))) {
      if (keepOrder) {
         x <- factor(x,
            levels=unique(x));
      } else {
         x <- factor(x,
            levels=sortFunc(unique(x)));
      }
   }
   soL <- strsplit(x=levels(x),
      split=split,
      fixed=fixed,
      perl=perl,
      useBytes=useBytes);
   so1 <- jamba::rbindList(soL);

   ## Note: the setdiff() is there to remove "" values
   so1levels <- setdiff(unique(unlist(apply(so1, 2, unique))), "");
   soSplitL <- strsplit(as.character(x),
      split=split,
      fixed=fixed,
      perl=perl,
      useBytes=useBytes);
   soLordered <- lapply(soSplitL, function(i){
      factor(i,
         levels=so1levels);
   });
   return(soLordered);
}
