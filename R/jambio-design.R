
#' Define experimental contrasts from sample groups
#'
#' Define experimental contrasts from sample groups
#'
#' This function is intended to define statistical contrasts
#' that compare one factor at a time. For two-factor designs,
#' it will create two-way contrasts, defined as the contrast
#' of pairwise contrasts.
#'
#' Input can be a character vector of group names, where by
#' default each factor is separated by an underscore `"_"`.
#' An example might be:
#'
#' `iFactors <- c("Control_Wildtype", "Control_Knockout",
#' "Treated_Wildtype", "Treated_Knockout")`
#'
#' In that case, there are two factors. The first factor
#' contains factor levels `c("Control", "Treated")`, and
#' the second factor contains factor levels
#' `c("Wildtype", "Knockout")`.
#'
#' Input can also be a `data.frame` (or compatible table-like
#' object including `data.table` and `tibble`). Each column
#' is considered a factor. From the example above, we can
#' create a `data.frame` using `jamba::rbindList()`,
#' see the Examples for more detail.
#'
#' `jamba::rbindList(strsplit(iFactors, "_"))`
#'
#' Lastly, if the input is a named vector, or a `data.frame`
#' with rownames,
#'
#'
#' This function will change any `"-"` in a factor name to
#' `"."` prior to detecting valid contrasts. Note that
#' `groups2contrasts()` does not call `base::make.names()`
#' because that function too aggressively converts characters
#' to `"."`. If data must be compliant with the rules used
#' by `base::make.names()`, run that function prior to calling
#' `groups2contrasts()`.
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

   if (jamba::igrepHas("factor|character", class(iFactors))) {
      #####################################################
      ## Vector input
      ##
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "splitting vector into groups");
      }
      if (length(names(iFactors)) == 0) {
         if (length(iSamples) == 0) {
            ## Create iSamples
            iSamples <- jamba::makeNames(rep("sample", length(iFactors)));
            names(iFactors) <- iSamples;
         } else if (length(iSamples) != length(iFactors)) {
            stop(paste0("length(iSamples) must be equal length(iFactors) ",
               "when there are no names(iFactors)."));
         }
         names(iFactors) <- iSamples;
      } else if (length(iSamples) == 0) {
         iSamples <- names(iFactors);
      } else {
         if (!any(iSamples %in% names(iFactors)) && length(iSamples) == length(iFactors)) {
            ## Use iSamples as-is
            names(iFactors) <- iSamples;
         } else if (!all(iSamples %in% names(iFactors))) {
            stop(paste0("iSamples is present in some not not all names(iFactors). ",
               "iSamples must either: all be present in names(iFactors); or ",
               "present in none of names(iFactors) and length(iSamples) == length(iFactors)."))
         } else {
            ## Re-order iFactors to match iSamples
            iFactors <- iFactors[match(iSamples, names(iFactors))];
         }
      }
      if (igrepHas("factor", class(iFactors))) {
         ## Convert factor to a data.frame where each column
         ## is a factor with ordered levels that match the order
         ## the factor levels appear in the original factor.
         iFactorsL <- strsplitOrdered(iFactors, factorSep);
         names(iFactorsL) <- names(iFactors);
         iFactorsLevels <- levels(iFactorsL[[1]]);
         iFactors <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            jamba::rbindList(
               strsplit(as.character(iFactors),
                  factorSep)));
         rownames(iFactors) <- names(iFactorsL);
         for (i in seq_len(ncol(iFactors))) {
            iFactors[,i] <- factor(iFactors[,i],
               levels=intersect(iFactorsLevels, iFactors[,i]));
         }
      } else {
         ## Convert to data.frame
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
      } else {
         colnames(iFactors) <- jamba::makeNames(
            rep("factor",
               length.out=ncol(iFactors)),
            renameOnes=TRUE);
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
   } else if (jamba::igrepHas("data.frame|dataframe", class(iFactors))) {
      #####################################################
      ## data.frame input
      ##
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "using existing data.frame");
      }
      if (length(rownames(iFactors)) == 0) {
         if (length(iSamples) == 0) {
            ## Create iSamples
            iSamples <- jamba::makeNames(rep("sample", nrow(iFactors)));
         } else if (length(iSamples) == nrow(iFactors)) {
            # use iSamples as-is
         } else {
            stop(paste0("iFactors has no rownames, and ",
               "length(iSamples) != nrow(iFactors). ",
               "Please make length(iSamples) == nrow(iFactor)"));
         }
      } else {
         if (length(iSamples) == 0) {
            iSamples <- rownames(iFactors);
         } else {
            if (!any(iSamples %in% iFactors) && length(iSamples) == nrow(iFactors)) {
               ## use iSamples as-is
            } else if (!all(iSamples %in% rownames(iFactors))) {
               stop(paste0("iSamples is not present in all rownames(iFactors). ",
                  "Either: all iSamples must be present in rownames(iFactors); or ",
                  "no iSamples are present in rownames(iFactors) and ",
                  "length(iSamples) == nrow(iFactors)."));
            } else {
               ## Subset or re-order iFactors using matching iSamples
               iFactors <- iFactors[match(iSamples, rownames(iFactors)),,drop=FALSE];
               if (verbose) {
                  jamba::printDebug("groups2contrasts(): ",
                     "Specifying iFactors[iSamples,]");
                  print(head(iFactors));
               }
            }
         }
         if (verbose) {
            printDebug("groups2contrasts(): ",
               "head(iFactors):");
            print(head(iFactors, 100));
         }
      }
      if (length(groupColumns) == 0) {
         if (length(colnames(iFactors)) == 0) {
            ## Create colnames
            groupColumns <- jamba::makeNames(
               renameOnes=TRUE,
               rep("factor",
                  length.out=ncol(iFactors)));
            colnames(iFactors) <- groupColumns;
         } else {
            groupColumns <- colnames(iFactors);
         }
      } else {
         if (!all(groupColumns %in% colnames(iFactors))) {
            stop(paste0("Not all groupColumns are in colnames(iFactors), please remedy."));
         }
         ## Use iFactors as-is
         #iFactors <- iFactors[,groupColumns,drop=FALSE];
      }
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "Specifying iFactors[,groupColumns,drop=FALSE]");
         jamba::printDebug("groups2contrasts(): ",
            "groupColumns:",
            groupColumns);
      }

      ## Use mixedSortDF() to sort the data.frame,
      ## which will honor factor level orders if present.
      ## To influence this sorting, use factors with ordered levels
      ## instead of character columns.
      iFactors <- jamba::mixedSortDF(iFactors,
         byCols=groupColumns);
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "iFactors:");
         print(head(iFactors));
      }

      ## rowGroups is the unique set of group names, used to keep the original order
      #rowGroups <- jamba::pasteByRowOrdered(iFactors[,groupColumns,drop=FALSE],
      #   sep=factorSep);
      ## Unclear whether to re-order columns to match groupColumns, for now we do not
      rowGroups <- jamba::pasteByRowOrdered(iFactors,
         sep=factorSep);
      if (length(rownames(iFactors)) == 0) {
         iFactors_names <- jamba::makeNames(rowGroups,
            suffix="_rep");
         rownames(iFactors) <- iFactors_names;
      } else {
         iFactors_names <- rownames(iFactors);
      }
      ## Assume for now sample rows and group columns
      sample2group <- split(iFactors_names, rowGroups);
      if (length(iDesign) == 0) {
         iDesign <- list2im(sample2group)[iFactors_names,as.character(unique(rowGroups)),drop=FALSE];
         if (all(iSamples %in% iFactors_names)) {
            iDesign <- iDesign[match(iSamples, iFactors_names),,drop=FALSE];
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
         colnames(iFactors) <- jamba::makeNames(
            rep("groupFactor",
               length.out=ncol(iFactors)),
            renameOnes=TRUE,
            suffix="_");
      }
      rownames(iFactors) <- unname(jamba::pasteByRow(iFactors, sep=factorSep));
      printDebug("iFactors:");print(iFactors);
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
   ## Check to make sure no factor levels contain "-"
   for (i in colnames(iFactors)) {
      if (igrepHas("-", iFactors[,i])) {
         iFactors[,i] <- gsub("-", ".", iFactors[,i]);
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
         jamba::printDebug("groups2contrasts(): ",
            "iDFcomponents:\n",
            iDFcomponents, sep="\n");
         jamba::printDebug("groups2contrasts(): ",
            "unique(iDFcomponents):\n",
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
      if (verbose) {
         jamba::printDebug("groups2contrasts(): ",
            "length(iContrastNamesInt):",
            length(iContrastNamesInt));
      }
      ## If length==0 then there are no valid interaction contrasts
      if (length(iContrastNamesInt) > 0 &&
         igrepHas("[(]", rownames(iContrastNamesInt[[1]]))) {
         return(iContrastNamesInt);
      }
      if (length(iContrastNamesInt) > 0 &&
         ncol(iContrastNamesInt) > 1 &&
         any(is.na(iContrastNamesInt[,1]))) {
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
         jamba::printDebug("      head(iContrastNames):");
         print(head(iContrastNames));
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
#' @param ignore.case logical passed to `jamba::provigrep()` indicating
#'    whether to ignore case-sensitive matching.
#' @param boundary logical indicating whether to require a word
#'    boundary at either the start or end of the control terms.
#'    When TRUE, it uses `perl=TRUE` by default, and allows either
#'    perl boundary or an underscore `"_"`.
#' @param perl logical indicating whether to use Perl regular
#'    expression pattern matching.
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
 controlTerms=c("WT|wildtype",
    "(^|[-_ ])(NT|NTC)($|[-_ ]|[0-9])",
    "ETOH",
    "control|ctrl|ctl",
    "Vehicle|veh",
    "none|empty|blank",
    "scramble",
    "ttx",
    "PBS",
    "knockout",
    "mutant"),
 sortFunc=jamba::mixedSort,
 preControlTerms=NULL,
 postControlTerms=NULL,
 ignore.case=TRUE,
 boundary=TRUE,
 perl=boundary,
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
      controlTerms <- unique(c(preControlTerms,
         controlTerms,
         postControlTerms));
      if (any(boundary)) {
         # Require regular expression boundary
         controlTerms1 <- unlist(lapply(controlTerms, function(i){
            paste0("(_|\\b)(", i, ")|(", i, ")(_|\\b)")
         }))
         if (any(!boundary)) {
            controlTerms <- c(controlTerms1,
               controlTerms);
         } else {
            controlTerms <- controlTerms1;
         }
      }
      xU <- jamba::provigrep(c(controlTerms, "."),
         sortFunc=sortFunc,
         perl=perl,
         ignore.case=ignore.case,
         x);
      xOrder <- order(match(x, xU));
      x <- x[xOrder];
      #attr(x, "controlTerms") <- controlTerms;
      x;
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

#' Curate vector into a data.frame
#'
#' Curate vector into a data.frame
#'
#' This function is intended to curate a vector into a data.frame
#' with specifically assigned colnames. It is intended to be a more
#' generic method of curation annotations than splitting a characteer
#' string by some delimiter, for example where the order of annotations
#' may differ entry to entry, but where there are known patterns
#' that are sufficient to describe an annotation column.
#'
#' That said, if annotations can be reliably split using a delimiter,
#' that method is often a better choice. In that case, this function
#' may be useful to make input data fit the expected format.
#'
#' For example from `c("Sample1_WT_LPS_1hour", "Sample2_KO_LPS_2hours")`
#' we can tell whether a sample is `KO` or `WT` by looking for that
#' substring.
#'
#' The `curationL` is a list with the following properties:
#'
#' * `names(curationL)` represent colnames to create in the output
#' data.frame.
#' * each list element contains a list of two-element vectors
#' * each two-element vector contains a substitution pattern and
#' substitution replacement
#'
#' When `matchWholeString=TRUE` the substitution patterns are extended
#' to match the whole string, using parentheses around the main pattern.
#' For example if the pattern is "KO" and replacement is "KO", then the
#' pattern is extended to "^.*KO.*$", so the entire string will be
#' replaced with "KO".
#'
#' Typically, `curationL` is derived from YAML formatted files, and
#' loaded into a list with this type of setup:
#'
#' `curationL <- yaml::yaml.load_file("curation.yaml")`.
#'
#' The generic YAML format is as follows:
#'
#' NewColname_1:
#' - - patternA
#'   - replacementA
#' - - patternB
#'   - replacementB
#' NewColname_2:
#' - - patternC
#'   - replacementC
#'
#' A specific example:
#'
#' Treatment:
#' - - LPS
#'   - LPS
#' - - Control|cntrl|ctrl
#'   - Control
#' Genotype:
#' - - WT|wildtype
#'   - WT
#' - - KO|knockout|knock
#'   - KO
#'
#' @param x character vector as input
#' @param curationL list containing curation rules, as described above, or
#'    a character vector of yaml files, which will be imported into
#'    a list format using `yaml::yaml.load_file()`.
#' @param matchWholeString logical indicating whether to match the whole
#'    string for each entry in `x`. If `matchWholeString=TRUE` then
#'    the substitution patterns are all extended where needed, in order
#'    to expand the pattern to match the whole string.
#' @param trimWhitespace logical indicating whether to trim leading and
#'    trailing whitespace characters from `x`.
#' @param whitespace character vector containing whitespace characters.
#' @param expandWhitespace logical indicating whether substitution patterns
#'    should be modified so any whitespace characters in the pattern will
#'    match the defined `whitespace` characters. For example when
#'    `expandWhitespace=TRUE`, the pattern `"_KO_"` will be modified to
#'    `"[ _]+KO[ _]+"` so the pattern will match `" KO "` and `"_KO_"`.
#' @param previous optional data.frame whose colnames may be present as
#'    names in `curationL`, or single vector with `length(previous)=length(x)`.
#'    If `previous` is supplied as a data.frame, and the curation colname
#'    is present in `colnames(previous)`, then unmatched substutition
#'    patterns will retain the data in the relevant column of `previous`.
#'    This mechanism allows editing single values in an existing column,
#'    based upon pattern matching in another column.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @family jam design functions
#'
#' @examples
#' set.seed(123);
#' x <- paste(
#'    paste0("file",
#'       sapply(1:5, function(i) {
#'          paste(sample(LETTERS, 5), collapse="")
#'       })),
#'    rep(c("WT", "Mut"), each=3),
#'    rep(c("Veh","EtOH"), 3),
#'    sep="_");
#' x;
#'
#' curationYaml <- c(
#' "Genotype:
#' - - WT|wildtype
#'   - WT
#' - - Mut|mutant
#'   - Mut
#' Treatment:
#' - - Veh|EtOH
#'   - \\1
#' File:
#' - - file([A-Z]+)
#'   - \\1
#' FileStem:
#' - - file([A-Z]+)
#'   - \\2");
#' # print the curation.yaml to show its structure
#' cat(curationYaml)
#' curationL <- yaml::yaml.load(curationYaml);
#' curateVtoDF(x, curationL);
#'
#' @export
curateVtoDF <- function
(x,
 curationL=NULL,
 matchWholeString=TRUE,
 trimWhitespace=TRUE,
 whitespace="_ ",
 expandWhitespace=TRUE,
 previous=NULL,
 verbose=TRUE,
 ...)
{
   ## Purpose is to take one vector x and apply a list of curation
   ## rules to the values, thereby creating a data.frame whose columns
   ## match the names of the curationL list.
   ##
   ## matchWholeString=TRUE will apply the regular expression pattern
   ## match to the whole string, even for patterns that do not
   ## specify leading or trailing match '^' and '$' respectively.
   ##
   ## trimWhitespace=TRUE will trim leading and trailing whitespace
   ## characters, including underscore '_' and dash '-'.
   ##
   ## expandWhitespace=TRUE will expand pattern matching to include
   ## any characters in whitespace, for example to match underscore '_'
   ## and space ' ' in the curation, which is helpful when file names
   ## or column headers may have been edited to replace spaces with
   ## non-whitespace characters.
   ##
   ## previous is a data.frame (or named list) with nrow=length(x)
   ## that contains values to be used when the grep patterns do not
   ## match any entries in x. It is mainly intended to be used by
   ## curateDFtoDF() for more complex column curation.
   ##
   ## Note that x is used to create rownames using jamba::makeNames(x) which
   ## ensures that rownames are unique. In this case, the resulting
   ## rownames will not match the input vector x.
   ##

   ## Make sure input whitespace does not already have square brackets
   whitespace <- gsub("^[[]|[]][+]*$", "", whitespace);
   ## Generate a pattern to remove leading and trailing whitespace
   trimRegexp <- paste0("^[", whitespace, "]+|[",
      whitespace, "]+$");

   ## First check if input is a list, or data.frame, and convert as needed
   if (igrepHas("data.frame", class(curationL))) {
      curationL <- curationDFtoL(curationL,
         whitespace=whitespace);
   }
   if (igrepHas("character", class(curationL)) &&
         all(file.exists(curationL))) {
      curationL <- do.call(c, lapply(curationL, yaml::yaml.load_file));
   }
   if (!igrepHas("list", class(curationL))) {
      stop("curateVtoDF() requires curationL in list format.");
   }

   expandGrep <- function
   (iGrep1)
   {
      ## Purpose is to expand grep pattern to include the whole string
      ## Todo: adjust iGrep2 replacement to increment the numeric
      ## references as needed.
      if (igrepHas("[^[][\\^]|^\\^", iGrep1)) {
         ## Do not add the leading ^
         if (igrepHas("[$]", iGrep1)) {
            ## Do not add trailing $
         } else {
            iGrep1 <- paste0("(", iGrep1, ").*$");
         }
      } else {
         if (igrepHas("[$]", iGrep1)) {
            ## Do not add trailing $
            iGrep1 <- paste0("^.*(", iGrep1, ")");
         } else {
            iGrep1 <- paste0("^.*(", iGrep1, ").*$");
         }
      }
      iGrep1;
   }
   ## Assignment here to force evaluation in this environment
   previous <- previous;

   curationValuesL <- lapply(nameVectorN(curationL), function(iName){
      if (verbose) {
         printDebug("curateVtoDF(): ",
            "Creating column:",
            iName);
      }
      iGrepL <- curationL[[iName]];
      x1 <- as.character(x);
      ## Todo: make it only replace matching entries, otherwise
      ## substitute NA for non-matched patterns.
      x1which <- integer(0);
      ## Iterate each from,to substitution for this column
      for (iGrepLx in iGrepL) {
         iGrep1 <- iGrepLx[[1]];
         if (expandWhitespace) {
            iGrep1 <- escapeWhitespaceRegexp(iGrep1,
               whitespace=whitespace);
         }
         iGrep2 <- iGrepLx[[2]];
         if (matchWholeString) {
            ## Todo: adjust iGrep2 so the numeric references are
            ## also changed where needed
            iGrep1 <- expandGrep(iGrep1);
         }
         ## Non-destructive replacement only using values
         ## matching the pattern
         x1which1 <- grep(iGrep1, x1);
         x1which <- unique(c(x1which, x1which1));
         if (length(x1which) > 0) {
            x1[x1which] <- gsub(iGrep1,
               iGrep2,
               x1[x1which]);
         }
      }
      if (trimWhitespace && igrepHas(trimRegexp, x1)) {
         if (verbose) {
            printDebug("curateVtoDF(): ",
               "   Trimming leading/trailing whitespace.");
         }
         x1 <- gsub(trimRegexp, "", x1);
      }
      ## Check for previous data, only when the grep pattern
      ## has not matched every entry in x
      if (length(previous) > 0 && length(x1which) < length(x)) {
         if (is.vector(previous)) {
            previousV <- previous;
            previousV[x1which] <- x1[x1which];
            ## Update the previous data
            ## Must check to see if the scope inside lapply() works as expected
            previous <- previousV;
            assign("previous", previous, pos=-1);
         } else if (iName %in% colnames(previous)) {
            previousV <- previous[[iName]];
            previousV[x1which] <- x1[x1which];
            ## Update the previous data.frame
            ## Must check to see if the scope inside lapply() works as expected
            previous[[iName]] <- previousV;
            assign("previous", previous, pos=-1);
         } else {
            ## Cannot do the thing
            previousV <- x1;
         }
         x1 <- previousV;
         if (verbose) {
            printDebug("curateVtoDF(): ",
               "   Curated a subset of values, and used previous data otherwise.");
         }
      }
      x1;
   });
   iDF <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      do.call(cbind, curationValuesL));
   rownames(iDF) <- jamba::makeNames(x);
   iDF;
}

#' Curate data.frame into a data.frame
#'
#' Curate data.frame into a data.frame
#'
#' This function takes a data.frame as input, where one or more
#' columns are expected to be used in data curation to create
#' another data.frame. This situation is useful when the final
#' desired data.frame depends upon values in more than one
#' column of the input data.frame.
#'
#' Specifically, this function is a wrapper around `curateVtoDF()`.
#'
#' Typically, `curationL2` is derived from YAML formatted files, and
#' loaded into a list with this type of setup:
#'
#' `curationL2 <- yaml::yaml.load_file("curation.yaml")`.
#'
#' The structure of curationL2:
#'
#' * `curationL2` is a list object, whose `names(curationL2)` are values
#' in `colnames(x)` and represent column of data used as input.
#' * each list element in `curationL2` is also a list, whose
#' `names` represent colnames to create or update in the output
#' `data.frame`.
#' * these lists contain character vectors `length=2` containing
#' a regular expression substitution pattern (see `base::gsub`),
#' and a replacement pattern.
#'
#' The list is processed in order, and names can be repeated as
#' necessary to apply the proper substitution patterns in the
#' order required. New columns created during the curation may also
#' be used in later curation steps.
#'
#' Example curation.yaml YAML format:
#'
#' From_ColnameA:
#'   To_ColnameC:
#'   - - patternA
#'     - replacementA
#'   - - patternB
#'     - replacementB
#'   To_ColnameD:
#'   - - patternC
#'     - replacementC
#'   - - patternD
#'     - replacementD
#' From_ColnameB:
#'   To_ColnameE:
#'   - - patternE
#'     - replacementE
#'   - - patternF
#'     - replacementF
#'
#' When the rule creates a colname already present in colnames(x),
#' then only values specifically matched by the substitution patterns
#' are modified. For example, this technique can be used to modify
#' the group assignment of a Sample_ID:
#'
#' Sample_ID:
#'   Group:
#'   - - Sample1234
#'     - WildType
#'
#' The rules above will match `"Sample1234"` in the `"Sample_ID"` column
#' of x, and assign `"WildType"` to the `"Group"` column only for
#' matching entries.
#'
#' In addition to values in `colnames(x)`, the "from" value may
#' also be `"rownames"` which will cause the curation rules to
#' act upon values in `rownames(x)` instead of values in a specific
#' column of `x`.
#'
#' Note that if a "to" column does not already exist, then all values
#' in the "from" column which do not match any substitution pattern
#' will be used to fill the remainder of the "to" column.
#' Once the "to" column exists, then only entries with a matching
#' substitution pattern are replaced using the replacement pattern.
#'
#' For example, for NanoString data, the column `"CartridgeWell"` can be
#' derived from `rownames(x)`, after which the new column `"CartridgeWell"`
#' can be used in subsequent curation steps.
#'
#' Additional notes:
#'
#' * The substitution pattern is automatically expanded to include the
#' whole input string, if not already present. For example supplying `"WT"`
#' will match `"^.*(WT).*$"`. However if the substitution pattern is
#' `"^.*(WT).*$"` then it will not be expanded.
#' * When the substitution pattern is expanded, the string is also enclosed
#' in parentheses `"()"` which means the replacement can use `"\\1"` to
#' use the successfully matched pattern as the output string. For example
#' if `"WT"` and `"Mutant"` are always valid genotypes, then it would
#' be sufficient to define substitution pattern `"WT|Mutant"` and
#' replacement pattern `"\\1"`.
#' * When the substitution pattern is expanded, and the string is enclosed
#' in parentheses, any parentheses in the substitution pattern are therefore
#' one level deeper, for example `"file([A-Z]+)"` will be expanded to
#' `"^.*(file([A-Z]+)).*$"`. See the example below, where the replacement
#' pattern uses `"\\2"` to use only the internal parentheses.
#'
#'
#' @param x data.frame
#' @param curationL2 list with curation rules as described above, or
#'    a character vector of yaml files, which will be imported into
#'    a list format using `yaml::yaml.load_file()`.
#' @param matchWholeString,trimWhitespace,whitespace,expandWhitespace
#'    arguments passed to `curateVtoDF()`.
#' @param keepAllColnames logical indicating whether to keep all colnames
#'    from `x` in addition to those created during curation.
#'    `keepAllColnames=FALSE` will only keep colnames specifically
#'    described in the `curationL2` list, while `keepAllColnames=TRUE`
#'    will keep all original colnames, and any colnames added during
#'    the curation steps.
#' @param verbose logical indicating whether to print verbose output
#' @param ... additional arguments are passed to `curateVtoDF()`
#'
#' @family jam design functions
#'
#' @examples
#' set.seed(123);
#' df <- data.frame(filename=paste(
#'    paste0("file",
#'       sapply(1:5, function(i) {
#'          paste(sample(LETTERS, 5), collapse="")
#'       })),
#'    rep(c("WT", "Mut"), each=3),
#'    rep(c("Veh","EtOH"), 3),
#'    sep="_"));
#' df;
#'
#' # Note a couple ways of accomplishing similar results:
#' # Genotype matches "WT|wildtype" and replaces with "WT",
#' # then matches "Mut|mutant" and replaces with "Mut"
#' #
#' # Treatment matches "Veh|EtOH" and simply replaces with
#' # whatever was matched
#' curationYaml <- c(
#' "filename:
#'   Genotype:
#'   - - WT|wildtype
#'     - WT
#'   - - Mut|mutant
#'     - Mut
#'   Treatment:
#'   - - Veh|EtOH
#'     - \\1
#'   File:
#'   - - file([A-Z]+)
#'     - \\1
#'   FileStem:
#'   - - file([A-Z]+)
#'     - \\2");
#' # print the curation.yaml to show its structure
#' cat(curationYaml)
#' curationL <- yaml::yaml.load(curationYaml);
#' curateDFtoDF(df, curationL);
#'
#' @export
curateDFtoDF <- function
(x,
 curationL2=NULL,
 matchWholeString=TRUE,
 trimWhitespace=TRUE,
 whitespace="_ ",
 expandWhitespace=TRUE,
 keepAllColnames=TRUE,
 verbose=TRUE,
 ...)
{
   ## Purpose is to provide a wrapper around curateVtoDF() when
   ## the input data is a data.frame.
   ##
   ## In this case, curationL is a list of curationL lists, named
   ## by the colname in x to use. Each column name is used in order,
   ## to subsequent call to curateVtoDF()
   ##
   ## One special case, "rownames" is allowed as a list name, which
   ## refers to the rownames(x) and not a formal column of x.
   ##
   if (!igrepHas("data.*frame|matrix|tibble|tbl", class(x))) {
      stop("curateDFtoDF() requires x as a data.frame or compatible object.");
   }
   ## Allow for curationL2 to be a vector of yaml files
   if (igrepHas("character", class(curationL2)) &&
         all(file.exists(curationL2))) {
      curationL2 <- do.call(c, lapply(curationL2, yaml::yaml.load_file));
   }

   for (i in names(curationL2)) {
      if (verbose) {
         jamba::printDebug("curateDFtoDF(): ",
            "Applying curation to column:",
            i);
      }
      if (i %in% "rownames") {
         iV <- rownames(x);
      } else if (i %in% colnames(x)) {
         iV <- x[[i]];
      } else {
         if (verbose) {
            jamba::printDebug("curateDFtoDF(): ",
               "Skipping rules for:",
               i);
         }
         next;
      }
      if (any(names(curationL2[[i]]) %in% colnames(x))) {
         ## If colnames already exist, send it as previous data
         ## to be kept until curated
         keepNames <- intersect(names(curationL2[[i]]),
            colnames(x));
         previous <- x[,keepNames,drop=FALSE];
      } else {
         previous <- NULL;
      }
      iDF <- curateVtoDF(x=iV,
         curationL=curationL2[[i]],
         matchWholeString=matchWholeString,
         trimWhitespace=trimWhitespace,
         whitespace=whitespace,
         expandWhitespace=expandWhitespace,
         previous=previous,
         verbose=verbose,
         ...);
      ## Add data.frame to the input data, which allows the new
      ## data to be available for use by subsequent rounds of curation
      x[,colnames(iDF)] <- iDF;
   }
   ## If not returning all colnames, return only those specifically curated
   if (!keepAllColnames) {
      curatedNames <- intersect(colnames(x),
         unique(unlist(lapply(curationL2, names))));
      x <- x[,curatedNames,drop=FALSE];
   } else {
      curatedNames2 <- unique(unlist(lapply(curationL2, names)));
      curatedNames <- c(intersect(colnames(x), curatedNames2),
         setdiff(colnames(x), curatedNames2));
      x <- x[,curatedNames,drop=FALSE];
   }
   return(x);
}

#' Escape whitespace in regular expression patterns
#'
#' Escape whitespace in regular expression patterns
#'
#' This function is intended to test for unescaped whitespace characters
#' in a regular expression pattern match string, and replace them with
#' escaped whitespace characters, possibly expanding the allowed
#' whitespace characters in the process.
#'
#' @param x character vector containing a regular expression pattern,
#'    or a character value that should be converted to a regular expression
#'    pattern.
#' @param whitespace character vector `length=1` containing whitespace
#'    characters, for example `" _"` will define space `" "` and
#'    underscore `"_"` both as whitespace characters.
#' @param maxN integer value for maximum iterations of the substitution,
#'    unfortunately the regular expression logic only matches once per
#'    iteration per string.
#' @param ... additional arguments are ignored.
#'
#' @family jam string functions
#'
#' @examples
#' x <- c("one two three", "one[ ]two three", "one[ 12] two three");
#' escapeWhitespaceRegexp(x);
#' # side-by-side summary of input and output
#' data.frame(input=paste0('"', x, '"'),
#'    output=paste0('"', escapeWhitespaceRegexp(x), '"'))
#'
#' @export
escapeWhitespaceRegexp <- function
(x,
 whitespace="_ ",
 maxN=20,
 ...)
{
   ## Purpose is to test for unescaped whitespace in a regular expression
   ## pattern, and replace any unescaped whitespace with a suitable
   ## regular expression pattern.
   ##
   ## It checks for whitespace which does not follow an opening [ bracket,
   ## and does not check for closing ] bracket since the absence of ]
   ## closing bracket would already cause an error with gsub().
   ##
   ## maxN is the maximum iterations when matching patterns, since it
   ## can only typically match one pattern per iteration as written.
   ## Once an iteration does not result in a change to x then iterations
   ## are stopped.
   ##
   ## Allowed conditions:
   ## start-of-line or closing bracket ],
   ## followed by no opening bracket [,
   ## followed by whitespace
   ##
   ## Test with x <- c("one two three", "one[ ]two three", "one[ 12] two three");
   ##
   ## Make sure input whitespace does not already have square brackets
   whitespace <- gsub("^[[]|[]][+]*$",
      "",
      whitespace);

   ## Define the patterns allowed
   whitespaceN3 <- paste0("((^|[]])[^[",
      whitespace,
      "]*)([",
      whitespace,
      "]+)");
   ## Define the substitution
   whitespaceT3 <- paste0("\\1[",
      whitespace,
      "]+");

   ## Iterate multiple times as needed, in order to escape each
   ## occurrence of whitespace
   for (i in 1:maxN) {
      iWhich <- igrep(whitespaceN3, x);
      xNew <- gsub(whitespaceN3, whitespaceT3, x[iWhich]);
      if (all(xNew == x[iWhich])) {
         break;
      }
      x[iWhich] <- xNew;
   }
   x;
}
