
#' Prepare Sashimi plot
#'
#' Prepare Sashimi plot
#'
#' @export
sashimiPrep <- function
(inputGenes,
exonsGR,
expandUpDown=c(0,0),
WHICHEXONS=NULL,
YMAX=NULL,
compressGaps=TRUE,
GAPWIDTH=100,
readCountCex=1.5,
readCountFont=1,
readCountSrt=0,
maxLwd=5,
readColor="#BB000099",
readCountPadding="         ",
bigWigFiles=NULL,
covL=NULL,
spliceCountsDF=NULL,
colorFunc2=function(n,...){head(rainbowJam(max(c(n, 5))), n)},
cex.main=1.7,
font.main=2,
alternateSpliceUpDown=FALSE,
srtSplice=0,
arcScale=20,
arcScaleFactor=100,
useAspect=2,
transformY=c("none", "sqrt", "log2"),
unionCoverages=FALSE,
unionLabel="unionCoverage",
ylimFactor=1.5,
spliceGRs=NULL,
spliceBuffer=3,
exactJunctionsOnly=TRUE,
absValues=TRUE,
countAsPercent=FALSE,
geneSymbolColname=c("geneSymbol","geneSymbolChrom","gene_id","gene_name"),
flipStrand=TRUE,
ignoreSpliceStrand=TRUE,
spliceOverlapType=c("within","any","start","end"),
reverseStrandCoords=TRUE,
showCoords=TRUE,
labelExons=FALSE,
labelExonName="exonLabel",
doPar=TRUE,
parLend="round",
doTitle=TRUE,
doMargins=TRUE,
doPlot=TRUE,
minReads=0,
VERBOSE=TRUE,
exonLabelColname="gene_nameExon",
verbose=VERBOSE,
...)
{
   ## Purpose is to mimic the sashimi plot from the MISO package,
   ## using bigWig coverage data and GRanges exons as input.
   ##
   ## unionCoverages=TRUE will load coverages from each file, then create one
   ## union coverage using their sum.
   ##
   ## transformY will transform the y-axis, either to normal scale, log2-scale,
   ## or square root scaled
   ##
   ## expandUpDown is an integer vector extended to length=2, which
   ## will expand the exon region upstream and downstream by the
   ## specified widths, respectively.
   ##
   ## TODO: make the arcs have higher offset based upon the visible distance
   ## on the plot.  I.e. if one exon is huge, and an arc skips that exon,
   ## make that arc have higher offset so it doesn't appear to be flat.
   ##
   ## doPlot=TRUE will draw the sashimiPlot
   ## doPlot=FALSE will not draw the sashimiPlot, but will still return
   ## data necessary to create sashimiPlots outside this function, e.g.
   ## to test vectorized drawing.
   ##
   ## ignoreSpliceStrand==TRUE will search for, and display, splice junctions
   ## from either strand. This option was necessitated by STAR junction output
   ## which sometimes incorrectly assigns to the wrong strand.
   ##
   ## spliceOverlapType is used as the type for subsetByOverlaps() when
   ## choosing which splice junctions to display. "within" will ensure all
   ## junctions are within the viewing range, but "any" may allow novel
   ## junctions which span beyond the upstream or downstream annotated
   ## exons.
   ##
   ## bigWigFiles is a vector of bigWig file paths, named by sample;
   ## if supplied, then covL is not used.
   ## covL is a list of "CoverageSignals" objects, analogous to a list
   ## of NumericLists, named by sample.
   ##
   if (suppressPackageStartupMessages(!require(GenomicRanges))) {
      stop("sashimiPlot() requires the GenomicRanges package.");
   }
   spliceOverlapType <- match.arg(spliceOverlapType);

   transformY <- match.arg(transformY);
   if (doPar) {
      parMfrow <- par("mfrow");
   }

   ## WHICHEXONS is a list named by gene, with integer values of exons
   ## to include.
   if (!is.null(WHICHEXONS)) {
      if (!class(WHICHEXONS) %in% c("list")) {
         WHICHEXONS <- nameVector(rep(list(WHICHEXONS),
            length.out=length(inputGenes)),
            inputGenes);
      } else {
         if (is.null(names(WHICHEXONS))) {
            names(WHICHEXONS) <- inputGenes;
         }
      }
   }
   geneSymbolColname <- head(intersect(geneSymbolColname, colnames(values(exonsGR))), 1);
   if (verbose) {
      printDebug("sashimiPlot() :",
         "geneSymbolColname:",
         geneSymbolColname);
   }

   ## Determine the most appropriate label column for inputExonsGR
   exonLabelColname <- intersect(exonLabelColname, colnames(values(exonsGR)));
   if (length(exonLabelColname) == 0) {
      exonLabelColname <- head(provigrep(c("^gene.*exon$", "exon$", "^exon", "^name$"),
         colnames(values(exonsGR))), 1);
      if (length(exonLabelColname) == 0) {
         stop("exonLabelColname is required to be in colnames(exonsGR) but was not.");
      }
   }

   ## Subset exonsGR for efficiency of downstream operations
   exonsGR <- exonsGR[values(exonsGR)[,geneSymbolColname] %in% inputGenes];
   if (verbose) {
      printDebug("length(exonsGR):", formatInt(length(exonsGR)));
   }

   ## Iterate each gene
   inputGenesV <- nameVector(seq_along(inputGenes), inputGenes);
   allGenePlots <- lapply(inputGenesV, function(inputGene1){
      retVals <- list();
      inputGene <- inputGenes[[inputGene1]];
      if (verbose) {
         printDebug("sashimiPrep(): ",
            "inputGene:",
            inputGene);
      }

      ## Subset exons for the gene of interest
      whichExonRows <- which(tolower(values(exonsGR)[,iCol]) %in% inputGene);
      inputExonsGR <- exonsGR[whichExonRows];
      if (length(names(inputExonsGR)) == 0) {
         names(inputExonsGR) <- makeNames(values(inputExonsGR)[,geneSymbolColname]);
      }
      if (verbose) {
         printDebug("sashimiPrep(): ",
            "inputExonsGR (length=",
            length(inputExonsGR),
            "):");
         print(inputExonsGR);
      }

      ## Optionally expand the exon range upstream and downstream, in
      ## case read coverage extends past the annotated exons
      if (length(expandUpDown) > 0 && any(expandUpDown) > 0) {
         if (verbose) {
            printDebug("sashimiPrep(): ",
               "Expanding exons with expandUpDown:",
               expandUpDown);
         }
         inputExonsGR <- expandExonsGR(inputExonsGR,
            expandUpDown=expandUpDown,
            inputExonsGRlabel=inputExonsGRlabel);
      }

      ## Optionally try to read splice junctions directly from a GRanges object
      if (!is.null(spliceGRs)) {
         verbose1 <- verbose;
         spliceCountsDF <- rbindList(lapply(names(spliceGRs), function(spliceGRname){
            if (verbose) {
               printDebug("sashimiPrep(): ",
                  "Loading splice junctions for:",
                  spliceGRname);
            }
            spliceGR <- spliceGRs[[spliceGRname]];

            ## Pull out only the splice junctions within the span of exons
            spliceGRgene <- subsetByOverlaps(spliceGR,
               range(inputExonsGR),
               ignore.strand=ignoreSpliceStrand,
               type=spliceOverlapType);
            iStrand <- as.character(strand(range(inputExonsGR)));
            if (verbose) {
               printDebug("sashimiPrept(): ",
                  "length(spliceGRgene):",
                  length(spliceGRgene));
            }

            if (length(spliceGRgene) > 0) {
               names(spliceGRgene) <- paste0("splice", padInteger(1:length(spliceGRgene)));
               spliceCountsDF1 <- spliceGR2junctionDF(spliceGRgene,
                  inputExonsGR=inputExonsGR,
                  spliceBuffer=spliceBuffer);

               spliceCountsDF1[,"gene"] <- inputGene;
               spliceCountsDF1[,"file"] <- spliceGRname;
               spliceCountsDF1;
            } else {
               NULL;
            }
         }));
         if (verbose) {
            printDebug("sashimiPrep(): ",
               "head(spliceCountsDF):");
            print(head(spliceCountsDF));
         }
         if (length(spliceCountsDF) > 0 &&
               nrow(spliceCountsDF) > 0) {
            if (verbose) {
               printDebug("sashimiPrep(): ",
                  "spliceCountsDF:");
               print(head(spliceCountsDF, 10));
            }
         }
         ## 17aug2016 added minimum read count filtering with minReads
         ## The code below removes any entries which are not "attached" to
         ## an annotated exon, making them effectively invisible.
         ## TODO: allow un-fixed junctions to draw an arc between exons,
         ## or off the edge of the figure.
         ##
         ## 20oct2017: coordFrom, coordTo can be used to display splice junctions
         ##
         if (length(spliceCountsDF) > 0) {
            if (exactJunctionsOnly) {
               if (verbose) {
                  printDebug("sashimiPrep(): ",
                     "filtering for",
                     " exact exon boundaries.");
               }
               spliceCountsDF <- subset(spliceCountsDF,
                  !exonFrom %in% c(NA, "") &
                  !exonTo %in% c(NA, "") &
                  !exonFrom == exonTo);
            }
            if (length(minReads) > 0 && minReads > 0) {
               if (verbose) {
                  printDebug("sashimiPrep(): ",
                     "filtering minReads:",
                     minReads);
               }
               spliceCountsDF <- subset(spliceCountsDF,
                  readCount >= minReads);
            }
         }
         retVals$spliceCountsDF <- spliceCountsDF;
      }
      ## End splice count data.frame

      ## TODO: consider removing WHICHEXONS and assume the input
      ## to this function has already decided which exons to display
      ##
      ## Get the subset of exons to use, if specified
      whichExonsSet <- WHICHEXONS[[inputGene]];
      if (!class(whichExonsSet) %in% c("list")) {
         whichExonsSet <- list(whichExonsSet);
      }
      if (is.null(names(whichExonsSet))) {
         names(whichExonsSet) <- makeNames(rep("whichExonSet",
            length(whichExonsSet)));
      }

         ## Iterate each subset of exons for this gene
         allExonPlots <- lapply(whichExonsSet, function(whichExons){
            if (length(whichExons) == 0 || whichExons %in% c(NA, "NA", "", 0, "0")) {
               whichExons <- seq_along(inputExonsGR);
            }

            ## Allow matching whichExons using the exon label
            if (igrepHas("character", whichExons) && length(exonLabelColname) > 0) {
               printDebug("Matching whichExons to exon labels.");
               whichExons <- rmNA(nameVector(match(whichExons,
                  values(inputExonsGR)[,exonLabelColname]),
                  whichExons));
            } else {
               l1 <- length(inputExonsGR);
               whichExons <- intersect(whichExons, seq_len(l1));
               if (length(exonLabelColname) > 0) {
                  names(whichExons) <- values(inputExonsGR)[whichExons,exonLabelColname]
               }
            }
            if (verbose) {
               printDebug("sashimiPrep(): ",
                  "whichExons (length=",
                  length(whichExons),
                  "):");
               print(whichExons);
            }
            ## Extract exons and sort by chromosome position
            inputExonsGR <- sort(inputExonsGR[whichExons]);

            ## Define color function
            colorFuncX <- function(n,...){colorFunc2(l1)};

            if (igrepHas("_", seqlevels(inputExonsGR))) {
               seqlevels(inputExonsGR) <- unique(unvigrep("_", seqlevels(inputExonsGR)));
            }
            seqlevels(inputExonsGR) <- seqlevelsInUse(inputExonsGR);
            if (compressGaps) {
               if (verbose) {
                  printDebug("sashimiPrep(): ",
                     "running compressGRgaps().");
               }
               ## 20oct2017: added upstream, upstreamGapWidth which
               ## effectively allows any reasonably close coordinates
               ## to be converted to compressed form, with anything outside
               ## the exon range being compressed about 10-fold in width.
               ## Otherwise splice junctions outside the viewing range
               ## become NULL coordinates
               newCoordsGR <- compressGRgaps(inputExonsGR,
                  gapWidth=GAPWIDTH,
                  upstream=1e7,
                  upstreamGapWidth=1e6,
                  downstream=1e7,
                  downstreamGapWidth=1e6);
               ## Get the function to convert reference coordinates
               ## to compressed coordinates.
               ref2compressed <- attr(newCoordsGR, "ref2compressed");
               lookupCoordDF <- attr(newCoordsGR, "lookupCoordDF");
            } else {
               newCoordsGR <- inputExonsGR;
               if (verbose) {
                  printDebug("sashimiPrep(): ",
                     "newCoordsGR (without compressing gaps):");
                  print(head(newCoordsGR, 5));
               }
               lookupCoordDF <- NULL;
               ref2compressed <- NULL;
            }
            ## Optionally propagate labels
            if (labelExons) {
               values(newCoordsGR)[,exonLabelColname] <- values(inputExonsGR)[,exonLabelColname]
            }

            ## Load all coverages upfront
            ## Updated 17aug2016 to use coverages only for the same strand as the exons
            ## We use the actual bigWigName SAMPLE_pos, but assign a name
            ## which removes the _pos, e.g. SAMPLE
            if (length(bigWigFiles) > 0) {
               exonStrand <- as.vector(head(strand(inputExonsGR), 1));
               if (igrepHas("_(pos|neg)", names(bigWigFiles))) {
                  if (exonStrand %in% "+") {
                     bigWigFilesUse <- vigrep("[._]pos", names(bigWigFiles));
                  } else {
                     bigWigFilesUse <- vigrep("[._]neg", names(bigWigFiles));
                  }
                  names(bigWigFilesUse) <- gsub("[._](pos|neg)", "", bigWigFilesUse);
               } else {
                  bigWigFilesUse <- nameVector(names(bigWigFiles));
               }
               if (verbose) {
                  printDebug("sashimiPlot(): ",
                     "bigWigFilesUse:",
                     bigWigFilesUse);
               }
               bigWigCovs <- lapply(bigWigFilesUse, function(bigWigName) {
                  if (verbose) {
                     printDebug("   bigWigName: ", bigWigName,
                        ", running bigWig2numericList()");
                     printDebug("bigWigFiles[bigWigName]:",
                        bigWigFiles[bigWigName]);
                     printDebug("sashimiPlot(): ",
                        "head(inputExonsGR):");
                     print(head(inputExonsGR));
                  }
                  strand(inputExonsGR) <- "*";
                  inputExonsCovNL <- bigWig2numericList(bigWigFiles[bigWigName],
                     inputExonsGR,
                     absValues=absValues,
                     verbose=verbose,
                     ...);
                  if (verbose) {
                     printDebug("length(inputExonsCovNL):",
                        formatInt(length(inputExonsCovNL)));
                  }
                  inputExonsCovNL;
               });
            } else {
               ## Use supplied coverages as a list named by sample, named by exon
               ## Note: if no covL is supplied, coverage is assumed to be zero,
               ## which can allow display of splice features without coverages.
               bigWigCovs <- lapply(covL, function(i){
                  iMatch <- match(names(inputExonsGR), names(i));
                  if (any(is.na(iMatch))) {
                     ## Fill with zeros
                     noMatch <- is.na(iMatch);
                     ix <- lapply(nameVector(width(inputExonsGR[noMatch]), names(inputExonsGR)[noMatch]), function(j){
                        rep(0, j);
                     });
                  } else {
                     ix <- NULL;
                  }
                  c(i[names(inputExonsGR)[!is.na(iMatch)]],
                     ix)[names(inputExonsGR)];
               });
            }
            if (unionCoverages && length(bigWigCovs) > 1) {
               bigWigCovs1 <- bigWigCovs[[1]];
               for (i in tail(seq_along(bigWigCovs), -1)) {
                  bigWigCovs1 <- bigWigCovs1 + bigWigCovs[[i]];
               }
               bigWigCovs <- list(unionCoverage=bigWigCovs1);
               names(bigWigCovs) <- unionLabel;
               rm(bigWigCovs1);
            }

            ## Establish the y-axis bounds
            if (is.null(YMAX)) {
               ymax <- 0;
               for (bigWigName in names(bigWigCovs)) {
                  inputExonsCovNL <- bigWigCovs[[bigWigName]];
                  ymax <- max(c(ymax, max(unlist(do.call(c, inputExonsCovNL)))));
                  #ymax <- max(c(ymax, max(abs(unlist(do.call(c, inputExonsCovNL))))));
                  if (verbose) {
                     printDebug("         ", bigWigName,
                        " ymax: ", ymax,
                        fgText=c("orange", "lightblue"));
                  }
               }
               ymax <- ymax * ylimFactor;
               #if (!is.null(spliceCountsDF)) {
               #   ymax <- ymax * ylimFactor;
               #}
               if (verbose) {
                  printDebug("      ymax: ", ymax,
                     fgText=c("orange", "lightgreen"));
               }
            } else {
               ymax <- YMAX;
            }
            ymin <- 0;
            arcExpandFactor <- 1.6 * arcScaleFactor/50;
            if (alternateSpliceUpDown) {
               #ymin <- -1 * (ymax / ylimFactor / 1.3);
               ymin <- -1 * (ymax / ylimFactor / arcExpandFactor);
            }
            if (transformY %in% "sqrt") {
               #yaxisDiff <- sqrt(ymax) / ylimFactor / 1.3;
               yaxisDiff <- sqrt(ymax) / ylimFactor / arcExpandFactor;
            } else if (transformY %in% "log2") {
               #yaxisDiff <- log2(ymax+1) / ylimFactor / 1.3;
               yaxisDiff <- log2(ymax+1) / ylimFactor / arcExpandFactor;
            } else {
               #yaxisDiff <- ymax / ylimFactor / 1.3;
               yaxisDiff <- ymax / ylimFactor / arcExpandFactor;
            }

            ## Prepare the plot space
            if (doPar) {
               par("mfrow"=rev(plotLayout(length(bigWigCovs))));
               if (!is.null(spliceCountsDF) && doMargins) {
                  if (doTitle) {
                     par("mar"=c(8,
                        par("mar")[c(2,3,4)]));
                  } else {
                     par("mar"=c(2,
                        par("mar")[c(2,3,4)]));
                  }
               }
               par("xpd"=FALSE);
            }
            ##########################################################
            ## Consider drawing coverages after splice junctions?
            allCoveragePlots <- lapply(nameVector(names(bigWigCovs)),
               function(bigWigName){
                  inputExonsCovNL <- bigWigCovs[[bigWigName]];
                  if (verbose) {
                     printDebug("sashimiPlot(): ",
                        " bigWigName:", bigWigName);
                     printDebug("   inputExonsCovNL:");
                     print(head(inputExonsCovNL, 5));
                     printDebug("   newCoordsGR sent to plotGRangesSimple():");
                     print(head(newCoordsGR, 5));
                  }
                  ## flipStrand=TRUE, reverseStrandCoords=TRUE,
                  pgrs1 <- plotGRangesSimple(GR=newCoordsGR,
                     NL=inputExonsCovNL,
                     plotType="coverage",
                     reverseStrandCoords=reverseStrandCoords,
                     flipStrand=flipStrand,
                     drawGrid=TRUE,
                     drawGaps=FALSE,
                     ylim=c(ymin,ymax),
                     showCoords=showCoords,
                     colorFunc=colorFunc2,
                     bty="n",
                     compressGaps=FALSE,
                     plotLabels=labelExons,
                     initialRadius=yaxisDiff/20,
                     transformY=transformY,
                     verbose=FALSE,
                     doPlot=doPlot,
                     ...);
                  if (verbose) {
                     printDebug("whichExons:", whichExons);
                     printDebug("length(newCoordsGR):", length(newCoordsGR));
                     printDebug("length(inputExonsCovNL):", length(inputExonsCovNL));
                  }
                  #if (doPar) {
                  #   par("lend"=parLend);
                  #}
                  inputGeneLabel <- gsub("[.]chr[0-9MTXYWZ]+[.](pos|neg)", "", inputGene);
                  if (doTitle && doPlot) {
                     title(main=paste(inputGeneLabel, "in", bigWigName),
                        font.main=font.main, cex.main=cex.main);
                  }

                  #######################################################
                  ## Draw arcs for splice junction read counts
                  spliceTexts <- NULL;
                  spliceCountsDF1 <- NULL;
                  if (verbose) {
                     printDebug("spliceCountsDF:");
                     print(head(spliceCountsDF));
                  }
                  if (!is.null(spliceCountsDF) && inputGene %in% spliceCountsDF[,"gene"]) {
                     cov1 <- rbindList(pgrs1[["cov1"]]);
                     if (verbose) {
                        printDebug("length(pgrs1[[cov1]]):", length(pgrs1[["cov1"]]));
                        printDebug("dim(cov1):", dim(cov1));
                        ch(cov1);
                     }
                     rownames(cov1) <- whichExons;
                     exonMaxSteps <- diff(range(as.numeric(rownames(cov1))));
                     if (unionCoverages) {
                        spliceCountsDF1 <- spliceCountsDF[spliceCountsDF[,"gene"] %in% inputGene,];
                        ## 20oct2017: Add some key that handles junctions
                        ## which do not snap to a known exon
                        ##
                        ## TODO: change naming convention for junctions not matching
                        ## an exon, to use closest exon name, and coordinate distance
                        ## from that exon. The goal is for someone to recognize
                        ## whether a junction was between exon 2 and 3, e.g.
                        ## 3a.-1700 would mean 1700 bases upstream exon "3a".
                        ## Another option is "2b.1700.3a" which would mean 1700
                        ## bases downstream "2b". In that case, "1700.1a" would mean
                        ## 1700 bases upstream the first exon "1a". Similarly,
                        ## "7f.1700" would mean 1700 bases downstream the last exon
                        ## "7f".
                        fromStr <- ifelse(!is.na(spliceCountsDF1[,"exonFrom"]),
                           spliceCountsDF1[,"exonFrom"],
                           paste0("c", spliceCountsDF1[,"coordFrom"]));
                        toStr <- ifelse(!is.na(spliceCountsDF1[,"exonTo"]),
                           spliceCountsDF1[,"exonTo"],
                           paste0("c", spliceCountsDF1[,"coordTo"]));
                        spliceCountsDF1[,"exonFromTo"] <- paste0(fromStr, "_", toStr);
                        #spliceCountsDF1[,"exonFromTo"] <-
                        #   pasteByRow(spliceCountsDF1[,c("exonFrom", "exonTo")], sep="_");
                        spliceCountsDF1readCts <- shrinkMatrix(spliceCountsDF1[,"readCount"],
                           groupBy=spliceCountsDF1[,"exonFromTo"],
                           shrinkFunc=sum);
                        ## Apply the counts to the original data.frame,
                        ## taking only the first row, keeping all the
                        ## annotation columns as-is
                        spliceCountsDF1 <- spliceCountsDF1[match(spliceCountsDF1readCts[,"groupBy"],
                           spliceCountsDF1[,"exonFromTo"]),];
                        spliceCountsDF1[,"readCount"] <- spliceCountsDF1readCts[,"x"];
                        spliceCountsDF1[,"file"] <- bigWigName;
                        maxCts <- max(spliceCountsDF1[,"readCount"]);
                     } else {
                        spliceSubset <- which(spliceCountsDF[,"file"] %in% bigWigName &
                              spliceCountsDF[,"gene"] %in% inputGene);
                        spliceCountsDF1 <- spliceCountsDF[spliceSubset,,drop=FALSE];
                        maxCts <- max(spliceCountsDF[spliceCountsDF[,"gene"] %in% inputGene,"readCount"]);
                     }
                     ## 20oct2017 convert reference coordinates to compressed
                     if (compressGaps) {
                        spliceCountsDF1[,"coordFrom"] <- round(ref2compressed(spliceCountsDF1[,"coordFrom"]));
                        spliceCountsDF1[,"coordTo"] <- round(ref2compressed(spliceCountsDF1[,"coordTo"]));
                        if (verbose) {
                           printDebug("sashimiPlot(): ",
                              "lookupCoordDF:");
                           ch(lookupCoordDF, 100);
                        }
                     }

                     spliceCountsDF1 <- mixedSortDF(spliceCountsDF1,
                        byCols=match(c("coordFrom", "coordTo"),
                           colnames(spliceCountsDF1)));
                     if (verbose) {
                        printDebug("sashimiPlot(): ",
                           "head(spliceCountsDF1):");
                        print(head(spliceCountsDF1), 100);
                     }
                     ## 20oct2017: Add some key that handles junctions
                     ## which do not snap to a known exon
                     fromStr <- ifelse(!is.na(spliceCountsDF1[,"exonFrom"]),
                        spliceCountsDF1[,"exonFrom"],
                        spliceCountsDF1[,"coordFrom"]);
                     toStr <- ifelse(!is.na(spliceCountsDF1[,"exonTo"]),
                        spliceCountsDF1[,"exonTo"],
                        spliceCountsDF1[,"coordTo"]);
                     spliceCountsDF1[,"exonFromTo"] <- paste0(fromStr, "_", toStr);
                     if (doPlot) {
                        par("xpd"=TRUE);
                     }

                     ## Sometimes there are still multiple entries per junction,
                     ## so we shrink them here for now
                     ## 20oct2017 added exonFromTo for non-canonical junctions
                     if (verbose) {
                        printDebug("head(spliceCountsDF1) before compression:");
                        ch(spliceCountsDF1);
                     }
                     spliceCountsDF1 <- shrinkDataFrame(spliceCountsDF1,
                        #groupBy=c("exonFrom", "exonTo", "gene", "file"),
                        groupBy=c("exonFromTo", "gene", "file"),
                        numShrinkFunc=list(readCount=sum, exonFromDist=mean, exonToDist=mean));

                     ## 06oct2017: added exonDiff so arcs could be ordered by the number of
                     ## exon hops, which should allow their left and right edges to be stacked
                     ## properly, with respect to other splice junction arcs
                     ##
                     ## 23oct2017: changed exonDiff to represent the coordinate
                     ## distance, instead of exon jumps, since not all junctions
                     ## will snap to an annotated exon
                     maxExonNum <- max(c(spliceCountsDF1$exonFrom,
                        spliceCountsDF1$exonTo), na.rm=TRUE);
                     #spliceCountsDF1[,"exonDiff"] <- abs(
                     #   rmNA(as.numeric(spliceCountsDF1$exonTo), naValue=maxExonNum+1) -
                     #   rmNA(as.numeric(spliceCountsDF1$exonFrom), naValue=0) );
                     spliceCountsDF1[,"exonDiff"] <- abs(spliceCountsDF1$coordTo -
                           spliceCountsDF1$coordFrom);
                     if (alternateSpliceUpDown) {
                        ## 20oct2017: changed to use coordFrom,coordTo
                        ## 23oct2017: changed to use coordFrom,exonDiff
                        spliceCountsDF1 <- mixedSortDF(spliceCountsDF1,
                           byCols=match(c("coordFrom","exonDiff"),
                              colnames(spliceCountsDF1))*c(1,+1));
                        #spliceCountsDF1 <- mixedSortDF(spliceCountsDF1, byCols=match(c("exonFrom","exonDiff"),
                        #   colnames(spliceCountsDF1))*c(1,+1));
                     } else {
                        ## 20oct2017: changed to use coordFrom,coordTo
                        ## 23oct2017: changed to use coordFrom,exonDiff
                        spliceCountsDF1 <- mixedSortDF(spliceCountsDF1,
                           byCols=match(c("coordFrom","exonDiff"),
                              colnames(spliceCountsDF1)));
                        #spliceCountsDF1 <- mixedSortDF(spliceCountsDF1,
                        #   byCols=match(c("exonFrom","exonDiff"),
                        #   colnames(spliceCountsDF1)));
                     }
                     ## 20oct2017: changed to use coordFrom
                     spliceCountsDF1[,"fromY"] <- unlist(lapply(split(spliceCountsDF1,
                        spliceCountsDF1$coordFrom), function(iDF1){
                           cumsum(c(0, iDF1$readCount))[seq_len(nrow(iDF1))];
                        }));
                     ## Now do the same for the exonTo
                     if (alternateSpliceUpDown) {
                        spliceCountsDF1 <- mixedSortDF(spliceCountsDF1,
                           byCols=match(c("coordTo","exonDiff"),
                              colnames(spliceCountsDF1))*c(1,+1));
                     } else {
                        spliceCountsDF1 <- mixedSortDF(spliceCountsDF1,
                           byCols=match(c("coordTo","exonDiff"),
                              colnames(spliceCountsDF1)));
                     }
                     ## 18dec2017: change the stacking order of splice junctions
                     ## so the shorter distance comes last.
                     spliceCountsDF1[,"toY"] <- unlist(lapply(split(spliceCountsDF1,
                        #spliceCountsDF1$exonTo), function(iDF1){
                        spliceCountsDF1$coordTo), function(iDF1){
                           cumsum(c(0, iDF1$readCount))[seq_len(nrow(iDF1))];
                        }));
                     if (verbose) {
                        printDebug("spliceCountsDF1:");
                        print(spliceCountsDF1);
                     }

                     ## Iterate each splice junction and draw an arc with text label
                     spliceTexts <- lapply(1:nrow(spliceCountsDF1), function(spliceRow) {
                        exonFrom <- as.character(spliceCountsDF1[spliceRow,"exonFrom"]);
                        exonTo <- as.character(spliceCountsDF1[spliceRow,"exonTo"]);
                        readCts <- spliceCountsDF1[spliceRow,"readCount"];
                        exonDiff <- abs(as.numeric(exonTo) - as.numeric(exonFrom));

                        ## New code for arc segments
                        if (transformY %in% "none") {
                           transFunc <- c;
                        } else if (transformY %in% "sqrt") {
                           transFunc <- sqrt;
                        } else if (transformY %in% "log2") {
                           transFunc <- function(x){log2(1+x)};
                        }
                        #if (exonFrom %in% rownames(cov1) && exonTo %in% rownames(cov1)) {
                        flipLabel <- (alternateSpliceUpDown && spliceCountsDF1[spliceRow,"exonDiff"]<=1);
                        if (as.character(strand(inputExonsGR)[1]) %in% "-") {
                           ## Negative strand
                           flipLabel <- flipLabel * -1;
                        }
                        ## 20oct2017: changed to use coordFrom, coordTo
                        x0 <- spliceCountsDF1[spliceRow,"coordFrom"];
                        x1 <- spliceCountsDF1[spliceRow,"coordTo"];
                        ## 20oct2017: changed to draw arc based upon horizontal
                        ## segment, then linearly scale from fromY to toY.
                        ## This change avoid the case where the segment has
                        ## a high slope, and would otherwise cause the arc
                        ## to be spin nearly sideways.
                        as1 <- arcSegments(x0=x0,
                           x1=x1,
                           y0=transFunc(spliceCountsDF1[spliceRow,"fromY"]),
                           y1=transFunc(spliceCountsDF1[spliceRow,"fromY"]),
                           doPlot=FALSE,
                           useAspect=useAspect,
                           verbose=FALSE,
                           scaling1=0.45 * sign(flipLabel-0.5))[[1]];
                        as1$y <- seq(from=0,
                           to=(spliceCountsDF1[spliceRow,"toY"]-
                                 spliceCountsDF1[spliceRow,"fromY"]),
                           length.out=length(as1$y)) + as1$y;
                        as2 <- as1;
                        as2$y <- as2$y + spliceCountsDF1[spliceRow,"readCount"];
                        ## Define a suitable midpoint for labeling
                        arcX12 <- mean(as2$x);
                        if (flipLabel) {
                           arcY12 <- as1$y[round(length(as1$y)/2)];
                        } else {
                           arcY12 <- as2$y[round(length(as2$y)/2)];
                        }
                        arcY12 <- (as1$y[round(length(as1$y)/2)] +
                              as2$y[round(length(as2$y)/2)])/2;
                        ## Convert the two arcs into a polygon
                        as1x1 <- c(head(as1$x, 1)-0.0001,
                           as1$x,
                           tail(as1$x, 1)+0.0001);
                        as2x1 <- c(head(as2$x, 1)-0.0001,
                           as2$x,
                           tail(as2$x, 1)+0.0001);
                        as1y1 <- c(head(as1$y, 1),
                           as1$y,
                           tail(as1$y, 1));
                        as2y1 <- c(head(as2$y, 1),
                           as2$y,
                           tail(as2$y, 1));
                        as1x <- c(as1x1, rev(as2x1), head(as1x1, 1));
                        as1y <- c(as1y1, rev(as2y1), head(as1y1, 1));
                        ## Draw the polygon
                        if (doPlot) {
                           polygon(x=as1x,
                              y=as1y,
                              border=alpha2col(readColor, alpha=1),
                              col=readColor,
                              lwd=1);
                        }
                        ## Adjust labels either with leading spaces or trailing spaces,
                        ## depending if the arcs are above or below the line
                        readCountPadding <- "";
                        readCtLabels <- ifelse(flipLabel,
                           paste0(formatInt(readCts), readCountPadding),
                           paste0(readCountPadding, formatInt(readCts)));
                        readCtAdj <- c(0.5,ifelse(flipLabel, 1.5, -0.5));
                        readCtAdj <- rep(0.5, length(readCtAdj));
                        arcLwd <- maxLwd;
                        list(x=arcX12,
                           y=arcY12,
                           xpoly=as1x,
                           ypoly=as1y,
                           polyColor=readColor,
                           labels=readCtLabels,
                           adj=readCtAdj,
                           col=makeColorDarker(readColor, fixAlpha=1, darkFactor=5),
                           cex=readCountCex,
                           font=readCountFont,
                           srt=srtSplice,
                           arcLwd=arcLwd,
                           readCts=readCts,
                           maxCts=maxCts,
                           maxLwd=maxLwd);
                     });
                     #if (verbose) {
                     #   printDebug("spliceTexts:");
                     #   print(head(spliceTexts));
                     #}
                     if (doPlot) {
                        for (spliceText in spliceTexts) {
                           if (!is.null(spliceText)) {
                              text(x=spliceText$x,
                                 y=spliceText$y,
                                 labels=spliceText$labels, #adj=c(0.5,0.5),#adj=spliceText$adj,
                                 col=spliceText$col,
                                 cex=spliceText$cex,
                                 font=spliceText$font,
                                 srt=spliceText$srt);
                              srt=#srtSplice);
                                 spliceTextDF <- data.frame(x=spliceText$x,
                                    y=spliceText$y,
                                    labels=spliceText$labels,
                                    adj=spliceText$adj,
                                    col=spliceText$col,
                                    cex=spliceText$cex,
                                    font=spliceText$font,
                                    srt=spliceText$srt);
                              #srt=srtSplice);
                              if (verbose) {
                                 printDebug("spliceTextDF:");
                                 ch(spliceTextDF, maxRows=100);
                              }
                           }
                        }
                        par("xpd"=FALSE);
                     }
                  }
                  list(inputExonsGR=inputExonsGR,
                     pgrs1=pgrs1,
                     spliceTexts=spliceTexts,
                     spliceCountsDF1=spliceCountsDF1,
                     lookupCoordDF=lookupCoordDF,
                     ref2compressed=ref2compressed);
               });
         });
      });
   if (doPar) {
      par("mfrow"=parMfrow);
   }
   invisible(allGenePlots);
}

#' Expand exon GRanges upstream and downstream
#'
#' Expand exon GRanges upstream and downstream
#'
#' This function is called by `sashimiPrep()`.
#'
#' @export
expandExonsGR <- function
(inputExonsGR,
expandUpDown=c(0,0),
expandLabels=c("upstream","downstream"),
inputExonsGRlabel="exonLabel",
...)
{
   ## Helper function to add upstream and downstream ranges to
   ## the exons being displayed
   expandUpDown <- rep(c(rep(expandUpDown, length.out=2), 0),
      length.out=2);
   if (expandUpDown[1] > 0) {
      if ("-" %in% as.character(strand(inputExonsGR))) {
         upGR <- tail(sort(inputExonsGR), 1);
      } else {
         upGR <- head(sort(inputExonsGR), 1);
      }
      upGR <- flank(upGR, width=expandUpDown[1], start=TRUE);
      if (inputExonsGRlabel %in% names(values(upGR))) {
         values(upGR)[,inputExonsGRlabel] <- expandLabels[1];
      }
      names(upGR) <- paste0(names(upGR), "_up");
   } else {
      upGR <- inputExonsGR[0];
   }
   if (expandUpDown[2] > 0) {
      if ("-" %in% as.character(strand(inputExonsGR))) {
         dnGR <- head(sort(inputExonsGR), 1);
      } else {
         dnGR <- tail(sort(inputExonsGR), 1);
      }
      dnGR <- flank(dnGR, width=expandUpDown[2], start=FALSE);
      if (inputExonsGRlabel %in% names(values(dnGR))) {
         values(dnGR)[,inputExonsGRlabel] <- expandLabels[2];
      }
      names(dnGR) <- paste0(names(dnGR), "_dn");
   } else {
      dnGR <- inputExonsGR[0];
   }
   return(c(upGR, dnGR));
}


#' Find closest exon to splice junction ends
#'
#' Find closest exon to splice junction ends
#'
#' @export
closestExonToJunctions <- function
(spliceGRgene,
 exonsGR,
 flipNegativeStrand=TRUE,
 reportActualCoords=FALSE,
 verbose=TRUE,
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
   updateColnames <- c("distFrom", "distTo", "nameFrom", "nameTo", "genesDiffer", "genesMatch", "tooFarFrom", "tooFarTo", "tooFar");
   if (any(updateColnames %in% colnames(values(spliceGRgene)))) {
      values(spliceGRgene) <- values(spliceGRgene)[,setdiff(colnames(spliceGRgene), updateColnames),drop=FALSE];
   }

   ## Distance from splice start to exon end
   ## ignore.strand=TRUE on resize makes the method search the left side of a splice site with the right side of an exon.
   if (verbose) {
      printDebug("closestExonToJunctions(): ",
         "Finding closest exons for splice starts.");
   }
   spliceStartExonEndD1 <- as.data.frame(
      distanceToNearest(
         resize(spliceGRgene,
            width=1,
            fix="start",
            ignore.strand=TRUE),
         resize(exonsGR,
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
      start(
         resize(spliceGRgene,
            width=1,
            fix="start",
            ignore.strand=TRUE)[spliceStartExonEndD1[,"queryHits"]]) -
      start(
         resize(exonsGR,
            width=1,
            fix="end",
            ignore.strand=TRUE)[spliceStartExonEndD1[,"subjectHits"]]));

   ## TODO: add the actual end coordinate of the exon matched
   ## end(exonsGR[spliceStartExonEndD1[,"subjectHits"]])
   spliceStartExonEndDactual1end <- end(exonsGR[spliceStartExonEndD1[,"subjectHits"]]);
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
      distanceToNearest(
         resize(spliceGRgene,
            width=1,
            fix="end",
            ignore.strand=TRUE),
         resize(exonsGR,
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
      start(
         resize(spliceGRgene,
            width=1,
            fix="end",
            ignore.strand=TRUE)[spliceEndExonStartD1[,"queryHits"]]) -
      start(
         resize(exonsGR,
            width=1,
            fix="start",
            ignore.strand=TRUE)[spliceEndExonStartD1[,"subjectHits"]]));

   ## TODO: add the actual start coordinate of the exon matched
   ## start(exonsGR[spliceEndExonStartD1[,"subjectHits"]])
   spliceEndExonStartDactual1start <- start(exonsGR[spliceEndExonStartD1[,"subjectHits"]]);
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
   values(spliceGRgene)[spliceStartExonEndD1[,"queryHits"],"nameFrom"] <- names(exonsGR[spliceStartExonEndD1[,"subjectHits"]]);
   values(spliceGRgene)[spliceEndExonStartD1[,"queryHits"],"nameTo"] <- names(exonsGR[spliceEndExonStartD1[,"subjectHits"]]);

   ## Update the stranded distance on the start side (from) and end side (to) of the splice junction
   values(spliceGRgene)[spliceStartExonEndD1[,"queryHits"],"distFrom"] <- spliceStartExonEndD1[,"strandedDistance"];
   values(spliceGRgene)[spliceEndExonStartD1[,"queryHits"],"distTo"] <- spliceEndExonStartD1[,"strandedDistance"];

   ## TODO: report the actual coordinate boundary being matched
   if (reportActualCoords) {
      if (verbose) {
         printDebug("closestExonToJunctions(): ",
            "Reporting actual coords.");
      }
      values(spliceGRgene)[spliceStartExonEndD1[,"queryHits"],"coordFrom"] <- spliceStartExonEndDactual1end;
      values(spliceGRgene)[spliceEndExonStartD1[,"queryHits"],"coordTo"] <- spliceEndExonStartDactual1start;
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
      switchCols1 <- intersect(fromToCols, colnames(values(spliceGRgene)));
      switchCols2 <- as.vector(matrix(nrow=2, switchCols1)[2:1,])
      negStrand <- (as.vector(strand(spliceGRgene)) %in% "-");
      if (any(negStrand)) {
         values(spliceGRgene)[negStrand,switchCols1] <-
            values(spliceGRgene)[negStrand,switchCols2];
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
#' @export
spliceGR2junctionDF <- function
(spliceGRgene,
exonsGR,
spliceBuffer=3,
spliceGRname="score",
geneExonSep="(:|_exon)",
useOnlyValidEntries=FALSE,
renameTooFar=TRUE,
scoreColname="score",
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

   ## Vectorized logic:
   ## - find closest start/end for each splice start/end
   ## - subset for splice start/end having same gene
   ##
   ## Call a method which encapsulates distanceToNearest(), calculates stranded distance,
   ## and knows to match splice start with exon end, etc.
   #spliceGRgeneVals <- closestExonToJunctions(spliceGRgene=spliceGRgene[,scoreColname], exonsGR=exonsGR,
   #   flipNegativeStrand=flipNegativeStrand, ...);
   spliceGRgene <- closestExonToJunctions(spliceGRgene=spliceGRgene[,scoreColname],
      exonsGR=exonsGR,
      flipNegativeStrand=flipNegativeStrand,
      verbose=verbose,
      ...)$spliceGRgene;
   #spliceGRgene <- spliceGRgeneVals$spliceGRgene;

   ## Check when genes match for the two junction sites
   ## Note: we changed from genesDiffer, because in cases where no gene is returned, the two
   ## sides of the junction would be equal, but still not represent the desired outcome.
   genesMatch <- (!is.na(values(spliceGRgene)[,"nameFrom"]) & !is.na(values(spliceGRgene)[,"nameTo"]) &
         (gsub(paste0(geneExonSep, ".*$"), "", values(spliceGRgene)[,"nameFrom"]) ==
               gsub(paste0(geneExonSep, ".*$"), "", values(spliceGRgene)[,"nameTo"]) ) );
   values(spliceGRgene)[,"genesMatch"] <- genesMatch;
   numGenesMatch <- sum(values(spliceGRgene)[,"genesMatch"]);
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
   tooFarFrom <- (abs(values(spliceGRgene)[,"distFrom"]) > spliceBuffer |
         is.na(values(spliceGRgene)[,"distFrom"]));
   tooFarTo <- (abs(values(spliceGRgene)[,"distTo"]) > spliceBuffer |
         is.na(values(spliceGRgene)[,"distTo"]));
   tooFar <- (tooFarFrom | tooFarTo);
   values(spliceGRgene)[,"tooFarFrom"] <- tooFarFrom;
   values(spliceGRgene)[,"tooFarTo"] <- tooFarTo;
   values(spliceGRgene)[,"tooFar"] <- tooFar;
   numTooFar <- sum(values(spliceGRgene)[,"tooFar"]);
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
         values(spliceGRgene)[tooFarFrom,"nameFrom"] <- paste(
            values(spliceGRgene)[tooFarFrom,"nameFrom"],
            values(spliceGRgene)[tooFarFrom,"distFrom"],
            sep=".");
      }
      if (any(tooFarTo)) {
         values(spliceGRgene)[tooFarTo,"nameTo"] <- paste(
            values(spliceGRgene)[tooFarTo,"nameTo"],
            values(spliceGRgene)[tooFarTo,"distTo"],
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
      toUse <- which(values(spliceGRgene)[,"genesMatch"] & !values(spliceGRgene)[,"tooFar"]);
   } else {
      ## Note that we still only return entries for which there is some nearby exon
      ## TODO: allow for un-named entries to have a temporary name for the purpose of
      ## returning the values.
      if (verbose) {
         printDebug("spliceGR2junctionDF(): ",
            "All entries having any nearest gene will be used, without regard to matching genes or distance from exon.");
      }
      toUse <- which(!is.na(values(spliceGRgene)[,"nameFrom"]) & !is.na(values(spliceGRgene)[,"nameTo"]));
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

   spliceColnames <- c("seqnames", "start", "end", "nameFrom", "nameTo", scoreColname, "strand");
   spliceCountsDF <- as.data.frame(spliceGRgene[toUse])[,spliceColnames,drop=FALSE];
   spliceCountsDF[,"nameFromTo"] <- pasteByRow(spliceCountsDF[,c("nameFrom","nameTo"),drop=FALSE],
      sep=" ",
      na.rm=TRUE);

   if (verbose) {
      printDebug("spliceGR2junctionDF(): ",
         "Shrinking matrix to combine counts spanning the same junctions.");
   }
   spliceCountsDFshrunk <- shrinkMatrix(spliceCountsDF[,scoreColname],
      groupBy=spliceCountsDF[,"nameFromTo"],
      shrinkFunc=sum);
   spliceCountsDFshrunk <- cbind(spliceCountsDFshrunk,
      rbindList(strsplit(spliceCountsDFshrunk[,"groupBy"], " ")));
   colnames(spliceCountsDFshrunk) <- c("nameFromTo", spliceGRname, "nameFrom", "nameTo");

   ## Now add the start and end coordinates, so the results can be plotted
   if (verbose) {
      printDebug("spliceGR2junctionDF(): ",
         "Adding ref, strand, start, end.");
   }
   matchFromTo <- match(spliceCountsDFshrunk[,"nameFromTo"], spliceCountsDF[,"nameFromTo"]);
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

   ## TODO: flip the from/to for negative stranded genes
   if (1 == 2 && flipNegativeStrand) {
      negStrand <- (spliceCountsDFshrunk[,"strand"] %in% "-");
      ## Exchange values in the exonFrom and exonTo column for negative strand entries
      ## using direct assignment so we avoid using temporary vectors
      spliceCountsDFshrunk[negStrand,c("exonFrom", "exonTo")] <- spliceCountsDFshrunk[negStrand,c("exonTo", "exonFrom")];
      spliceCountsDFshrunk[negStrand,c("nameFrom", "nameTo")] <- spliceCountsDFshrunk[negStrand,c("nameTo", "nameFrom")];
      spliceCountsDFshrunk[negStrand,c("geneSymbolFrom", "geneSymbolTo")] <-
         spliceCountsDFshrunk[negStrand,c("geneSymbolTo", "geneSymbolFrom")];

   }

   ## Create junctionID only after we decided how to orient exonFrom and exonTo
   ## for negative strand entries, avoiding re-creating the name for those rows
   spliceCountsDFshrunk[,"junctionID"] <- pasteByRow(
      spliceCountsDFshrunk[,c("exonFrom","exonTo")],
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

#' Compress GRanges gaps to a fixed width
#'
#' Compress GRanges gaps to a fixed width
#'
#' @export
compressGRgaps <- function
(GR,
 gapWidth=NULL,
 keepValues=FALSE,
 upstream=10000,
 upstreamGapWidth=gapWidth*3,
 downstream=10000,
 downstreamGapWidth=gapWidth*3,
 refGR=NULL,
 refGRnames=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to take a GRanges object, usually a set of transcript exons
   ## and compress coordinates so the gaps between segments are each a fixed-width.
   ##
   ## keepValues=TRUE will try to keep all the values(GR) annotation columns
   ##
   ## It also returns a conversion lookup function in attribute "ref2compressed"
   ##
   ## gapWidth=NULL will set gapWidth to half the largest contiguous range,
   ## or a minimum gapWidth of 10.
   ##
   ## refGR is an optional GRanges object used to define a reference set of
   ## gaps, therefore ignoring gaps in GR and using refGR in its place.
   ##
   ## refGRnames is an optional vector of GR names, which is used to subset
   ## refGR if supplied, or GR if refGR is not supplied.
   ##
   ## TODO: allow fixed feature width as well as fixed gap width, for the
   ## scenario where some features are substantially wider than others and
   ## would otherwise cause the small features to be too small to visualize
   ## effectively.
   ##
   ## TODO: fix issue when given features with no gap
   ##

   ## First check for refGR
   if (length(refGR) > 0 && igrepHas("GRanges", class(refGR))) {
      if (verbose) {
         printDebug("compressGRgaps(): ",
            "Using refGR to define gaps.");
      }
      if (length(refGRnames) > 0) {
         if (length(names(refGR)) > 0) {
            refGR <- refGR[names(refGR) %in% refGRnames];
         } else {
            stop("refGRnames is supplied but names(refGR) is empty.");
         }
      }
      cgL <- compressGRgaps(GR=refGR,
         gapWidth=gapWidth,
         keepValues=keepValues,
         upstream=upstream,
         upstreamGapWidth=upstreamGapWidth,
         downstream=downstream,
         downstreamGapWidth=downstreamGapWidth,
         refGR=NULL,
         refGRnames=NULL,
         verbose=verbose,
         ...);
      ref2compressed <- attr(cgL, "ref2compressed");
      gapWidth <- attr(cgL, "gapWidth");
      GR <- ref2compressedGR(GR,
         ref2compressed=ref2compressed);
      attr(GR, "gapWidth") <- gapWidth;
      return(GR);
   }

   if (length(gapWidth) == 0) {
      gapWidth <- max(c(10,
         round(max(width(reduce(GR))) / 2)));
      if (verbose) {
         printDebug("compressGRgaps(): ",
            "determined gapWidth:",
            format(gapWidth,
               scientific=FALSE,
               big.mark=","));
      }
   }
   if (length(upstreamGapWidth) == 0) {
      upstreamGapWidth <- gapWidth * 3;
   }
   if (length(downstreamGapWidth) == 0) {
      downstreamGapWidth <- gapWidth * 3;
   }

   ## Define the exon-exon distance, to see if they are adjacent or not
   #disGR <- disjoin(GR);
   disGR <- reduce(GR);
   if (is.null(names(disGR))) {
      names(disGR) <- makeNames(rep("disGR", length(disGR)));
   }
   if (length(disGR) > 1) {
      exonGaps <- sapply(head(seq_along(disGR), -1), function(i1){
         i2 <- i1 + 1;
         distance(disGR[i1], disGR[i2]) > 0;
      });
      names(exonGaps) <- head(names(disGR), -1);
   } else {
      exonGaps <- NULL;
   }

   ## Now reset coordinates using fixed gap width
   #newWidths <- head(intercalate(width(disGR),
   #   (gapWidth)*exonGaps), -1);
   newWidths <- head(intercalate(width(disGR),
      (gapWidth)*exonGaps),
      length(disGR)*2-1);

   ##
   newCoords <- cumsum(newWidths);
   newCoordsM <- matrix(c(0, newCoords), ncol=2, byrow=TRUE);
   newCoordsM[,1] <- newCoordsM[,1] + 1;

   lookupCoordDF <- mixedSortDF(unique(data.frame(
      refCoord=c(start(disGR),end(disGR)),
      coord=c(newCoordsM[,1],newCoordsM[,2])
   )));
   if (upstream > 0) {
      lookupCoordDF <- rbind(
         data.frame(refCoord=head(lookupCoordDF$refCoord,1)-upstream,
            coord=head(lookupCoordDF$coord,1)-upstreamGapWidth),
         lookupCoordDF);
   }
   if (downstream > 0) {
      lookupCoordDF <- rbind(
         lookupCoordDF,
         data.frame(refCoord=tail(lookupCoordDF$refCoord,1)+downstream,
            coord=tail(lookupCoordDF$coord,1)+downstreamGapWidth)
      );
   }
   ref2compressed <- approxfun(x=lookupCoordDF[,1],
      y=lookupCoordDF[,2],
      method="linear");
   newCoordsGR <- GR;
   values(newCoordsGR)[,"refStart"] <- start(newCoordsGR);
   values(newCoordsGR)[,"refEnd"] <- end(newCoordsGR);
   ## Re-define the start and end, but since they must be
   ## valid start-end ranges we must set them both at once
   ranges(newCoordsGR) <- IRanges(start=ref2compressed(start(GR)),
      end=ref2compressed(end(GR)))

   ## Add some attributes for later ease of use
   attr(newCoordsGR, "ref2compressed") <- ref2compressed;
   attr(newCoordsGR, "lookupCoordDF") <- lookupCoordDF;
   attr(newCoordsGR, "gapWidth") <- gapWidth;

   ## Convert reference coordinates to feature-based coordinates,
   ## e.g. to coordinates relative to the length of a transcript
   if (1 == 2) {
      newCoordsGRo <- order(GR);
      newCoordsGR <- newCoordsGR[newCoordsGRo];
      txCoords <- matrix(ncol=2, data=c(rep(1, length(newCoordsGR)), cumsum(width(newCoordsGR))));
      txCoords[2:length(newCoordsGR)] <- head(txCoords[,2], -1) + 1;
      if ("-" %in% strand(newCoordsGR[1])) {
         colnames(txCoords) <- c("featureEnd","featureStart");
      } else {
         colnames(txCoords) <- c("featureStart","featureEnd");
      }
      values(newCoordsGR)[,"featureStart"] <- txCoords[,"featureStart"];
      values(newCoordsGR)[,"featureEnd"] <- txCoords[,"featureEnd"];
      txCoords <- txCoords[,c("featureStart","featureEnd")];
   }

   newCoordsGR;
}

#' Create a ref2compressed function to compress GR gaps
#'
#' Create a ref2compressed function to compress GR gaps
#'
#' This function takes a set of GRanges which are to be maintained
#' with fixed aspect ratio, and it defines a function to compress
#' coordinates of the gaps between GRanges features.
#'
#' @export
make_ref2compressed <- function
(GR,
gapWidth=NULL,
keepValues=FALSE,
upstream=50000,
upstreamGapWidth=gapWidth*3,
downstream=50000,
downstreamGapWidth=gapWidth*3,
refGR=NULL,
refGRnames=NULL,
verbose=FALSE,
...)
{
   ## Purpose is to use a GRanges object to train a coordinate
   ## compression function, which assigns the gaps to a fixed
   ## width, but allows GRanges features to keep their original
   ## width.
   if (length(gapWidth) == 0) {
      gapWidth <- max(c(10,
         round(max(width(reduce(GR))) / 2)));
      if (verbose) {
         printDebug("compressGRgaps(): ",
            "determined gapWidth:",
            format(gapWidth,
               scientific=FALSE,
               big.mark=","));
      }
   }
   if (length(upstreamGapWidth) == 0) {
      upstreamGapWidth <- gapWidth * 3;
   }
   if (length(downstreamGapWidth) == 0) {
      downstreamGapWidth <- gapWidth * 3;
   }

   ## Define the exon-exon distance, to see if they are adjacent or not
   #disGR <- disjoin(GR);
   disGR <- reduce(GR);
   if (is.null(names(disGR))) {
      names(disGR) <- makeNames(rep("disGR", length(disGR)));
   }
   if (length(disGR) > 1) {
      exonGaps <- sapply(head(seq_along(disGR), -1), function(i1){
         i2 <- i1 + 1;
         distance(disGR[i1], disGR[i2]) > 0;
      });
      names(exonGaps) <- head(names(disGR), -1);
   } else {
      exonGaps <- NULL;
   }

   ## Now reset coordinates using fixed gap width
   #newWidths <- head(intercalate(width(disGR),
   #   (gapWidth)*exonGaps), -1);
   newWidths <- head(intercalate(width(disGR),
      (gapWidth)*exonGaps),
      length(disGR)*2-1);

   ##
   newCoords <- cumsum(newWidths);
   newCoordsM <- matrix(c(0, newCoords), ncol=2, byrow=TRUE);
   newCoordsM[,1] <- newCoordsM[,1] + 1;

   lookupCoordDF <- mixedSortDF(unique(data.frame(
      refCoord=c(start(disGR),end(disGR)),
      coord=c(newCoordsM[,1],newCoordsM[,2])
   )));
   if (upstream > 0) {
      lookupCoordDF <- rbind(
         data.frame(refCoord=head(lookupCoordDF$refCoord,1)-upstream,
            coord=head(lookupCoordDF$coord,1)-upstreamGapWidth),
         lookupCoordDF);
   }
   if (downstream > 0) {
      lookupCoordDF <- rbind(
         lookupCoordDF,
         data.frame(refCoord=tail(lookupCoordDF$refCoord,1)+downstream,
            coord=tail(lookupCoordDF$coord,1)+downstreamGapWidth)
      );
   }
   ## TODO: expand the range of coordinates so approxfun() will not fail
   ref2compressed <- approxfun(x=lookupCoordDF[,1],
      y=lookupCoordDF[,2],
      method="linear");

   ## Add some attributes for later ease of use
   ref2compressedGR <- function(GR) {
      if (all(c("refStart","refEnd") %in% colnames(values(GR)))) {
         ## Re-compress the original reference coordinates
         start(GR) <- values(GR)[,"refStart"];
         end(GR) <- values(GR)[,"refEnd"];
      } else {
         values(GR)[,"refStart"] <- start(GR);
         values(GR)[,"refEnd"] <- end(GR);
      }
      ranges(GR) <- IRanges(
         start=ref2compressed(start(GR)),
         end=ref2compressed(end(GR))
      );
      return(GR);
   }
   retVals <- list();
   retVals$GR <- ref2compressedGR;
   retVals$x <- ref2compressed;

   attr(retVals, "lookupCoordDF") <- lookupCoordDF;
   attr(retVals, "gapWidth") <- gapWidth;
   ## GR might be used to allow adding a feature then refreshing the function
   attr(retVals, "GR") <- GR;

   ## TODO: Convert reference coordinates to feature-based coordinates,
   ## e.g. to coordinates relative to the length of a transcript

   ## Note the format to convert any GRanges coordinates
   ## newCoordsGR <- ref2compressedGR(GR,
   ##    ref2compressed=ref2compressed);

   return(retVals);
}

#' Simplify XY coordinates to minimal line segments
#'
#' Simplify XY coordinates to minimal line segments
#'
#' This function takes a numeric matrix of x,y coordinates
#' and returns the minimal matrix of x,y coordinates that represents
#' the same line segments. It is intended in cases where there
#' is a long repeated line segment that could be represented by
#' far fewer points.
#'
#' @param xy numeric matrix with two columns representing x,y coordinates.
#' @param minN integer value to define the minimal number of repeated
#'    values before the compression is used. By definition, three consecutive
#'    points must have the same slope in order for compression to be
#'    effective, otherwise the original coordinates will be returned.
#' @param restrictDegrees numeric vector of degrees to restrict the
#'    simplification. For exmample `restrictDegrees=c(0,180)` will only
#'    simplify horizontal lines.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' xy <- cbind(
#'    x=c(1,1:15,1),
#'    y=c(0,0,0,0,1,2,3,4,5,5,5,5,6,0,0));
#' par("mfrow"=c(1,1));
#' plot(xy, cex=2);
#' points(simplifyXY(xy), pch=20, col="orange", cex=3);
#' points(simplifyXY(xy, restrictDegrees=c(0,180)), pch=17, col="purple");
#' legend("topleft", pch=c(1, 20, 17),# box.col="black", border="black",
#'    pt.cex=c(2, 3, 1),
#'    col=c("black", "orange", "purple"),
#'    legend=c("all points", "simplified points", "only simplified horizontal"));
#'
#' @export
simplifyXY <- function
(xy,
 minN=3,
 restrictDegrees=NULL,
 ...)
{
   ## Purpose is to take a matrix of x,y coordinates and return
   ## the minimal x,y coordinates that describes the same line segments.
   ## It represents only the first and last point for each consecutive
   ## line segments sharing the same angle.
   xDiff <- xy[,1] - c(head(xy[,1], 1), head(xy[,1], -1));
   yDiff <- xy[,2] - c(head(xy[,2], 1), head(xy[,2], -1));

   ## First compress simple repeated xy coordinates
   if (any(tail(xDiff, -1) == 0 & tail(yDiff, -1) == 0)) {
      repXY <- which(tail(xDiff, -1) == 0 & tail(yDiff, -1) == 0) + 1;
      xy <- xy[-repXY,,drop=FALSE];
      xDiff <- xDiff[-repXY];
      yDiff <- yDiff[-repXY];
   }

   ## Next determine the angle from point to point,
   ## two non-repeated points sharing the same angle must be
   ## on the same line, so we only need the first and last
   ## points to define the segment.
   xyAngle <- atan2(x=xDiff,
      y=yDiff);
   if (length(restrictDegrees) > 0) {
      xyAngle <- rad2deg(xyAngle);
   }
   yRle <- Rle(xyAngle);

   if (any(runLength(yRle) >= minN)) {
      rl <- runLength(yRle);
      rv <- runValue(yRle);
      rlDF <- data.frame(start=cumsum(c(1, head(rl, -1))),
         end=cumsum(rl),
         rl=rl,
         rv=rv);
      if (length(restrictDegrees) > 0) {
         ## && any(rlDF$rv %in% restrictDegrees)
         #
         yWhich <- (rlDF$rv %in% restrictDegrees & rlDF$rl > 1);
         yKeep <- unique(c(1,
            cumsum(
               rep(ifelse(yWhich, rlDF$rl, 1),
                  ifelse(yWhich, 1, rlDF$rl)))
            ));
      } else {
         yKeep <- unique(c(1, rlDF$end));
      }
      xyNew <- xy[yKeep,,drop=FALSE];
      if (1 == 2) {
         rlDF1 <- data.frame(idx=c(
            cumsum(c(1, head(rl, -1))),
            cumsum(rl)),
            rl=rep(rl, 2),
            rv=rep(rv, 2)
         )
         rlDF2 <- unique(mixedSortDF(
            data.frame(idx=c(rlDF[,1], rlDF[,2]),
               rl=rep(rl, 2),
               rv=rep(rv, 2)),
            byCols=1));
         xyNew <- xy[rlDF2$idx,,drop=FALSE];
      }
      xyNew;
   } else {
      xy;
   }
}

#' Compress genome coordinates of a matrix of polygons
#'
#' Compress genome coordinates of a matrix of polygons
#'
#' This function takes a two-column numeric matrix of polygons
#' where the x coordinate is the genomic position, and y coordinate
#' is the coverage. It uses `ref2compressed$x` to convert coordinates
#' to compressed coordinates.
#'
#' For regions that have been compresssed, it then compresses the
#' y-coordinate information to roughly one y value per integer in
#' compressed coordinate space, using the runmax across the window of
#' coverages compressed to this value.
#'
#' @export
compressPolygonM <- function
(polyM,
 ref2compressed,
 ...)
{
   ## Purpose is to compress coordinates in coverage polygons
   ## polyM <- covPolyML[[1]]
   polyMlengths <- diff(unique(c(0, which(is.na(polyM[,1])), nrow(polyM))));
   polyDF <- as.data.frame(polyM);
   polyDF$newX <- ref2compressed$x(polyDF$x);
   polyDF$ratio <- c(NA, diff(polyDF$newX) / diff(polyDF$x));
   polyMrepN <- rep(seq_along(polyMlengths), polyMlengths);
   polyDF$n <- polyMrepN;
   iMedianRatios <- shrinkMatrix(polyDF$ratio,
      groupBy=polyDF$n,
      shrinkFunc=function(x){median(x, na.rm=TRUE)})$x;
   polyDF$medianRatio <- rep(iMedianRatios,
      polyMlengths);
   polyDFL <- split(polyDF, polyMrepN);
   whichComp <- which(iMedianRatios < 0.5);
   whichNorm <- which(iMedianRatios >= 0.5);
   polyDFLnew <- lapply(nameVector(whichComp), function(k){
      iDF <- polyDFL[[k]];
      baseline <- iDF[1,"y"];
      iDF <- iDF[!is.na(iDF[,1]),,drop=FALSE];
      iDFu <- iDF[match(unique(iDF$x), iDF$x),,drop=FALSE];

      iRange <- range(iDF$x, na.rm=TRUE);
      iRangeNew <- range(iDF$newX, na.rm=TRUE);
      iN <- ceiling(diff(iRange)*head(iDF$medianRatio, 1));
      iMedRatio <- head(iDF$medianRatio, 1);
      iMultiple <- ceiling(1/iMedRatio);
      iSeqNew <- seq(from=iRangeNew[1],
         to=iRangeNew[2],
         by=iMedRatio);
      iSeqSub <- seq(from=iRangeNew[1],
         to=iRangeNew[2],
         length.out=iN);
      iSeqSub[iSeqSub > max(iSeqNew)] <- max(iSeqNew);

      ## use approx() to fill in holes
      iYnew <- approx(x=iDFu$newX, y=iDFu$y, xout=iSeqNew)$y;
      iYrunmax <- caTools::runmax(iYnew, k=iMultiple, endrule="keep");
      iYsub <- approx(x=iSeqNew,
         y=iYrunmax,
         xout=iSeqSub)$y;
      if (1 == 2) {
         plot(x=iDFu$newX, y=iDFu$y, cex=2);
         polygon(x=c(iDFu$newX[1], iDFu$newX, tail(iDFu$newX, 1)),
            y=c(0,iDFu$y, 0), col="grey85");
         points(x=iSeqNew, y=iYnew, pch=20, cex=2, col="red");
         points(x=iSeqNew, y=iYrunmax, pch=17, cex=1, col="orange");
         points(x=iSeqSub, y=iYsub, pch=17, cex=2, col="purple");
      }
      if (head(iYsub, 1) != baseline) {
         iYsub <- c(baseline, iYsub);
         iSeqSub <- c(head(iSeqSub, 1), iSeqSub);
      }
      if (tail(iYsub, 1) != baseline) {
         iYsub <- c(iYsub, baseline);
         iSeqSub <- c(iSeqSub, tail(iSeqSub, 1));
      }
      iM <- cbind(x=c(iSeqSub, NA),
         y=c(iYsub, NA));
   });
   newPolyDFL <- list();
   newPolyDFL[names(polyDFL)[whichNorm]] <- lapply(nameVector(whichNorm), function(k){
      as.matrix(renameColumn(polyDFL[[k]][,c("newX","y")],
         from="newX", to="x"));
   });
   newPolyDFL[names(polyDFLnew)] <- polyDFLnew;
   newPolyDFL <- newPolyDFL[mixedSort(names(newPolyDFL))];
   newPolyM <- rbindList(newPolyDFL);
   #plot(newPolyM, pch=".");
   #polygon(newPolyM, col=rainbowJam(46));
   return(newPolyM);
}

#' Convert exon coverage to polygons
#'
#' Convert exon coverage to polygons
#'
#' @export
exoncov2polygon <- function
(gr,
 covNames=NULL,
 baseline=0,
 gapWidth=250,
 doPlot=FALSE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to take exon coverage in the form of NumericList
   ## associated with GRanges, and produce a list of polygons
   ## with baseline=0
   ##
   ## Workflow:
   ## - input GRanges with coverages in the values columns
   ##   stored as NumericList
   ## create polygons

   if (!igrepHas("granges", class(gr))) {
      stop("Input gr should be a GRanges object.");
   }
   if (length(covNames) == 0) {
      covNames <- vigrep("pos|neg", colnames(values(gr)));
   }
   if (length(covNames) == 0) {
      stop("covNames must be colnames(values(gr)).");
   }

   ## Define compressed coordinate space
   ref2c <- make_ref2compressed(
      subset(gr, feature_type %in% "exon"),
      gapWidth=gapWidth);

   ## Extend baseline to length of gr, so the baseline applies
   ## to each exon
   baseline <- rep(baseline, length.out=length(gr));

   ## List of lists of polygons
   covPolyL <- lapply(nameVector(covNames), function(iName){
      polyL <- lapply(names(gr), function(iGRname){
         yVals1 <- unlist(values(gr[iGRname])[[iName]]) + baseline;
         ## Note: we define xVals1 width using yVals1 width,
         ## because this GRanges might have compressed coordinates
         ## and therefore the width(gr) is not an accurate measure
         ## of the actual width of coverage data
         xVals1 <- seq(from=start(gr[iGRname]),
            to=end(gr[iGRname]),
            length.out=length(yVals1));
         xVals <- c(rep(head(xVals1, 1), 2) - 0.5,
            xVals1,
            rep(tail(xVals1, 1), 2) + 0.5,
            head(xVals1, 1) - 0.5);
         yVals <- c(0,
            head(yVals1, 1),
            yVals1,
            tail(yVals1, 1),
            0,
            0);
         ## TODO: compress coordinates where y-value doesn't change
         xy <- cbind(x=xVals, y=yVals);
         ## Simplify the coordinates, typically removing up to 90% rows
         xy <- simplifyXY(xy,
            restrictDegrees=c(0,180));
      });
   });
   ## List of coordinate matrices suitable for vectorized polygon()
   covPolyML <- lapply(covPolyL, function(iL){
      tail(rbindList(lapply(iL, function(iM){
         rbind(cbind(x=NA, y=NA),
            iM)
      })), -1);
   });

   ## Compress polygon
   compCovPolyML <- lapply(covPolyML, function(iL){
      iML <- compressPolygonM(iL, ref2compressed=ref2c);
   });

   ## Optionally plot coverage
   if (doPlot) {
      par("mfrow"=c(2,1));
      polyCol <- colorjam::rainbowJam(length(covPolyL[[1]]));
      plot(rbindList(covPolyML),
         pch=".",
         xaxt="n",
         col="transparent");
      polygon(rbindList(covPolyML),
         border=polyCol,
         col=polyCol);
      plot(rbindList(compCovPolyML),
         pch=".",
         xaxt="n",
         col="transparent");
      polygon(compCovPolyML[[1]],
         border=polyCol,
         col=polyCol);
      polygon(compCovPolyML[[2]],
         border=polyCol,
         col=polyCol);
      par("mfrow"=c(1,1));
   }
   retVals <- list(covPolyML=covPolyML,
      compCovPolyML=compCovPolyML);
   return(retVals);
}
