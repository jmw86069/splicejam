

#' Find closest exon to splice junction ends
#'
#' Find closest exon to splice junction ends
#'
#' This function is used to annotate splice junction GRanges entries
#' based upon the closest compatible stranded exon boundary.
#'
#' This function should usually be called by `spliceGR2junctionDF()`
#' and not called directly.
#'
#' @param spliceGRgene GRanges representing splice junctions
#' @param exonsGR GRanges representing flattened exons per gene, as
#'    is produced by `flattenExonsByGene()`.
#' @param flipNegativeStrand logical indicating whether to flip
#'    the orientation of features on the negative strand.
#' @param sampleColname character value matching the colname that
#'    defines distinct sample identifier, for which junctions will
#'    be kept separate.
#' @param reportActualCoords logical indicating whether to report
#'    genomic coordinates or transcriptome coordinates. (Work in
#'    progress.)
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
closestExonToJunctions <- function
(spliceGRgene,
 exonsGR,
 flipNegativeStrand=TRUE,
 sampleColname="sample_id",
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
      if (verbose) {
         printDebug("closestExonToJunctions(): ",
            "switchCols1:",
            switchCols1);
         printDebug("closestExonToJunctions(): ",
            "switchCols2:",
            switchCols2);
      }
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
#' This function takes a GRanges object representing multiple splice
#' junction ranges, with associated scores, and returns a data.frame
#' summary of junctions with annotated boundaries using a set
#' of gene exon models. Junctions whose ends are within `spliceBuffer`
#' distance are combined, and the scores are summed.
#'
#' By default, junctions not within `spliceBuffer` of a compatible exon
#' boundary are named by the nearest exon boundary, and the distance
#' upstream or downstream from the boundary.
#'
#' Multiple samples can be processed together, and the results will
#' be aggregated within each sample, using `sampleColname`. The results
#' in that case may be cast to wide format using `nameFromTo` as the
#' row identifier, `score` as the value column, and `sampleColname` as
#' the new column headers.
#'
#' @param spliceGRgene GRanges object containing splice junctions, where
#'    the `scoreColname` contains numeric scores.
#' @param exonsGR GRanges object containing flattened exons by gene,
#'    as is provided by `flattenExonsByGene()`.
#' @param spliceBuffer integer distance allowed from a compatible exon
#'    boundary, for a junction read to be snapped to that boundary.
#' @param useOnlyValidEntries logical indicating whether to remove
#'    junctions that do not align with a compatible exon boundary.
#' @param renameTooFar logical indicating whether junctions are
#'    named by the nearest exon boundary and the distance to that
#'    boundary.
#' @param scoreColname,sampleColname colnames in `values(spliceGRgene)`
#'    to define the score, and `sample_id`.
#' @param flipNegativeStrand logical indicating whether to flip the
#'    orientation of negative strand features when matching exon
#'    boundaries. This argument is passed to `closestExonToJunctions()`.
#' @param returnGRanges logical indicating whether to return GRanges,
#'    or by default, `data.frame`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
spliceGR2junctionDF <- function
(spliceGRgene,
 exonsGR,
 spliceBuffer=3,
 geneExonSep="(:|_exon)",
 useOnlyValidEntries=FALSE,
 renameTooFar=TRUE,
 scoreColname="score",
 sampleColname="sample_id",
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
   sampleColname <- intersect(sampleColname, colnames(values(spliceGRgene)));

   ## Vectorized logic:
   ## - find closest start/end for each splice start/end
   ## - subset for splice start/end having same gene
   ##
   ## Call a method which encapsulates distanceToNearest(), calculates stranded distance,
   ## and knows to match splice start with exon end, etc.
   #spliceGRgeneVals <- closestExonToJunctions(spliceGRgene=spliceGRgene[,scoreColname], exonsGR=exonsGR,
   #   flipNegativeStrand=flipNegativeStrand, ...);
   spliceGRgene <- closestExonToJunctions(spliceGRgene=spliceGRgene[,c(scoreColname,sampleColname)],
      exonsGR=exonsGR,
      flipNegativeStrand=flipNegativeStrand,
      sampleColname=sampleColname,
      verbose=verbose)$spliceGRgene;
      #...)$spliceGRgene;
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

   spliceColnames <- c("seqnames", "start", "end", "nameFrom", "nameTo",
      scoreColname, "strand", sampleColname);
   spliceCountsDF <- as.data.frame(spliceGRgene[toUse])[,spliceColnames,drop=FALSE];
   spliceCountsDF[,"nameFromTo"] <- pasteByRow(spliceCountsDF[,c("nameFrom","nameTo"),drop=FALSE],
      sep=" ",
      na.rm=TRUE);
   spliceCountsDF[,"nameFromToSample"] <- pasteByRow(
      spliceCountsDF[,c("nameFromTo", sampleColname),drop=FALSE],
      sep=":!:");

   if (verbose) {
      printDebug("spliceGR2junctionDF(): ",
         "Shrinking matrix to combine counts spanning the same junctions.");
   }
   spliceCountsDFshrunk <- renameColumn(
      shrinkMatrix(spliceCountsDF[,scoreColname,drop=FALSE],
         groupBy=spliceCountsDF[,"nameFromToSample"],
         shrinkFunc=sum),
      from="groupBy",
      to="nameFromToSample");
   if (length(sampleColname) > 0) {
      spliceCountsDFshrunk[,c("nameFromTo",sampleColname)] <- rbindList(
         strsplit(spliceCountsDFshrunk[,"nameFromToSample"], ":!:"));
   } else {
      spliceCountsDFshrunk[,"nameFromTo"] <- spliceCountsDFshrunk[,"nameFromToSample"];
   }
   spliceCountsDFshrunk[,c("nameFrom", "nameTo")] <- rbindList(
      strsplit(spliceCountsDFshrunk[,"nameFromTo"], " "));
   #spliceCountsDFshrunk <- cbind(spliceCountsDFshrunk,
   #   rbindList(strsplit(spliceCountsDFshrunk[,"groupBy"], " ")));
   #colnames(spliceCountsDFshrunk) <- c("nameFromTo", spliceGRname, "nameFrom", "nameTo");

   ## Now add the start and end coordinates, so the results can be plotted
   if (verbose) {
      printDebug("spliceGR2junctionDF(): ",
         "Adding ref, strand, start, end.");
   }
   matchFromTo <- match(spliceCountsDFshrunk[,"nameFromToSample"],
      spliceCountsDF[,"nameFromToSample"]);
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

   ## Create junctionID only after we decided how to orient exonFrom and exonTo
   ## for negative strand entries, avoiding re-creating the name for those rows
   spliceCountsDFshrunk[,"junctionID"] <- pasteByRow(
      spliceCountsDFshrunk[,c("exonFrom","exonTo",sampleColname),drop=FALSE],
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
#' This function is deprecated, in factor of `make_ref2compressed()`.
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
#' @family jam GRanges functions
#' @family jam RNA-seq functions
#'
#' @return list with `trans_grc` which is class `"trans"` suitable
#'    for use in ggplot2 functions; `transform` a function that converts
#'    chromosome coordinates to compressed coordinates; `inverse` a function
#'    that converts compressed coordinates to chromosome coordinates;
#'    `scale_x_grc` a function used similar to `ggplot2::scale_x_continuous()`
#'    during ggplot2 creation; `gr` a function that compresses coordinates
#'    in a `GRanges` object; `grl` a function that compresses coordinates in
#'    a `GRangesList` object. Attributes `"lookupCoordDF"` is a two-column
#'    data.frame with chromosome coordinates and compressed coordinates,
#'    which is used to create the other transformation functions via
#'    `stats::approx()`; `"gapWidth"` the gap width used, since it can
#'    be programmatically defined; `"gr"` the `GRanges` input data used
#'    to train the transformation.
#'
#' @family jam GRanges functions
#'
#' @param gr GRanges object containing regions not to compress. Regions
#'    which are unstranded gaps are compressed to fixed width.
#' @param gapWidth integer value used for fixed gap width, or when
#'    NULL the gap width is defined as 3 times the median feature width.
#' @param keepValues logical indicating whether to keep feature values
#'    in the GRanges data.
#' @param upstream,downstream,upstreamGapWidth,downstreamGapWidth used
#'    to define the compression of coordinates upstream and downstream
#'    the supplied GRanges. In reality, the upstream range and upstream
#'    gap width defines a multiplier, and all upstream coordinates are
#'    compressed through zero. Similarly, all downstream coordinates
#'    are compressed to 10 billion, which is roughly 3 times the size
#'    of the human genome.
#' @param nBreaks the default number of x-axis coordinate breaks used
#'    in ggplot labeling.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @seealso `grl2df()`
#'
#' @export
make_ref2compressed <- function
(gr,
 gapWidth=200,
 keepValues=FALSE,
 upstream=50000,
 upstreamGapWidth=gapWidth*3,
 downstream=50000,
 downstreamGapWidth=gapWidth*3,
 nBreaks=7,
 verbose=FALSE,
...)
{
   ## Purpose is to use a GRanges object to train a coordinate
   ## compression function, which assigns the gaps to a fixed
   ## width, but allows GRanges features to keep their original
   ## width.
   if (length(gapWidth) == 0) {
      gapWidth <- max(c(10,
         round(median(width(reduce(gr))) / 2)));
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
   #disGR <- disjoin(gr);
   disGR <- reduce(gr);
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
   newWidths <- suppressWarnings(
      head(intercalate(width(disGR),
         (gapWidth)*exonGaps),
         length(disGR)*2-1)
      );

   ##
   if (verbose) {
      printDebug("compressGRgaps(): ",
         "newCoords");
   }
   newCoords <- cumsum(newWidths);
   newCoordsM <- matrix(c(0, newCoords), ncol=2, byrow=TRUE);
   newCoordsM[,1] <- newCoordsM[,1] + 1;

   lookupCoordDF <- jamba::mixedSortDF(unique(data.frame(
      refCoord=c(start(disGR),end(disGR)),
      coord=c(newCoordsM[,1],newCoordsM[,2])
   )));
   if (upstream > 0) {
      ## Add one upstream point
      refCoord1 <- head(lookupCoordDF$refCoord,1);
      coord1 <- head(lookupCoordDF$coord,1);
      upExtended1 <- (refCoord1 - upstream);
      upComp1 <- (coord1 - upstreamGapWidth);
      lookupCoordDF <- rbind(
         data.frame(refCoord=upExtended1,
            coord=upComp1),
         lookupCoordDF);

      ## Extend the upstream gap to 1
      minRefCoord <- 1;
      if (upExtended1 > minRefCoord) {
         upExtended1 <- minRefCoord;
         newRefDiff <- (refCoord1 - minRefCoord);
         upComp1 <- coord1 - (upstreamGapWidth * newRefDiff) / upstream;
         lookupCoordDF <- rbind(
            data.frame(refCoord=upExtended1,
               coord=upComp1),
            lookupCoordDF);
      }
   }
   if (downstream > 0) {
      ## Add one downstream point
      refCoord2 <- tail(lookupCoordDF$refCoord,1);
      coord2 <- tail(lookupCoordDF$coord,1);
      downExtended2 <- (refCoord2 + downstream);
      downComp2 <- (coord2 + downstreamGapWidth);
      lookupCoordDF <- rbind(
         lookupCoordDF,
         data.frame(refCoord=downExtended2,
            coord=downComp2)
      );
      ## Extend the upstream gap to 10 Gb
      maxRefCoord <- 3e10;
      if (downExtended2 < maxRefCoord) {
         downExtended2 <- maxRefCoord;
         newRefDiff <- (maxRefCoord - refCoord2);
         downComp2 <- coord2 + (downstreamGapWidth * newRefDiff) / downstream;
         lookupCoordDF <- rbind(
            lookupCoordDF,
            data.frame(refCoord=downExtended2,
               coord=downComp2)
         );
      }
   }
   ## TODO: expand the range of coordinates so approxfun() will not fail
   if (verbose) {
      printDebug("compressGRgaps(): ",
         "approxfun");
   }
   ref2compressed <- approxfun(x=lookupCoordDF[,1],
      y=lookupCoordDF[,2],
      method="linear");
   compressed2ref <- approxfun(x=lookupCoordDF[,2],
      y=lookupCoordDF[,1],
      method="linear");

   ## Add a convenience function for GRanges objects
   if (verbose) {
      printDebug("compressGRgaps(): ",
         "ref2compressedGR");
   }
   ref2compressedGR <- function(gr) {
      if (length(names(gr)) > 0) {
         GRnames <- names(gr);
      }
      if (all(c("refStart","refEnd") %in% colnames(values(GR)))) {
         ## Re-compress the original reference coordinates
         ranges(gr) <- IRanges(
            start=values(gr)[,"refStart"],
            end=values(gr)[,"refEnd"]
         );
      } else {
         values(gr)[,"refStart"] <- start(gr);
         values(gr)[,"refEnd"] <- end(gr);
      }
      ranges(gr) <- IRanges(
         start=ref2compressed(start(gr)),
         end=ref2compressed(end(gr))
      );
      if (length(GRnames) > 0) {
         names(gr) <- GRnames;
      }
      return(gr);
   }
   ## Add a convenience function for GRangesList objects
   ref2compressedGRL <- function(grl) {
      if (length(names(grl)) == 0) {
         names(grl) <- makeNames(rep("grl", length(grl)));
      }
      grlc1 <- ref2compressedGR(grl@unlistData);
      grlc <- GenomicRanges::split(grlc1,
         factor(
            rep(names(grl), elementNROWS(grl)),
            levels=names(grl)));
      grlc;
   }
   if (verbose) {
      printDebug("compressGRgaps(): ",
         "retVals");
   }

   ## Custom breaks function
   breaks_gr <- function(x, limits=NULL, n=nBreaks, xFixed=lookupCoordDF[,1], verbose=FALSE, ...) {
      xvals <- unique(sort(xFixed));
      xvals <- xvals[xvals >= min(x) & xvals <= max(x)];
      if (verbose) {
         printDebug("breaks_gr(): ",
            "x:", x);
         printDebug("breaks_gr(): ",
            "limits:", limits);
         printDebug("breaks_gr(): ",
            "n:", n);
      }
      if (n > length(xvals)) {
         return(xvals);
      }
      idx1 <- round(seq.int(from=1, to=length(xvals), length.out=n));
      return(xvals[idx1]);
   }
   minor_breaks <- function(b, limits=NULL, n=2, xFixed=lookupCoordDF[,1], verbose=FALSE, ...) {
      if (verbose) {
         printDebug("minor_breaks(): ",
            "x:", x);
         printDebug("minor_breaks(): ",
            "limits:", limits);
         printDebug("minor_breaks(): ",
            "n:", n);
      }
      nUse <- length(b) * n;
      ref2compressed(breaks_gr(compressed2ref(b),
         limits=compressed2ref(limits),
         n=nUse,
         xFixed=xFixed,
         verbose=verbose));
   }
   labels <- function(x, n=10) {
      scales::comma_format()(breaks_gr(x, n=n))
   }
   retVals <- list();

   ## Make the custom trans function
   trans_grc <- scales::trans_new(name="compressed_gr",
      transform=ref2compressed,
      inverse=compressed2ref,
      breaks=breaks_gr,
      minor_breaks=minor_breaks,
      format=scales::comma_format(),
      domain=range(lookupCoordDF[,1]));

   ## Make the actual scale_x_gr_compressed() function
   scale_x_grc <- function(..., trans=trans_grc){
      scale_x_continuous(name="grc",
         breaks=waiver(),#trans_grc$breaks,
         minor_breaks=waiver(),#trans_grc$minor_breaks,
         labels=function(...){scales::comma(...)},#trans_grc$labels,
         trans=trans,
         ...)
   }

   ## Functions required by scales::trans_new()
   retVals$transform <- ref2compressed;
   retVals$inverse <- compressed2ref;

   ## Convenience functions for GenomicRanges objects
   retVals$gr <- ref2compressedGR;
   retVals$grl <- ref2compressedGRL;

   #retVals$breaks <- breaks_gr;
   #retVals$minor_breaks <- minor_breaks;
   #retVals$format <- scales::comma_format;
   #retVals$labels <- labels;

   ## ggplot2 functions to compress the x-axis
   retVals$scale_x_grc <- scale_x_grc;
   retVals$trans_grc <- trans_grc;

   ## TODO: custom breaks() function that prioritizes choosing labels
   ## at borders of the compressed ranges, then sensible labels between them


   attr(retVals, "lookupCoordDF") <- lookupCoordDF;
   attr(retVals, "gapWidth") <- gapWidth;
   ## GR might be used to allow adding a feature then refreshing the function
   attr(retVals, "gr") <- gr;

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
#' @family jam spatial functions
#'
#' @param xy numeric matrix of two columns x,y
#' @param minN integer minimum number of consecutive points to
#'    cause compression to occur. It requires at least 3 points
#'    inherent to the algorithm.
#' @param restrictDegrees numeric vector of allowed angles in degrees
#'    (range 0 to 360) of angles allowed to be compressed, when supplied
#'    all other angles will not be compressed.
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
#' is the coverage. It uses `ref2compressed$transform` to convert
#' coordinates to compressed coordinates, as output from
#' `make_ref2compressed()`.
#'
#' For regions that have been compresssed, it then compresses the
#' y-coordinate information to roughly one y value per integer in
#' compressed coordinate space, using the runmax across the window of
#' coverages compressed to this value.
#'
#' @family jam spatial functions
#' @family jam RNA-seq functions
#'
#' @return data.frame with the same colnames, with reduced rows
#'    for polygons where the coordinate compression defined in
#'    `ref2c` is above `minRatio` ratio. The compressed coverage roughly
#'    represents the max coverage value for each point.
#'
#' @param polyM data.frame including `x`,`y` coordinates of polygons,
#'    `cov` and `gr` to group each polygon.
#' @param ref2c list output from `make_ref2compressed()`.
#' @param minRatio the minimum ratio of coordinate compression required
#'    before polygon resolution is downsampled.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
compressPolygonM <- function
(polyM,
 ref2c,
 minRatio=3,
 verbose=FALSE,
 ...)
{
   ## Purpose is to compress coordinates in coverage polygons
   ##
   ## General strategy is to determine the relative scaling of
   ## each polygon, and downsample polygons proportional to the
   ## amount of compression on the x-axis.
   if (any(is.na(polyM[,1]))) {
      data_style <- "base";
      polyMlengths <- diff(unique(c(0, which(is.na(polyM[,1])), nrow(polyM))));
      polyMratios <- shrinkMatrix(polyDF$ratio,
         groupBy=polyDF$n,
         shrinkFunc=function(x){median(x, na.rm=TRUE)})$x;
   } else {
      data_style <- "fortify";
      idrows <- pasteByRow(polyM[,c("cov","gr")]);
      idrows <- factor(idrows, levels=unique(idrows));
      polyML <- split(polyM, idrows);
      polyMratios <- unlist(lapply(polyML, function(i){
         xr1 <- range(i[,1]);
         xc1 <- ref2c$transform(xr1);
         xv1 <- rev(sort(c(diff(xr1)+1, diff(xc1)+1)));
         xv1[1] / xv1[2];
      }));
      polyMlengths <- sdim(polyML)[,1];
   }
   if (verbose) {
      printDebug("compressPolygonM(): ",
         "polyMratios:", head(polyMratios, 20));
   }

   ## data.frame describing the compression and polygon
   polyDF <- as.data.frame(polyM);
   polyDF$newX <- ref2c$transform(polyDF$x);
   polyDF$ratio <- c(NA, diff(polyDF$newX) / diff(polyDF$x));
   polyMrepN <- rep(seq_along(polyMlengths), polyMlengths);
   polyDF$n <- polyMrepN;

   polyDF$medianRatio <- rep(polyMratios,
      polyMlengths);

   polyDFL <- split(polyDF, polyMrepN);
   whichComp <- which(polyMratios >= minRatio);
   whichNorm <- which(polyMratios < minRatio);
   if (length(whichComp) == 0) {
      return(polyM);
   }
   ## Compress those polygons whose coordinates get compressed
   polyDFLnew <- lapply(nameVector(whichComp), function(k){
      iDF <- polyDFL[[k]];
      baseline <- iDF[1,"y"];
      iDF <- iDF[!is.na(iDF[,1]),,drop=FALSE];
      iDFu <- iDF[match(unique(iDF$x), iDF$x),,drop=FALSE];

      iRange <- range(iDF$x, na.rm=TRUE);
      iRangeNew <- range(iDF$newX, na.rm=TRUE);
      iMedRatio <- head(iDF$medianRatio, 1);
      iN <- ceiling(diff(iRange)/iMedRatio);
      iMultiple <- ceiling(iMedRatio);
      iSeqNew <- seq(from=iRange[1],
         to=iRange[2],
         by=iMedRatio);
      iSeqSub <- seq(from=iRange[1],
         to=iRange[2],
         length.out=iN);
      #iSeqNew <- seq(from=iRangeNew[1],
      #   to=iRangeNew[2],
      #   by=iMedRatio);
      #iSeqSub <- seq(from=iRangeNew[1],
      #   to=iRangeNew[2],
      #   length.out=iN);
      iSeqSub[iSeqSub > max(iSeqNew)] <- max(iSeqNew);

      ## use approx() to fill in holes
      #iYnew <- approx(x=iDFu$newX,
      iYnew <- approx(x=iDFu$x,
         y=iDFu$y,
         xout=iSeqNew)$y;
      ## Take running max value, then use approx
      iYrunmax <- caTools::runmax(iYnew,
         k=iMultiple,
         endrule="keep");
      if (1 == 2) {
         printDebug("head(iDFu):");print(head(iDFu));
         printDebug("head(iYnew):", head(iYnew));
         printDebug("head(iMultiple):", head(iMultiple));
         printDebug("head(iN):", head(iN));
         printDebug("head(iSeqNew):", head(iSeqNew));
         printDebug("head(iYrunmax):", head(iYrunmax));
         printDebug("head(iSeqSub):", head(iSeqSub));
      }
      iYsub <- approx(x=iSeqNew,
         y=iYrunmax,
         xout=iSeqSub)$y;
      if (head(iYsub, 1) != baseline) {
         iYsub <- c(baseline, iYsub);
         iSeqSub <- c(head(iSeqSub, 1), iSeqSub);
      }
      if (tail(iYsub, 1) != baseline) {
         iYsub <- c(iYsub, baseline);
         iSeqSub <- c(iSeqSub, tail(iSeqSub, 1));
      }
      iM <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         x=iSeqSub,
         y=iYsub,
         gr=head(iDF$gr,1),
         cov=head(iDF$cov, 1));
   });
   newPolyDFL <- list();
   newPolyDFL[names(polyDFL)[whichNorm]] <- lapply(nameVector(whichNorm), function(k){
      polyDFL[[k]][,c("x","y","cov","gr"),drop=FALSE];
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
#' This function is a workhorse function that converts a GRanges
#' object containing column values with NumericList coverage data,
#' into a full data.frame sufficient to define ggplot2 and other
#' coverage polygon plots.
#'
#' An interesting argument is `baseline` which allows each exon
#' in the `gr` GRanges object to be offset from zero, in order to
#' make certain features visually easier to distinguish.
#'
#' This function also calls `simplifyXY()` which reduces the
#' stored polygon detail for regions whose coordinates are compressed
#' on the x-axis, taking roughly the max value for each point.
#'
#' The default output is roughly similar to `broom::tidy()` in
#' that it converts a custom R object into a tidy data.frame
#' suitable for use by ggplot2 and other tidy workflows.
#'
#' The function `getGRcoverageFromBw()` takes a set of bigWig files
#' and returns a GRanges object whose columns contain NumericList data,
#' which is the intended input for `exoncov2polygon()`.
#'
#' @family jam GRanges functions
#' @famly jam RNA-seq functions
#'
#' @param gr GRanges where `colnames(values(gr))` is present in `covNames`,
#'    and contains data with class `NumericList`.
#' @param covNames character vector contained in `colnames(values(gr))`.
#' @param baseline numeric vector of length 0, 1 or `length(gr)`. If
#'    `baseline` has names matching `names(gr)` they will be used for
#'    each `gr` entry; if `baseline` is not named, it is extended
#'    to `length(gr)`. The `baseline` value is added to the coverage
#'    for each exon to offset the polygon as needed.
#' @param gapWidth numeric value sent to `make_ref2compressed()`.
#' @param coord_style character value to define the output style:
#'    `"base"` returns a matrix with polygons separated by a row of `NA`
#'    values; `"fortify"` returns a `data.frame` intended for ggplot2,
#'    with columns `"cov"` and `"gr"` indicating the values in `covNames` and
#'    `names(gr)` used to separate each polygon.
#' @param ref2c optional list containing output from `make_ref2compressed()`,
#'    used to compress the GRanges coordinates.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
exoncov2polygon <- function
(gr,
 covNames=NULL,
 sample_id=NULL,
 baseline=NULL,
 gapWidth=250,
 doPlot=FALSE,
 coord_style=c("fortify", "base", "list", "all"),
 ref2c=NULL,
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
   coord_style <- match.arg(coord_style);

   if (!igrepHas("granges", class(gr))) {
      stop("Input gr should be a GRanges object.");
   }
   if (length(covNames) == 0) {
      covNames <- provigrep(c("pos|[+]$", "neg|[-]$"),
         colnames(values(gr)));
   }
   if (length(covNames) == 0) {
      stop("covNames must be colnames(values(gr)).");
   }
   retVals <- list();

   ## Define compressed coordinate space
   #ref2c <- make_ref2compressed(
   #   subset(gr, feature_type %in% "exon"),
   #   gapWidth=gapWidth);

   ## Extend baseline to length of gr, so the baseline applies
   ## to each exon
   baselineV <- nameVector(
      rep(0,
         length.out=length(gr)),
      names(gr));
   if (length(baseline) > 0) {
      if (length(names(baseline)) > 0) {
         baselineV[names(baseline)] <- baseline;
      } else {
         baselineV[] <- baseline;
      }
   }

   ## List of lists of polygons
   if (verbose) {
      printDebug("exoncov2polygon(): ",
         "covNames:",
         covNames);
   }
   covPolyL <- lapply(nameVector(covNames), function(iName){
      polyL <- lapply(nameVectorN(gr), function(iGRname){
         yVals1a <- unlist(values(gr[iGRname])[[iName]]);
         iBase <- baselineV[iGRname];
         yVals1 <- yVals1a + iBase;
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
         yVals <- c(iBase,
            head(yVals1, 1),
            yVals1,
            tail(yVals1, 1),
            iBase,
            iBase);
         if (length(xVals) != length(yVals)) {
            printDebug("length(xVals):", length(xVals));
            printDebug("length(yVals):", length(yVals));
         }
         ## TODO: compress coordinates where y-value doesn't change
         if (length(xVals) > 0 && length(yVals) > 0) {
            xy <- cbind(x=xVals, y=yVals);
            ## Simplify the coordinates, typically removing up to 90% rows
            xy <- simplifyXY(xy,
               restrictDegrees=c(0,180));
         } else {
            xy <- cbind(x=0, y=0)[0,,drop=FALSE];
         }
      });
   });
   if ("list" %in% coord_style) {
      return(covPolyL);
   }
   retVals$covPolyL <- covPolyL;

   if ("base" %in% coord_style) {
      ## coordinate matrix suitable for vectorized polygon() by separating
      ## each polygon with a row of NA values
      covPolyML <- lapply(covPolyL, function(iL){
         tail(rbindList(lapply(iL, function(iM){
            rbind(cbind(x=NA, y=NA),
               iM)
         })), -1);
      });
      return(covPolyML);
   }
   if ("fortify" %in% coord_style) {
      covPolyML <- lapply(nameVectorN(covPolyL), function(iN){
         iL <- covPolyL[[iN]];
         rbindList(lapply(nameVectorN(iL), function(iN2){
            iM <- iL[[iN2]];
            data.frame(check.names=FALSE,
               stringsAsFactors=FALSE,
               iM,
               cov=iN,
               gr=iN2);
         }));
      });
      retVals$covPolyML <- covPolyML;

      ## Compress polygon
      if (length(ref2c) > 0) {
         compCovPolyML <- lapply(covPolyML, function(iL){
            iML <- compressPolygonM(iL,
               ref2c=ref2c,
               coord_style=coord_style);
         });
         retVals$compCovPolyML <- compCovPolyML;
         #return(compCovPolyML);
         covPolyML <- compCovPolyML;
      }

      covPolyDF <- rbindList(covPolyML);
      covPolyDF$cov <- factor(covPolyDF$cov,
         levels=unique(c(covNames, covPolyDF$cov)));
      covPolyDF$gr <- factor(covPolyDF$gr,
         levels=unique(c(names(gr),
            covPolyDF$gr)));
      ## Optionally add sample_id
      if (length(sample_id) > 0) {
         cov2sample <- nameVector(sample_id, covNames);
         covPolyDF$sample_id <- cov2sample[covPolyDF$cov];
      }
      return(covPolyDF);
   }

      ## Compress polygon
      if (length(ref2c) > 0) {
         compCovPolyML <- lapply(covPolyML, function(iL){
            iML <- compressPolygonM(iL,
               ref2compressed=ref2c,
               coord_style=coord_style);
         });
         retVals$compCovPolyML <- compCovPolyML;
      }

   ## Optionally plot coverage
   if (doPlot) {
      if ("base" %in% coord_style) {
         if (exists(compCovPolyML)) {
            par("mfrow"=c(2,1));
         }
         polyCol <- colorjam::rainbowJam(length(covPolyL[[1]]));
         plot(rbindList(covPolyML),
            pch=".",
            xaxt="n",
            col="transparent");
         polygon(rbindList(covPolyML),
            border=polyCol,
            col=polyCol);
         if (exists(compCovPolyML)) {
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
      }
   }
   return(retVals);
}

#' Get coverage for GRanges from bigWig files
#'
#' Get coverage for GRanges from bigWig files
#'
#' This function takes a GRanges object to define genomic regions,
#' for which coverage data is loaded for each `bwUrls` input file.
#'
#' @return DataFrame object, whose colnames are defined using
#'    `makeNames(basename(bwUrls))`. Each column is type `NumericList`,
#'    which is a list of numeric coverage values.
#'
#' @family jam GRanges functions
#' @family RNA-seq functions
#'
#' @param gr GRanges object
#' @param bwUrls character vector of full file paths or web URLs
#'    to bigWig files, suitable for use by `rtracklayer::import()`.
#' @param addGaps logical indicating whether gaps between GRanges
#'    should be added to the query. Gaps are determined using
#'    `getGRgaps()`.
#' @param feature_type_colname,gap_feature_type,default_feature_type
#'    When `addGaps=TRUE` a
#'    new column named using `feature_type_colname` is added to `values(gr)`,
#'    whose value for gap regions is `gap_feature_type`. When
#'    `feature_type_colname` is already present in `gr` it is not modified,
#'    otherwise the column is created with value `default_feature_type`.
#'    By default, this function adds a column `"feature_type"` with
#'    value `"gap"`.
#' @param ... additional arguments are ignored.
#'
#' @export
getGRcoverageFromBw <- function
(gr,
 bwUrls,
 addGaps=FALSE,
 gap_feature_type="gap",
 default_feature_type="exon",
 feature_type_colname="feature_type",
 verbose=FALSE,
 ...)
{
   ## Purpose is to get coverage from bigWig files,
   ## returning GRanges with columns containing NumericList
   ## get coverage across the regions of interest.
   ##
   ## TODO: consider splitting bwUrls into bwUrlsPos, bwUrlsNeg
   ## in order to allow strand-specificity
   if (!igrepHas("GRanges", class(gr))) {
      stop("gr must be a GRanges object.");
   }
   if (length(names(bwUrls)) == 0) {
      names(bwUrls) <- makeNames(
         gsub("[.](bw|bigWig)$",
            "",
            ignore.case=TRUE,
            basename(bwUrls)));
   }
   default_feature_type <- head(c(default_feature_type, "exon"), 1);
   gap_feature_type <- head(c(gap_feature_type, "gap"), 1);
   if (addGaps) {
      if (verbose) {
         printDebug("getGRcoverageFromBw(): ",
            "addGRgaps()");
      }
      newValues <- list(feature_type=gap_feature_type);
      names(newValues)[1] <- feature_type_colname;
      if (!feature_type_colname %in% colnames(values(gr))) {
         values(gr)[[feature_type_colname]] <- default_feature_type;
      }
      gr <- addGRgaps(gr,
         newValues=newValues,
         ...);
   }
   covL <- lapply(nameVectorN(bwUrls), function(iBw){
      bwUrl <- bwUrls[[iBw]];
      if (verbose) {
         printDebug("getGRcoverageFromBw(): ",
            "Importing bwUrl:",
            bwUrl);
      }
      cov1 <- rtracklayer::import(bwUrl,
         selection=rtracklayer::BigWigSelection(gr),
         as="NumericList");
   });
   #grNew <- gr[,0];
   values(gr)[,names(covL)] <-S4Vectors::DataFrame(covL);
   #if (addGaps) {
   #   values(grNew)[,feature_type_colname] <- values(gr)[,feature_type_colname];
   #}
   return(gr);
}

#' Combine GRanges coverage replicates
#'
#' Combine GRanges coverage replicates
#'
#' This function takes a GRanges object as output from
#' `getGRcoverageFromBw()` and combines the coverages into
#' one coverage per strand for each `covName` (equivalent
#' to `sample_id`). Each coverage value is multiplied by
#' its `scaleFactors` value, then the sum is returned for
#' each strand, for each `covName` (`sample_id`).
#'
#' The strand is inferred by the
#' presence of negative values, where any negative value
#' indicates the column is negative strand.
#'
#' @return GRanges object whose colnames contain the
#'    `covName` (`sample_id`) for each observed strand,
#'    with coverage combined taking the sum of individual
#'    coverages after multiplying each by `scaleFactors`.
#'
#' @family jam GRanges functions
#' @family jam RNA-seq functions
#'
#' @param gr GRanges object containing coverage data in columns
#'    containing NumericList class data.
#' @param covNames character vector of `colnames(values(gr))`
#'    that contain coverage in NumericList format.
#' @param covName character vector with length equal to
#'    `length(covNames)` representing the `sample_id` for each
#'    `covNames` entry.
#' @param strands character vector, or NULL, indicating the strand
#'    for which the coverage data was obtained. When `NULL` the strand
#'    is inferred by the presence of any negative values.
#' @param scaleFactors numeric vector length equal to `length(covNames)`
#'    of values to multiply by each coverage result, in order to
#'    normalized coverages to each other.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
combineGRcoverage <- function
(gr,
 covNames=NULL,
 covName=NULL,
 strands=NULL,
 scaleFactors=1,
 verbose=FALSE,
 ...)
{
   ## Purpose is to combine replicate coverage values
   ##
   ## covName is used to group together replicates of a sample_id
   ## covNames is a unique name per coverage
   if (any(!covNames %in% colnames(values(gr)))) {
      stop("combineGRcoverage() requires covNames be present in colnames(values(gr)).");
   }
   if (length(scaleFactors) == 0) {
      scaleFactors <- 1;
   }
   if (length(covName) == 0) {
      covName <- "cov";
   }
   covName <- rep(covName, length.out=length(covNames));
   names(covName) <- covNames;
   if (length(scaleFactors) != length(covNames)) {
      scaleFactors <- rep(scaleFactors, length.out=length(covNames));
   }
   names(scaleFactors) <- covNames;
   if (length(strands) == 0) {
      strands <- factor(sapply(seq_along(covNames), function(i){
         iCov <- covNames[[i]];
         ifelse(any(any(values(gr)[[iCov]] * scaleFactors[i] < 0) &
               all(values(gr)[[iCov]] * scaleFactors[i] <= 0)),
            "-",
            "+")
      }), levels=c("+", "-"));
   }
   if (!"factor" %in% class(strands)) {
      strands <- factor(strands,
         levels=jamba::provigrep(c("[+]|pos|plus", "[-]|neg|minus", "."),
            unique(strands)));
   }
   covNamesL <- split(covNames, paste0(covName, strands));
   covNameL <- split(covName, paste0(covName, strands));
   for (iCovN in names(covNamesL)) {
      iCovV <- covNamesL[[iCovN]];
      if (verbose) {
         printDebug("combineGRcoverage(): ",
            "iCovN:", iCovN);
         printDebug("combineGRcoverage(): ",
            "iCovV:", iCovV);
      }
      iCovX <- Reduce("+",
         lapply(iCovV, function(i){
            values(gr)[[i]] * scaleFactors[i]
         }));
      if (verbose) {
         printDebug("combineGRcoverage(): ",
            "iCovN:",
            iCovN);
      }
      values(gr)[[iCovN]] <- iCovX;
   }
   gr <- gr[,!colnames(values(gr)) %in% covNames];
   attr(gr, "covNames") <- names(covNamesL);
   attr(gr, "covName") <- cPasteUnique(covNameL);
   return(gr);
}

#' Prepare Sashimi plot data
#'
#' Prepare Sashimi plot data
#'
#' This function is the workhorse function used to produce
#' Sashimi plots, and is intended to be a convenient wrapper
#' function for several other individual functions.
#'
#' At a minimum, a Sashimi plot requires three things:
#'
#' 1. The gene of interest, with corresponding exon model.
#' 2. RNA-seq coverage data.
#' 3. Splice junction data.
#'
#' There is some required pre-processing before running
#' `prepareSashimi()`:
#'
#' * Prepare flattened exons by gene using `flattenExonsByGene()`
#' and corresponding data, including `exonsByGene`, `cdsByGene`,
#' and `tx2geneDF`. Verify the gene exon model data using
#' `gene2gg()`.
#' * Find file paths, or web URLs, for a set of bigWig coverage
#' files, representing RNA-seq coverage for each strand, for
#' the samples of interest. Test the coverage data using
#' `getGRcoverageFromBw()` for a small set of GRanges data.
#' * Find file paths, or web URLs, for a set of BED6 or BED12
#' format files, note that it cannot currently use bigBed format
#' due to limitations in the `rtracklayer` package.
#' Test the splice junction data using `rtracklayer::import()`
#' for a small range of GRanges features, then send the data
#' to `spliceGR2junctionDF()` to prepare a data.frame summary.
#'
#' The basic input for coverage and junction data is a data.frame,
#' which defines each file path or url, the type of data
#' `"bw"` or `"junction"`, and the biological sample `"sample_id"`.
#' Any file path compatible with `rtracklayer::import()` will
#' work, including web URLs and local files. When using a web URL
#' you may need to use `"https://"` format to force the use
#' of secure web requests, but this requirement varies by country.
#'
#' @return list containing `ggSashimi` a ggplot2 graphical object
#'    containing a full Sashimi plot; `ggCov` the RNA-seq coverage
#'    subset of the Sashimi plot; `ggJunc` the splice junction
#'    subsset of the Sashimi plot; `ref2c` the output of
#'    `make_ref2compressed()` used for ggplot2 coordinate
#'    visualization; `covDF`, `juncDF` data.frame objects
#'    with the raw data used to create ggplot2 objects;
#'    `covGR`, `juncGR` the GRanges objects used to create
#'    the data.frames; `gr` the GRanges object representing the
#'    exons for the gene of interest; `juncLabelDF` the data.frame
#'    containing exon label coordinates used to add labels to
#'    the splice junction arcs.
#'
#' @param flatExonsByGene GRangesList named by gene, whose GRanges
#'    elements are flattened, disjoint, non-overlapping genomic ranges
#'    per gene.
#' @param filesDF data.frame with columns `url`, `sample_id`, `type`,
#'    where: `url` is any valid file path or URL compatible with
#'    `base::read.table()`; `sample_id` is an identified representing
#'    a biological sample, used to group common files together;
#'    `type` is one of `"bw"` for bigWig coverage, `"junction"` for
#'    BED12 format splice junctions.
#' @param gene character string of the gene to prepare, which must be
#'    present in `names(flatExonsByGene)`.
#' @param gapWidth numeric value of the fixed width to use for
#'    gaps (introns) between exon features. If `NULL` then
#'    `getGRgaps()` will use the default based upon the median exon
#'    width.
#' @param addGaps logical indicating whether to include gap regions
#'    in the coverage plot, for example including introns or intergenic
#'    regions. When `compressGR=TRUE` then gaps regions are
#'    down-sampled using running maximum signal with roughly the same
#'    x-axis resolution as uncompressed regions.
#' @param baseline numeric vector named by `names(flatExonsByGene)`
#'    where baseline is used to adjust the y-axis baseline position
#'    above or below zero.
#' @param compressGR logical indicating whether to compress GRanges
#'    coordinates in the output data, where gaps/introns are set
#'    to a fixed width. When `ref2c` is not supplied, and
#'    `compressGR=TRUE`, then `ref2c` is created using
#'    `make_ref2compressed()`.
#' @param ref2c list object output from `make_ref2compressed()` used
#'    to compress axis coordinates, to compress polygon coverage
#'    data in compressed regions, and to adjust splice junction arcs
#'    using compressed coordinates.
#' @param gapname the default feature_type value to use for gaps when
#'    `addGaps=TRUE`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passed to `make_ref2compressed()`,
#'    `getGRcoverageFromBw()`, `exoncov2polygon()`.
#'
#' @family RNA-seq functions
#' @family jam plot functions
#'
#' @examples
#' # First assemble a data.frame of coverage bigWig and junction BED files
#' if (1 == 2) {
#' ## Steps below will capture all replicate coverage and junction
#' ## files for the Farris et all dataset, and although are relatively
#' ## efficient (5-10 seconds per gene figure for 8 samples), are
#' ## too slow to be allowed in CRAN package docs.
#' ##
#' baseurl <- "https://orio.niehs.nih.gov/ucscview/farrisHub/mm10/";
#' bedext <- ".STAR_mm10.combinedJunctions.bed";
#' bwext <- c("492_1.sickle.merged.cutadapt.STAR_mm10.pos.bw",
#'    "492_1.sickle.merged.cutadapt.STAR_mm10.neg.bw");
#' c1 <- c("CA1", "CA2");
#' r1 <- c("CB", "DE");
#' bedsamples <- paste0(rep(c1, each=2), "_", r1);
#' bedurls <- paste0(baseurl,
#'    bedsamples,
#'    bedext);
#' bedurls;
#' bwsamples <- paste0(rep(c1, each=4), rep(r1, each=2));
#' bwurls <- paste0(baseurl,
#'    bwsamples,
#'    bwext);
#' bwurls;
#' filesDF <- data.frame(url=c(bedurls, bwurls),
#'    type=rep(c("junction", "bw"), c(4,8)),
#'    sample_id=gsub("_", "", c(bedsamples, bwsamples)));
#' filesDF;
#'
#' ## Next assemble exons by gene from a GTF file
#' if (!require(rtracklater)) {
#'    stop("rtracklayer is required.");
#' }
#' vM12gtf <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gtf.gz";
#' vM12gtfBase <- basename(vM12gtf);
#' if (!file.exists(vM12gtfBase)) {
#'    curl::curl_download(url=vM12gtf,
#'       destfile=vM12gtfBase);
#' }
#' localDb <- file.path(".", "vM12gtf.txdb");
#' if (!file.exists(localDb)) {
#'    vM12txdb <- GenomicFeatures::makeTxDbFromGFF(vM12gtfBase);
#'    AnnotationDbi::saveDb(x=vM12txdb, file=localDb);
#' } else {
#'    vM12txdb <- AnnotationDbi::loadDb(file=localDb);
#' }
#' exonsByTx <- GenomicFeatures::exonsBy(vM12txdb,
#'    by="tx",
#'    use.names=TRUE);
#' values(exonsByTx@unlistData)[,"transcript_id"] <- rep(names(exonsByTx),
#'    elementNROWS(exonsByTx));
#'
#' ## Use the GTF file to define a tx2geneDF data.frame
#' tx2geneDF <- makeTx2geneFromGtf(GTF=vM12gtfBase);
#'
#' ## Next flatten exons per gene, optionally using detectedTx
#' ## to limit the exon models to observed transcripts.
#' flatExonsByGene <- flattenExonsByGene(exonsByTx=exonsByTx,
#' #   detectedTx=detectedTx,
#'    tx2geneDF=tx2geneDF);
#'
#' ## Finally, you can choose a gene, and any sample_id from
#' ## filesDF$sample_id to prepare Sashimi plot data
#' shGria1 <- prepareSashimi(flatExonsByGene=flatExonsByGene,
#'    gene="Gria1",
#'    sample_id=c("CA2CB", "CA2DE"),
#'    filesDF=filesDF);
#'
#' ## Print the structure of data returned
#' jamba::sdim(shGria1);
#'
#' ## Of initial interest is the ggplot2 sashimi plot
#' print(shGria1$ggSashimi);
#'
#' ## You can just plot coverage
#' print(shGria1$ggCov);
#'
#' ## You can just plot junction reads
#' print(shGria1$ggJunc);
#'
#' ## For more intricate manipulations, use the data.frames
#' head(shGria1$covDF);
#' head(shGria1$juncDF);
#'
#' ## For kicks, the ref2c data
#' sdim(shGria1$ref2c)
#' ## the shGria1$ref2c$trans_grc is usable in ggplot like this
#' ## + trans_new(x=shGria1$ref2c$trans_grc)
#' ## or
#' ## + scale_x_continuous(trans=shGria1$ref2c$trans_grc)
#'
#' }
#'
#' @export
prepareSashimi <- function
(flatExonsByGene=NULL,
 filesDF=NULL,
 gene,
 sample_id=NULL,
 minJunctionScore=10,
 gapWidth=200,
 addGaps=TRUE,
 baseline=0,
 compressGR=TRUE,
 ref2c=NULL,
 gap_feature_type="intron",
 default_feature_type="exon",
 feature_type_colname="feature_type",
 exon_label_type=c("none", "repel", "mark"),
 junc_label_type=c("repel", "mark", "none"),
 return_data=c("ggCov", "ggJunc", "ggSashimi", "covDF", "juncDF", "ref2c", "all"),
 junc_color=alpha2col("goldenrod3", 0.7),
 junc_fill=alpha2col("goldenrod1", 0.4),
 coord_method=c("coord", "scale", "none"),
 scoreFactor=1,
 verbose=FALSE,
 ...)
{
   ## Purpose it to wrapper several functions used to prepare various
   ## types of data for Sashimi plots
   ##
   ## TODO: Allow filtering junction scores based upon the exon coverage
   ## so the threshold is proportional to the coverage.
   junc_label_type <- match.arg(junc_label_type);
   exon_label_type <- match.arg(exon_label_type);
   coord_method <- match.arg(coord_method);
   retVals <- list();
   #if (length(return_data) > 1 || "all" %in% return_data) {
   #}

   ## Validate filesDF
   if (length(filesDF) == 0 ||
         !jamba::igrepHas("tibble|tbl|data.*frame", class(filesDF))) {
      stop("filesDF must be a data.frame or equivalent");
   }
   if (!all(c("url","sample_id","type") %in% colnames(filesDF))) {
      stop("filesDF must contain colnames 'url', 'sample_id', and 'type'.");
   }
   ## validate other input
   if (!igrepHas("GRangesList", class(flatExonsByGene))) {
      stop("flatExonsByGene must be GRangesList")
   }
   if (length(sample_id) == 0) {
      sample_id <- unique(filesDF$sample_id);
   }

   ############################################
   ## Extract exons
   gr <- flatExonsByGene[gene]@unlistData;
   if (any(c("all") %in% return_data)) {
      retVals$gr <- gr;
   }

   ## Compress GRanges coordinates
   if (length(ref2c) == 0) {
      if (compressGR) {
         if (verbose) {
            printDebug("prepareSashimi(): ",
               "running make_ref2compressed.");
         }
         ref2c <- make_ref2compressed(gr=gr,
            gapWidth=gapWidth,
            ...);
      } else {
         ref2c <- NULL;
         coord_method <- "none";
      }
   }
   if (any(c("all", "ref2c") %in% return_data)) {
      retVals$ref2c <- ref2c;
   }

   ############################################
   ## Get coverage data
   bwFilesDF <- filesDF[filesDF$type %in% "bw" &
         filesDF$sample_id %in% sample_id,,drop=FALSE];
   bwUrls <- nameVector(bwFilesDF[,c("url","url")]);
   bwSamples <- nameVector(bwFilesDF[,c("sample_id","url")]);
   bwUrlsL <- split(bwUrls, unname(bwSamples));
   if ("scale_factor" %in% colnames(bwFilesDF)) {
      bwScaleFactors <- rmNA(naValue=1,
         bwFilesDF$scale_factor);
   } else {
      bwScaleFactors <- rep(1, length(bwUrls));
   }
   names(bwScaleFactors) <- names(bwUrls);
   #if (any(bwScaleFactors != 1)) {
      printDebug("prepareSashimi(): ",
         "bwScaleFactors:",
         bwScaleFactors);
   #}
   if (verbose) {
      printDebug("prepareSashimi(): ",
         "bwUrls:");
      if (length(bwUrls) > 0) {
         print(data.frame(bwUrls));
      } else {
         printDebug("No coverage files", fgText="red");
      }
      printDebug("GRanges:");
      print(gr);
   }
   if (length(bwUrls) > 0) {
      covGR <- getGRcoverageFromBw(gr=gr,
         bwUrls=bwUrls,
         addGaps=addGaps,
         gap_feature_type=gap_feature_type,
         default_feature_type=default_feature_type,
         feature_type_colname=feature_type_colname,
         verbose=verbose,
         ...);
      ## Combine coverage per strand
      covGR2 <- combineGRcoverage(covGR,
         covName=bwSamples,
         scaleFactors=bwScaleFactors,
         covNames=names(bwUrls));
      #retVals$covGR <- covGR2;
      ## Obtain the new set of covNames
      covNames <- attr(covGR2, "covNames");
      covName <- attr(covGR2, "covName");
      if (any(c("all", "covDF") %in% return_data)) {
         retVals$covGR <- covGR;
         retVals$covGR2 <- covGR2;
      }

      ## Create polygon data.frame
      covDF <- exoncov2polygon(covGR2,
         ref2c=ref2c,
         covNames=covNames,
         sample_id=covName,
         coord_style="fortify",
         ...);
      if (any(c("all", "covDF") %in% return_data)) {
         retVals$covDF <- covDF;
      }
      ## Add feature_type data
      if (feature_type_colname %in% colnames(values(covGR2)) &&
            !feature_type_colname %in% colnames(covDF)) {
         covDF[,feature_type_colname] <- values(covGR2)[match(
            as.character(covDF$gr),
            names(covGR2)),feature_type_colname];
      }

      ## ggplot2 for exon coverage data
      if (any(c("all","ggCov","ggSashimmi") %in% return_data)) {
         ggCov <- ggplot(covDF, aes(x=x, y=y, group=gr, fill=gr)) +
            geom_shape(show.legend=FALSE) +
            theme_jam() +
            scale_fill_jam() +
            facet_grid(~sample_id, scales="free_y");
         if ("scale" %in% coord_method) {
            ggCov <- ggCov +
               scale_x_continuous(trans=ref2c$trans_grc);
         } else if ("coord" %in% coord_method) {
            ggCov <- ggCov +
               coord_trans(x=ref2c$trans_grc);
         }
         if (any(c("all", "ggCov") %in% return_data)) {
            retVals$ggCov <- ggCov;
         }
         ########################################
         ## Optional exon labels
         covDFsub <- (as.character(covDF$gr) %in% names(gr));
         covDFlab <- covDF[covDFsub,,drop=FALSE];

         exonLabelDF1 <- shrinkMatrix(covDFlab[,c("x","y")],
            groupBy=pasteByRowOrdered(covDFlab[,c("gr", "sample_id")], sep=":!:"),
            shrinkFunc=function(x){mean(range(x))});
         exonLabelDF1[,c("gr","sample_id")] <- rbindList(
            strsplit(as.character(exonLabelDF1$groupBy), ":!:"));
         exonLabelDF1$gr <- factor(exonLabelDF1$gr, levels=unique(exonLabelDF1$gr));
         exonLabelDF <- renameColumn(exonLabelDF1,
            from="groupBy",
            to="gr_sample");
         retVals$exonLabelDF <- exonLabelDF;
         if ("mark" %in% exon_label_type) {
            ggExonLabels <- ggforce::geom_mark_hull(data=exonLabelDF,
               aes(x=x, y=y, label=gr, group=gr_sample),
               fill="transparent",
               label.fontsize=10,
               label.fill=alpha2col("white", 0.3),
               colour="transparent",
               concavity=1,
               label.buffer=unit(2, 'mm'),
               con.cap=0,
               expand=unit(1, "mm"),
               con.border="one",
               con.colour="navy",
               label.colour="navy");
            ggCov <- ggCov + ggExonLabels;
         } else if ("repel" %in% junc_label_type) {
            yMax <- max(exonLabelDF$y);
            yUnit <- 10^floor(log10(yMax));
            yMaxUse <- floor(yMax/yUnit)*yUnit;
            ggExonLabels <- ggrepel::geom_text_repel(data=exonLabelDF,
               inherit.aes=FALSE,
               aes(x=x, y=y,
                  group=gr_sample,
                  fill="transparent",
                  label=gr),
               #nudge_y=yMaxUse-juncLabelDF$y, vjust=1, ## Used to fix labels at certain height
               angle=90,
               vjust=1,
               direction="y",
               point.padding=0
            );
            ggCov <- ggCov + ggExonLabels;
         } else {
            ggExonLabels <- NULL;
         }
      }
   }

   ############################################
   ## Load Junctions
   juncFilesDF <- filesDF[filesDF$type %in% "junction" &
         filesDF$sample_id %in% sample_id,c("url","sample_id"),drop=FALSE];
   juncUrls <- nameVector(juncFilesDF);
   juncSamples <- nameVector(juncFilesDF[,c("sample_id","sample_id")]);
   juncUrlsL <- split(juncUrls, juncSamples);
   if ("score_factor" %in% colnames(juncFilesDF)) {
      juncScaleFactors <- rmNA(naValue=1,
         juncFilesDF$scale_factor);
   } else {
      juncScaleFactors <- rep(1, length(juncUrls));
   }
   names(juncScaleFactors) <- names(juncUrls);
   #juncUrls <- nameVector(
   #   filesDF[filesDF$type %in% "junction" & filesDF$sample_id %in% sample_id,c("url","sample_id"),drop=FALSE]);
   if (verbose) {
      printDebug("prepareSashimi(): ",
         "juncUrls:");
      if (length(juncUrls) > 0) {
         print(data.frame(juncUrls));
      } else {
         printDebug("No junction files", fgText="red");
      }
   }
   if (length(juncUrls) > 0) {
      juncBedGR <- GRangesList(lapply(nameVectorN(juncUrls), function(iBedName){
         iBed <- juncUrls[[iBedName]];
         if (verbose) {
            printDebug("prepareSashimi(): ",
               "Importing bed:",
               iBed);
         }
         bed1 <- rtracklayer::import(iBed,
            which=range(gr));
         values(bed1)$score <- as.numeric(values(bed1)$name) * juncScaleFactors[iBedName];
         values(bed1)[,c("juncNames")] <- iBedName;
         values(bed1)[,c("sample_id")] <- juncSamples[iBedName];
         bed1;
      }))@unlistData;
      ## Create junction summary data.frame
      if (verbose) {
         printDebug("prepareSashimi(): ",
            "running spliceGR2junctionDF()");
      }
      juncDF1 <- spliceGR2junctionDF(spliceGRgene=juncBedGR,
         exonsGR=gr,
         sampleColname="sample_id");
      juncGR <- as(renameColumn(juncDF1, from="ref", to="seqnames"), "GRanges");
      names(juncGR) <- makeNames(values(juncGR)[,"nameFromToSample"]);
      ## Convert junctions to polygons usable by geom_diagonal_wide()
      #      juncPolyDF <- grl2df(setNames(GRangesList(subset(juncGR, score > minJunctionScore)), sample_id),
      if (length(baseline) == 0) {
         baseline <- 0;
      }
      juncDF <- grl2df(split(juncGR, values(juncGR)[["sample_id"]]),
         shape="junction",
         ref2c=ref2c,
         scoreFactor=juncScaleFactors,
         scoreArcFactor=0.3,
         baseline=baseline,
         verbose=verbose,
         doStackJunctions=TRUE);
      if (!"sample_id" %in% colnames(juncDF)) {
         juncDF <- renameColumn(juncDF,
            from="grlNames",
            to="sample_id");
      }
      juncDF <- juncDF[,!colnames(juncDF) %in% c("grlNames"),drop=FALSE];

      if ("Gria1_exon13 Gria1_exon15" %in% values(juncGR)[["nameFromTo"]]) {
         printDebug("juncDF1:");
         print(subset(juncGR, nameFromTo %in% "Gria1_exon13 Gria1_exon15"));
      }
      if ("Gria1_exon13 Gria1_exon15" %in% juncDF1$nameFromTo) {
         printDebug("juncDF1:");
         print(subset(juncDF1, nameFromTo %in% "Gria1_exon13 Gria1_exon15"));
      }
      if ("Gria1_exon13 Gria1_exon15" %in% juncDF$nameFromTo) {
         printDebug("juncDF:");
         print(subset(juncDF, nameFromTo %in% "Gria1_exon13 Gria1_exon15"));
      }
      if (any(c("all", "juncDF") %in% return_data)) {
         retVals$juncDF1 <- juncDF1;
         retVals$juncDF <- juncDF;
      }

      ## define junction label positions
      #juncLabelDF1 <- subset(mutate(juncCoordDF, id_name=makeNames(id)), grepl("_v1_v3$", id_name));
      juncLabelDF1 <- subset(plyr::mutate(juncDF, id_name=makeNames(id)),
         grepl("_v1_v[23]$", id_name));
      juncLabelDF <- renameColumn(
         shrinkMatrix(juncLabelDF1[,c("x","y","score"),drop=FALSE],
            groupBy=juncLabelDF1[,"nameFromToSample"]),
         from="groupBy",
         to="nameFromToSample");
      juncLabelDF[,c("nameFromTo","sample_id")] <- rbindList(
         strsplit(juncLabelDF[,"nameFromToSample"], ":!:"));
      juncLabelDF[,c("nameFrom", "nameTo")] <- rbindList(
         strsplit(juncLabelDF[,"nameFromTo"], " "));

      if (any(c("all", "juncLabelDF") %in% return_data)) {
         retVals$juncLabelDF <- juncLabelDF;
      }

      ## ggplot2 for junction data
      if (any(c("all","ggJunc","ggSashimi") %in% return_data)) {
         ggJunc <- ggplot(juncDF) +
            ggforce::geom_diagonal_wide(aes(x=x, y=y, group=id),
               color=junc_color,
               fill=junc_fill,
               strength=0.4);
         ggJunc <- ggJunc +
            facet_grid(~sample_id, scales="free_y");
         if ("mark" %in% junc_label_type) {
            ggJuncLabels <- ggforce::geom_mark_hull(data=juncLabelDF,
               aes(x=x, y=y, group=nameFromToSample, label=round(score)),
               fill="transparent",
               label.fontsize=10,
               label.fill=alpha2col("white", 0.3),
               colour="transparent",
               concavity=1,
               label.buffer=unit(2, 'mm'),
               con.cap=0,
               expand=unit(1, "mm"),
               con.border="one",
               con.colour="navy",
               label.colour="navy");
            ggJunc <- ggJunc + ggJuncLabels;
         } else if ("repel" %in% junc_label_type) {
            yMax <- max(juncLabelDF$y);
            yUnit <- 10^floor(log10(yMax));
            yMaxUse <- floor(yMax/yUnit)*yUnit;yMaxUse;
            ggJuncLabels <- ggrepel::geom_text_repel(data=juncLabelDF,
               inherit.aes=FALSE,
               aes(x=x, y=y,
                  group=nameFromToSample,
                  fill="transparent",
                  label=scales::comma(round(score))),
               #nudge_y=yMaxUse-juncLabelDF$y, vjust=1, ## Used to fix labels at certain height
               angle=90, vjust=0.5,
               direction="y",
               point.padding=0
               );
            ggJunc <- ggJunc + ggJuncLabels;
         } else {
            ggJuncLabels <- NULL;
         }
         ggJunc <- ggJunc +
            theme_jam() + scale_fill_jam();
         if ("scale" %in% coord_method) {
            ggJunc <- ggJunc +
               scale_x_continuous(trans=ref2c$trans_grc);
         } else if ("coord" %in% coord_method) {
            ggJunc <- ggJunc +
               coord_trans(x=ref2c$trans_grc);
         }
         if (any(c("all", "ggJunc") %in% return_data)) {
            retVals$ggJunc <- ggJunc;
         }
      }
   }

   #########################################
   ## ggplot2 for Sashimi plot
   if (any(c("all","ggSashimi") %in% return_data)) {
      if (length(bwUrls) > 0) {
         ggSashimi <- ggCov;
         if (length(juncUrls) > 0) {
            ggSashimi <- ggSashimi +
               geom_diagonal_wide(data=juncDF,
                  aes(x=x, y=y, group=id),
                  color=junc_color,
                  fill=junc_fill,
                  strength=0.4);
            if (length(ggJuncLabels) > 0) {
               ggSashimi <- ggSashimi + ggJuncLabels;
            }
         }
      } else {
         if (length(juncUrls) > 0) {
            ggSashimi <- ggJunc;
         } else {
            ggSashimi <- NULL;
         }
      }
      retVals$ggSashimi <- ggSashimi;
   }
   return(retVals);
}

