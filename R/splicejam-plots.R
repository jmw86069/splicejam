
#' Create an interactive 3-D BGA plotly visualization
#'
#' Create an interactive 3-D BGA plotly visualization
#'
#' This function takes the output from `made4::bga()` (class `"bga"`)
#' and creates a 3D plotly visualization which shows the samples,
#' sample centroids, and optionally a subset of genes in the same
#' BGA component space.
#'
#' For general guidance, when the data contains four or more sample
#' groups, this function can be called to produce a 3-D plot.
#' BGA itself produces `n-1` dimensions, therefore when there are
#' fewer than four groups, a 3-D plot cannot be produced.
#'
#' When analyzing data with fewer than four groups, it is sometimes
#' helpful to cluster samples individually. This approach is roughly
#' equivalent to performing PCA or COA directly, but is sometimes
#' convenient in being consistent with calling `made4::bga()`, as
#' it also applies the same scaling method used in `made4:bga()`
#' prior to calling PCA or COA. In that case, the groups can still
#' be displayed by `bgaPlotly3d()` by defining the `sampleGroups`
#' argument. The result will display typical PCA or COA coordinates,
#' but visually displays the centroid of each sample group.
#'
#' This function is especially useful when there are multiple
#' sample groups that belong to a "supergroup" -- a higher level
#' of grouping. Consider an experiment with two cell lines treated
#' over multiple time points. The samples would be assigned to
#' a group based upon cell_line and time, then each group would
#' be assigned to a supergroup based upon cell_line. The output
#' plot will display each sample group, group centroid, then
#' connect group centroids in order by the supergroup.
#' The order of groups within each supergroup is determined
#' either by factor levels, or by the order returned by
#' `jamba::mixedSort()` which performs a proper alphanumeric sort.
#'
#' ## Background on `bga` object content
#'
#' Output from `made4::bga()` is an object with `class` containing `"bga"`.
#' It is a list with three elements:
#'
#' 1. `"ord"`: ordination, the initial clustering results from either `"coa"`
#' (default) or `"pca"` clustering across all individual samples. This data
#' is fed into the second phase, group-based clustering.
#' 2. `"bet"`: the between group analysis result, derived from `"ord"`.
#' 3. `"fac"`: `factor` vector of groups across the `colnames(x)` of
#' input data.
#'
#' The `"bet"` contents are described below, assuming the results of
#' `made::bga()` are stored as `bgaInfo`:
#'
#' * gene (row) coordinates are stored:
#'
#'    * `bgaInfo$bet$co`: unscaled coordinates for each component axis
#'    * `bgaInfo$bet$c1`: scaled coordinates for each component axis
#'
#' * group centroid coordinates are stored:
#'
#'    * `bgaInfo$bet$li`: unscaled group centroid coordinates
#'    * `bgaInfo$bet$ls`: scaled group centroid coordinates
#'
#' * sample centroid coordinates are stored:
#'
#'    * `bgaInfo$bet$ls`: scaled coordinates for each component axis
#'
#' @return plotly object sufficient to render an HTML page containing
#' a 3-D BGA plot.
#'
#' @param bgaInfo object of class `c("coa", "bga")` created by `made4::bga()`.
#' @param axes `integer` vector length 3 indicating the BGA axes to use
#'    for the 3-D visualization. It is sometimes helpful to use
#'    `axes=c(2,3,4)` especially when the first component suggests
#'    a large experimental batch effect, for example when
#'    combining data from multiple GEO studies, often the first
#'    BGA component is indicative of the different experimental
#'    platforms.
#' @param superGroups `character` or factor vector with length
#'    `length(bgaInfo$fac)`, in the same order. A super group consists
#'    of multiple sample groups, as defined by values in `bgaInfo$fac`.
#' @param superGroupAlpha `numeric` value between 0 and 1 indicating the
#'    alpha transparency to use for the vector that connects groups
#'    within each super group.
#' @param superGroupLwd `numeric` line width when connecting sample group
#'    centroids by `superGroups`.
#' @param arrowSmoothFactor `numeric` value typically ranging from 1 to 5,
#'    affecting the amount of curve smoothing when generating a
#'    3-D spline curve to connect more than two super groups.
#'    When there are only two group members in one super group,
#'    the two groups will always use a straight line. For three or
#'    more groups, a 3-D spline is used, in order to help visualize
#'    trends by means of a curve.
#' @param colorSub `character` vector of R colors, whose names contain values in
#'    `bgaInfo$fac`. These colors are used for the samples, sample
#'    centroids, and are shaded via a gradient when connecting multiple
#'    centroids using `superGroups`.
#' @param drawVectors `character` vector of one or more types of vectors:
#'    * `"centroids"` draws a vector from the origin to each group centroid.
#'    * `"genes"` draws a vector from the origin to each gene.
#' @param drawSampleLabels `logical` indicating whether to label each
#'    individual sample point, as opposed to labeling only the
#'    sample group centroid.
#' @param ellipseType `character` value indicating how to draw sample
#'    centroid ellipses, where `"none"` draws no ellipse; `"alphahull"`
#'    draws a wireframe; and `"ellipsoid"` draws
#'    a partially transparent shell. Note that plotly does not handle
#'    transparency well in 3-D space, sometimes fully obscuring
#'    background components.
#' @param ellipseOpacity `numeric` value between 0 and 1 indicating the
#'    transparency of the ellipse, when `ellipseType` is not `"none"`;
#'    0 is fully transparent, and 1 is non-transparent.
#' @param useScaledCoords `logical` indicating whether to use the BGA
#'    scaled coordinates, or the raw coordinates. By default the raw
#'    coordinates are used in order to preserve the relative contribution
#'    of each BGA component to the overall variability (or inertial
#'    moment.)
#' @param geneColor `character` color used when displaying gene points.
#' @param geneScaleFactor `numeric` value used to expand the position of
#'    genes relative to the origin, intended to help visualize the vector
#'    position of genes relative to the origin.
#' @param geneLwd,geneAlpha `numeric` used to customize the line width
#'    and alpha transparency of gene points, respectively.
#' @param geneLabels `character` vector of gene labels to use instead
#'    of rownames. The `names(geneLabels)` are expected to match
#'    `rownames(bgaInfo$bet$co)`. This option is intended to allow
#'    display of gene symbols, instead of an assay identifier like
#'    probe set name.
#' @param maxGenes `integer` maximum number of genes to display.
#' @param centroidLwd,centroidAlpha `numeric` used to customize the
#'    line width and alpha transparency of vectors drawn to sample
#'    centroids, respectively.
#' @param textposition_centroid,textposition_sample `character` vector
#'    containing values used as the `plotly` argument `"textposition"`
#'    to position text labels for centroids, and samples. Each
#'    vector is extended to the number of labels, which allows
#'    each label to be individually positioned.
#' @param sceneX,sceneY,sceneZ `numeric` used to define the camera
#'    position in 3-D space as defined by plotly. Frankly, the ability
#'    to customize the starting position is fairly confusing, but these
#'    values are faithfully passed along to the corresponding plotly
#'    call.
#' @param main `character` string used as a title on the plotly
#'    visualization.
#' @param sampleGroups `character` or `factor` vector, whose
#'    names are expected to contain values in
#'    `rownames(bgaInfo$ord$ord$tab)`, which should also be `colnames(x)`
#'    of the input data to `made4::bga(x, ...)`. The values in
#'    `sampleGroups` are used to re-group samples, which effectively
#'    calculates new sample group centroids used in this
#'    visualization. It is intended to be used when the `made4::bga()`
#'    is used without grouping, for example when each sample replicate
#'    is defined as its own group. The end result is that centroids
#'    and ellipsoids are calculated for each group during visualization,
#'    to help organize the points displayed.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional parameters are ignored
#'
#' @family jam plot functions
#' @family jam spatial functions
#'
#' @export
bgaPlotly3d <- function
(bgaInfo,
 axes=c(1,2,3),
 superGroups=NULL,
 superGroupAlpha=0.3,
 superGroupLwd=25,
 arrowSmoothFactor=8,
 colorSub=NULL,
 drawVectors=c("none","centroids","genes"),
 highlightGenes=NULL,
 drawSampleLabels=TRUE,
 ellipseType=c("ellipsoid","alphahull","none"),
 ellipseOpacity=0.2,
 useScaledCoords=FALSE,
 geneColor="#555555",
 geneScaleFactor=1,
 geneLwd=2,
 geneAlpha=0.7,
 geneLabels=NULL,
 maxGenes=50,
 centroidLwd=15,
 centroidAlpha=0.3,
 textposition_centroid="top center",
 textposition_sample="middle right",
 sceneX=-1.25,
 sceneY=1.25,
 sceneZ=1.25,
 main="",
 plot_bgcolor="#eeeeee",
 paper_bgcolor="#dddddd",
 sampleGroups=NULL,
 debug=FALSE,
 verbose=TRUE,
 ...)
{
   ## Purpose is to take BGA data and create a 3D plotly visualization
   ##
   ## superGroups is a named vector whose names match names(bgaInfo$fac)
   ## and whose values are super groups, which are defined as sets of
   ## sample groups as defined by bgaInfo$fac.
   ##
   if ("none" %in% drawVectors) {
      drawVectors <- "none";
   }
   if (!jamba::check_pkg_installed("made4")) {
      stop("bgaPlotly3d() requires the made4 package.");
   }
   if (!jamba::check_pkg_installed("plotly")) {
      stop("bgaPlotly3d() requires the plotly package.");
   }
   if (jamba::igrepHas("list", class(bgaInfo)) &&
         "bgaInfo" %in% names(bgaInfo)) {
      bgaInfo <- bgaInfo$bgaInfo;
   }

   ellipseType <- match.arg(ellipseType);
   ## bgaInfo$fac    sample group factor
   ##
   ## bgaInfo$bet$li sample centroid coordinates
   ##
   ## bgaInfo$bet$ls sample coordinates
   ##
   ## bgaInfo$bet$co gene coordinates
   ##    rownames(bgaInfo$bet$co)[made4:::genes(bgaInfo$bet$co, n=10)]
   ##    will return top 10 genes
   ##

   ## Step 1:
   ## - draw samples as small spheres
   ##   - (optional) label samples
   ## - add lines connecting them to the sample centroid
   ##   - (optional) label the sample centroid
   ##
   ## Step 2: (optional)
   ## - create and draw an ellipsoid around each sample centroid
   ##
   ## Step 3: (optional)
   ## - group and order sample centroids by supergroup info
   ## - fit a 3D spline connecting each supergroup
   ## - draw the supergroup spline with arrow at the end
   ##
   ## Step 4:
   ## - add gene points using top N genes
   ##   - (optional) draw line from origin to each gene

   ## Pull out coordinates for convenience
   ## pre-sample coordinates (un-scaled)
   samplecoords_unscaled <- bgaInfo$bet$ls;
   ## group coordinates (scaled)
   groupcoords_scaled <- bgaInfo$bet$l1;
   ## group coordinates (un-scaled)
   groupcoords_unscaled <- bgaInfo$bet$li;
   ## Calculate the scaling factor adjustment
   scale_factor1 <- head(groupcoords_scaled / groupcoords_unscaled, 1);
   scale_factor <- scale_factor1[rep(1, nrow(samplecoords_unscaled)),,drop=FALSE];
   samplecoords_scaled <- samplecoords_unscaled * scale_factor;
   rownames(samplecoords_scaled) <- rownames(samplecoords_unscaled);

   ## Confirm the sample grouping
   if (length(sampleGroups) == 0 || identical(sampleGroups, bgaInfo$fac)) {
      sampleGroups <- bgaInfo$fac;
      names(sampleGroups) <- rownames(bgaInfo$ord$ord$tab);
   } else {
      ## Use different sample groupings, re-calculate sample centroids
      if (!is.factor(sampleGroups)) {
         sampleGroups <- factor(sampleGroups,
            levels=unique(as.character(sampleGroups)));
      }
      if (length(names(sampleGroups)) == 0) {
         names(sampleGroups) <- rownames(bgaInfo$ord$ord$tab);
      } else {
         sampleGroups <- sampleGroups[rownames(bgaInfo$ord$ord$tab)];
      }
      if (length(sampleGroups) != nrow(bgaInfo$ord$ord$tab)) {
         stop("sampleGroups does not match dimensions in bgaInfo, nrow(bgaInfo$ord$ord$tab)");
      }
      ## re-calculate sample centroids and update bgaInfo
      if (verbose) {
         printDebug("bgaPlotly3d(): ",
            "Re-calculating sample centroids using sampleGroups:");
         print(sampleGroups);
      }
      samplecoords_unscaled1 <- samplecoords_unscaled;
      groupcoords_scaled1 <- groupcoords_scaled;
      groupcoords_unscaled1 <- groupcoords_unscaled;
      groupcoords_unscaled <- shrinkMatrix(samplecoords_unscaled,
         groupBy=sampleGroups,
         shrinkFunc=mean,
         returnClass="matrix");
      colnames(groupcoords_unscaled) <- colnames(groupcoords_unscaled1);
      scale_factor <- scale_factor1[rep(1, nrow(groupcoords_unscaled)),,drop=FALSE];
      groupcoords_scaled <- groupcoords_unscaled * scale_factor;
      rownames(groupcoords_scaled) <- rownames(groupcoords_unscaled);
      if (verbose) {
         printDebug("bgaPlotly3d(): ",
            "groupcoords_unscaled:");
         print(groupcoords_unscaled);
         printDebug("bgaPlotly3d(): ",
            "groupcoords_scaled:");
         print(groupcoords_scaled);
      }
      bgaInfo$bet$li <- groupcoords_unscaled;
      bgaInfo$bet$l1 <- groupcoords_scaled;
      bgaInfo$fac <- sampleGroups;
   }

   ## Define some colors
   sampleColors <- jamba::nameVector(
      colorjam::group2colors(as.character(sampleGroups),
      colorSub=colorSub),
      names(sampleGroups));
   groupNameColors <- jamba::nameVector(
      sampleColors[match(unique(sampleGroups), sampleGroups)],
      unique(sampleGroups));
   if (verbose) {
      jamba::printDebug("bgaPlotly3d(): ",
         "sampleGroups:");
      print(sampleGroups);
      if (length(superGroups) > 0) {
         jamba::printDebug("bgaPlotly3d(): ",
            "superGroups:");
         print(superGroups);
      }
      jamba::printDebug("bgaPlotly3d(): ",
         "sampleColors:",
         names(sampleColors),
         fgText=list("orange","dodgerblue",sampleColors));
      jamba::printDebug("bgaPlotly3d(): ",
         "groupNameColors:",
         names(groupNameColors),
         fgText=list("orange","dodgerblue",groupNameColors));
   }


   ############################################################
   ## Define axis columns
   ## Sample coordinate colname
   Xs <- colnames(samplecoords_unscaled)[axes[1]];
   Ys <- colnames(samplecoords_unscaled)[axes[2]];
   Zs <- colnames(samplecoords_unscaled)[axes[3]];
   if (useScaledCoords) {
      ## Sample centroid coordinate colname
      Xsc <- colnames(groupcoords_scaled)[axes[1]];
      Ysc <- colnames(groupcoords_scaled)[axes[2]];
      Zsc <- colnames(groupcoords_scaled)[axes[3]];
      ## Gene coordinate colname
      Xg <- colnames(bgaInfo$bet$c1)[axes[1]];
      Yg <- colnames(bgaInfo$bet$c1)[axes[2]];
      Zg <- colnames(bgaInfo$bet$c1)[axes[3]];
   } else {
      ## Sample centroid coordinate colname
      Xsc <- colnames(groupcoords_unscaled)[axes[1]];
      Ysc <- colnames(groupcoords_unscaled)[axes[2]];
      Zsc <- colnames(groupcoords_unscaled)[axes[3]];
      ## Gene coordinate colname
      Xg <- colnames(bgaInfo$bet$co)[axes[1]];
      Yg <- colnames(bgaInfo$bet$co)[axes[2]];
      Zg <- colnames(bgaInfo$bet$co)[axes[3]];
   }
   axesVs <- jamba::nameVector(1:3, c(Xs,Ys,Zs));
   axesVsc <- jamba::nameVector(1:3, c(Xsc,Ysc,Zsc));
   axesVg <- jamba::nameVector(1:3, c(Xg,Yg,Zg));
   if (verbose) {
      jamba::printDebug("bgaPlotly3d(): ",
         "axesVs:",
         names(axesVs));
      jamba::printDebug("bgaPlotly3d(): ",
         "axesVsc:",
         names(axesVsc));
      jamba::printDebug("bgaPlotly3d(): ",
         "axesVg:",
         names(axesVg));
   }


   ############################################################
   ## Sample segments to centroid
   if (useScaledCoords) {
      dfSamplesDF <- data.frame(
         Label=rownames(samplecoords_scaled),
         samplecoords_scaled[,names(axesVs)],
         groupcoords_scaled[as.character(sampleGroups),names(axesVsc)]);
   } else {
      dfSamplesDF <- data.frame(
         Label=rownames(samplecoords_unscaled),
         samplecoords_unscaled[,names(axesVs)],
         groupcoords_unscaled[as.character(sampleGroups),names(axesVsc)]);
   }
   if (verbose) {
      printDebug("bgaPlotly3d(): ",
         "head(dfSamplesDF):");
      print(head(dfSamplesDF, 12));
      printDebug("bgaPlotly3d(): ",
         "head(samplecoords_unscaled):");
      print(head(samplecoords_unscaled, 12));
   }
   dfSamples <- rownames(dfSamplesDF);
   dfSampleLinesDF <- dfWide2segments(
      dfSamplesDF,
      axes1=names(axesVs),
      axes2=names(axesVsc));
   dfSampleLinesDF$Name <- rep(names(sampleGroups), each=3);
   dfSampleLinesDF$groupName <- rep(sampleGroups, each=3);
   dfSampleLinesDF$Symbol <- rep(c("circle","circle-open","x"),
      length(sampleGroups));
   dfSampleLinesDF$size <- rep(c(20,0,0),
      length(sampleGroups));
   dfSampleLinesDF$color <- sampleColors[dfSampleLinesDF$Name];
   i1 <- jamba::igrep("^circle$", dfSampleLinesDF$Symbol);
   i2 <- jamba::igrep("^circle-open$", dfSampleLinesDF$Symbol);
   i3 <- jamba::igrep("^x$", dfSampleLinesDF$Symbol);
   dfSampleLinesDF[,"Label"] <- "";
   if (drawSampleLabels) {
      dfSampleLinesDF[i1,"Label"] <- dfSampleLinesDF$Name[i1];
   } else {
      dfSampleLinesDF[i1,"Label"] <- "";
   }
   dfSampleLinesDF[i2,"Label"] <- as.character(dfSampleLinesDF$groupName)[i2];

   i2keep <- match(unique(sampleGroups), dfSampleLinesDF[i2,"Label"]);
   dfSampleLinesDF[i2,"Label"] <- "";
   dfSampleLinesDF[i2[i2keep],"Label"] <- as.character(dfSampleLinesDF$groupName)[i2[i2keep]];
   dfSampleLinesDF[,"textposition"] <- textposition_sample;
   if (drawSampleLabels) {
      dfSampleLinesDF[i1,"textposition"] <- textposition_sample;
   }
   dfSampleLinesDF[i2[i2keep],"textposition"] <- textposition_centroid;
   if (verbose) {
      jamba::printDebug("bgaPlotly3d(): ",
         "head(dfSampleLinesDF, 20):");
      print(head(dfSampleLinesDF, 20));
   }


   ############################################################
   ## Optionally add vectors from origin to each centroid
   if (jamba::igrepHas("sample|centroid", drawVectors)) {
      axesVscO <- paste0("origin_", names(axesVsc));
      if (useScaledCoords) {
         dfCVDF <- data.frame(
            Label=rownames(groupcoords_scaled),
            jamba::renameColumn(
               groupcoords_scaled[sampleGroups,names(axesVsc)]*0,
               from=names(axesVsc),
               to=paste0("origin_", names(axesVsc))),
            groupcoords_scaled[sampleGroups,names(axesVsc)]);
      } else {
         dfCVDF <- data.frame(
            Label=rownames(groupcoords_unscaled),
            jamba::renameColumn(groupcoords_unscaled[sampleGroups,names(axesVsc)]*0,
               from=names(axesVsc),
               to=paste0("origin_", names(axesVsc))),
            groupcoords_unscaled[sampleGroups,names(axesVsc)]);
      }
      #dfCV <- rownames(dfCVDF);
      dfCVLDF <- dfWide2segments(dfCVDF,
         axes1=names(axesVsc), axes2=axesVscO);
      dfCVLDF$Name <- rep(names(sampleGroups), each=3);
      dfCVLDF$groupName <- rep(sampleGroups, each=3);
      dfCVLDF$Symbol <- rep(c("circle","circle-open","x"), length(sampleGroups));
      dfCVLDF$size <- rep(c(2,0,0), length(sampleGroups));
      dfCVLDF$color <- sampleColors[dfCVLDF$Name];
      i1 <- jamba::igrep("^circle$", dfCVLDF$Symbol);
      i2 <- jamba::igrep("^circle-open$", dfCVLDF$Symbol);
      i3 <- jamba::igrep("^x$", dfCVLDF$Symbol);
      dfCVLDF[,"Label"] <- "";
      dfCVLDF[i1,"Label"] <- dfCVLDF$Name[i1];
      dfCVLDF[i2,"Label"] <- as.character(dfCVLDF$groupName)[i2];
      i2keep <- match(unique(sampleGroups), dfCVLDF[i2,"Label"]);
      dfCVLDF[i2,"Label"] <- "";
      dfCVLDF[i2[i2keep],"Label"] <- as.character(dfCVLDF$groupName)[i2[i2keep]];
      dfCVLDF[,"textposition"] <- rep(textposition_centroid,
         length.out=nrow(textposition_centroid));
      dfCVLDF[i1,"textposition"] <- rep(textposition_centroid,
         length.out=length(i1));
      dfCVLDF[i2[i2keep],"textposition"] <- rep(textposition_centroid,
         length.out=length(i2[i2keep]));
   }


   ############################################################
   ## Optionally add vectors from origin to genes
   if (jamba::igrepHas("gene|row", drawVectors)) {
      axesVgO <- paste0("origin_", names(axesVg));
      if (useScaledCoords) {
         iGenes <- rownames(bgaInfo$bet$c1);
         dfVgDF <- data.frame(
            Label=iGenes,
            jamba::renameColumn(
               bgaInfo$bet$c1[iGenes,names(axesVg)]*0,
               from=names(axesVg),
               to=paste0("origin_", names(axesVg))),
            bgaInfo$bet$c1[iGenes,names(axesVg)]*geneScaleFactor);
      } else {
         iGenes <- rownames(bgaInfo$bet$co);
         dfVgDF <- data.frame(
            Label=iGenes,
            jamba::renameColumn(
               bgaInfo$bet$co[iGenes,names(axesVg)]*0,
               from=names(axesVg),
               to=paste0("origin_", names(axesVg))),
            bgaInfo$bet$co[iGenes,names(axesVg)]*geneScaleFactor);
      }
      ## Determine distance from origin
      if (length(highlightGenes) > 0) {
         if (verbose) {
            jamba::printDebug("bgaPlotly3d(): ",
               "Using ",
               length(highlightGenes),
               " highlightGenes");
         }
         iGenes <- rownames(dfVgDF)[tolower(rownames(dfVgDF)) %in% tolower(highlightGenes)];
         maxGenes <- length(iGenes);
      } else {
         geneDist <- apply(abs(dfVgDF[,names(axesVg)]), 1, function(i){
            splicejam:::geomean(i, offset=1);
         });
         iGenes <- head(iGenes[order(-geneDist)], maxGenes);
         if (verbose) {
            jamba::printDebug("bgaPlotly3d(): ",
               "Using ",
               length(iGenes),
               " based upon distance from origin.");
         }
      }
      dfVgDF <- dfVgDF[rownames(dfVgDF) %in% iGenes,,drop=FALSE];

      #dfCV <- rownames(dfCVDF);
      dfVgLDF <- dfWide2segments(dfVgDF,
         axes1=names(axesVg), axes2=axesVgO);
      dfVgLDF$Name <- rep(dfVgDF$Label, each=3);
      dfVgLDF$groupName <- rep(dfVgDF$Label, each=3);
      if (length(geneLabels) > 0 && names(geneLabels) > 0) {
         geneMatch <- match(dfVgLDF$groupName, names(geneLabels));
         print(table(is.na(geneMatch)));
         dfVgLDF[,"groupName"] <- geneLabels[dfVgLDF$Name];
         dfVgLDF[,"Name"] <- dfVgLDF[,"groupName"];
         printDebug("head(dfVgLDF):");
         print(head(dfVgLDF, 20));
      }
      dfVgLDF$Symbol <- rep(c("circle","circle-open","x"), nrow(dfVgDF));
      dfVgLDF$size <- rep(c(2,0,0), nrow(dfVgDF));
      if (length(geneColor) > 1 && length(names(geneColor)) > 0) {
         dfVgLDF$color <- jamba::rmNA(geneColor[dfVgLDF$Name],
            naValue="#555555");
      } else {
         dfVgLDF$color <- rep(geneColor, length.out=nrow(dfVgLDF));
      }
      i1 <- jamba::igrep("^circle$", dfVgLDF$Symbol);
      i2 <- jamba::igrep("^circle-open$", dfVgLDF$Symbol);
      i3 <- jamba::igrep("^x$", dfVgLDF$Symbol);
      dfVgLDF[,"Label"] <- "";
      dfVgLDF[i1,"Label"] <- dfVgLDF$Name[i1];
      dfVgLDF[i2,"Label"] <- as.character(dfVgLDF$groupName)[i2];
      if (length(geneLabels) > 0) {
         i2keep <- match(unique(geneLabels[iGenes]), dfVgLDF[i2,"Label"]);
      } else {
         i2keep <- match(unique(iGenes), dfVgLDF[i2,"Label"]);
      }
      dfVgLDF[i2[i2keep],"Label"] <- as.character(dfVgLDF$groupName)[i2[i2keep]];
      dfVgLDF[i2,"Label"] <- "";
      dfVgLDF[,"textposition"] <- "top center";
      # dfVgLDF[i1,"textposition"] <- "top center";
      # dfVgLDF[i2[i2keep],"textposition"] <- "top center";
      if (verbose) {
         jamba::printDebug("bgaPlotly3d(): ",
            "head(dfVgLDF, 10):");
         print(head(dfVgLDF, 10));
         jamba::printDebug("bgaPlotly3d(): ",
            "head(dfVgDF, 10):");
         print(head(dfVgDF, 10));
      }
   }

   ############################################################
   ## First create plotly with sample points, lines, labels
   if (verbose) {
      jamba::printDebug("bgaPlotly3d(): ",
         "Creating plotly object with samples and sample centroids.");
   }
   if (debug == 1) {
      return(dfSampleLinesDF);
   }
   p9 <- plotly::plot_ly(
      data=dfSampleLinesDF,
      name="Samples",
      type="scatter3d",
      x=as.formula(paste0("~", Xs)),
      y=as.formula(paste0("~", Ys)),
      z=as.formula(paste0("~", Zs)),
      mode="markers+lines",
      text=dfSampleLinesDF$Label,
      # textposition=dfSampleLinesDF$textposition,
      hoverinfo="text",
      line=list(
         width=5,
         color=dfSampleLinesDF$color
      ),
      hoverlabel=list(
         bgcolor=dfSampleLinesDF$color,
         bordercolor=jamba::setTextContrastColor(dfSampleLinesDF$color),
         font=list(
            color=jamba::setTextContrastColor(dfSampleLinesDF$color)
         )
      ),
      marker=list(
         size=dfSampleLinesDF$size,
         color=dfSampleLinesDF$color
      ),
      textfont=list(
         color=dfSampleLinesDF$color,
         size=(40-dfSampleLinesDF$size)/2
      )
   );
   if (debug == 2) {
      return(p9);
   }
   p10 <- p9;


   ####################################################
   ## Optionally add centroid lines to the origin
   if (jamba::igrepHas("sample|centroid", drawVectors)) {
      if (verbose) {
         jamba::printDebug("bgaPlotly3d(): ",
            "Adding vectors for sample centroids.");
         print(head(dfCVLDF, 20));
      }
      p10 <- p10 %>% plotly::add_trace(
         name="CentroidLines",
         type="scatter3d",
         x=dfCVLDF[[Xsc]],
         y=dfCVLDF[[Ysc]],
         z=dfCVLDF[[Zsc]],
         mode="lines",
         # text=~dfSampleLinesDF$Label,
         # textposition="top center",
         line=list(width=centroidLwd,
            color=dfCVLDF$color),
         opacity=centroidAlpha
         #hoverlabel=list(
         #   bgcolor=dfSampleLinesDF$color,
         #   bordercolor=jamba::setTextContrastColor(dfSampleLinesDF$color),
         #   font=list(
         #      color=jamba::setTextContrastColor(dfSampleLinesDF$color))),
         #marker=list(
         #   size=dfSampleLinesDF$size,
         #   color=dfSampleLinesDF$color),
         #textfont=list(
         #   color=dfSampleLinesDF$color,
         #   size=(40-dfSampleLinesDF$size)/2)
      );
   }


   ####################################################
   ## Optionally add gene lines to the origin
   if (jamba::igrepHas("gene|row", drawVectors)) {
      if (verbose) {
         jamba::printDebug("bgaPlotly3d(): ",
            "Adding vectors for gene rows.");
         print(head(dfVgLDF, 20));
      }
      p10 <- p10 %>% plotly::add_trace(
         type="scatter3d",
         x=dfVgLDF[[Xg]],
         y=dfVgLDF[[Yg]],
         z=dfVgLDF[[Zg]],
         mode="markers+lines", #"markers+lines+text"
         text=~dfVgLDF$Label,
         line=list(width=geneLwd,
            color=dfVgLDF$color),
         opacity=geneAlpha,
         hoverlabel=list(
            bgcolor=dfVgLDF$color,
            bordercolor=jamba::setTextContrastColor(dfVgLDF$color),
            font=list(
               color=jamba::setTextContrastColor(dfVgLDF$color))),
         marker=list(
            size=dfVgLDF$size,
            color=dfVgLDF$color),
         textfont=list(
            color=dfVgLDF$color,
            size=(40-dfVgLDF$size)/4)
      );
   }


   ####################################################
   ## Create ellipses as 3d surfaces
   if (!ellipseType %in% "none") {
      sampleGroupsL <- split(names(sampleGroups), sampleGroups);
      if (verbose) {
         jamba::printDebug("bgaPlotly3d(): ",
            "Creating sample centroid ellipsoids.");
      }
      ###################################
      ## Iterate each group
      p11 <- p10;
      for (iGroup in names(sampleGroupsL)) {
         if (verbose) {
            jamba::printDebug("   bgaPlotly3d(): ",
               "Creating ellipsoid for group:",
               iGroup,
               ".");
         }
         iGroupSamples <- sampleGroupsL[[iGroup]];
         dfSamplesDFsub <- dfSamplesDF[iGroupSamples,,drop=FALSE];
         bgaS2coords <- dfSamplesDFsub[,c(Xs,Ys,Zs),drop=FALSE];
         ## Expand the points with radial jitter
         ## to help the ellipsoid calculation
         if (nrow(bgaS2coords) <= 150) {
            if (verbose) {
               jamba::printDebug("   bgaPlotly3d(): ",
                  "jitter_norm()");
            }
            bgaS2coords <- dfSamplesDF[rep(iGroupSamples, each=100),c(Xs,Ys,Zs),drop=FALSE];
            for (k in 1:3) {
               ## Use the range along each axis to impart variance relative to the visual noise space
               newC <- jitter_norm(bgaS2coords[,k],
                  amount=diff(range(dfSamplesDF[,c(Xs,Ys,Zs)[k]]))/50);
               bgaS2coords[,k] <- newC;
            }
         }
         if (ellipseType %in% "ellipsoid") {
            ## Return a mesh3d object
            ## formerly rgl::ellipse3d()
            #iEllipse <- ellipse3d.default(cov(bgaS2coords)*1,
            iEllipse <- rgl::ellipse3d(cov(bgaS2coords)*1,
               centre=colMeans(dfSamplesDFsub[,c(Xsc,Ysc,Zsc),drop=FALSE]),
               subdivide=3,
               smooth=TRUE,
               level=0.95);
         } else {
            iEllipse <- list(vb=t(bgaS2coords));
         }
         dfSampleLinesDFsub <- dfSampleLinesDF[match(iGroupSamples, dfSampleLinesDF$Name),,drop=FALSE];
         iGroupColor1 <- jamba::nameVector(
            dfSampleLinesDFsub[,"color"],
            iGroupSamples);
         iGroupColor <- head(iGroupColor1, 1);
         if (verbose) {
            printDebug("   bgaPlotly3d(): ",
               "hull name:",
               paste0(iGroup,
                  " group hull (", closestRcolor(iGroupColor), ")"));
         }
         iName <- head(paste0(iGroup,
            " group hull (",
            colorjam::closestRcolor(iGroupColor),
            ")"), 1);
         iEllipseDF <- data.frame(x=iEllipse$vb[1,],
            y=iEllipse$vb[2,],
            z=iEllipse$vb[3,],
            name=factor(rep(iGroup, length(iEllipse$vb[2,])),
               levels=levels(sampleGroups)),
            color=rep(iGroupColor,
               length(iEllipse$vb[2,]))
         );
         ## Add to the plotly object
         p11 <- p11 %>% plotly::add_trace(
            data=iEllipseDF,
            name=iGroup,
            type="mesh3d",
            showlegend=TRUE,
            inherit=FALSE,
            legendgroup="ellipsoids",
            hoverinfo="text",
            text=iGroup,
            hoverlabel=list(
               bgcolor=iGroupColor,
               bordercolor=jamba::setTextContrastColor(iGroupColor),
               font=list(
                  color=jamba::setTextContrastColor(iGroupColor))),
            flatshading=FALSE,
            #vertexcolor=iEllipseDF$color,
            vertexcolor=~color,
            opacity=ellipseOpacity,
            x=~x,
            y=~y,
            z=~z,
            alphahull=0);
      }
      p10 <- p11;
   }


   ############################################################
   ## superGroups
   if (length(superGroups) > 0) {
      if (verbose) {
         printDebug("bgaPlotly3d(): ",
            "Processing supergroups.");
      }
      if (length(names(superGroups)) == 0 &&
            length(superGroups) == length(sampleGroups)) {
         names(superGroups) <- names(sampleGroups);
      }
      if (all(names(superGroups) %in% names(sampleGroups))) {
         if (useScaledCoords) {
            superGroupsDF <- data.frame(superGroups=superGroups,
               sampleGroups=sampleGroups[names(superGroups)],
               color=sampleColors[names(superGroups)],
               groupcoords_scaled[as.character(sampleGroups[names(superGroups)]),c(Xsc,Ysc,Zsc)]);
         } else {
            superGroupsDF <- data.frame(superGroups=superGroups,
               sampleGroups=sampleGroups[names(superGroups)],
               color=sampleColors[names(superGroups)],
               groupcoords_unscaled[as.character(sampleGroups[names(superGroups)]),c(Xsc,Ysc,Zsc)]);
         }
         superGroupDF <- unique(superGroupsDF);
         superGroupDF <- jamba::mixedSortDF(superGroupDF);
         if (verbose) {
            printDebug("bgaPlotly3d(): ",
               "head(superGroupDF, 30):");
            print(head(superGroupDF, 30));
         }
         for (iDF in split(superGroupDF, superGroupDF$superGroups)) {
            if (nrow(iDF) <= 1) {
               next;
            }
            iSuperGroup <- iDF[1,"superGroups"];
            if (verbose) {
               printDebug("   bgaPlotly3d(): ",
                  "head(iDF, 30):");
               print(head(iDF, 30));
            }
            ## Iterate each set of centroids by making a 3d spline
            ## between centroids
            #arrowSmoothFactor <- 4;
            iDFspline <- spline3d(x=iDF[,Xsc],
               y=iDF[,Ysc],
               z=iDF[,Zsc],
               lengthFactor=arrowSmoothFactor);
            colnames(iDFspline) <- c(Xsc, Ysc, Zsc);
            ## Smoothing
            iDFnew <- iDF[rep(1:nrow(iDF), each=arrowSmoothFactor),];
            iDFnew[,c(Xsc, Ysc, Zsc)] <- iDFspline;
            iDFnew$newColor <- colorRampPalette(iDF$color, alpha=TRUE)(nrow(iDFnew))
            iDFnew <- jamba::renameColumn(iDFnew,
               from=c(Xsc, Ysc, Zsc),
               to=c(Xs, Ys, Zs));
            ## Now add to plotly object
            p10 <- p10 %>% plotly::add_trace(
               data=iDFnew,
               name=iSuperGroup,
               legendgroup="supergroups",
               inherit=FALSE,
               type="scatter3d",
               text=iSuperGroup,
               mode="lines",
               hoverinfo="text",
               line=list(color=iDFnew$newColor,
                  width=superGroupLwd),
               opacity=superGroupAlpha,
               x=as.formula(paste0("~",Xs)),
               y=as.formula(paste0("~",Ys)),
               z=as.formula(paste0("~",Zs)));
         }
      } else {
         if (verbose) {
            printDebug("bgaPlotly3d(): ",
               "names(superGroups) were not all present in names(bgaInfo$fac).");
         }
      }
   }

   ## Optionally adjust the 3D orientation
   scene <- list(
      camera=list(
         up=list(
            x=0,
            y=1,
            z=0),
         eye=list(
            x=sceneX,
            y=sceneY,
            z=sceneZ),
         xaxis=list(
            tickwidth=2,
            gridwidth=4),
         zaxis=list(
            gridwidth=5),
         yaxis=list(
            gridwidth=5)
         ),
         aspectmode="data"
      );
   p10 <- p10 %>% plotly::layout(title=main,
      dragmode="orbit",
      scene=scene,
      paper_bgcolor=paper_bgcolor,
      plot_bgcolor=plot_bgcolor
      );

   return(p10);
}

#' Convert data.frame to plotly line segment format
#'
#' Convert data.frame to plotly line segment format
#'
#' Purpose is to take a data.frame containing (x,y) or (x,y,z) coordinates
#' for two points defined on one row, and create a data.frame with 3 rows
#' where blank coordinates are used to separate each line segment.
#' The main benefit in using this function is that it enables vectorized
#' line segment drawing when used with a plotting function that recognizes
#' line breaks via empty rows, for example plotly.
#'
#' @return data.frame where every two rows represents a line segment
#' to be drawn by plotly, followed by a blank line indicating the line
#' ends.
#'
#' @param x numeric matrix of coordinates, where x,y,z coordinates are
#'    present twice per row, and whose colnames are defined by
#'    `axes1` and `axes2`.
#' @param axes1 vector of colnames(x) corresponding to the x,y,z coordinates
#'    of the first point in a line segment. It can be of any length, but
#'    is intended for length 2 or 3.
#' @param axes2 vector of colnames(x) corresponding to the x,y,z coordinates
#'    of the second point in a line segment. It must be the same length as
#'    `axes1`.
#' @param ... additional parameters are ignored.
#'
#' @family jam spatial functions
#'
#' @export
dfWide2segments <- function
(x,
 axes1,
 axes2,
 ...)
{
   ## Purpose is to take a data.frame containing (x,y) or (x,y,z) coordinates
   ## for two points defined on one row, and create a data.frame with 3 rows
   ## where blank coordinates are used to separate each line segment.

   ## First make sure axes1 and axes2 are colnames(x), but
   ## allow for re-using the same colname
   axes1 <- axes1[axes1 %in% colnames(x)];
   axes2 <- axes2[axes2 %in% colnames(x)];
   if (length(axes1) != length(axes2)) {
      stop("dfWide2segments() requires axes1 and axes2 have identical length.");
   }
   axes1v <- jamba::nameVector(seq_along(axes1),
      axes1);
   segmentsDF <- as.data.frame(do.call(cbind,
      lapply(axes1v, function(i){
         intercalate(x[,axes1[i]],
            x[,axes2[i]],
            NA)
      })
   ));
   return(segmentsDF);
}

#' Intercalate list values into a vector
#'
#' Intercalate list values into a vector
#'
#' This function takes a list of vectors, and intercalates them
#' by inserting alternating values from each vector to make one
#' combined vector. It is analogous to shuffling a deck of cards,
#' except that there could be more than two stacks of cards shuffled
#' together.
#'
#' @return vector of values derived from each input list `...`, with
#' one element per list in order. In the event any list is shorter than
#' the others, its values are recycled so each list is equal length.
#'
#' @param ... parameters are expected as multiple vectors.
#'
#' @family jam list functions
#'
#' @examples
#' A <- LETTERS[1:10];
#' B <- letters[1:10];
#' intercalate(A, B)
#'
#' @export
intercalate <- function
(...)
{
   ## Purpose is to take a list of vectors, and intercalate their values, e.g.
   ## list1 <- paste("name1", letters[1:10], sep="");
   ## list2 <- paste("name2", letters[1:10], sep="");
   ## intercalate(list1, list2);
   ## name1a, name2a, name1b, name2b, name1c, name2c, etc.
   ##
   ## The special case where there are two lists, and the first has
   ## one element more than the second, then the second will only have
   ## its values in between the first, e.g.
   ## A B A B A B A
   ##
   ## Note: rmNULL() will remove empty lists
   aList <- jamba::rmNULL(list(...));
   if (length(aList) == 1 && class(aList[[1]]) %in% "list") {
      aList <- aList[[1]];
   }
   ## do.call will automatically repeat any vector to fill each row
   ## up to the maximum number of columns.
   if (length(unique(lengths(aList))) > 1) {
      ## Unequal lengths, to avoid warning should we expand them?
   }
   aMatrix <- do.call(rbind, aList);
   newVector <- as.vector(aMatrix);

   ## The special case where intercalating two vectors,
   ## where the second vector has one fewer entry, we
   ## will not repeat the last entry.
   ## E.g.
   ## c("A","A","A")
   ## c("B","B")
   ##
   ## desired output is
   ## c("A","B","A","B","A")
   if (length(aList) == 2 && length(aList[[1]]) == (length(aList[[2]]) + 1)) {
      newVector <- head(newVector, -1);
   }
   return(newVector);
}

#' Calculate a spline curve fit in 3-D
#'
#' Calculate a spline curve fit in 3-D
#'
#' This function takes 3-dimensional points as input, and fits
#' a 3-D spline curve to connect the points. It essentially calls
#' `stats::spline()` on one dimension at a time, relative to the
#' row rank.
#'
#' @param x numeric vector or numeric `matrix` or `data.frame` with `colnames(x)`
#'    containing `c("x","y","z")`.
#' @param y,z optional numeric vectors, used only when x is also a numeric
#'    vector.
#' @param lengthFactor numeric multiplier used to control the number
#'    of intermediate points during smoothing between each provided point.
#' @param verbose logical indicating whether to print verbose output.
#'
#' @return
#' numeric matrix with `colnames(x)=c("x","y","z")`.
#'
#' @family jam spatial functions
#'
#' @export
spline3d <- function
(x,
 y=NULL,
 z=NULL,
 lengthFactor=4,
 verbose=FALSE,
 ...)
{
   ## Purpose is to provide 3D spline()
   if (jamba::igrepHas("matrix|data.frame", class(x))) {
      coordsT1 <- cbind(t=1:nrow(x), x);
      colnames(coordsT1)[1:4] <- c("t","x","y","z");
   } else {
      coordsT1 <- cbind(t=1:length(x), x=x, y=y, z=z);
   }
   ## Define the number of points to include
   ts <- seq(from=min(coordsT1[,"t"]),
      to=max(coordsT1[,"t"]),
      length.out=nrow(coordsT1)*lengthFactor);
   d2 <- apply(coordsT1[,-1], 2, function(u) {
      spline(coordsT1[,"t"], u, xout=ts)$y;
   });
   if (verbose) {
      printDebug("spline3d(): ",
         "head(d2, 20):")
      print(head(d2, 20));
   }
   return(d2);
}

#' Convert data.frame to categorical colors
#'
#' Convert data.frame to categorical colors
#'
#' This function is a temporary placeholder function which simply
#' provides categorical colors for values in each column of the input
#' `data.frame`. If colors are provided using the named color vector
#' `colorSub`, they are used to substitute values in each column,
#' otherwise `colorjam::group2colors()` is used to create categorical
#' colors, which itself uses `colorjam::rainbowJam()`.
#'
#' Note that providing `colorSub` helps keep colors consistent, otherwise
#' each column is independently colorized.
#'
#' @return data.frame of the same dimensions as the input `x` data.frame,
#' where values have been substituted with R colors.
#'
#' @param x data.frame
#' @param colorSub vector of colors, whose names are intended to match
#'    values in the input data.frame `x`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional parameters are passed to `colorjam::group2color()`.
#'
#' @family jam color functions
#'
#' @examples
#' colorSub1 <- colorjam::group2colors(LETTERS[1:6]);
#' df <- data.frame(one=LETTERS[1:6],
#'    two=rep(LETTERS[3:4], each=3),
#'    three=rep(LETTERS[5:6], 3));
#' dfColors <- df2colorSub(df, colorSub1);
#' jamba::imageByColors(dfColors, cellnote=df);
#'
#' # Or in one step
#' jamba::imageByColors(df2colorSub(df), cellnote=df);
#'
#' @export
df2colorSub <- function
(x,
 colorSub=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is a simple intermediate replacement for df2groupColors()
   ## for the special case where all values per data.frame column are
   ## categorical. Colors are assigned by matching with names(colorSub)
   ## where possible, then group2colors() otherwise.
   x2 <- as.data.frame(do.call(cbind, lapply(seq_len(ncol(x)), function(i) {
      j <- x[[i]];
      colorjam::group2colors(j, colorSub=colorSub, ...);
   })));
   rownames(x2) <- rownames(x);
   colnames(x2) <- colnames(x);
   x2;
}

#' Apply jitter using normal distribution
#'
#' Apply jitter using normal distribution
#'
#' This function applies a jitter (noise) to points by adding
#' random values from a normal distribution, where the argument `sd`
#' is used to apply the jitter magnitude.
#'
#' The default jitter is defined as 1/10 the median difference between
#' unique, finite input values, which can be scaled using the argument
#' `factor`. Note the use of "unique" input values, which ensures the
#' presence of duplicate values does not skew the jitter toward zero.
#' When applied to three dimensions, this results in jitter consistently
#' scaled relative to the range of values in each dimension. That is,
#' points will appear to have a radial jitter.
#'
#' @family jam plot functions
#'
#' @param x numeric vector
#' @param factor numeric value to define the magnitude of jitter,
#'    multiplied by the default jitter which is 1/10 the median
#'    difference between unique values in `x`.
#' @param amount optional numeric value indicating a fixed `sd`
#'    standard deviation passed to `stats::rnorm()`.
#' @param ... additional arguments are ignored.
#'
#' @export
jitter_norm <- function
(x,
 factor=1,
 amount=NULL,
 ...)
{
   ## Determine the amount of jitter relative to the median distance
   ## between points
   if (length(amount) == 0) {
      x <- x[is.finite(x)];
      z <- diff(range(x));
      if (z == 0) {
         z <- abs(range(x)[1]);
      }
      if (z == 0) {
         z <- 1;
      }
      xx <- unique(sort(round(x, 3 - floor(log10(z)))));
      if (length(xx) == 1) {
         dAmount <- xx / 10;
      } else {
         dAmount <- median(diff(xx));
      }
      if (length(dAmount) == 0 || dAmount == 0) {
         dAmount <- z / 10;
      }
      amount <- (factor / 5) * dAmount;
   }

   ## Choose distances from normal distribution
   xDists <- rnorm(length(x), sd=amount * 0.6);
   xValue <- x + xDists;
   return(xValue);
}
