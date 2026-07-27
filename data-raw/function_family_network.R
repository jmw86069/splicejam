#
# R code to load pkgdown family names,
# parse the package family tags,
# then create igraph network to visualize.


# get order currently in '_pkgdown.yml'
pkgdownorder <- gsub("^.+[(].|.[)].*$", "",
   jamba::vigrep("has.concept", readLines("_pkgdown.yml")))
pkgdownorder


# split into list by family
fn_list <- split(fn_family[,2], fn_family[,1])

# check missing names
setdiff(names(fn_list), pkgdownorder)
# internal
setdiff(pkgdownorder, names(fn_list))
# splicejam data

fn_order <- fn_list[pkgdownorder];
names(fn_order) <- pkgdownorder;
fn_order



## Use roxygen2 to parse the package for family assignments
blocks <- roxygen2::parse_package(".")
# parse block function, data, or rdname
block_names <- sapply(seq_along(blocks), function(bnum) {
   b <- blocks[[bnum]];
   use_name <- b$object$topic;
   if (length(use_name) == 0) {
      use_name <- roxygen2::block_get_tag(b, "rdname")$val
   }
  use_name
})
# names(blocks) <- block_names;

family_dfs <- (lapply(seq_along(blocks), function(bnum) {
   bname <- block_names[[bnum]];
   jamba::printDebug(bnum);
   b <- blocks[[bnum]];
   btags <- roxygen2::block_get_tags(b, tags="family")
   bfamily <- unique(sapply(btags, function(i){
      gsub("\n.*", "", i$val)
   }))
   if (length(bfamily) == 0) {
     data.frame(family=character(0), bname=character(0))
   } else {
      data.frame(
         family=bfamily,
         name=rep(bname, length.out=length(bfamily)))
   }
}))
family_df <- jamba::mixedSortDF(unique(jamba::rbindList(family_dfs)))
family_df[, 1] <- factor(family_df[, 1],
   levels=unique(c(pkgdownorder, family_df[, 1])))
jamba::mixedSortDF(family_df)

family_el <- as.matrix(family_df)
g <- igraph::graph_from_edgelist(family_el, directed=FALSE)
igraph::V(g)$size <- 5;
igraph::V(g)$nodeType <- ifelse(igraph::V(g)$name %in% family_el[, 2], "Gene", "Set")
igraph::V(g)$size <- c(Set=8, Gene=5)[igraph::V(g)$nodeType];
igraph::V(g)$color <- colorjam::group2colors(igraph::V(g)$nodeType)
igraph::V(g)$frame.color <- jamba::makeColorDarker(igraph::V(g)$color)
igraph::V(g)$font <- 2;
g <- multienrichjam::set_igraph_layout(g=g,
   layout=multienrichjam::layout_with_qfrf(repulse=3.1))


options("jam.shadow.n"=16, "jam.shadow.r"=0.1)


multienrichjam::jam_igraph(g,
   nodegroups=g_nodesets,
   edge_bundling="none",
   mark.alpha=0.05,
   mark.expand=3,
   use_shadowText=TRUE,
   node_factor_l=list(nodeType=c(Gene=1, Set=1)),
   label_factor_l=list(nodeType=c(Gene=1.2, Set=2)),
   mark.groups=unname(use_groups_all))


## Alternative using Sugiyama hierarchical layout
# igraph::plot.igraph(g, layout=igraph::layout_with_sugiyama(g)$layout[,2:1])
sugi_layout <- igraph::layout_with_sugiyama(g)$layout[,2:1];
sugi_layout[, 1] <- scale(sugi_layout[, 1]) * 5 + rnorm(nrow(sugi_layout)) / 100
sugi_layout[, 2] <- scale(sugi_layout[, 2]) * 10
igraph::V(g)$label <- gsub(" functions$", "", igraph::V(g)$name)

g_nodesets <- multienrichjam::get_cnet_nodeset(g)
use_groups_all <- lapply(jamba::nameVector(igraph::V(g)$name[igraph::V(g)$nodeType %in% "Set"]), function(v){
   c(v, igraph::neighbors(g, v=v)$name)
})


multienrichjam::jam_igraph(g, layout=sugi_layout,
   nodegroups=g_nodesets,
   mark.alpha=0.05,
   use_shadowText=TRUE,
   node_factor_l=list(nodeType=c(Gene=0.4, Set=1)),
   label_factor_l=list(nodeType=c(Gene=0.8, Set=1)),
   mark.groups=unname(use_groups_all))


################################################
#
# New function categories, see below:
#
################################################
## Splicejam core
# sashimiDataConstants()
# splicejamFigure()
# launchSashimiApp()


## Detected transcripts
# defineDetectedTx
# detectedTxInfo


## Data import functions
# getGRcoverageFromBw
# import_juncs_from_bed
# psl2df


## Splicejam data
# sjenvtest
# test_cov_gr
# test_cov_wide_gr
# test_exon_gr
# test_exon_wide_gr
# test_junc_gr
# test_junc_wide_gr



## Plot utility functions
# bgaPlotly3d()
# jitter_norm
# spline3d


## Design functions
# curateDFtoDF
# curateVtoDF
# group2contrasts
# runDiffSplice
# sortSamples


## GenomicRanges functions
# addGRgaps
# addGRLgaps
# getGRgaps
# getGRLgaps
# annotateGRfromGR
# annotateGRLfromGRL
# assignGRLexonNames
# findOverlapsGRL
# jam_isDisjoint
# flattenExonsBy
# sortGRL


## GTF functions
# describeGtfAttrNames
# getGtfAttrs
# makeTx2geneFromGtf
# readGtf


## ALE and codon functions
# ale2violin
# getFirstStrandedFromGRL
# tx2ale
# dna2codon
# codonUsage2df
# jamCai


## Shiny prep functions
# sashimiAppConstants
# sashimiAppServer
# sashimiAppUI


## Sashimi prep functions
# prepareSashimi
# plotSashimi
# gene2gg
# grl2df
# exoncov2polygon
# make_ref2compressed


## Splicejam ggplot2 customizations
# geom_diagonal_wide_arc
# stat_diagonal_wide_arc
# stat_unpack_polygon
# StatDiagonalWideArc
# to_basic.GeomShape


## Internal utility functions
# stackJunctions
# factor2label
# df2colorSub
# spliceGR2junctionDF
# closestExonToJunctions
# combineGRcoverage
# intercalate
# list2im
# shrinkMatrix
# geomean
# jamGeomean
# internal_junc_score
# compressPolygonM
# dfWide2segments
# simplifyXY
# escapeWhitespaceRegexp
# strsplitOrdered


## Candidate to make internal-only
# getFirstStrandedFromGRL
# factor2label
# stackJunctions
