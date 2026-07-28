
# link plotly sjfig panels

#' Add links to Splicejam Figure for Plotly
#' 
#' Add links to Splicejam Figure for Plotly, experimental
#' approach to add some type of link_id to be used
#' to link multiple panels.
#' 
#' @keywords internal
#' @noRd
#' @examples
#' data(sjenvtest)
#' sjfig <- splicejamFigure(sjenv=sjenvtest, gene="Gria1", use_memoise=TRUE)
#' sjfig2 <- link_sjfig_plotly(sjfig)
link_sjfig_plotly <- function
(sjfig,
 ...)
{
   # Grab data from combined panels
   idf1 <- jamba::rbindList(lapply(sjfig$cp_sashimi_list, function(gg){
      gg$data
   }))
   # idf1 <- sjfig$cp_sashimi$data;
  
   # adjust coverage coordinates
   idf1c <- subset(idf1, type %in% "coverage" & !grepl("^gap", feature))
   idf1c_xrall <- sapply(split(idf1c$x, factor(idf1c$name)), function(i){
      irange <- range(unlist(i));
      addto <- c(0, 0);
      jamba::cPaste(irange + addto, sep=":")
   })
   idf1c_xr <- idf1c_xrall[!duplicated(idf1c_xrall)]
  
   # adjust coverage gap coordinates
   idf1cg <- subset(idf1, type %in% "coverage" & grepl("^gap", feature))
   idf1cg_xrall <- sapply(split(idf1cg$x, factor(idf1cg$name)), function(i){
      irange <- range(unlist(i));
      # addto <- c(0.5, -1.5)
      addto <- c(-1, 2);
      # addto <- c(0, 0);
      jamba::cPaste(irange + addto, sep=":")
   })
   idf1cg_xr <- idf1cg_xrall[!duplicated(idf1cg_xrall)]
   idf1c_xr[names(idf1cg_xr)] <- idf1cg_xr;
   
   # adjust junction coordinates
   idf1j <- subset(idf1, type %in% "junction")
   idf1j_xr <- sapply(split(idf1j$x, factor(idf1j$name)), function(i){
      irange <- range(unlist(i));
      # addto <- c(0.5, -1.5)
      addto <- c(1.5, -1.5)
      jamba::cPaste(irange + addto, sep=":")
   })
   
   # get gene model coordinates
   idf3 <- sjfig$cp_gene$data;
   idf3gaps_xr <- sapply(split(idf3gaps$x, factor(idf3gaps$gr_name)), function(i){
      jamba::cPaste(range(i) + c(-0, 0), sep=":")
   })
   idf3gaps_xr2 <- rev(rev(idf3gaps_xr)[!duplicated(rev(idf3gaps_xr))])
   
   # merge one data.frame
   u1c <- sort((idf1c_xr));
   u1cV1 <- u1c
   u1j <- sort((idf1j_xr));
   u3 <- sort((idf3gaps_xr2));
   umerge <- jamba::mergeAllXY(
      data.frame(xrange=u1c, coverage=names(u1c)),
      data.frame(xrange=u1j, junction=names(u1j)),
      data.frame(xrange=u3, gene=names(u3)));
   # umerge
   
   # rescue gap with no matching gene region,
   # which are within 2 nt of the next gene feature
   ufix <- which(!is.na(umerge$coverage) & is.na(umerge$gene));
   ufix <- ufix[ufix < nrow(umerge)];
   um1 <- as.numeric(gsub(":.+", "", umerge$xrange))
   um2 <- as.numeric(gsub("^.+:", "", umerge$xrange))
   if (length(ufix) > 0) {
      ufp1 <- abs(um1[ufix] - um1[ufix + 1]) <= 2;
      ufp2 <- abs(um2[ufix] - um2[ufix + 1]) <= 2;
      ufp3 <- !is.na(umerge$gene[ufix + 1]);
      ufp12 <- ufp1 & ufp2 & ufp3;
      if (any(ufp12)) {
         um1[ufix[ufp12]] <- um1[ufix[ufp12] + 1]
         um2[ufix[ufp12]] <- um2[ufix[ufp12] + 1]
         um12 <- paste0(um1[ufix[ufp12]], ":", um2[ufix[ufp12]])
         umerge$xrange[ufix[ufp12]] <- um12;
         names(um12) <- umerge$coverage[ufix[ufp12]];
         u1c[names(um12)] <- um12;
      }
   }
   umerge2 <- (jamba::mergeAllXY(
      data.frame(xrange=u1c, coverage=names(u1c)),
      data.frame(xrange=u1j, junction=names(u1j)),
      data.frame(xrange=u3, gene=names(u3))))
   u1c_changed <- which(u1c != u1cV1)
   if (length(u1c_changed) > 0) {
      #
      cg1 <- which(idf1cg_xrall %in% u1cV1[u1c_changed]);
      cgmatch <- match(idf1cg_xrall[cg1], u1cV1);
      idf1cg_xrall[cg1] <- u1c[cgmatch];
   }

   # now add to all ggplot entries
   idf1change <- c(idf1c_xrall, idf1cg_xrall, idf1j_xr);
   # table(idf1$name %in% names(idf1change))
   # table(subset(idf1, !name %in% names(idf1change))$type)
   # head(subset(idf1, !name %in% names(idf1change) & type %in% "junction"))
   # head(subset(idf1, !name %in% names(idf1change) & type %in% "coverage"))
   idf1$link_id <- idf1change[as.character(idf1$name)];
   # propagate to each cp plot
   for (i in seq_along(sjfig$cp_sashimi_list)) {
      idf5 <- sjfig$cp_sashimi_list[[i]]$data;
      idf5$link_id <- idf1change[as.character(idf5$name)];
      sjfig$cp_sashimi_list[[i]]$data <- idf5;
   }

   # now add to gene model
   for (j in seq_along(sjfig$cp_genes)) {
      idf7 <- sjfig$cp_genes[[j]]$data;
      idf7$link_id <- idf3gaps_xr[as.character(idf7$gr_name)];
      sjfig$cp_genes[[j]]$data <- idf7;

   }
   # idf3gaps_xr
  
   # re-assemble panels
   # convert spacers to plotly if needed
   sjfig$cp_spacers <- lapply(sjfig$cp_spacers, function(cp_spacer){
      if (!inherits(cp_spacer, "plotly")) {
         plotly::ggplotly(cp_spacer + ggplot2::theme_void())
      } else {
         cp_spacer
      }
   })
   cp_list <- c(sjfig$cp_sashimi_list,
      sjfig$cp_spacers,
      sjfig$cp_genes)
   plotlys <- lapply(cp_list, function(icp){
      jam_ggplotly(icp)
   })
   #
   sjfig$cp <- plotly::subplot(plotlys,
      shareX=TRUE,
      shareY=FALSE,
      nrows=sjfig$layout_nrow) %>%
      plotly::layout(
         margin=list(t=60, b=50, l=80, r=30));
   
   return(invisible(sjfig));
}
