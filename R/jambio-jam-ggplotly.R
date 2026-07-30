
#' Internal function to convert multi-plot ggplot to plotly panels
#' @keywords internal
#' @noRd
jam_ggplotly <- function(p, ...)
{
   # check p$data for list columns 'x', 'y'
   if (is.list(p$data$x)) {
      # jamba::printDebug("jam_ggplotly(): ",
      #    "head(p$data):");
      # print(head(p$data, 10));# debug
      # jamba::printDebug("jam_ggplotly(): ",
      #    "sclass(p$data):");
      # print(jamba::sclass(p$data));# debug
      # jamba::printDebug("jam_ggplotly(): ",
      #    "jamba::tcount(p$data$name):");
      # print(head(jamba::tcount(p$data$name), 20));# debug
      name_ct <- lengths(p$data$x);

      # insert 0 at beginning and end
      multixs <- which(name_ct > 1);
      for (multix in multixs) {
         p$data$x[[multix]] <- c(head(p$data$x[[multix]], 1),
            p$data$x[[multix]],
            tail(p$data$x[[multix]], 1));
         p$data$y[[multix]] <- c(0,
            p$data$y[[multix]],
            0);
      }
      name_ct <- lengths(p$data$x);

      use_rows <- rep(seq_along(name_ct), name_ct);
      new_data <- p$data[use_rows, , drop=FALSE];
      new_data$x <- unname(unlist(p$data$x))
      new_data$y <- unname(unlist(p$data$y))
     
      ## Insert 'text' column for plotly label
      # junction_label: (none)
      # junction: feature, score
      # coverage: x, y, feature
      use_text <- ifelse(new_data$type %in% "coverage",
         paste0("<b>Coverage:</b>\n",
            "Feature: '",
            new_data$feature, "'\n",
            "Score: ",
            jamba::formatInt(round(new_data$y)), "\n",
            "Coordinate: ",
            jamba::formatInt(round(new_data$x)), "\n"),
         ifelse(new_data$type %in% "junction",
            paste0("<b>Junction:</b>\n",
               "From: '", new_data$nameFrom, "'\n",
               "To:     '", new_data$nameTo, "'\n",
               new_data$feature, "\n",
               "Score: ",
               jamba::formatInt(round(new_data$score)), "\n"),
            ""))
      new_data$text <- use_text;
      
      # new_data <- plotly::highlight_key(new_data,
      #    key=~feature);
      p$data <- new_data;
      # add text to aesthetics
      p <- p + ggplot2::aes(text=text);
   }
   # cp <- plotly::ggplotly(p,
   #    tooltip=c("text"));
   cp <- plotly::layout(
      plotly::ggplotly(p,
         tooltip=c("text")),
      showlegend=FALSE)
   #
   cp$x$data <- lapply(cp$x$data, function(xdata){
      if (any(grepl("points|fills", xdata$hoveron))) {
         xdata$hoveron <- "fills+points";
      }
      xdata
   })
   #
   # cp <- plotly::highlight(cp,
   #    "plotly_hover",
   #    opacityDim=0.8,
   #    selected=plotly::attrs_selected(
   #       line=list(color="#444444")));
   cp;
}
