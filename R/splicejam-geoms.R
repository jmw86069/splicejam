#' splicejam extensions to ggforce
#'
#' ggforce makes heavy use of the ggproto class system to extend the
#' functionality of ggplot2. In general the actual classes should be of little
#' interest to users as the standard ggplot2 api of using geom_* and stat_*
#' functions for building up the plot is encouraged.
#'
#' @name splicejam-extensions
#' @rdname splicejam-extensions
#'
#' @family jam ggplot2 functions
#'
NULL

#' Draw an area defined by an upper and lower diagonal into an arc
#'
#' The `geom_diagonal_wide_arc()` function draws a *thick* diagonal, that is, a
#' polygon confined between a lower and upper [diagonal][geom_diagonal]. As with
#' the diagonal functions in `ggforce`, the wide diagonal variant is horizontal.
#' This function joins two adjacent diagonals into one wider arc.
#'
#' @section Aesthetics:
#' geom_diagonal_wide_arc understand the following aesthetics
#' (required aesthetics are in bold):
#'
#' - **x**
#' - **y**
#' - **group**
#' - color
#' - size
#' - linetype
#' - alpha
#' - lineend
#'
#' @inheritParams ggforce::geom_shape
#' @inheritParams ggplot2::stat_identity
#'
#' @param n The number of points to create for each of the bounding diagonals
#'
#' @param strength The proportion to move the control point along the x-axis
#' towards the other end of the bezier curve
#'
#' @name geom_diagonal_wide_arc
#' @rdname geom_diagonal_wide_arc
#'
#' @family jam ggplot2 functions
#'
#' @examples
#' data <- data.frame(
#'   x = c(1, 2, 2, 1, 2, 3, 3, 2),
#'   y = c(1, 2, 3, 2, 2, 1, 2, 3),
#'   group = c(1, 1, 1, 1, 1, 1, 1, 1)
#' )
#'
#' ggplot(data) +
#'   geom_diagonal_wide_arc(aes(x, y, group=group))
#'
#' # The strength control the steepness
#' ggplot(data, aes(x, y, group = group)) +
#'   geom_diagonal_wide_arc(strength=0.75, alpha=0.5, fill='red') +
#'   geom_diagonal_wide_arc(strength=0.25, alpha=0.5, fill='blue')
#'
#' # The diagonal_wide_arc geom uses geom_shape under the hood, so corner rounding
#' # etc are all there
#' ggplot(data) +
#'   geom_diagonal_wide_arc(aes(x, y, group=group), radius=unit(5, 'mm'))
#'
NULL

#' @rdname splicejam-extensions
#' @format NULL
#' @usage NULL
#'
#' @family jam ggplot2 functions
#'
#' @export
StatDiagonalWideArc <- ggproto('StatDiagonalWideArc', Stat,
   setup_data = function(data, params) {
      if (any(!table(data$group) %in% c(8))) {
         stop('Each group must consist of 8 points')
      }
      data
   },
   compute_panel = function
   (data,
    scales,
    strength=0.5,
    n=100)
   {
      # Keep original order of groups
      data$group_factor <- factor(data$group,
         levels=(unique(data$group)));

      data <- data[order(data$group_factor, data$x, data$y),]
      group_data <- data$group;
      new_group <- rep(
         jamba::makeNames(
            rep(unique(group_data), each=2)),
         each=2);
      new_group <- factor(new_group,
         levels=unique(new_group));
      lower_logic <- c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE);
      lower <- data[lower_logic, ]
      lower$group <- new_group;
      upper <- data[!lower_logic, ]
      upper$group <- new_group;
      lower <- ggforce:::add_controls(lower,
         strength);
      upper <- ggforce:::add_controls(upper[rev(seq_len(nrow(upper))), ],
         strength)
      lower <- ggforce::StatBezier$compute_layer(lower, list(n=n))
      upper <- ggforce::StatBezier$compute_layer(upper, list(n=n))
      lower$group <- as.integer(gsub("_v[12]$", "", levels(new_group)[lower$group]));
      upper$group <- as.integer(gsub("_v[12]$", "", levels(new_group)[upper$group]));
      lower <- lower[order(lower$group, lower$x),];
      upper <- upper[order(upper$group, -upper$x),];
      diagonals <- rbind(lower, upper);
      #diagonals$index <- NULL
      diagonals[order(diagonals$group_factor),,drop=FALSE];
   },
   required_aes = c('x', 'y', 'group'),
   extra_params = c('na.rm', 'n', 'strength')
)

#' @rdname geom_diagonal_wide_arc
#'
#' @family jam ggplot2 functions
#'
#' @export
stat_diagonal_wide_arc <- function
(mapping=NULL,
 data=NULL,
 geom='shape',
 position='identity',
 n=100,
 strength=0.5,
 na.rm=FALSE,
 show.legend=NA,
 inherit.aes=TRUE,
 ...)
{
   layer(
      stat = StatDiagonalWideArc,
      data = data,
      mapping = mapping,
      geom = geom,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
         na.rm = na.rm,
         n = n,
         strength = strength,
         ...)
   )
}

#' @rdname geom_diagonal_wide_arc
#'
#' @family jam ggplot2 functions
#'
#' @export
geom_diagonal_wide_arc <- function
(mapping = NULL,
 data = NULL,
 stat = 'diagonal_wide_arc',
 position = 'identity',
 n = 100,
 na.rm = FALSE,
 strength = 0.5,
 show.legend = NA,
 inherit.aes = TRUE,
 ...)
{
   layer(
      data = data,
      mapping = mapping,
      stat = stat,
      #geom = ggforce::GeomShape,
      geom = GeomPolygon,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
         na.rm = na.rm,
         n = n,
         strength = strength,
         ...)
   )
}
