
#' Convert list to incidence matrix
#'
#' Convert list to incidence matrix
#'
#' This function uses some very efficient methods from the `arules`
#' package to convert a list of items into an incidence matrix.
#' The conversion creates a `transactions` object, which takes
#' advantage of the efficient code, then extracts the matrix from
#' that object. Unfortunately, that object class does not permit
#' counting each item, though it may be possible through another
#' method.
#'
#' The incidence matrix will have rownames matching list elements,
#' and colnames match the list names.
#'
#' In future the uniqueness step may be carried out using
#' `S4Vectors` list methods, which are also highly efficient
#' specifically in handling very large lists.
#'
#' @param x list of character vectors.
#' @param makeUnique logical indicating whether to enforce uniqueness on
#'    each vector in the list. If `makeUnique=FALSE` and there are
#'    duplicated values, then the resulting incidence matrix will be
#'    updated to reflect the count of each element.
#' @param verbose logical indicating whether to print verbose output.
#'
#'
#' @export
list2im <- function
(x,
 makeUnique=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to convert a list of vectors into an incident matrix
   ## using the arules package
   if (!suppressPackageStartupMessages(require(arules))) {
      stop("The arules package is required by list2im().");
   }
   if (makeUnique) {
      x <- lapply(x, unique);
   } else {
      xCt <- lapply(x, tcount);
      x <- lapply(x, unique);
   }
   if (verbose) {
      printDebug("list2im(): ",
         "Converting to '", "transactions", "' object.");
   }
   xT <- as(x, "transactions");

   if (verbose) {
      printDebug("list2im(): ",
         "Converting to '", "matrix", "' object.");
   }
   xM <- as.matrix(xT@data*1);
   colnames(xM) <- xT@itemsetInfo[,1];
   rownames(xM) <- xT@itemInfo[,1];
   if (!makeUnique && any(unlist(xCt) > 1)) {
      if (verbose) {
         printDebug("list2im(): ",
            "Applying item counts to the incidence matrix for ",
            format(big.mark=",", length(xCt)), " items.");
      }
      for (i in seq_along(xCt)) {
         if (any(xCt[[i]] > 1)) {
            xM[names(xCt[[i]]),i] <- xCt[[i]];
         }
      }
   }
   return(xM);
}
