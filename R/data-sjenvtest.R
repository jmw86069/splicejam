
#' Splicejam Environment Test Data
#' 
#' Splicejam Environment Test Data, minimal subset of
#' Farris et al bulk RNA-seq data from mouse hippocampus.
#' It contains only 'CellType' values CA1 and CA2, with
#' both 'Compartment' values CB (cell body) and
#' DE (dendrites).
#' 
#' @family splicejam data
#' 
#' @format `environment` with components suitable as
#'    input to `splicejamFigure()` and other functions
#'    which accept environment as input, typically
#'    using argument 'sjenv'.
#' \describe{
#'   \item{tx2geneDF}{`data.frame` with colnames including
#'     the minimum required: 'gene_name', 'transcript_id'.}
#'   \item{detectedGenes}{`character` vector of gene symbols.}
#'   \item{detectedTx}{`character` vector of 'transcript_id'.}
#'   \item{flatExonsByGene}{`GRangesList` named by 'gene_name'}
#'   \item{flatExonsByTx}{`GRangesList` named by 'transcript_id'}
#'   \item{txdb}{`TxDb` object used to create various exon
#'     `GRanges` and `GRangesList` intermediate objects.
#'     Not essential for processing, but often useful.}
#'   \item{filesDF}{`data.frame` with colnames including
#'     the minimum required: 'sample_id', 'url', 'type'}
#' }
"sjenvtest"
