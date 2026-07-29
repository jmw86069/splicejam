# Guidance for implementing TxDb and org.Hs.eg.db as input

## Goals

The current package is called splicejam. The goal is to expand splicejam
to support the use of Bioconductor ‘TxDb’ annotation packages as input,
instead of requiring a GTF file.

However, using TxDb input also requires using an organism gene
annotation package, for example human uses ‘org.Hs.eg.db’. This gene
annotation is necessary to convert the ‘gene_id’ used in ‘TxDb’ to the
‘gene_name’ used in splicejam. In this case, the SYMBOL from
‘org.Hs.eg.db’ will be stored as ‘gene_name’, for consistency with other
splicejam processes.

## Current Approach

The current workflow in
[`sashimiDataConstants()`](https://jmw86069.github.io/splicejam/reference/sashimiDataConstants.md)
works best with GTF file input.

- It then creates `data.frame` called ‘tx2geneDF’ which contains
  colnames ‘transcript_id’, ‘gene_name’, ‘gene_id’.
- It also uses the GTF file to create ‘txdb’ which is a `TxDb` object,
  using the txdbmaker package.
- The ‘txdb’ is used to create ‘exonsByTx’ using
  [`GenomicFeatures::exonsBy()`](https://rdrr.io/pkg/GenomicFeatures/man/transcriptsBy.html),
  and ‘cdsByTx’ also using
  [`GenomicFeatures::cdsBy()`](https://rdrr.io/pkg/GenomicFeatures/man/transcriptsBy.html).
  Those objects are `GRangesList` objects named by ‘transcript_id’.
- The exonsByTx and cdsByTx are used to create flatExonsByTx using
  [`flattenExonsBy()`](https://jmw86069.github.io/splicejam/reference/flattenExonsBy.md)
  in this splicejam package.

## Desired workflow

- Write a new function
  [`splicejamDataFromTxDb()`](https://jmw86069.github.io/splicejam/reference/splicejamDataFromTxDb.md)
  which requires two key arguments:

  1.  ‘txdb’ - specifically a Bioconductor package like
      ‘TxDb.Hsapiens.UCSC.mm10.knownGene’.
  2.  ‘ann_lib’ - a Bioconductor annotation package like ‘org.Mm.eg.db’.
      This package should be the same organism used in ‘txdb’, however
      we assume the user will provide the correct input and we will skip
      any validation for now. This package should be provided as a
      `character` string to be used with
      [`genejam::freshenGenes()`](https://jmw86069.github.io/genejam/reference/freshenGenes.html)
      later.

- It should have optional arguments

  - ‘detectedTx’ which would be a `character` vector of ‘transcript_id’
    entries to use in the
    [`flattenExonsBy()`](https://jmw86069.github.io/splicejam/reference/flattenExonsBy.md)
    step later.
  - ‘detectedGenes’ which would be a `character` vector of ‘gene_name’
    values to limit the number of entries retained in ‘flatExonsByGene’.
  - Internally, the ‘txdb’ data will contain all ‘transcript_id’ and all
    ‘gene_name’ entries. When ‘detectedGenes’ and ‘detectedTx’ are not
    provided, they will be derived from the full set of ‘transcript_id’
    and ‘gene_name’ entries available.
  - When ‘detectedTx’ is provided, it should also limit the available
    ‘gene_name’ entries via the tx2geneDF `data.frame` which contains
    the ‘transcript_id’ and ‘gene_name’ association. The tx2geneDF
    `data.frame` should also be subset to retain only those rows
    containing ‘detectedTx’ in the ‘transcript_id’ column.
  - When ‘detectedGenes’ is provided, also subset ‘tx2geneDF’ by
    matching ‘detectedGenes’ with the ‘gene_name’ column of ‘tx2geneDF’.
  - Make sure ‘detectedTx’ only contains entries in ‘tx2geneDF’.
  - Make sure ‘detectedGenes’ only contains entries in ‘tx2geneDF’.

- The TxDb should be used to prepare ‘exonsByTx’ and ‘cdsByTx’.

- I think ‘transcript_id’ can be associated with ‘gene_id’ using
  `GenomicFeatures::transcriptsBy(txdb, by='gene', use.names=TRUE)`. It
  will return a `GRangesList` named by ‘gene_id’, where each ‘gene_id’
  contains a `GRanges` object named by ‘transcript_id’. There might be a
  more direct approach.

- The TxDb should be used to prepare ‘exonsByGene’ specifically so the
  `names(exonsByGene)` will provide the ENTREZID entries which are
  available in the txdb data. I assume it may not include all possible
  ENTREZID values. The ‘ENTREZID’ should be stored as ‘gene_id’ in the
  ‘tx2geneDF’ `data.frame` later.

  - The ENTREZID entries, called ‘gene_id’ in splicejam, should be used
    to determine SYMBOL for each, and the SYMBOL should be stored as
    ‘gene_name’ in ‘tx2geneDF’. The ‘gene_name’ will become the primary
    name in the ‘flatExonsByGene’ data later. The genejam package,
    available by Github ‘jmw86069/genejam’, uses this approach to
    convert ENTREZID to SYMBOL. Notice the ann_lib is

  &nbsp;

      use_df <- data.frame(ENTREZID=gene_id_values);
      gene_df <- genejam::freshenGenes(use_df, ann_lib='org.Mm.eg.db', empty_rule=='original');

  It returns a `data.frame` with new column ‘SYMBOL’. When any ENTREZID
  entry is not found, it is instructed to return the original entry,
  which is the ENTREZID.

  - Now the data should be available to create `data.frame` ‘tx2geneDF’
    with colnames ‘transcript_id’, ‘gene_id’, ‘gene_name’.
  - This step should probably be performed in a stand-alone `function`
    that takes ‘txdb’, ‘ann_lib’ as input, and returns `data.frame`
    ‘tx2geneDF’ as output. Name it
    [`makeTx2geneFromTxdb()`](https://jmw86069.github.io/splicejam/reference/makeTx2geneFromTxdb.md)

- After ‘tx2geneDF’ is available, it should be subset by using
  ‘detectedTx’ in the ‘transcript_id’ column, and by using
  ‘detectedGenes’ in the ‘gene_name’ column. Afterward, ‘detectedTx’
  should only contain entries in ‘tx2geneDF’ column ‘transcript_id’, and
  ‘detectedGenes’ should only contain entries in the ‘tx2geneDF’ column
  ‘gene_name’.

- Next, ‘exonsByTx’ and ‘cdsByTx’ should be subset using ‘detectedTx’ to
  match `names(exonsByTx)` and `names(cdsByTx)`.

- Then ‘flatExonsByTx’ should be derived.

- Then ‘flatExonsByGene’ should be derived.

- Then these data objects should be stored into a new `environment`.

  - ‘flatExonsByGene’
  - ‘flatExonsByTx’
  - ‘tx2geneDF’
  - ‘detectedTx’
  - ‘detectedGenes’

- Other optional arguments:

  - ‘filesDF’ is an optional `data.frame` with colnames ‘sample_id’,
    ‘url’, ‘type’, and may contain other colnames.
  - ‘color_sub’ is an optional `character` vector of R colors whose
    names should all match ’filesDF\$sample_id' if 'filesDF' is
    provided. If 'filesDF' is provided, but 'color_sub' is not provided,
    'color_sub' should be derived using Github package
    'jmw86069/colorjam' like this: \`color_sub \<-
    colorjam::group2colors(unique(filesDF\$sample_id));\`

- Note that the Github packages ‘colorjam’ and ‘genejam’ are already
  installed. Also installed: ‘TxDb.Mmusculus.UCSC.knownGene’ and
  ‘org.Mm.eg.db’.

## Testing

- Using ‘mm10’ mouse genome as described above, you can re-use the same
  ‘filesDF’ in the test data ‘sjenvtest’ since it also uses the ‘mm10’
  mouse genome. You can use ‘Gria1’ as the test gene.
