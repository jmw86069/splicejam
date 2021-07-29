# shiny todo optimization/performance tuning ideas:

* UNNECESSARY: consider how to keep sashimi plot updates to one pass,
preventing sample_id and gene changes from being out of
sync, and triggering a partial update to the plot before
both are ready.

   * `get_sashimi_data()` might be that step

* UNNECESSARY: move memoise function definitions outside the reactive({}) blocks

   * e.g. prepareSashimi_m
   * also check `use_memoise` and create non-memoise version

* UNNECESSARY: move gene_exon_structure `gene2gg()` into separate reactive function
that depends upon:

```
show_gene_model_d() # if FALSE then return NULL
flatExonsByGene1 <- get_flat_gene_exons_plot();
gene <- get_active_gene();
get_gene_coords()
exonLabelSize=14 + exon_label_size_d() + font_sizing_d()
```

* consider moving plotly sashimi creation into separate function

```
gg_ly <- get_plotly_sashimi() # reactive
```

* consider moving ggplot2 sashimi creation into separate function

```
gg_ly <- getgg_plot2_sashimi() # reactive
```

* COMPLETE: allow gene/transcript/exon model panel height to be adjusted

* TODO: when there are 2+ panel columns, display gene model below each column

* TODO: when "show detected transcripts" is not checked, re-calculate flat exons for that gene

* organize dependency tree and see if there are shortcuts in the process

* File issue with ggrepel, ask if labels outside plot range can be hidden,
for example when using `coord_cartesian()` to limit visual range.
