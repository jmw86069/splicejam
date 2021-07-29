
observe

   * gene search, coord, junction minimum inputs
   * enable "Update" button

observe

   * input$gene
   * if "blank" then
   
      * disable slider bars
      * disable "Update" button

   * if not "blank" then

      * get_flat_gene_exons()
      * update coord slider bar
      * update exon name slider bar
      * update coord label

get_default_gene()

   * depends
      * get_gene_choices()


get_flat_gene_exons()

   * input$gene
   * returns exons for gene query

get_active_gene()

   * input$calc_gene_params
   * create list of gene query variables
   * disable "Update" button
   * return list

get_flat_gene_exons_plot()

   * depends
      * get_active_gene()
   * returns list(flatExonsByGene, gene_vals)

get_gene_coords()

   * depends
      * get_flat_gene_exons_plot()
   * return NULL if gene="blank"
   * return gene_coords from slider bar logic

get_sashimi_data()

   * depends
      * get_flat_gene_exons_plot()
      * get_sample_id_dt()
   * creates prepareSashimi_m() if it does not exist
   * isolate(input$gene)
   * sashimi_data <- prepareSashimi_m()
   * repeats some steps if some data is NULL
   * progress bar update

output$sashimiplot_output <- renderUI()

   * depends
      * get_sashimi_data()
      * share_y_axis_d()
      * do_plotly_d()
      * show_gene_model_d()
         * show_tx_model_d()
         * show_detected_tx_d()
         * exon_label_size_d()
         * font_sizing_d()
      * panel_height_d()
      * font_sizing_d()
