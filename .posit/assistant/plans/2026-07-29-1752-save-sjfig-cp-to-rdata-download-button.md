# Plan: Save sjfig_cp to RData — Download Button Feature

## Feasibility Assessment

**Yes, this is fully feasible without re-running `splicejamFigure()`.**

`get_sashimi_plot()` is a `shiny::reactive()` (server lines 555–682). Shiny reactives cache their return value and only re-execute when their reactive dependencies change. A `downloadHandler` that calls `get_sashimi_plot()` inside a Shiny session will receive the already-computed `sjfig$cp` — no repeat of the expensive `splicejamFigure()` call, as long as the user has not changed any plot inputs since the last render.

---

## What is `sjfig_cp`?

`get_sashimi_plot()` calls `splicejamFigure()` internally and returns `sjfig$cp` (server line 681), where `cp` is the **composed plot** — either:

- A `patchwork`/`cowplot` assembled ggplot2 figure (non-plotly mode), or
- A `plotly` subplot object (plotly mode)

The local variable `sjfig_cp <- get_sashimi_plot()` at line 870 is just an alias for this return value.

> **Optional scope note**: `get_sashimi_plot()` currently returns *only* `sjfig$cp`, not the full `sjfig` list (which also contains individual `gg_sashimi`, `gg_gene` components, `sashimi_data`, timing, etc.). If saving the full intermediate data is desired in the future, `get_sashimi_plot()` could be refactored to return the full list — but that is **out of scope for this feature**.

---

## Changes Required

### 1. `R/splicejam-shiny-server.R` — Add `downloadHandler`

Add a new `output$download_sjfig_cp` handler after the `get_sashimi_plot()` reactive definition (after line ~682). The handler:

1. Calls `sjfig_cp <- get_sashimi_plot()` to retrieve the cached reactive value.
2. Checks for `NULL` (no plot computed yet); the button should be disabled in the UI when there's no plot — but a NULL guard is still needed for safety.
3. Saves with `save(sjfig_cp, file = file)` so that `load("filename.RData")` restores the object under the name `sjfig_cp`.
4. Generates a descriptive filename: `sashimi_<gene>_<timestamp>.RData`.

```r
output$download_sjfig_cp <- shiny::downloadHandler(
   filename = function() {
      gene <- tryCatch(shiny::isolate(input$gene), error = function(e) "unknown")
      paste0("sashimi_", gene, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData")
   },
   content = function(file) {
      sjfig_cp <- get_sashimi_plot()
      if (length(sjfig_cp) == 0) {
         # Nothing to save; write an empty placeholder so the download still completes
         save(list = character(0), file = file)
         return(invisible(NULL))
      }
      save(sjfig_cp, file = file)
   }
)
```

### 2. `R/splicejam-shiny-ui.R` — Add `downloadButton`

Add a `shiny::downloadButton` near the existing "Update Sashimi Plots" `actionButton` (lines 358–361). Placing it immediately after keeps related actions together.

```r
shiny::actionButton(
   inputId = "calc_gene_params",
   label = "Update Sashimi Plots"),
shiny::downloadButton(
   outputId = "download_sjfig_cp",
   label = "Save Figure Data (.RData)")
```

---

## Behaviour Summary

| Scenario | Result |
|---|---|
| Plot has been rendered | `get_sashimi_plot()` returns cached `sjfig_cp`; file saved immediately |
| No plot yet (gene = "blank") | `sjfig_cp` is NULL; empty RData written; user sees a zero-byte or stub file |
| User changes inputs but hasn't clicked "Update" | Reactive is invalidated; `get_sashimi_plot()` re-runs `splicejamFigure()` on download click |

> The third scenario is worth documenting for the user: if plot inputs have been changed but "Update Sashimi Plots" has not been clicked, the download will trigger a fresh `splicejamFigure()` run. Clicking "Update" first is the recommended workflow.

---

## Files Modified

| File | Change |
|---|---|
| `R/splicejam-shiny-server.R` | Add `output$download_sjfig_cp` downloadHandler (~5 lines) |
| `R/splicejam-shiny-ui.R` | Add `downloadButton("download_sjfig_cp", ...)` after "Update Sashimi Plots" button |

---

## Out of Scope (for now)

- Saving the full `sjfig` object (all components + data): requires refactoring `get_sashimi_plot()` to return the full list.
- Exporting the figure as a PNG/PDF: a separate feature noted in your question as "two new features" — to be addressed separately.
- Disabling the download button when no plot is available: possible with `shinyjs::toggleState`, noted as a nice-to-have enhancement.
