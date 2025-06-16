make_sub_heatmap <- function(
    res,
    highlight = NULL) {
    tryCatch(
        expr = {
            ht <- Heatmap(t(res[highlight, ]),
                top_annotation = HeatmapAnnotation(df = annotations_boolean[highlight, ], show_legend = FALSE),
                name = " ",
                show_row_names = TRUE, show_column_names = TRUE,
                show_row_dend = FALSE, show_column_dend = FALSE,
                column_names_gp = grid::gpar(fontsize = 14),
                row_names_gp = grid::gpar(
                    fontsize = 14,
                    heatmap_legend_param = list(
                        title_gp = grid::gpar(fontsize = 14), # Title font size
                        labels_gp = grid::gpar(fontsize = 14) # Labels font size
                    )
                )
            )
            ht <- draw(ht, merge_legend = TRUE)
            ht
        },
        error = function(e) {
            print("Error in make_sub_heatmap:")
            print(e)
        }
    )
}
