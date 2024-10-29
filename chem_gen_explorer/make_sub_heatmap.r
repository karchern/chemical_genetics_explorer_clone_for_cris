make_sub_heatmap <- function(
    res,
    highlight = NULL) {
    tryCatch(
        expr = {
            ht <- Heatmap(res[highlight, ],
                right_annotation = rowAnnotation(df = annotations_boolean[highlight, ], col = annot_columns_of_interest_colors),
                name = " ",
                show_row_names = TRUE, show_column_names = TRUE,
                show_row_dend = FALSE, show_column_dend = FALSE,
                column_names_gp = grid::gpar(fontsize = 14),
                row_names_gp = grid::gpar(fontsize = 14)
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
