library(tidyverse)
library(GGally)
library(here)
library(shiny)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(circlize)
library(GetoptLong)
library(shinydashboard)
library(DT)

source("app_utils.r")

env <- new.env() # This is crucial for this app, so don't remove!

fitness_data <- read_csv("Buni_compiled.csv") %>%
    select(Name, scaledLFC, contrast, FDR) %>%
    rename(
        gene = Name,
        log2fc_scaled = scaledLFC,
        condition = contrast,
        fdr = FDR
    ) %>%
    mutate(log2fc_scaled_empbayes = log2fc_scaled * (1 - fdr)) %>%
    pivot_wider(names_from = condition, id_cols = gene, values_from = log2fc_scaled_empbayes, values_fill = 0) %>%
    column_to_rownames("gene") %>%
    as.matrix()

correlation_matrix <- cor(t(fitness_data), method = "pearson")

annotations <- read_csv("essentiality_table_all_libraries_240818.csv")
annotations <- annotations[annotations$locus_tag %in% rownames(fitness_data), ] %>%
    select(
        all_of(annot_columns_of_interest) # comes from chemical_genetics_utils.r
    ) %>%
    distinct()
annotations <- annotations[match(rownames(fitness_data), annotations$locus_tag), ]
annotations_boolean <- annotations %>%
    mutate(across(-c(locus_tag), make_annot_col_bool))
annotations <- annotations %>%
    column_to_rownames("locus_tag") %>%
    as.data.frame()
annotations_boolean <- annotations_boolean %>%
    column_to_rownames("locus_tag") %>%
    as.data.frame()

stopifnot(all(dim(fitness_data)[1] == dim(correlation_matrix)[1]))
stopifnot(all(rownames(fitness_data) == colnames(correlation_matrix)))
stopifnot(all(rownames(fitness_data) == rownames(annotations)))
stopifnot(all(rownames(fitness_data) == rownames(annotations_boolean)))

make_heatmap <- function() {
    l <- rep(TRUE, dim(correlation_matrix)[1])
    env$row_index <- which(l)

    ht <- Heatmap(correlation_matrix,
        name = "effect size",
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE
    )
    ht <- draw(ht, merge_legend = TRUE)
    ht
}

make_sub_heatmap <- function(res, highlight = NULL) {
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
        }
    )
}

make_pw_scatter <- function(res, highlight = NULL) {
    tryCatch(
        expr = pairs(t(res[highlight, ])),
        error = function(e) {
        }
    )
}
body <- dashboardBody(
    fluidRow(
        column(
            width = 12,
            box(
                title = "Gene correlations over conditions", width = NULL, solidHeader = TRUE, status = "primary",
                originalHeatmapOutput("ht", height = 1500, width = 1500, containment = TRUE)
            )
        )
        # fluidRow(
        #     column(
        #         width = 8,
        #         height = 8,
        #         box(
        #             title = "Pairwise scatter plots", width = NULL, solidHeader = TRUE, status = "primary",
        #             plotOutput("pairwise_scatters")
        #         )
        #         # box(
        #         #     title = "Output", width = NULL, solidHeader = TRUE, status = "primary",
        #         #     HeatmapInfoOutput("ht", title = NULL)
        #         # ),
        #         # box(
        #         #     title = "Note", width = NULL, solidHeader = TRUE, status = "primary",
        #         #     htmlOutput("note")
        #         # ),
        #     ),
        #     # column(
        #     #     width = 4,
        #     #     box(
        #     #         title = "MA-plot", width = NULL, solidHeader = TRUE, status = "primary",
        #     #         plotOutput("ma_plot")
        #     #     ),
        #     #     box(
        #     #         title = "Volcanno plot", width = NULL, solidHeader = TRUE, status = "primary",
        #     #         plotOutput("volcanno_plot")
        #     #     ),
        #     #     box(
        #     #         title = "Result table of the selected genes", width = NULL, solidHeader = TRUE, status = "primary",
        #     #         DTOutput("res_table")
        #     #     )
        #     # ),
        #     tags$style("
        #         .content-wrapper, .right-side {
        #             overflow-x: auto;
        #         }
        #         .content {
        #             min-width:1500px;
        #         }
        #     ")
        # )
    ),
    fluidRow(
        column(
            width = 12,
            id = "column2",
            box(
                title = "Fitness data of subselection", width = NULL, solidHeader = TRUE, status = "primary",
                plotOutput("sub_heatmap", height = "400px")
            ),
            # column(
            #     width = 4,
            #     box(
            #         title = "MA-plot", width = NULL, solidHeader = TRUE, status = "primary",
            #         plotOutput("ma_plot")
            #     ),
            #     box(
            #         title = "Volcanno plot", width = NULL, solidHeader = TRUE, status = "primary",
            #         plotOutput("volcanno_plot")
            #     ),
            # ),
            tags$style("
            .content-wrapper, .right-side {
                overflow-x: auto;
            }
            .content {
                min-width:1500px;
            }
        ")
        )
    ),
    fluidRow(
        column(
            width = 6, # Adjust the width as needed (total should be 12)
            box(
                title = "Pairwise scatter of effect sizes", width = NULL, solidHeader = TRUE, status = "primary",
                plotOutput("pairwise_scatters", height = "600px")
            )
        )
    ),
    fluidRow(
        column(
            width = 12,
            box(
                title = "Annotations of highlighted genes", width = NULL, solidHeader = TRUE, status = "primary",
                DTOutput("res_table")
            )
        )
    )
)

brush_action <- function(df, input, output, session) {
    row_index <- unique(unlist(df$row_index))
    selected <- env$row_index[row_index]

    output[["sub_heatmap"]] <- renderPlot({
        make_sub_heatmap(fitness_data, selected)
    })

    output[["pairwise_scatters"]] <- renderPlot({
        make_pw_scatter(fitness_data, selected)
    })

    output[["res_table"]] <- renderDT(
        datatable(annotations[selected, ], rownames = TRUE)
    )

    #     output[["note"]] <- renderUI({
    #         if (!is.null(df)) {
    #             HTML(qq("<p>Row indices captured in <b>Output</b> only correspond to the matrix of the differential genes. To get the row indices in the original matrix, you need to perform:</p>
    # <pre>
    # l = res$padj <= @{input$fdr} &
    #     res$baseMean >= @{input$base_mean} &
    #     abs(res$log2FoldChange) >= @{input$log2fc}
    # l[is.na(l)] = FALSE
    # which(l)[row_index]
    # </pre>
    # <p>where <code>res</code> is the complete data frame from DESeq2 analysis and <code>row_index</code> is the <code>row_index</code> column captured from the code in <b>Output</b>.</p>"))
    #         }
    #     })
}

ui <- dashboardPage(
    dashboardHeader(title = "Chemical genetics data explorer"),
    dashboardSidebar(
        # selectInput("fdr", label = "Cutoff for FDRs:", c("0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05)),
        # numericInput("base_mean", label = "Minimal base mean:", value = 0),
        # numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 0),
        # actionButton("filter", label = "Generate heatmap")
        sliderInput("heatmaprangelow", min = 0, max = 50, value = 0, label = "lower percentage of input genes to show"),
        sliderInput("heatmaprangehigh", min = 50, max = 100, value = 100, label = "higher percentage of input genes to show"),
        textAreaInput("genes_to_viz", label = "Comma-separated list of genes"),
        actionButton("viz_specified_genes", label = "Subset heatmap!")
    ),
    body
)

server <- function(input, output, session) {
    ht <- make_heatmap()
    makeInteractiveComplexHeatmap(input, output, session, ht, "ht",
        brush_action = brush_action
    )
    observeEvent(input$viz_specified_genes,
        {
            selected <- prep_char_selection(input$genes_to_viz)
            output[["pairwise_scatters"]] <- renderPlot({
                make_pw_scatter(fitness_data, selected)
            })
            output[["sub_heatmap"]] <- renderPlot({
                make_sub_heatmap(fitness_data, selected)
            })
            output[["res_table"]] <- renderDT(
                datatable(annotations[selected, ], rownames = TRUE)
            )
        },
        ignoreNULL = FALSE
    )
}

shinyApp(ui, server)
