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
source("make_sub_heatmap.r")
source("make_pw_scatter.r")

env <- new.env() # This is crucial for this app, so don't remove!

make_heatmap <- function(correlation_matrix = NULL, dendrogram) {
    env$row_index <- 1:dim(correlation_matrix)[1]

    ht <- Heatmap(correlation_matrix,
        name = "effect size",
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE
    )
    ht <- draw(ht, merge_legend = TRUE)
    ht
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
    ),
    fluidRow(
        column(
            width = 4,
            id = "column2",
            box(
                title = "Fitness data of subselection", width = NULL, solidHeader = TRUE, status = "primary",
                plotOutput("sub_heatmap", width = "600px", height = "1500px")
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
        ),
        column(
            width = 8,
            box(
                title = "Annotations of highlighted genes", width = NULL, solidHeader = TRUE, status = "primary",
                DTOutput("res_table")
            ),
            fluidRow(
                column(
                    width = 12, # Adjust the width as needed (total should be 12)
                    box(
                        title = "Pairwise scatter of effect sizes", width = NULL, solidHeader = TRUE, status = "primary",
                        plotOutput("pairwise_scatters", height = "900px")
                    )
                )
            )
        )
    ),
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
        # sliderInput("range", min = 0, max = 100, value = c(0, 100), label = "percentages to subselect things to."),
        radioButtons("range",
            label = "Which part of the correlation matrix to show?",
            choices = list("1st third" = "0_33", "2nd third" = "33_66", "3rd third" = "66_100", "Everything" = "0_100"),
        ),
        # sliderInput("heatmaprangehigh", min = 50, max = 100, value = 100, label = "higher percentage of input genes to show"),
        textAreaInput("genes_to_viz", label = "Comma-separated list of genes"),
        actionButton("viz_specified_genes", label = "Subset heatmap!")
    ),
    body
)

server <- function(input, output, session) {
    # loads fitness_data and annotations; computes correlation_matrix, ensures data integrity.
    assign_list_entries_to_global_env(load_all(
        fitness_data_path = "Buni_compiled.csv",
        annotations_path = "essentiality_table_all_libraries_240818.csv",
        subset_perc_low = 0,
        subset_perc_high = 100
    ))

    ht <- make_heatmap(
        correlation_matrix
    )
    makeInteractiveComplexHeatmap(input, output, session, ht, "ht",
        brush_action = brush_action
    )

    observeEvent(input$range, {
        subset_perc_low <- as.numeric(str_split(input$range, "_")[[1]][1])
        subset_perc_high <- as.numeric(str_split(input$range, "_")[[1]][2])

        assign_list_entries_to_global_env(load_all(
            fitness_data_path = "Buni_compiled.csv",
            annotations_path = "essentiality_table_all_libraries_240818.csv",
            subset_perc_low = subset_perc_low,
            subset_perc_high = subset_perc_high
        ))

        ht <- make_heatmap(
            correlation_matrix
        )
        makeInteractiveComplexHeatmap(input, output, session, ht, "ht",
            brush_action = brush_action
        )
    })

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
