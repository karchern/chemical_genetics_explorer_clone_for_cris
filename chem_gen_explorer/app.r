library(tidyverse)
library(GGally)
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

options(shiny.maxRequestSize = 250 * 1024^2)

make_heatmap <- function(correlation_matrix = NULL, dendrogram) {
    env$row_index <- 1:dim(correlation_matrix)[1]
    env$column_index <- 1:dim(correlation_matrix)[2]


    ht <- Heatmap(correlation_matrix,
        name = "effect size",
        cluster_rows = dendrogram,
        cluster_columns = dendrogram,
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
                title = "Current selection of genes", width = NULL, solidHeader = TRUE, status = "info",
                textOutput("current_selection_textbox")
            ),
            box(
                title = "Annotations of highlighted genes", width = NULL, solidHeader = TRUE, status = "info",
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
    coluumn_index <- unique(unlist(df$column_index))
    inters <- union(row_index, coluumn_index)
    selected <- env$row_index[inters]

    output[["sub_heatmap"]] <- renderPlot({
        make_sub_heatmap(fitness_data, selected)
    })

    output[["pairwise_scatters"]] <- renderPlot({
        make_pw_scatter(fitness_data, selected)
    })

    output[["res_table"]] <- renderDT(
        datatable(annotations[selected, ], rownames = TRUE)
    )

    output[["current_selection_textbox"]] <- renderText({
        str_c(rownames(fitness_data)[selected], collapse = ",")
    })
}

effect_size_names <- c("scaledLFC", "EB")

ui <- dashboardPage(
    dashboardHeader(title = "Chemical genetics data explorer"),
    dashboardSidebar(
        actionButton("load_default_input_data", label = HTML("Load default input data<br>(for debugging only!)")),
        radioButtons("range",
            label = "Which part of the correlation matrix to show?",
            choices = list(
                "Everything" = "0_100",
                "1st 20%" = "0_20",
                "2nd 20%" = "19_40",
                "3rd 20%" = "39_60",
                "4th 20%" = "59_80",
                "5th 20%" = "79_100"
            ),
        ),
        fileInput("fitness_data", "Upload fitness table"),
        fileInput("annotation_data", "Upload annotation table"),
        fileInput("precomputed_dendrogram", "Precomputed dendrogram (optional)", accept = c(".rds", ".RData", ".RDS", ".rdata")),
        checkboxGroupInput(
            "annot_columns_of_interest",
            label = "Select annotation columns of interest:",
            choices = c(
                "locus_tag" = "locus_tag",
                "gene_name" = "gene_name",
                "product" = "product",
                "uniProt" = "uniProt",
                "eggNOG_OGs_Bacteria" = "eggNOG_OGs_Bacteria",
                "eggNOG_OGs_Bacteroidia" = "eggNOG_OGs_Bacteroidia",
                "EggNOG_PFAMs" = "EggNOG_PFAMs",
                "eggNOG_OGs" = "eggNOG_OGs",
                "EggNOG_KEGG_Pathway" = "EggNOG_KEGG_Pathway",
                "EggNOG_KEGG_Module" = "EggNOG_KEGG_Module",
                "EggNOG_KEGG_ko" = "EggNOG_KEGG_ko",
                "EggNOG_GOs" = "EggNOG_GOs"
            ),
            selected = c(
                "locus_tag",
                "gene_name",
                "product",
                "uniProt",
                "EggNOG_PFAMs",
                "EggNOG_KEGG_Pathway",
                "EggNOG_KEGG_Module",
                "EggNOG_KEGG_ko"
            ),
            inline = TRUE
        ),
        radioButtons("effect_size_name", "Column name holding effect\nsize to use for correlation (if EB, make sure you have a\ncolumn holding q-values called FDR)", effect_size_names),
        actionButton("load_input_data", label = "Load input data"),
        textAreaInput("genes_to_viz", label = "Comma-separated list of genes"),
        actionButton("trigger_genes_to_viz", label = "Subset heatmap!"),
        textInput("gene_you_want_to_zoom_in_on", label = "gene of interest (correlated genes will be highlighted)"),
        numericInput("number_of_genes_you_want_to_zoom_in_around", "number of genes to correlate", value = 5, min = 0, max = 1000),
        radioButtons(
            "type_of_correlation",
            label = "Select correlation type:",
            choices = list(
                "Most correlated genes" = "correlated",
                "Most anti-correlated genes" = "anti_correlated"
            ),
            selected = "correlated"
        ),
        actionButton("trigger_gene_you_want_to_zoom_in_on", label = "Zoom in on gene!")
    ),
    body
)

server <- function(input, output, session) {
    # loads fitness_data and annotations; computes correlation_matrix, ensures data integrity.
    observeEvent(input$load_input_data, {
        subset_perc_low <- as.numeric(str_split(input$range, "_")[[1]][1])
        subset_perc_high <- as.numeric(str_split(input$range, "_")[[1]][2])
        tryCatch(
            assign_list_entries_to_global_env(load_all(
                # fitness_data_path = "Buni_compiled.csv",
                fitness_data_path = input$fitness_data$datapath,
                # annotations_path = "essentiality_table_all_libraries_240818.csv",
                annotations_path = input$annotation_data$datapath,
                dendrogram_path = input$precomputed_dendrogram$datapath,
                subset_perc_low = subset_perc_low,
                subset_perc_high = subset_perc_high,
                load_what = input$effect_size_name,
                annot_columns_of_interest = input$annot_columns_of_interest
            )),
            error = function(e) {
                print("Loading of data failed, make sure to set both fitness and annotation data...")
            }
        )

        ht <- make_heatmap(
            correlation_matrix,
            hclust_object
        )
        makeInteractiveComplexHeatmap(input, output, session, ht, "ht",
            brush_action = brush_action
        )
    })

    observeEvent(input$load_default_input_data, {
        subset_perc_low <- as.numeric(str_split(input$range, "_")[[1]][1])
        subset_perc_high <- as.numeric(str_split(input$range, "_")[[1]][2])
        tryCatch(
            assign_list_entries_to_global_env(load_all(
                fitness_data_path = "Buni_compiled.csv",
                # fitness_data_path = input$fitness_data$datapath,
                annotations_path = "essentiality_table_all_libraries_240818.csv",
                # annotations_path = input$annotation_data$datapath,
                dendrogram_path = input$precomputed_dendrogram$datapath,
                subset_perc_low = subset_perc_low,
                subset_perc_high = subset_perc_high,
                load_what = input$effect_size_name,
                annot_columns_of_interest = input$annot_columns_of_interest
            )),
            error = function(e) {
                print("Loading of data failed, make sure to set both fitness and annotation data...")
            }
        )

        ht <- make_heatmap(
            correlation_matrix,
            hclust_object
        )
        makeInteractiveComplexHeatmap(input, output, session, ht, "ht",
            brush_action = brush_action
        )
    })

    observeEvent(input$range, {
        if (is.null(input$fitness_data)) {
            return()
        }

        subset_perc_low <- as.numeric(str_split(input$range, "_")[[1]][1])
        subset_perc_high <- as.numeric(str_split(input$range, "_")[[1]][2])

        tryCatch(
            assign_list_entries_to_global_env(load_all(
                # fitness_data_path = "Buni_compiled.csv",
                fitness_data_path = input$fitness_data$datapath,
                # annotations_path = "essentiality_table_all_libraries_240818.csv",
                annotations_path = input$annotation_data$datapath,
                dendrogram_path = input$precomputed_dendrogram$datapath,
                subset_perc_low = subset_perc_low,
                subset_perc_high = subset_perc_high,
                load_what = input$effect_size_name,
                annot_columns_of_interest = input$annot_columns_of_interest
            )),
            error = function(e) {
                print("Loading of data failed, make sure to set both fitness and annotation data...")
            }
        )

        ht <- make_heatmap(
            correlation_matrix,
            hclust_object
        )
        makeInteractiveComplexHeatmap(input, output, session, ht, "ht",
            brush_action = brush_action
        )
    })

    observeEvent(input$trigger_genes_to_viz,
        {
            if (is.null(input$fitness_data)) {
                return()
            }
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
            output[["current_selection_textbox"]] <- renderText({
                str_c(rownames(fitness_data)[selected], collapse = ",")
            })
        },
        ignoreNULL = FALSE
    )

    observeEvent(input$trigger_gene_you_want_to_zoom_in_on,
        {
            gene_index_of_interest <- which(rownames(correlation_matrix) == input$gene_you_want_to_zoom_in_on)
            if (input$type_of_correlation == "correlated") {
                selected <- order(correlation_matrix[gene_index_of_interest, ], decreasing = TRUE)[1:input$number_of_genes_you_want_to_zoom_in_around]
            } else if (input$type_of_correlation == "anti_correlated") {
                selected <- c(gene_index_of_interest, order(correlation_matrix[gene_index_of_interest, ], decreasing = FALSE)[1:(input$number_of_genes_you_want_to_zoom_in_around - 1)])
            } else {
                stop("Unknown type of correlation selected.")
            }
            output[["pairwise_scatters"]] <- renderPlot({
                make_pw_scatter(fitness_data, selected)
            })
            output[["sub_heatmap"]] <- renderPlot({
                make_sub_heatmap(fitness_data, selected)
            })
            output[["res_table"]] <- renderDT(
                datatable(annotations[selected, ], rownames = TRUE)
            )
            output[["current_selection_textbox"]] <- renderText({
                str_c(rownames(fitness_data)[selected], collapse = ",")
            })
        },
        ignoreNULL = FALSE,
        ignoreInit = TRUE
    )
}

shinyApp(ui, server)
