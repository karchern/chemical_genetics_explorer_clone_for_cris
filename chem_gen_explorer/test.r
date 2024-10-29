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

# loads fitness_data and annotations; computes correlation_matrix, ensures data integrity
all_data <- load_all(
    fitness_data_path = "Buni_compiled.csv",
    annotations_path = "essentiality_table_all_libraries_240818.csv",
    subset_perc_low = 0,
    subset_perc_high = 100
)
fitness_data <- all_data$fitness_data
annotations <- all_data$annotations
annotations_boolean <- all_data$annotations_boolean
correlation_matrix <- all_data$correlation_matrix

res <- fitness_data
highlight <- c(1, 2, 3)
ht <- Heatmap(res[highlight, ],
    right_annotation = rowAnnotation(df = annotations_boolean[highlight, ], col = annot_columns_of_interest_colors),
    name = " ",
    show_row_names = TRUE, show_column_names = TRUE,
    show_row_dend = FALSE, show_column_dend = FALSE,
    column_names_gp = grid::gpar(fontsize = 14),
    row_names_gp = grid::gpar(fontsize = 14)
)
