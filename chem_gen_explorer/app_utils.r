library(tidyverse)

annot_columns_of_interest <- c(
    "locus_tag",
    "gene_name",
    "product",
    "uniProt",
    # "eggNOG_OGs_Bacteria",
    # "eggNOG_OGs_Bacteroidia",
    "EggNOG_PFAMs",
    # "eggNOG_OGs",
    "EggNOG_KEGG_Pathway",
    "EggNOG_KEGG_Module",
    "EggNOG_KEGG_ko"
    # "EggNOG_GOs"
)

annot_columns_of_interest_colors <- map(
    annot_columns_of_interest,
    \(x) {
        # return(c("Not annotated" = "blue", "Annotated" = "green"))
        return(c("No annot." = "#CCCCCC", "Annot" = "darkgreen"))
    }
)
names(annot_columns_of_interest_colors) <- annot_columns_of_interest


load_fitness_data <- function(path) {
    data <- read_csv(path) %>%
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
    return(data)
}

load_annotations <- function(path, fitness_data = NULL) {
    annotations <- read_csv(path)
    annotations <- annotations[annotations$locus_tag %in% rownames(fitness_data), ] %>%
        select(
            all_of(annot_columns_of_interest) # comes from chemical_genetics_utils.r
        ) %>%
        distinct()
    annotations <- annotations[match(rownames(fitness_data), annotations$locus_tag), ]
    annotations <- annotations %>%
        column_to_rownames("locus_tag") %>%
        as.data.frame()
    return(annotations)
}

make_annot_col_bool <- function(co) {
    factor(case_when(
        is.na(co) ~ "No annot.",
        co == "-" ~ "No annot.",
        co == "no" ~ "No annot.",
        co == "No" ~ "No annot.",
        co == "" ~ "No annot.",
        co == " " ~ "No annot.",
        co == "blank" ~ "No annot.",
        co == "Blank" ~ "No annot.",
        TRUE ~ "Annot"
    ), levels = c("No annot.", "Annot"))
}

prep_char_selection <- function(highlight) {
    if (is.character(highlight)) {
        highlight <- str_replace_all(highlight, " ", "")
        highlight <- str_replace_all(highlight, "\n", "")
        # also remove white spaces..
        highlight <- str_split(highlight, ",")[[1]]
    }
    return(highlight)
}

ensure_data_concordance <- function(
    fitness_data = NULL,
    correlation_matrix = NULL,
    annotations = NULL) {
    stopifnot(all(dim(fitness_data)[1] == dim(correlation_matrix)[1]))
    stopifnot(all(rownames(fitness_data) == colnames(correlation_matrix)))
    stopifnot(all(rownames(fitness_data) == rownames(annotations)))
}

get_indices_from_heatmap_range <- function(num_rows, subset_perc_low, subset_perc_high) {
    i_low <- round(num_rows * subset_perc_low / 100)
    i_high <- round(num_rows * subset_perc_high / 100)
    return(list(i_low = i_low, i_high = i_high))
}


load_all <- function(
    fitness_data_path = NULL,
    annotations_path = NULL,
    subset_perc_low = NULL,
    subset_perc_high = NULL,
    cor_meth = "pearson") {
    fitness_data <- load_fitness_data(fitness_data_path)
    fi <- get_indices_from_heatmap_range(nrow(fitness_data), subset_perc_low, subset_perc_high)
    fitness_data <- fitness_data[fi$i_low:fi$i_high, ]
    annotations <- load_annotations(annotations_path, fitness_data)
    annotations_boolean <- annotations %>%
        mutate_all(make_annot_col_bool)
    correlation_matrix <- cor(t(fitness_data), method = cor_meth)
    hclusto <- hclust(dist(correlation_matrix))
    ensure_data_concordance()
    return(
        list(
            "fitness_data" = fitness_data,
            "annotations" = annotations,
            "annotations_boolean" = annotations_boolean,
            "correlation_matrix" = correlation_matrix,
            "hclust_object" = hclusto
        )
    )
}

assign_list_entries_to_global_env <- function(li) {
    for (entry_i in 1:length(li)[1]) {
        na <- names(li)[entry_i]
        entry <- li[[entry_i]]
        assign(na, entry, .GlobalEnv)
    }
}
