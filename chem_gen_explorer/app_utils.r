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
