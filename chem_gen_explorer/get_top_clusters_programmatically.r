library(tidyverse)
library(ComplexHeatmap)

# Load
cor_matrix <- read_csv("/Users/karcher/Correlation_matrixes/Buni_Correlationmatrix.csv") %>%
    as.data.frame() %>%
    column_to_rownames("...1") %>%
    as.matrix()

check_clusters <- function(cor_matrix = NULL,
                           distances_of_cor_matrix = NULL,
                           n_top_joins = NULL,
                           num_clusters_to_break_at = NULL,
                           min_mean_cluster_similarity = NULL,
                           nn_tmp = NULL) {
    hc_a <- hclust(distances_of_cor_matrix)
    if (n_top_joins > length(hc_a$height)) {
        quit("n_top_joins is too large, set to a smaller value")
    }
    cut_res <- cutree(hc_a, h = hc_a$height[n_top_joins])
    cluster_counts <- table(cut_res)
    mean_cluster_similarities <- numeric(length(cluster_counts))
    cluster_members <- list()
    for (cluster_index in 1:length(cluster_counts)) {
        cluster_id <- names(cluster_counts)[cluster_index]
        cluster_size <- cluster_counts[cluster_index]
        mean_cluster_similarities[cluster_index] <- mean(cor_matrix[cut_res == cluster_id, cut_res == cluster_id])
        cluster_members[[cluster_id]] <- rownames(cor_matrix)[cut_res == cluster_id]
    }
    if (nn_tmp %% 10 == 0) {
        print(str_c("Heuristic still running, having ", sum(cluster_counts >= 2), " clusters with mean similarity of ", round(mean(mean_cluster_similarities[cluster_counts >= 2]), 3)))
    }

    if (sum(cluster_counts >= 2) == num_clusters_to_break_at || (mean(mean_cluster_similarities[cluster_counts >= 2]) < 0.75)) {
        return(
            tibble(
                cluster_id = names(cluster_counts),
                cluster_size = cluster_counts,
                mean_cluster_similarity = mean_cluster_similarities,
                cluster_members = cluster_members,
                cluster_members_comma_separated = map_chr(cluster_members, ~ str_c(., collapse = ","))
            ) %>%
                filter(cluster_size > 1)
        )
    } else {
        n_top_joins <- n_top_joins + 1
    }
}

get_top_clusters <- function(cor_matrix = NULL,
                             distances_of_cor_matrix = NULL,
                             n_top_joins = 100, # This is the number of top joins to start with,
                             num_clusters_to_break_at = 150, # If 150 clusters with at least 2 members are reached, break, or...
                             min_mean_cluster_similarity = 0.75 # If the mean similarity of clusters with at least 2 members is below this, break as well
) {
    print("Starting heuristic to get informative clusters...")
    nn_tmp <- 0
    while (TRUE) {
        # start with n_top_joins <- 50
        res <- check_clusters(
            cor_matrix,
            distances_of_cor_matrix,
            n_top_joins,
            num_clusters_to_break_at = num_clusters_to_break_at,
            min_mean_cluster_similarity = min_mean_cluster_similarity,
            nn_tmp = nn_tmp
        )
        if (!is.numeric(res)) {
            break
        } else {
            n_top_joins <- res
            nn_tmp <- nn_tmp + 1
        }
    }
    return(
        res %>%
            # arrange by mean_cluster_similarity, but also favour larger clusters
            # arrange(desc(mean_cluster_similarity)) %>%
            mutate(order_metric = mean_cluster_similarity * (cluster_size^(1 / 5))) %>%
            arrange(desc(order_metric)) %>%
            select(-order_metric) %>%
            filter(cluster_size > 1)
    )
}


# @ Eva: Exchange the following line with the path to your correlation matrix
# This will run for a few minutes
top_clusters <- get_top_clusters(
    cor_matrix,
    dist(cor_matrix)
    n_top_joins = 100 # Sensible default
)
