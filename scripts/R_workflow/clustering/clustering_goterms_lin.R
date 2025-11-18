##########################################################
# CLUSTERING 

# Load libraries
suppressMessages({
    library("clusterProfiler")
    library("AnnotationDbi")
    library("GO.db")
    library("org.Hs.eg.db")
    library("dplyr")
    library("ggplot2")
    library("enrichplot")
    library("GOSemSim")
    library("dotenv")
})

# --- get arguments ---
args <- commandArgs(trailingOnly=TRUE)
GO_result <- args[1]  # This should be an object of class enrichResult from clusterProfiler
output_file <- args[2]  # e.g., "GO_lin_clustering.tsv"

# --- 1. Calculate semantic similarity and clustering ---
calculate_sim_clustering <- function(GO_result) {
    GO_result_sim <- pairwise_termsim(GO_result, measure = "Lin")
    print("Calculated pairwise semantic similarity between GO terms...")
    GO_result_matrix <- GO_result_sim@termsim
    GO_dist_matrix <- 1- GO_result_matrix
    GO_hclustered <- hclust(as.dist(GO_dist_matrix), method = "average")


    # --- 2. Convert clustering to hierarchical table structure ---
    go_hcluster_table_maker <- function(GO_hclustered) {
    n_leaves <- length(GO_hclustered$labels)

    # Empty table structure
    table <- data.frame(
        id = integer(0),
        name = character(0),
        parent = integer(0),
        height = numeric(0),
        stringsAsFactors = FALSE
    )

    # Add leaves (original GO term labels)
    leaves <- data.frame(
        id = 1:n_leaves,
        name = GO_hclustered$labels,
        parent = NA_integer_,
        height = 0,
        stringsAsFactors = FALSE
    )
    table <- rbind(table, leaves)

    # Add internal cluster nodes
    for (i in seq_len(nrow(GO_hclustered$merge))) {
        left <- GO_hclustered$merge[i, 1]
        right <- GO_hclustered$merge[i, 2]
        
        left_id <- if (left < 0) -left else n_leaves + left
        right_id <- if (right < 0) -right else n_leaves + right
        parent_id <- n_leaves + i
        
        # Add parent node
        parent_row <- data.frame(
        id = parent_id,
        name = "",
        parent = NA_integer_,
        height = GO_hclustered$height[i],
        stringsAsFactors = FALSE
        )
        table <- rbind(table, parent_row)
        
        # Assign relationships
        table$parent[table$id == left_id] <- parent_id
        table$parent[table$id == right_id] <- parent_id
    }

    return(table)
    }

    go_hcluster_table <- go_hcluster_table_maker(GO_hclustered)
    cat("Performed hierarchical clustering of GO terms...\n")

    # --- 3. Attach p.adjust to leaf nodes ---
    # Assume GO_result@result has columns ID (GO ID) and Description
    go_hcluster_table$go_id <- ifelse(
        go_hcluster_table$name != "",
        GO_result@result$ID[match(go_hcluster_table$name, GO_result@result$Description)],
        NA_character_
    )

    # Map p.adjust using go_id
    leaf_padj <- setNames(GO_result@result$p.adjust, GO_result@result$ID)
    go_hcluster_table$padjust <- ifelse(
        !is.na(go_hcluster_table$go_id),
        leaf_padj[go_hcluster_table$go_id],
        NA_real_
    )

    # --- 4. Helper: Find the most specific common ancestor term ---
    find_common_ancestor <- function(go_terms) {
        if (length(go_terms) == 0) return(NA_character_)
        
        ancestors_list <- lapply(go_terms, function(go_id) {
            anc <- GOBPANCESTOR[[go_id]]
            if (!is.null(anc)) unique(c(anc, go_id)) else go_id
        })
        
        common_ancestors <- Reduce(intersect, ancestors_list)
        common_ancestors <- common_ancestors[!is.na(common_ancestors)]
        if (length(common_ancestors) == 0) return(NA_character_)
        
        depths <- sapply(common_ancestors, function(a) {
            if (!is.null(GOBPANCESTOR[[a]])) length(GOBPANCESTOR[[a]]) else 0
        })
        
        chosen <- common_ancestors[which.max(depths)]
        
        label <- tryCatch(Term(GOTERM[[chosen]]), error = function(e) NA_character_)
        return(label)
    }

    # --- 5. Helper: Get all leaf nodes under a given cluster ---
    get_leaves <- function(cluster_id, table) {
        idx <- which(!is.na(table$parent) & table$parent == cluster_id)
        children <- table$id[idx]
        print(children)
        leaves <- c()
        for (c in children) {
            nm <- table$go_id[table$id == c]
            print(nm)
            if (length(nm) == 1 && !is.na(nm) && nm != "") {
            leaves <- c(leaves, c)
            } else if (length(nm) == 0 || is.na(nm) || nm == "") {
            leaves <- c(leaves, get_leaves(c, table))
            } 
        }
        unique(leaves)
    }

    # --- 6. Assign labels and average p.adjust to internal nodes ---
    go_hcluster_table$label <- ifelse(go_hcluster_table$name != "", go_hcluster_table$name, NA_character_)
    internal_ids <- go_hcluster_table$id[go_hcluster_table$name == ""]

    for (cid in internal_ids) {
        leaf_ids <- get_leaves(cid, go_hcluster_table)
        leaf_terms <- go_hcluster_table$go_id[go_hcluster_table$id %in% leaf_ids]

        # Assign common ancestor label
        label <- find_common_ancestor(leaf_terms)
        go_hcluster_table$label[go_hcluster_table$id == cid] <- label

        # Compute average p.adjust
        go_hcluster_table$padjust[go_hcluster_table$id == cid] <- 
            mean(go_hcluster_table$padjust[go_hcluster_table$id %in% leaf_ids], na.rm = TRUE)
    }

    # --- 7. Add virtual root node ---
    root_id <- max(go_hcluster_table$id) + 1
    root_node <- data.frame(
        id = as.integer(root_id),
        name = "",
        parent = as.integer(NA),
        height = max(go_hcluster_table$height) + 1,
        label = "biological_process",
        go_id = "",
        padjust = mean(go_hcluster_table$padjust, na.rm = TRUE),
        stringsAsFactors = FALSE
    )

    go_hcluster_table <- rbind(go_hcluster_table, root_node)

    # Assign root as parent to top-level clusters
    top_level <- go_hcluster_table$id[is.na(go_hcluster_table$parent) & go_hcluster_table$id != root_id]
    go_hcluster_table$parent[go_hcluster_table$id %in% top_level] <- root_id

    cat("Cluster hierarchy with averaged p.adjust values created successfully!\n")
    return(go_hcluster_table)
}

# --- Run the clustering ---
clustering_results <- calculate_sim_clustering(GO_result)

# save clustering results
write.table(clustering_results, file = output_file, sep = "\t", row.names=FALSE, quote=FALSE)