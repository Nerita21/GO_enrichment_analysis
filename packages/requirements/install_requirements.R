# Load required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install_if_missing <- function(pkgs, lib) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p)
    } else {
      message(p, " already installed, skipping.")
    }
  }
}
install_if_missing(c("systemfonts", "ggforce", "scatterpie", "dplyr", "ggplot2", "ggsci", "dotenv")
)

bioc_install_if_missing <- function(pkgs, lib) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      BiocManager::install(p, ask = FALSE, update = FALSE)
    } else {
      message(p, " already installed, skipping.")
    }
  }
}

bioc_install_if_missing(
  c("clusterProfiler", "AnnotationDbi", "Biostrings", "KEGGREST", "GO.db", "org.Hs.eg.db", "msigdb", "GOSemSim")
)