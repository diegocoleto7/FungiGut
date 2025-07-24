#!/usr/bin/env Rscript

library(tidyverse)
library(phyloseq)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_directory> <output_rds_file>")
}

input_dir <- args[1]
output_rds <- args[2]

sample_files <- list.files(path = input_dir,
                           pattern = "\\.txt$",
                           full.names = TRUE)

# Read and combine
all_samples_df <- map_dfr(sample_files, function(file_path) {
  sample_id <- tools::file_path_sans_ext(basename(file_path))

  # Read
  df <- read_delim(file_path,
                   delim = "\t",
                   comment = "@",
                   col_names = c("TaxID", "Level", "TaxPath", "TaxPathSN", "Percentage"),
                   show_col_types = FALSE)

  if (nrow(df) == 0 || all(is.na(df$Percentage))) {
    return(tibble(
      Sample = sample_id,
      Percentage = NA_real_,
      TaxID = NA_character_,
      Kingdom = NA_character_,
      Phylum = NA_character_,
      Class = NA_character_,
      Order = NA_character_,
      Family = NA_character_,
      Genus = NA_character_,
      Species = NA_character_
    ))
  }

  df %>%
    filter(!is.na(Percentage)) %>%
    mutate(Percentage = as.numeric(Percentage)) %>%
    separate(
      col = TaxPathSN,
      into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = "\\|",
      fill = "right",
      extra = "merge"
    ) %>%
    mutate(
      Sample = sample_id,
      TaxID = coalesce(Species, Genus, Family, Order, Class, Phylum, Kingdom)
    )
})
sample_ids <- unique(all_samples_df$Sample)

# Build OTU table
otu_partial <- all_samples_df %>%
  filter(!is.na(TaxID)) %>%
  select(Sample, TaxID, Percentage) %>%
  distinct(TaxID, Sample, .keep_all = TRUE) %>%
  pivot_wider(
    names_from = Sample,
    values_from = Percentage,
    values_fill = 0
  )

missing <- setdiff(sample_ids, colnames(otu_partial))
for (s in missing) {
  otu_partial[[s]] <- 0
}

otu_matrix <- otu_partial %>%
  select(TaxID, all_of(sort(sample_ids))) %>%
  column_to_rownames(var = "TaxID") %>%
  as.matrix()

otu_table_obj <- otu_table(otu_matrix, taxa_are_rows = TRUE)

# Build taxonomy table
tax_df <- all_samples_df %>%
  filter(!is.na(TaxID) & TaxID != "") %>%
  select(TaxID, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  distinct(TaxID, .keep_all = TRUE) %>%
  column_to_rownames(var = "TaxID") %>%
  as.matrix()

tax_table_obj <- tax_table(tax_df)

# Create phyloseq object
physeq <- phyloseq(otu_table_obj, tax_table_obj)
physeq_species <- tax_glom(physeq, taxrank = "Species")

# Save  RDS
saveRDS(physeq_species, file = output_rds)
