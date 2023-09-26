
# Create a MAE from PathoScope files
# Output is compatible with animalcules

# Setup -----
suppressPackageStartupMessages({
  library(tidyverse)
  library(animalcules)
})

# Step 0: USER-INPUT FILES -----

stem <- "/restricted/projectnb/infant-microbiome/work/aodom/Sepsis"
# Directory of PathoScope files
input_patho_files <- file.path(stem, "Patho_out")

# Read in metadata file (make sure to specify delimiter appropriately)
# Example:
metadata_table <- read.table(file.path(stem, "Data", "Cleaned_metadata_allsamp.tsv"),
                             sep = "\t",
                             header = TRUE,
                             row.names = 1,
                             stringsAsFactors = FALSE,
                             strip.white = TRUE)

# What directory to save MAE to?
# Example:
MAE_out <- file.path(stem, "Data", "raw_nofilter.RDS")

# Step 1: Create counts data from PathoScope files -----
count_table <- animalcules::read_pathoscope_data(input_dir = input_patho_files,
                                  pathoreport_file_suffix = "-sam-report.tsv")$countdat

# Step 1.5: Add "Xs" to columns
count_table %<>% rename_with(.fn = function(x) paste0("X", x))

# Step 2: Create MAE -----

# Choose only the samples in metadata that have counts data as well
sample_overlap <- intersect(colnames(count_table), rownames(metadata_table))
if (length(sample_overlap) < length(colnames(count_table))){
  print(paste("The following samples don't have metadata info:",
              paste(colnames(count_table)[which(!colnames(count_table) %in% sample_overlap)],
                    collapse = ",")))
  
  count_table <- count_table[,which(colnames(count_table) %in% sample_overlap)]
}
metadata_table <- metadata_table[match(colnames(count_table), rownames(metadata_table)),,drop = FALSE]

# print("read in done!")
# Test and fix the constant/zero row
row.remove.index <- c()
if (sum(base::rowSums(as.matrix(count_table)) == 0) > 0){
  row.remove.index <- which(base::rowSums(as.matrix(count_table)) == 0)
  count_table <- count_table[-row.remove.index,]
}

ids <- rownames(count_table)
tids <- unlist(lapply(ids, FUN = grep_tid))
if (sum(is.na(tids)) > 0){
  tid_remove <- which(is.na(tids))
  ids <- ids[-tid_remove]
  tids <- tids[-tid_remove]
  count_table <- count_table[-tid_remove,]
}

taxonLevels <- find_taxonomy(tids)
tax_table <- find_taxon_mat(ids, taxonLevels)


# replace spaces in tax name with underscore
tax_table <- as.data.frame(apply(tax_table,
                                 2,
                                 function(x)gsub('\\s+', '_',x)))

# create MAE object
se_mgx <-
  count_table %>%
  base::data.matrix() %>%
  S4Vectors::SimpleList() %>%
  magrittr::set_names("MGX")

se_colData <-
  metadata_table %>%
  S4Vectors::DataFrame()

se_rowData <-
  tax_table %>%
  base::data.frame() %>%
  dplyr::mutate_all(as.character) %>%
  dplyr::select(superkingdom, phylum, class, order, family, genus, species) %>%
  S4Vectors::DataFrame()

microbe_se <-
  SummarizedExperiment::SummarizedExperiment(assays = se_mgx,
                                             colData = se_colData,
                                             rowData = se_rowData)
mae_experiments <-
  S4Vectors::SimpleList(MicrobeGenetics = microbe_se)

MAE <-
  MultiAssayExperiment::MultiAssayExperiment(experiments = mae_experiments,
                                             colData = se_colData)

saveRDS(MAE, MAE_out)

