metadata <- read_csv("Data/Cleaned_Metadata_infants.csv")

# This code identifies the filepaths of fastqs in the infant-microbiome project
find_sample <- function(samplenum) {
  l1 <- list.files(path = "/restricted/projectnb/infant-microbiome/data",
                   pattern = paste(samplenum, "_*", sep = ""),
                   full.names = TRUE, 
                   include.dirs = FALSE, recursive = TRUE)
  ind1 <- stringr::str_detect(string = l1, pattern = paste0("/", samplenum))

  all_files <- l1[ind1]
  ind_unique <- sapply(stringr::str_split(all_files, "/"),
                       function(x) tail(x, n = 1)) |>
    duplicated()
  return(all_files[!ind_unique])
}

tmp <- metadata
out <- lapply(tmp$`Sample ID`, find_sample)
names(out) <- tmp$`Sample ID`

out

# -----------------------------------------------------------------------------
# Copy all files to new folder
copy_func <- function(file) {
  file.copy(from = file,
            to = "/restricted/projectnb/infant-microbiome/data/death_sepsis",
            overwrite = FALSE, recursive = FALSE)
}
sapply(na.omit(unname(unlist(out))), copy_func)

# Where are the duplicates
all_lengths <- sapply(out, function(x) length(x))
dup_names <- names(all_lengths)[all_lengths >2]

# Remove any existing cat files
lcat <- list.files(
  path = "/restricted/projectnb/infant-microbiome/data/death_sepsis",
  pattern = "cat", full.names = TRUE, 
  include.dirs = FALSE, recursive = FALSE)
file.remove(lcat)

# Get names of all duplicate files
mk_files_mat <- function(num) {
  l2 <- list.files(
    path = "/restricted/projectnb/infant-microbiome/data/death_sepsis",
    pattern = paste0("*R", num, "_001.fastq.gz"), full.names = FALSE, 
    include.dirs = FALSE, recursive = TRUE)
  nums <- sapply(stringr::str_split(l2, "_"), function(x) head(x, n = 1))
  out_mat <- matrix(l2[nums %in% dup_names], ncol = 2, byrow = TRUE)
  return(out_mat)
}

# These are all files that need to be concatenated
write.table(rbind(mk_files_mat(1), mk_files_mat(2)),
            "/restricted/projectnb/infant-microbiome/data/death_sepsis/sample_paths.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")





