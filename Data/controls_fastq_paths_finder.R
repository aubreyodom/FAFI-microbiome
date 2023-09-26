metadata <- read_csv("Data/Cleaned_Metadata_infants_controls.csv")

# This code identifies the filepaths of fastqs in the infant-microbiome project
find_sample <- function(samplenum) {
  l1 <- list.files(path = "/restricted/projectnb/infant-microbiome/data",
                   pattern = paste(samplenum, "_*", sep = ""),
                   full.names = TRUE, 
                   include.dirs = FALSE, recursive = TRUE)
  ind1 <- stringr::str_detect(string = l1, pattern = paste0("/", samplenum, "_"))
  
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
stem <- "/restricted/projectnb/infant-microbiome/data/death_sepsis_controls"

# Copy all files to new folder
copy_func <- function(file) {
  file.copy(from = file,
            to = stem,
            overwrite = FALSE, recursive = FALSE)
}

sapply(na.omit(unname(unlist(out))), copy_func)

