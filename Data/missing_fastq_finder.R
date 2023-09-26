all_files_patho <- c(list.files(path = "Patho_out",
                                pattern = "*.tsv", all.files = FALSE,
                                full.names = FALSE, recursive = TRUE,
                                ignore.case = FALSE, 
                                include.dirs = TRUE, no.. = FALSE))
patho_IDs <- sapply(strsplit(all_files_patho, split = "-sam"), function(x) x[[1]])
all_files_fastq <- read_table("Data/sample_fastq_paths.txt", col_names = FALSE) %>%
  unlist()
all_fastq_fixed <- all_files_fastq %>% strsplit(all_files_fastq, split = "/") %>%
  lapply(function(x) tail(x, n = 1)) %>% unlist()
all_IDs <- paste(sapply(strsplit(all_fastq_fixed, split = "_"), function(x) x[[1]]),
                 sapply(strsplit(all_fastq_fixed, split = "_"), function(x) x[[2]]), sep = "_")

# Which IDs are missing?
all_IDs[!(all_IDs %in% patho_IDs)]

# Indices of ID files that are not in the folder...
which(!(all_IDs %in% patho_IDs))

# Append .fastq
missing_fastq <- all_files_fastq[!(all_IDs %in% patho_IDs)]

write.table(missing_fastq, "Data/missing_fastq.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
