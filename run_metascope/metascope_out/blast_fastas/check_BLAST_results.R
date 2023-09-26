
# Script to analyze BLAST check results

# Setup
library(tidyverse)

# Obtain .csv paths
stem <- "~/infant-microbiome/work/aodom/Sepsis/run_metascope/metascope_out/blast_fastas"
all_inputs <- list.files(stem, pattern = "X3333_Weissella_ceti_ntBlastResults.csv", full.names = FALSE)

# Function to get statistics of matches
get_stats <- function(input) {
  short_name <- stringr::str_remove(input, "_ntBlastResults.csv")
  strep <- file.path(stem, input) %>%
  readr::read_csv(col_names = FALSE, show_col_types = FALSE) %>%
    tidyr::separate(X1, into = c(NA, "seq_num"), "_",
                    remove = TRUE, convert = TRUE, extra = "drop")
  
  sum_across <- strep %>% group_by(seq_num, X3) %>% summarise(total = n(), .groups = "drop") %>%
    arrange(seq_num, desc(total))
  
  get_maxes <- sum_across %>% group_by(seq_num) %>% filter(total == max(total))
  
  pct <- get_maxes %>% group_by(X3) %>%summarise(new_total = n(), .groups = "drop") %>%
    summarise(out = max(new_total) / sum(new_total) * 100, .groups = "drop") %>%
    pull(out) %>% round(2) %>% paste0('%')
  actual <- get_maxes %>% group_by(X3) %>% summarise(all_sum = sum(total)) %>% filter(all_sum == max(all_sum)) %>%
    pull(X3)
  message(short_name, ": likely ", actual, " ", pct)
  return(get_maxes)
}

# Vectorize functions
all_results <- plyr::alply(all_inputs, 1, get_stats)
names(all_results) <- all_inputs %>% stringr::str_remove("_ntBlastResults.csv")
