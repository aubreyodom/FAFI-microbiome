
# Script to analyze BLAST check results

# Setup
library(tidyverse)

# Obtain .csv paths
stem <- "/restricted/projectnb/singlecell/aodom/fafi_SRA/Sepsis/FAFI-microbiome/run_metascope/metascope_out/blast_fastas"
all_files <- list.files(stem, full.names = FALSE, include.dirs = FALSE)

ind1 <- stringr::str_detect(all_files, ".csv$")
ind2 <- stringr::str_detect(all_files, "Haemophilus_")

# Filter files matching the pattern *_Haemophilus*.csv
all_inputs <- all_files[ind1 & ind2]

# Function to get statistics of matches
get_stats <- function(input) {
  short_name <- stringr::str_remove(input, "_ntBlastResults.csv")
  strep <- file.path(stem, input) %>%
    readr::read_csv(col_names = FALSE, show_col_types = FALSE) %>%
    tidyr::separate(X1, into = c(NA, "seq_num"), "_",
                    remove = TRUE, convert = TRUE, extra = "drop")
  
  # All summed guesses, unfiltered
  sum_across <- strep %>% group_by(seq_num, X3) %>% summarise(total = n(), .groups = "drop") %>%
    arrange(seq_num, desc(total))
  
  # Get top hit
  get_maxes <- sum_across %>% group_by(seq_num) %>% filter(total == max(total))
  
  pct <- get_maxes %>% group_by(X3) %>%summarise(new_total = n(), .groups = "drop") %>%
    summarise(out = max(new_total) / sum(new_total) * 100, .groups = "drop") %>%
    pull(out) %>% round(2) %>% paste0('%')
  actual <- get_maxes %>% group_by(X3) %>% summarise(all_sum = sum(total)) %>% filter(all_sum == max(all_sum)) %>%
    pull(X3)
  message(short_name, ": likely ", actual, " ", pct)
  
  # Second most likely
  get_maxes_2 <- sum_across %>% group_by(seq_num) %>% dplyr::slice(2)
  likely_2 <- sort(table(get_maxes_2$X3), decreasing = TRUE)[1] |> names()
  
  ## total for each seq
  props <- sum_across |> group_by(seq_num) |> summarize(by_seq_denom = sum(total)) |>
    left_join(get_maxes_2, by = "seq_num") |>
    mutate(prop = total/by_seq_denom)
  pct_2 <- round(mean(props$prop)*100, 2) |> paste0('%')
  message(short_name, ": SECOND MOST likely ", likely_2, " ", pct_2)
  
  return(get_maxes)
}

# Vectorize functions
all_results <- plyr::alply(all_inputs, 1, get_stats)
names(all_results) <- all_inputs %>% stringr::str_remove("_ntBlastResults.csv")
