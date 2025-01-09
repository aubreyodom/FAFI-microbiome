# Which samples have generally most genus-level
counts_table["Alkalibacterium",] |> unlist() |> sort(decreasing = TRUE) |> extract(1:5)

X10192
X9308

# Most prevalent Alkalibacterium species
microbe_ind <- stringr::str_detect(rownames(counts_lcpm_species), "^Alkalibacterium")
prelim_microbe_counts <- counts_lcpm_species[microbe_ind, ]
test_microbe <- sort(rowSums(prelim_microbe_counts), decreasing = TRUE)

# Which samples have most of a given species
this_species <- "Alkalibacterium subtropicum"
counts_lcpm_species[this_species, ] |> unlist() |>
  sort(decreasing = TRUE) |> head()

# Investigate present Avibacterium species
counts_lcpm_species["Avibacterium",] |> unlist() |> sort(decreasing = TRUE) |> magrittr::extract(1:5)
microbe_ind <- stringr::str_detect(rownames(counts_lcpm_species), "^Avibacterium")
prelim_microbe_counts <- counts_lcpm_species[microbe_ind, ] 
unique(rownames(prelim_microbe_counts))
prelim_microbe_counts |> unlist() |> sort(decreasing = TRUE) |> magrittr::extract(1:5)

# Most prevalent Haemophilus species
counts_lcpm_species["Haemophilus",] |> unlist() |> sort(decreasing = TRUE) |> magrittr::extract(1:5)
microbe_ind <- stringr::str_detect(rownames(counts_lcpm_species), "^Haemophilus")
prelim_microbe_counts <- counts_lcpm_species[microbe_ind, ] 
rowSums(prelim_microbe_counts)
unique(rownames(prelim_microbe_counts))

result <- apply(prelim_microbe_counts, 1, function(row) {
  max_value <- max(row)                        # Find the maximum value
  max_column <- names(row)[which.max(row)]      # Find the column name
  c(Max_Value = max_value, Column_Name = max_column)
})

result