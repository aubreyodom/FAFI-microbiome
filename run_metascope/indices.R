library(MetaScope)

download_refseq("bacteria", reference = TRUE,
                representative = TRUE,
                compress = TRUE,
                patho_out = FALSE,
                out_dir = "/restricted/projectnb/pathoscope/reflib/2023_index_rep_bacteria",
                quiet = FALSE)