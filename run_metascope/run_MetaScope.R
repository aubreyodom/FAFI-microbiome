# Take in arguments from bash script
args <- commandArgs(trailingOnly = TRUE)

readPath1 <- args[1]
readPath2 <- args[2]
indexDir <- args[3]
expTag <- args[4] 
outDir <- args[5]
tmpDir <- args[6]
threads <- args[7]
targets <- stringr::str_split(args[8], ",")[[1]]
filters <- stringr::str_split(args[9], ",")[[1]]

# Time this!
now <- Sys.time()

# Load MetaScope
if (!requireNamespace("devtools", quietly = TRUE)) {
    devtools::install_github("compbiomed/metascope")
}
library(MetaScope)

# Align to targets
do_this <- function(x) stringr::str_replace_all(tolower(x), c(" " = "_"))
targets_ <- do_this(targets) 
filters_ <- do_this(filters)


target_map <- align_target_bowtie(read1 = readPath1,
                                  read2 = readPath2,
                                  lib_dir = indexDir,
                                  libs =  targets_,
                                  align_dir = tmpDir,
                                  align_file = expTag,
                                  overwrite = TRUE,
                                  threads = threads,
                                  quiet = FALSE,
                                  bowtie2_options = "--sensitive-local -k 100 --score-min L,20,1.0")

# Align to filters
#output <- paste(paste0(outDir, expTag), "filtered", sep = ".")
#final_map <- filter_host_bowtie(reads_bam = target_map,
#                                lib_dir = indexDir,
#                                libs = filters_,
#                                make_bam = FALSE,
#                                output = output,
#                                threads = threads,
#                                overwrite = TRUE,
#                                quiet = FALSE)

# MetaScope ID
output <- file.path(outDir, paste0(expTag, ".metascope_id.csv"))
metascope_id(target_map, input_type = "bam", aligner = "bowtie2",
             NCBI_key = "01d22876be34df5c28f4aedc479a2674c809",
             num_species_plot = 0,
             quiet = FALSE,
             out_dir = outDir)

file.copy(target_map, outDir)

message(capture.output(Sys.time() - now))
