#!/usr/bin/env Rscript
# Run Structure in parallel through an R interface
library(parallel)
library(starmie)
library(readr)
set.seed(666351)
# prelim setup
# number of Ks to try
structure_binary <- "/usr/local/bioinf/bin/structure"
# run structure in parallel
n.cores <- 20L
tryK <- 2L:15L
n.runs <- 1L:20L

mainparams <- "./structure_output/mainparams"
extraparams <- "./structure_output/extraparams"
extraparamsLinkage <- "./structure_output/extraparamsLinkage"
unlinked_file <- "cache/png_unlinked.str"
unlinked_map <- "cache/png_unlinked.map"
complete_file <- "cache/png_complete.str"
complete_map <- "cache/png_complete.map"

stopifnot(file.exists(unlinked_file))
stopifnot(file.exists(complete_file))
stopifnot(file.exists(mainparams))
stopifnot(file.exists(extraparams))
stopifnot(file.exists(extraparamsLinkage))

N <- strsplit(system(paste("wc -l", unlinked_file), intern = TRUE), " ")[[1]][1]
nsnps <- strsplit(system(paste("wc -l", unlinked_map), intern = TRUE), " ")[[1]][1]


strout_gen <- function(prefix, x,y) {
  paste0(prefix, stringr::str_pad(x, width = 2, pad = 0), 
    "_K_", stringr::str_pad(y, width = 2, pad = 0), ".out")
}
# output files (will be appended by out_f from structure)
out_files <- outer(n.runs, tryK, 
  function(x,y) strout_gen("str_run_unlinked_", x, y))

log_files <- gsub("out", "log", out_files)


# create function to run structure, we will curry it to setup our two models
run_structure <- function(out_file, log_file, input_file, mp, ep, n, l) {
  k <- as.integer(stringr::str_extract(out_file, "[0-9]{2}\\b"))
  paste(structure_binary, "-K", k, "-L", l,  "-N", n,
   "-i", input_file, "-m", mainparams, "-e", extraparams, 
    "-D", round(runif(1) * 1e8), "-o", out_file, "&>", log_file)
}

# unlinked model run ----------------------------------------------------------
unlinked_run <- function(out_file, log_file) {
  cmd_out <- run_structure(out_file, log_file, unlinked_file, 
    mainparams, extraparams, N, nsnps)
  system(cmd_out)
  cmd_out
}

# unlinked model run ----------------------------------------------------------
print("Running unlinked model")
print(paste("Number of samples", N))
print(paste("Number of markers", nsnps))
st_runs <- mcmapply(unlinked_run, out_files, log_files, 
                    mc.cores = n.cores, mc.set.seed = TRUE)

print("Ran STRUCTURE with the following input commands")
print(st_runs)

str_filelist <- mapply(loadStructure, paste0(out_files, "_f"), 
                       log_files, SIMPLIFY = FALSE)
# save parsed output
write_rds(str_filelist, file = "structure_output/unlinked_model.rds")

if(file.exists("structure_output/unlinked_model.rds")) {
  system("rm str_run_unlinked_*.out_f")
  system("rm str_run_unlinked_*.log")
  system("mv seed.txt structure_output/str_run_unlinked_seeds.txt")
}

# linkage model run ----------------------------------------------------------
out_files <- outer(n.runs, tryK, 
  function(x,y) strout_gen("str_run_linkagemod_", x, y))

log_files <- gsub("out", "log", out_files)

N <- strsplit(system(paste("wc -l", complete_file), intern = TRUE), " ")[[1]][1]
nsnps <- strsplit(system(paste("wc -l", complete_map), intern = TRUE), " ")[[1]][1]

linked_run <- function(out_file, log_file) {
    cmd_out <- run_structure(out_file, log_file, complete_file, 
    mainparams, extraparamsLinkage, N, nsnps)
    system(cmd_out)
    cmd_out
}

print("Running linkage model")
print(paste("Number of samples", N))
print(paste("Number of markers", nsnps))

st_runs2 <- mcmapply(linked_run, out_files, log_files, 
                    mc.cores = n.cores, mc.set.seed = TRUE)

print(st_runs2)

str2_filelist <- mapply(loadStructure, paste0(out_files, "_f"), 
                       log_files, SIMPLIFY = FALSE)

write_rds(str_filelist, file = "structure_output/linkage_model.rds")

if(file.exists("structure_output/linkage_model.rds")) {
  system("rm str_run_linkagemod_*.out_f")
  system("rm str_run_linkagemod_*.log")
  system("mv seed.txt structure_output/str_run_linkagemod_seeds.txt")
}


