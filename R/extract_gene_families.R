library(tidyverse)

library(wbData)

gids_230 <- wb_load_gene_ids(230)
gids_294 <- wb_load_gene_ids(294)

# from http://www.treefam.org/download
# note: 1.7 GB

# download.file("http://www.treefam.org/static/download/treefam_family_data.tar.gz")





# Extract only the TFxxx.nh.emf as they're typically the smallest file
# and they contain all the information we need

# untar("data/treefam_family_data.tar.gz",
#       files = "treefam_family_data/*.nh.emf")


all_fams <- list.files("treefam_family_data",
                       pattern = "\\.nh\\.emf$") |>
  str_split_i("\\.", 1)

length(all_fams)
head(all_fams)

gene_by_fam <- map_dfr(
  all_fams,
  \(fam){
    tibble(
      family = fam,
      seqname =  paste0("treefam_family_data/", fam, ".nh.emf") |>
        read_lines() |>
        str_subset(pattern = "^SEQ caenorhabditis_elegans") |>
        str_split_i(" ", 8)
    )
  },
  .progress = TRUE)

# qs::qsave(gene_by_fam, "outs/241127_genes_by_family_unannotated.qs")


gene_by_fam <- qs::qread("241127_gene_by_family.qs") |>
  mutate(gene_id_230 = wb_seq2id(seqname, gids_230, warn_missing = TRUE),
         gene_id_294 = wb_clean_gene_names(gene_id_230),
         gene_name = i2s(gene_id_294, gids_294, warn_missing = TRUE))


# several of the 28 correspond to transposons that have been "suppressed"
gene_by_fam |> filter(is.na(gene_name))

stopifnot(identical( is.na(gene_by_fam$gene_name),
                     is.na(gene_by_fam$gene_id_294) ))


# gene_by_fam |>
#   filter(!is.na(gene_name)) |>
#   select(family, gene_id = gene_id_294, gene_name) |>
#   qs::qsave("outs/241127_genes_by_family.qs")







# tests and explorations

# # from http://www.treefam.org/static/download/species_list.txt
# taxon_id <- 6239
# 
# 
# fams <- read_tsv("data/treefam9_family_annotation.tab")
# 
# fams
# 
# fams |>
#   head() |>
#   as.data.frame()
# 
# # not a perfect match
# table( sort(unique(data_avail$family)) %in% 
#          sort(unique(fams$TF_name)) )
# 
# table( sort(unique(fams$TF_name)) %in% 
#          sort(unique(data_avail$family)) )
# 
# 
# 
# 
# fam_data <- untar("data/treefam_family_data.tar.gz",
#                   list = TRUE)
# 
# fam_data |> length()
# 
# data_avail <- enframe(fam_data,
#         value = "path",
#         name = NULL) |>
#   filter(startsWith(path, "treefam_family_data/TF")) |>
#   separate_wider_regex(
#     cols = path,
#     patterns = c("^treefam_family_data/",
#                  family = "TF[0-9]+", "\\.",
#                  filetype = "[a-z._]+"),
#     cols_remove = FALSE)
# 
# 
# data_avail
# length(unique(data_avail$family))
# length(unique(data_avail$filetype))
# 
# stopifnot(
#   nrow(data_avail) ==
#     length(unique(data_avail$family)) *
#     length(unique(data_avail$filetype)))
# )
# 
# 
# 
# unique(data_avail$filetype)
# 
# 
# # TF315535 contains appg
# my_fam <- "TF315535"
# 
# 
# xx <- data_avail |>
#   filter(family == my_fam) |>
#   pull(path)














