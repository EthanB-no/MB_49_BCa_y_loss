library(DESeq2)
library(ggplot2)
library(dplyr)
library(purrr)

files <- list.files("./GSE229233_RAW/", pattern = "\\.txt", full.names = TRUE)

# read all tables into a named list
counts_list<- lapply(files, function(f) {
  read.delim(file = f, sep = "\t")
})

names(counts_list) <- sub("\\.txt$", "", files)

#inspect counts list 
length(counts_list)
head(counts_list[[1]])

lapply(counts_list, names)

which(!sapply(counts_list, function(df) "gene_name" %in% names(df)))
counts_list <- counts_list[-12]


# Extract gene_id and the readcount column
counts_df <- lapply(counts_list, function(df) {
  sample_col <- grep("_readcount$", names(df), value=TRUE)
  df %>% select(gene_name, all_of(sample_col))
})

# Collapse within each file: sum counts per gene_name
counts_df_summarized <- lapply(counts_df, function(df) {
  sample_col <- grep("_readcount$", names(df), value = TRUE)
  df %>%
    group_by(gene_name) %>%
    summarize(readcount = sum(.data[[sample_col]], na.rm = TRUE)) %>%
    rename_with(~ sub("_readcount$", "", sample_col), readcount)
})

# Merge them by gene_name
counts_table <- reduce(counts_df_summarized, full_join, by = "gene_name")

#check for uniqueness
is_unique <- !duplicated(counts_table$gene_name) & 
  !duplicated(counts_table$gene_name, fromLast = TRUE)

sum(is_unique)

write.csv(counts_table, file = "./11-10_MB_49_counts")

