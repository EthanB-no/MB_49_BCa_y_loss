library(dplyr)
library(DESeq2)
library(dplyr)
library(purrr)
library(sva)

MB_counts <- read.csv("./11-10_MB_49_counts")
BC_counts <- read.csv("./pc_counts_5-20.csv")

BC_counts <- BC_counts %>%
  dplyr::rename(gene_name = X)


#inner join on gene names
combined_counts <- inner_join(MB_counts, BC_counts, by = c("gene_name"))

combined_counts$X <- NULL
combined_counts[is.na(combined_counts)] <- 0
# assuming your first column is gene_name
sample_names <- colnames(combined_counts)
sample_names <- sample_names[!sample_names %in% c("gene_name")]



colData <- data.frame(
  names = sample_names,
  batch = case_when(
    grepl("^BC", sample_names) ~ "BC",
    sample_names == "Norm_bla" ~ "BC",  # exception
    TRUE ~ "MB"
  ),
  stringsAsFactors = FALSE
)



# ensure colData is properly aligned
colData <- colData[match(colnames(combined_counts[,-1]), colData$names), ]
rownames(colData) <- colData$names


colData <- colData[!is.na(colData$names), ]

count_matrix <- combined_counts %>%
  select(-gene_name) %>%
  as.matrix()
rownames(count_matrix) <- combined_counts$gene_name

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = colData,
  design = ~ batch  # for now, simple model
)

vsd <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")

mat_compare <- assay(vsd)

mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, batch = colData$batch)
assay(vsd) <- mat
counts_batch_corrected <- assay(vsd)


write.csv(as.data.frame(counts_batch_corrected), "./11-12_VST_counts_batch_corrected.csv", row.names = TRUE)
