
#exercise 1.1
library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)
BiocManager::install("DESeq2")
library(DESeq2)

setwd("~/qbio490/qbio_data_analysis_echo/analysis_data")

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")

GDCdownload(query)
sum_exp <- GDCprepare(query)

#exercise 1.2
bool_age_na = is.na(colData(sum_exp)$age_at_diagnosis)
bool_age_na

patients_data = colData(sum_exp)  # contains the clinical data
counts = assays(sum_exp)$"HTSeq - Counts" # contains the counts data
#4 patients have NA data

# counts = data.frame(counts)

patients_data = patients_data[!bool_age_na,]
counts = counts[, !bool_age_na]

patients_data$age_category = ifelse(patients_data$age_at_diagnosis < 50, "Young", "Old")

patients_data$age_category = factor(patients_data$age_category, levels = c("Young", "Old"))

#exercise 1.3
if (all(rownames(counts) == names(rowRanges(sum_exp)))){
  print("Changing row names")
  rownames(counts) = rowRanges(sum_exp)$external_gene_name
}

counts_row_sums = rowSums(counts)
low_counts_mask = counts_row_sums >= 10
sum(low_counts_mask)
counts = counts[low_counts_mask, ]

#exercise 2.1
dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = patients_data,
                             design = ~age_category)

dds_obj = DESeq(dds)
resultsNames(dds_obj)  # see what comparisons got run

# get the young vs. old comparison
results = results(dds_obj, format = "DataFrame", contrast = c("age_category", "young", "old"))

#exercise 2.2
my_df = data.frame(x = c('b', 'd', 'c', 'e', 'a'),
                   y = c(2,4,3,5,1))
order_indices = order(my_df$y)
# expect 5, 1, 3, 2 4
order_indices

my_df = my_df[order_indices, ]
my_df

#exercise 2.3
row_order = order(results$padj)
results = results[row_order,]
head(results)

# a) SCARNA5. It is highly expressed in older patients.
# b) The full name is Small Cajal Body-Specific RNA 5. It is a kind of snoRNA, primarily functioning in biogenesis of 
# snRNPs and modifying RNA polymerase.

#exercise 2.4
log2FoldChange_threshold = 1
padj_threshold = 0.05

log2FoldChange_greater = results$log2FoldChange > log2FoldChange_threshold
padj_lesser = results$padj < padj_threshold
results = log2FoldChange_greater & padj_lesser

# Exercise 2.5
fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

# fill in your plot code here!
# be sure to relabel the axes!
# note: you can perform the log transformation directly in the plot function
plot(x = results$log2FoldChange,
     y = -log10(results$padj),
     xlab = "Log 2 Fold Change (young/old)", # be sure the specify that it's young over old!
     ylab = "-log10 (adjusted) p-value",
     pch = 20) # smaller solid circles

# these lines put the lines on the plot
# abline() plots straight lines on an R plot.
# v argument is for a vertical line, h argument is for a horizontal line, col argument is color
abline(v=c(-log2(fc_threshold), log2(fc_threshold)), h= c(-log10(p_threshold)), col="green")

# Exercise 2.6
write.csv(x = results,
          file = "~/qbio490/qbio_data_analysis_echo/week6_DESeq2/results.csv",
          row.names = FALSE)








