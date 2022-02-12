library(BiocManager)
BiocManager::install("SummarizedExperiment") 
library(TCGAbiolinks)
library(SummarizedExperiment)
setwd("~/qbio490/qbio_data_analysis_echo/analysis_data")

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts")

GDCdownload(query)
sum_exp <- GDCprepare(query)

str(sum_exp)

#exercise 2.2
counts = assays(sum_exp)$"HTSeq - Counts"[1:5, 1:5]
#[1:5, 1:5] means first five rows and first five columns
head(rowData(sum_exp))
colData(sum_exp)[1:5, 25:29]
metadata(sum_exp)

#exercise 2.3
dim(colData(sum_exp))
dim(rowData(sum_exp))
dim(assays(sum_exp)$"HTSeq - Counts")

#exercise 2.4
str(colData(sum_exp))
head(colData(sum_exp))
colnames(colData(sum_exp))

#exercise 2.5
colData(sum_exp)$age_at_diagnosis

#exercise 2.6
colData(sum_exp)$age_at_diagnosis[1:10]

#exercise 2.7
colData(sum_exp)$age_at_diagnosis = (colData(sum_exp)$age_at_diagnosis)/365

#get rid of NA values by substituting with age at index
#colData(sum_exp)$age_at_diagnosis[is.na(colData(sum_exp)$age_at_diagnosis)] = colData(sum_exp)$age_at_index[is.na(colData(sum_exp)$age_at_diagnosis)]
colData(sum_exp)[, colData(sum_exp)$age_at_diagnosis < 50] = colData(sum_exp)$age_category["Young"]
colData(sum_exp)[, colData(sum_exp)$age_at_diagnosis >= 50] = colData(sum_exp)$age_category["Old"]

#exercise 2.8
colData(sum_exp)$age_category = ifelse(colData(sum_exp)$age_at_diagnosis < 50, "Young", "Old")

#exercise 2.9
head(rowData(sum_exp))
dim(rowData(sum_exp))

#exercise 2.10
"APC" %in% rowData(sum_exp)[,2]
#True
"EGFR" %in% rowData(sum_exp)[,2]
#True

#exercise 2.11
assays(sum_exp)$"HTSeq - Counts"[20:25,30:35]
#columns: patients and rows: genes

#exercise 2.12
geneA_mask = (rowData(sum_exp)$external_gene_name == "APC")
geneB_mask = (rowData(sum_exp)$external_gene_name == "EGFR")
sum(geneA_mask)
sum(geneB_mask)

ensembl_geneA = rowData(sum_exp)$ensembl_gene_id[geneA_mask == TRUE]
ensembl_geneB = rowData(sum_exp)$ensembl_gene_id[geneB_mask == TRUE]

#exercise 2.13
#ensembl gene id is a row in assays(sum_exp)$"HTSeq - Counts"

#exercise 2.14
min(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ])
max(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ])
summary(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB, ])

#exercise 2.15
plot((assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ]),
          (assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB, ]),
          xlab = "APC", 
          ylab = "EGFR"
)
#They are positively correlated

#exercise 2.16
bool_age_na = is.na(colData(sum_exp)$age_category)
num_na = sum(bool_age_na)
num_na

#exercise 2.17
age_cat_no_NAs = colData(sum_exp)$age_category[(colData(sum_exp)$age_category != bool_age_na)]

#exercise 2.18
length(age_cat_no_NAs)

dim(colData(sum_exp)) #gives number of rows then number of columns
dim(colData(sum_exp))[1] #gives number of rows
dim(colData(sum_exp))[2] #gives number of columns

num_na + length(age_cat_no_NAs) == dim(colData(sum_exp))[1]

#exercise 2.19
dim(assays(sum_exp)$"HTSeq - Counts")
#521 patients

#exercise 2.20
identical(rownames(colData(sum_exp)), colnames(assays(sum_exp)$"HTSeq - Counts"))
gene_counts = (assays(sum_exp)$"HTSeq - Counts")[ensembl_geneA, !bool_age_na]

#exercise 2.21
length(gene_counts) == length(age_cat_no_NAs)
#gene counts was created using the counts as a dataframe. The rows, again, are the genes. Therefore, in the data
#frame, it will look for all indices that have ensembl_geneA as the gene. Then, it will scan all the columns that 
#have a valid age. This should match up with age_cat_no_NAs because this vector is all the indices with no NA,
#and gene counts should also have no NA. 

#exercise 2.22
boxplot(gene_counts ~ age_cat_no_NAs, xlab = "Age category", ylab = "Gene counts for APC")
#Some observations: There doesn't seem to be a significant difference in gene counts for APC by age. However, there seems to be
#a larger range of gene counts for old patients. 

