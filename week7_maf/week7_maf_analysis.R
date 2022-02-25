
# exercise 1.1
BiocManager::install("maftools")
library(maftools)

library(TCGAbiolinks)
library(ggplot2)

setwd("~/qbio490/analysis_data")

# exercise 1.2
clinic <- data.table::fread("~/qbio490/analysis_data/clinic.csv",
                            data.table = F)
colnames(clinic)[colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

# exercise 1.3
length(colnames(clinic))
# length oc colnames(clinic) = 77

length(colnames(clinic) == "bcr_patient_barcode")
# length = 77
# data type of booleans
colnames(clinic) == "bcr_patient_barcode"
# all false no trues

# exercise 1.4
mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)
maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)

# exercise 1.5
maf_dataframe <- data.table::fread("~/qbio490/analysis_data/GDCdata/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv",
                                  data.table = F)
  
clinic <- data.table::fread("~/qbio490/analysis_data/clinic.csv",
                            data.table = F)
colnames(clinic)[colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)
maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)

# exercise 2.1
maf_object
str(maf_object)

maf_object@data
str(maf_object@data)

maf_object@clinical.data
str(maf_object@clinical.data)
# They share the Tumor Sample Barcode column

# exercise 3.1
oncoplot(maf = maf_object,
         top = 10)
oncoplot(maf = maf_object,
         top = 20)
ggsave("~/qbio490/qbio_data_analysis_echo/week7_maf/oncoplot.png")

# exercise 3.2
# choose gene TP53. TP53 provides instructions on creating the tumor suppressor protein p53. p53 regulates cell division
# and prevents cells from proliferating too quickly. If there is a mutation in TP53, it would cause disruptions in generating
# p53, thus leading to no checks on cell growth.

# exercise 3.3
clinic = maf_object@clinical.data
young_patient_ids = c(clinic$Tumor_Sample_Barcode[clinic$age_category == "young"])
young_maf = subsetMaf(maf = maf_object,
                      tsb = young_patient_ids)

old_maf = subsetMaf(maf = maf_object,
                    tsb = c(clinic$Tumor_Sample_Barcode[clinic$age_category == "old"]))

# exercise 3.4
coOncoplot(m1 = young_maf, 
           m2 = old_maf, 
           m1Name = "MAF data of young patients", 
           m2Name = "MAF data of old patients")

ggsave("~/qbio490/qbio_data_analysis_echo/week7_maf/coOncoplot_old_young.png")

# exercise 3.5
# There are many genes where the mutation rate is not significantly different. The biggest difference in mutation rates
# among young and old patients are APC (63% in young and 75% in old) and TTN (45% in young and 54% in old). All mutation rates are 
# within 12% of each other.

# exercise 3.6
lollipopPlot(maf_object, gene = "TP53")
ggsave("~/qbio490/qbio_data_analysis_echo/week7_maf/tp53_lollipop.png")

# exercise 3.7
lollipopPlot2(m1 = young_maf, 
              m2 = old_maf, 
              m1_name = "MAF data of young patients",
              m2_name = "MAF data of old patients",
              gene = "TP53")
# the gene is more commonly mutated in old patients. Most mutations for both young and old patients are at p53 itself. This 
# is most likely because p53 is the protein that most directly controls cell proliferation. Missense mutations are the most common.
# One intersting thing to note is that at one position in p53, there are many mutated samples among old patients, much more than
# any other position and much more than the highest sample size in the young patients.

ggsave("~/qbio490/qbio_data_analysis_echo/week7_maf/tp53_lollipop_compare.png")

# exercise 3.8
# There are 30 samples that have only a mutation in gene B but no mutation in gene A

# exercise 3.9
# b = 7
# c = 2
# d = 35
# e = 37
# f = 42

# exercise 3.10
geneA_maf <- subsetMaf(maf = maf_object,
                       genes = "TP53")
geneB_maf <- subsetMaf(maf = maf_object, 
                       genes = "KRAS")

# exercise 3.11
help("subsetMaf")
# subsetMaf makes a maf_object that only has MAF data for the desired gene
str(geneA_maf)

str(geneA_maf@data)
# no, they don't all only have one mutation. There are missenses, positions, etc. that vary.
length(geneA_maf@data$Tumor_Sample_Barcode)
# There are not the same number of samples. This is only for those that have a mutation in TP53.

# exercise 3.12
mut_bc_geneA = c(geneA_maf@data$Tumor_Sample_Barcode)
mut_bc_geneB = c(geneB_maf@data$Tumor_Sample_Barcode)

num_mut_geneA = length(mut_bc_geneA)
# 222 patients
num_mut_geneB = length(mut_bc_geneB)
#165 patients

mut_bc_geneAB = intersect(mut_bc_geneA, mut_bc_geneB)
num_mut_geneAB = length(mut_bc_geneAB)
# 78 patients have a mutation in both genes

# exercise 3.13
num_mut_geneA_only = num_mut_geneA - num_mut_geneAB
# 144 with mutation in only TP53
num_mut_geneB_only = num_mut_geneB - num_mut_geneAB
# 87 with mutation in only KRAS

# exercise 3.14
num_neither_mutation = length(maf_object@clinical.data$Tumor_Sample_Barcode) - num_mut_geneA_only - num_mut_geneB_only - num_mut_geneAB
# 231 with neither mutation 

contig_table <- matrix(c(num_mut_geneAB, 
                         num_mut_geneB_only,
                         num_mut_geneA_only,
                         num_neither_mutation), 
                       nrow=2)
contig_table

# exercise 3.15
fe_results <- fisher.test(contig_table)
fe_results
# There is a very high rate of co-mutation for TP53 and KRAS; the p-value is 0.05794. If the hypothesis test takes in p = 0.05, however,
# this will not be considered statistically significant because 0.05794 > 0.05. This tells us that there exists a relation
# in mutation status for KRAS and TP53. 

