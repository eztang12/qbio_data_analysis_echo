setwd("~/qbio490/qbio_data_analysis_echo/analysis_data")

library(maftools)
library(survival)
library(survminer)
library(SummarizedExperiment)
library(TCGAbiolinks)

# get gene counts for EFGR stratified by gender
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")

sum_exp <- GDCprepare(query)

patients_data = colData(sum_exp)
counts = assays(sum_exp)$"HTSeq - Counts"

EGFR_mask = (rowData(sum_exp)$external_gene_name == "EGFR")
ensembl_EFGR = rowData(sum_exp)$ensembl_gene_id[EGFR_mask == TRUE]

gene_counts = (assays(sum_exp)$"HTSeq - Counts")[ensembl_EFGR,]
summary(gene_counts)

# figure 1 - not much difference
boxplot(gene_counts ~ patients_data$gender, xlab = "Gender", ylab = "Gene counts for EFGR")


# get survival data for patients with relative low EGFR counts
miss_death <- which(is.na(patients_data$days_to_death))
patients_data$days_to_death[miss_death] <- patients_datadays_to_last_follow_up[miss_death]
patients_data$death_event <- ifelse(patients_data$vital_status == "Alive", 0, 1)
patients_data$EFGR_count = ifelse(gene_counts > 2905, "High", "Low")
low_counts = patients_data[patients_data$EFGR_count == "Low",]

surv_object <- Surv(time = low_counts$days_to_death, 
                    event = low_counts$death_event)

gender_fit <- surv_fit(surv_object ~ low_counts$gender, data = low_counts)

survplot = ggsurvplot(gender_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

# figure 2 - km plot
p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p

# do another survival plot, now only for those with significantly high EGFR counts
# patients_data$EFGR_count = ifelse(gene_counts > 3824, "High", "Low")

high_counts = patients_data[patients_data$EFGR_count == "High",]
new_surv_object <- Surv(time = high_counts$days_to_death, 
                        event = high_counts$death_event)
gender_counts_fit <- surv_fit(new_surv_object ~ high_counts$gender, data = high_counts)

egfr_survplot = ggsurvplot(gender_counts_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

# figure 3 - km plot for those with high counts, > third quartile
p2 = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p2




clin_query <- GDCquery(project = "TCGA-COAD", data.category = "Clinical", file.type = "xml")
# GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
names(clinic)[names(clinic)=="days_to_last_followup"] <- "days_to_last_follow_up"
colnames(clinic)[1] <- "Tumor_Sample_Barcode"
clinic

maf_dataframe <- data.table::fread("~/qbio490/analysis_data/GDCdata/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv",
                                   data.table = F)

mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)
maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)
clinic = maf_object@clinical.data
female_patient_ids = c(clinic$Tumor_Sample_Barcode[clinic$gender == "FEMALE"])
female_maf = subsetMaf(maf = maf_object,
                      tsb = female_patient_ids)

male_maf = subsetMaf(maf = maf_object,
                    tsb = c(clinic$Tumor_Sample_Barcode[clinic$gender == "MALE"]))

coOncoplot(m1 = female_maf, 
           m2 = male_maf, 
           m1Name = "MAF data of female patients", 
           m2Name = "MAF data of male patients")

lollipopPlot2(m1 = female_maf, 
              m2 = male_maf, 
              m1_name = "MAF data of female patients",
              m2_name = "MAF data of male patients",
              gene = "EGFR")







