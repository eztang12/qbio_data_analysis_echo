if(!require(BiocManager)) install.packages("BiocManager")

# the double colon syntax calls a function from a specific package
# this avoids loading the entire package
# in this case, we need to download TCGAbiolinks from Bioconductor using BiocManager
if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")

# this just loads a package
library(TCGAbiolinks)
getwd()
setwd("~/qbio490/qbio_data_analysis_echo/analysis_data/")
clin_query <- GDCquery(project = "TCGA-COAD", data.category = "Clinical", file.type = "xml")
# Only use this line ONCE! Comment out after you have downloaded the data. 
#GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
# Just adding an underscore between follow and up
names(clinic)[names(clinic)=="days_to_last_followup"] <- "days_to_last_follow_up"

str(clinic)
head(clinic)

#exercise 1.1
#str displays the structure of an object, similar to summary
#head displays the first lines

#exercise 1.2
colnames(clinic)
View(clinic$vital_status)

#exercise 2.1
plot(clinic$age_at_initial_pathologic_diagnosis, clinic$weight, xlab = "Age at initial pathologic diagnosis", ylab = "Weight")

#exercise 2.2
clinic$race_list <- as.character(clinic$race_list)
unique(clinic$race_list)
boxplot(clinic$age_at_initial_pathologic_diagnosis ~ clinic$race_list, xlab = "Race", ylab = "Age at initial pathologic diagnosis")

#exercise 2.3
miss_race <- which(clinic$race_list == "")
clinic$race_list[miss_race] = "No data"

#exercise 2.4
min(clinic$age_at_initial_pathologic_diagnosis)
max(clinic$age_at_initial_pathologic_diagnosis)
mean(clinic$age_at_initial_pathologic_diagnosis)
median(clinic$age_at_initial_pathologic_diagnosis)
summary(clinic$age_at_initial_pathologic_diagnosis)

#exercise 2.5
length(which(clinic$age_at_initial_pathologic_diagnosis < 50))
length(which(clinic$age_at_initial_pathologic_diagnosis >= 50))
young <- which(clinic$age_at_initial_pathologic_diagnosis < 50)
old <- which(clinic$age_at_initial_pathologic_diagnosis >= 50)

#exercise 2.6
young_patient_ids <- clinic$bcr_patient_barcode[young]
old_patient_ids <- clinic$bcr_patient_barcode[old]

#exercise 2.7
clinic$age_category <- ifelse(clinic$age_at_initial_pathologic_diagnosis < 50, "young", "old")

#exercise 2.8
clinic [1,1] #this is the top left entry of the dataframe. This is patient barcode TCGA-3L-AA1B
clinic[1,] #this is the entire first row. This is all the patient information for patient TCGA-3L-AA1B
clinic[2:5,] #this is the entire second to fifth rows. This is all the patient information for patients TCGA-4N-A93T,
# TCGA-4T-AA8H, TCGA-5M-AAT4, and TCGA-5M-AAT6
clinic[,3] #this is the tumor tissue site data for all patients in the dataset

#exercise 2.9
young_clinic <- data.frame(clinic[clinic$age_category == "young",])
old_clinic <- data.frame(clinic[clinic$age_category == "old",])

#exercise 2.10
young_clinic_one_line <- data.frame(clinic[clinic$age_at_initial_pathologic_diagnosis < 50,])
identical(dim(young_clinic), dim(young_clinic_one_line))

install.packages("survival")
install.packages("survminer")
library(survival)
library(survminer)

#exercise 3.1
miss_death <- which(is.na(clinic$days_to_death))
clinic$days_to_death[miss_death] <- clinic$days_to_last_follow_up[miss_death]

#exercise 3.2
View(clinic$vital_status)
clinic$death_event <- ifelse(clinic$vital_status == "Alive", 0, 1)

# We initialize a 'survival' object first, which contains the data we need.
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)

# We then create a fit object
age_fit <- surv_fit( surv_object ~ clinic$age_category, data = clinic )

#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot = ggsurvplot(age_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p

# save the plot as a png
# if you want to present the data, change the strata labels!
ggsave("../week4_clinical/kmplot_by_age.png", plot = p, width = 12, height = 9)


#One main takeaway is how age affects survival rates. While both survival rates go down over over time, older patients' survival rates
#decrease at a much higher rate than younger patients. Some questions it would pose is what about age makes this occur: for example,
#is it due to increased risk of other disease? Weakening of related body systems? Later diagnosis time? 
#One thing to improve may be specificity. The time axis has no specifics on what unit of time (although we know that it is by days
#from the clinical dataframe, viewers who only see the plot will not be able to know). Furthermore, having the legend include what
#p refers to would also be helpful.

