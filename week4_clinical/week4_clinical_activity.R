setwd("~/qbio490/qbio_data_analysis_echo/analysis_data/")
clinic_read_in <- read.csv("~/qbio490/qbio_data_analysis_echo/analysis_data/clinic.csv")

#Written activity
#1. A categorical variable is a variable that is qualitative. It assigns an observation to a particular group or category. 
#An example of a categorical variable is sex. A discrete variable is a quantitative, whole number variable. An example of 
#a discrete variable is number of days since cancer diagnosis. A continuous variable is a variable where the value has to be
#measured and is not a concrete whole number. An example of a continuous variable is height. 

#2
colnames(clinic_read_in)
is.na(clinic_read_in$number_of_first_degree_relatives_with_cancer_diagnosis)
sum(is.na(clinic_read_in$number_of_first_degree_relatives_with_cancer_diagnosis))
#variable: number of first degree relatives with cancer diagnosis

#3
#The variable is discrete because it counts number of people. There can only be a whole number of people as a possile
#observation, so the variable cannot be categorical (qualitative and grouping) or continuous (any real number, not just
#whole numbers)

#4
#https://ashpublications.org/blood/article/134/12/960/374916/Analysis-of-153-115-patients-with-hematological
#The study found that for those with blood cancer, their first-degree relatives (parent, child, 
#or sibling) were at a higher risk of also developing blood cancer. The study found that the risk of developing blood 
#cancer is partly inherited from DNA changes passed on from parents. 
#https://www-clinicalkey-com.libproxy2.usc.edu/#!/content/playContent/1-s2.0-S095980490900344X
#The study found that family history is a very significant risk factor for developing small-cell lung cancer.
#Heritable cancers are usually marked by early onset of cancer, and the study showed that patients had a 4.6x higher
#likelihood of developing cancer if they have a first-degree relative with small-cell lung cancer.

#5
is.na(clinic_read_in$race_list)
#Chosen variable: race. It's a categorical variable that can be measured through surveys and/or DNA ancestry tests.
#This is a categorical variable because it is a qualitative variable in which people can be grouped for analysis.

#6
#1) Patients from minority groups might have more first degree relatives with cancer diagnosis because they are more likely to be in
#environmental conditions that pose risks to their health, including proximity to pollution and lack of access to health care. 
#2) Patients with more first-degree relatives with cancer diagnosis might have a lower survival rate for colorectal cancer. 
#3) Patients from minority groups might have a lower survival rate for colorectal cancer. 

#7


#Coding
#1
#create a boxplot comparing number of first-degree relatives with cancer diagnosis per race
par(mar=c(15,2,1,1))
boxplot(clinic_read_in$number_of_first_degree_relatives_with_cancer_diagnosis ~ clinic_read_in$race_list, xlab = "Race", 
        ylab = "Number of first degree relatives with cancer diagnosis", las = 2)

library(survival)
library(survminer)

#predict NA days to death with days to last follow up
clinic_read_in$days_to_death[is.na(clinic_read_in$days_to_death)] <- clinic_read_in$days_to_last_follow_up[is.na(clinic_read_in$days_to_death)]
#define a column saying whether or not death occurred
clinic_read_in$death_event <- ifelse(clinic_read_in$vital_status == "Alive", 0, 1)
#create a survival object with days to death and whether or not death occurred
surv_object <- Surv(time = clinic_read_in$days_to_death, 
                    event = clinic_read_in$death_event)
#fit race data to survival object
race_fit <- surv_fit( surv_object ~ clinic_read_in$race_list, data = clinic )

#plot survival data for CRC according to race
survplot = ggsurvplot(race_fit, 
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

#Repeat but with first-degree cancer diagnosis data
relative_fit <- surv_fit( surv_object ~ clinic_read_in$number_of_first_degree_relatives_with_cancer_diagnosis, data = clinic )

#plot survival data for CRC according to race
survplot = ggsurvplot(relative_fit, 
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

