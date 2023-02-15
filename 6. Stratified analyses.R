###################  1. Age<70 vs age>=70 ###################################### 

# OS 
# 1. cr13
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr13.RData")
dput(names(stracr13))
stra<-subset(stracr13, select = c("id", "Diagnosis_year", "Age", "Gender", 
                                  "chemradther", "Stage_at_diagnosis","Location", "MSI_gent", "brafmut", "krasmut", 
                                  "cimphi"))


stracr13<-merge(stra, unicox_cr13, by = 'id')

summary(stracr13$Stage_at_diagnosis)

save(stracr13,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")

# Age <70
agey<-stracr13[stracr13$Age<70, ] ### 
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = agey, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=agey, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'Age<70'
  size = nrow(agey)
  
})


## for age > 79
ageo<-stracr13[stracr13$Age>=70, ] ### 382
####### OS
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = ageo, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=ageo, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'Age>=70'
  size = nrow(ageo)
  
})

# bind results
sagecr13os<-rbind(AUC_ally, AUC_allo)
save(sagecr13os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

## PFS
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")
stracr13<-stracr13[!is.na(stracr13$recurr_cp), ]

## for age <79
agey<-stracr13[stracr13$Age<70, ] 

covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = agey, y=TRUE, x = TRUE)})

overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=agey, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allypm<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DFS'
  Group<-'Age<70'
  size = nrow(agey)
  
})


# old people
ageo<-stracr13[stracr13$Age>=70, ] 

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = ageo, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=ageo, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allopm<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DFS'
  Group<-'Age>70'
  size = nrow(ageo)
})

sagecr13dfs<-rbind(AUC_allypm, AUC_allopm)
resave(sagecr13dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

################## sex ##########################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")

# female
female<-stracr13[stracr13$Gender=='Female', ] ### 
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = female, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=female, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'Female'
  size = nrow(female)
  
})


male<-stracr13[stracr13$Gender=='Male', ]

stracr13[stracr13$Gender=='Male', ]

### 382
####### OS
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = male, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=male, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'Male'
  size = nrow(male)
  
})

# bind results
sexcr13os<-rbind(AUC_ally, AUC_allo)
resave(sexcr13os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

########## PFS ##########
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")
stracr13<-stracr13[!is.na(stracr13$recurr_cp), ]

## for female
female<-stracr13[stracr13$Gender=='Female', ] 

covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = female, y=TRUE, x = TRUE)})

overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=female, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allypm<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DFS'
  Group<-'Female'
  size = nrow(female)
  
})


# male
male<-stracr13[stracr13$Gender=='Male', ] 

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = male, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=male, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allopm<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DFS'
  Group<-'Male'
  size = nrow(male)
})

sexcr13dfs<-rbind(AUC_allypm, AUC_allopm)
resave(sexcr13dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")


####### location proximal vs distal vs rectum ################### 
########## OS #########
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")

sg12<-subset(stracr13, Location=='Distal colon') 

covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg12, y=TRUE, x = TRUE)})

overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=sg12, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'Distal colon'
  size = nrow(sg12)
  
})


sg3<-subset(stracr13, Location=='Proximal colon')
####### OS
covariates<-c("Gundert_2019_CR13")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg3, y=TRUE, x = TRUE)})

overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=sg3, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'Proximal colon'
  size = nrow(sg3)
  
  
})

sg4<-subset(stracr13, Location=='Rectum')

covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg4, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=sg4, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_all4<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'Rectum'
  size = nrow(sg4)
  
  
})

locationcr13os<-rbind(AUC_ally, AUC_allo, AUC_all4)
resave(locationcr13os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

###################### DFS #################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")

stracr13<-stracr13[!is.na(stracr13$recurr_cp), ]

sg12<-subset(stracr13,  Location=='Distal colon') ### 
## PFS
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg12, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=sg12, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DFS'
  Group<-'Distal colon'
  size = nrow(sg12)
  
})

sg3<-subset(stracr13, Location=='Proximal colon') ### 

## PFS
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg3, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=sg3, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'Proximal colon'
  size = nrow(sg3)
  
  
})

################### location IV
sg4<-subset(stracr13, Location=='Rectum')

covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg4, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=sg4, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_all4<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'Rectum'
  size = nrow(sg4)
  
})

locationcr13dfs<-rbind(AUC_ally, AUC_allo, AUC_all4)
resave(locationcr13dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

####### 5. no treatment vs with treatment #####
### OS ##
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")

summary(stracr13$chemradther)
treat<-subset(stracr13, chemradther=='Yes') ### 
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = treat, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=treat, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'treat'
  size = nrow(treat)
  
})


### notreat
notreat<-subset(stracr13, chemradther=='No')### 1769
####### OS
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = notreat, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=notreat, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'notreat'
  size = nrow(notreat)
  
  
})

# bind results
treatcr13os<-rbind(AUC_ally, AUC_allo)
resave(treatcr13os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")


########### DFS #####
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")
stracr13<-stracr13[!is.na(stracr13$recurr_cp), ]
###  treat
treat<-subset(stracr13, chemradther=='Yes') ### 
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = treat, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=treat, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'treat'
  size = nrow(treat)
  
})

### notreat
notreat<-subset(stracr13, chemradther=='No') ### 
## PFS
covariates<-c("Gundert_2019_CR13")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = notreat, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=notreat, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'notreat'
  size = nrow(notreat)
  
  
})

# bind results
treatcr13dfs<-rbind(AUC_ally, AUC_allo)
resave(treatcr13dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

################## 6. MSI vs MSS ##########
### OS ##
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")
summary(stracr13$MSI_gent)
# MSI
MSI<-subset(stracr13, MSI_gent=='MSI') ### 223

covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = MSI, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=MSI, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'MSI'
  size = nrow(MSI)
  
})


## MSS
MSS<-subset(stracr13, MSI_gent=='MSS')### 1769
####### OS
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = MSS, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=MSS, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'MSS'
  size = nrow(MSS)
  
  
})

# bind results
msicr13os<-rbind(AUC_ally, AUC_allo)
resave(msicr13os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")


## DSS
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")
stracr13<-stracr13[!is.na(stracr13$recurr_cp), ]
# MSI
MSI<-subset(stracr13, MSI_gent=='MSI') ### 

covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = MSI, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=MSI, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'MSI'
  size = nrow(MSI)
  
})

## MSS
MSS<-subset(stracr13, MSI_gent=='MSS') ### 
## PFS
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = MSS, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=MSS, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'MSS'
  size = nrow(MSS)
  
})

msi13dfs<-rbind(AUC_ally, AUC_allo)
resave(msi13dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

############# 7. BRAF mutation ##########
### OS ##
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")
summary(stracr13$brafmut)
## mut
mut<-subset(stracr13, brafmut=='BRAF mut') ### 
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mut, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=mut, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'mut'
  size = nrow(mut)
  
})

## wt
wt<-subset(stracr13, brafmut=='BRAF wt')### 1769
####### OS
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = wt, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=wt, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'wt'
  size = nrow(wt)
  
})

brafcr13os<-rbind(AUC_ally, AUC_allo)
resave(brafcr13os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")


## DFS
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")
stracr13<-stracr13[!is.na(stracr13$recurr_cp), ]
## mut
mut<-subset(stracr13, brafmut=='BRAF mut') ### 
## PFS
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mut, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=mut, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'mut'
  size = nrow(mut)
  
})

## wt
wt<-subset(stracr13, brafmut=='BRAF wt') ### 
## PFS
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = wt, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=wt, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'wt'
  size = nrow(wt)
  
})

brafcr13dfs<-rbind(AUC_ally, AUC_allo)
resave(brafcr13dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")


########## 8. KRAS mutation ##########
### OS ##
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")
summary(stracr13$krasmut)
## mut
mut<-subset(stracr13, krasmut=='KRAS mut') ### 
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mut, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=mut, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'mut'
  size = nrow(mut)
  
})

## wt
wt<-subset(stracr13, krasmut=='KRAS wt')### 1769
####### OS
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = wt, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=wt, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'wt'
  size = nrow(wt)
  
  
})

krascr13os<-rbind(AUC_ally, AUC_allo)
resave(krascr13os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

## DFS
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")
stracr13<-stracr13[!is.na(stracr13$recurr_cp), ]
## mut
mut<-subset(stracr13, krasmut=='KRAS mut') ### 
## PFS
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mut, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=mut, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'mut'
  size = nrow(mut)
  
})


## wt
wt<-subset(stracr13, krasmut=='KRAS wt') ### 
## PFS
covariates<-c("Gundert_2019_CR13")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = wt, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=wt, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'wt'
  size = nrow(wt)
  
  
})

krascr13dfs<-rbind(AUC_ally, AUC_allo)
resave(krascr13dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

######### 9. CIMP mutation ##################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")
summary(stracr13$cimphi)
high<-subset(stracr13, cimphi=='CIMP-high') ### 
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = high, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=high, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'high'
  size = nrow(high)
  
})


######## low
low<-subset(stracr13, cimphi=='CIMP-low')### 1769
####### OS
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = low, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=low, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'low'
  size = nrow(low)
  
  
})

cimpcr13os<-rbind(AUC_ally, AUC_allo)
resave(cimpcr13os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

## DFS
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/stracr13.RData")
stracr13<-stracr13[!is.na(stracr13$recurr_cp), ]
high<-subset(stracr13, cimphi=='CIMP-high') ### 

covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = high, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=high, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'high'
  size = nrow(high)
})

######## low
low<-subset(stracr13, cimphi=='CIMP-low') ### 
####### PFS
covariates<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = low, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=low, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'low'
  size = nrow(low)
  
})

cimpcr13dfs<-rbind(AUC_ally, AUC_allo)
resave(cimpcr13dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

## bind all results together
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

os<-rbind(sagecr13os,sexcr13os, locationcr13os, treatcr13os, msicr13os, brafcr13os, krascr13os, cimpcr13os )
dfs<-rbind(sagecr13dfs,sexcr13dfs, locationcr13dfs, treatcr13dfs, msi13dfs, brafcr13dfs, krascr13dfs, cimpcr13dfs )

cr13<-cbind(os, dfs)

write.csv(cr13, 
          "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/tablesfigures/stra_cr13.csv")

######################### c14 ############################
###################  1. Age<70 vs age>=70 ###################################### 

# OS 
# 1. c14
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c14.RData")

stra<-subset(stradata, select = c("id", "Diagnosis_year", "Age", "Gender", 
                                  "chemradther", "Stage_at_diagnosis","Location", "MSI_gent", "brafmut", "krasmut", 
                                  "cimphi"))

strac14<-merge(stra, unicox_c14, by = 'id')

summary(strac14$Stage_at_diagnosis)

save(strac14,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")

dput(names(strac14))

# Age <70
agey<-strac14[strac14$Age<70, ] ### 
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = agey, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=agey, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'Age<70'
  size = nrow(agey)
  
})


## for age > 79
ageo<-strac14[strac14$Age>=70, ] ### 382
####### OS
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = ageo, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=ageo, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'Age>=70'
  size = nrow(ageo)
  
})

# bind results
sagec14os<-rbind(AUC_ally, AUC_allo)
resave(sagec14os, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

## PFS
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")
strac14<-strac14[!is.na(strac14$recurr_cp), ]

## for age <79
agey<-strac14[strac14$Age<70, ] 

covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = agey, y=TRUE, x = TRUE)})

overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=agey, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allypm<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DFS'
  Group<-'Age<70'
  size = nrow(agey)
  
})

# old people
ageo<-strac14[strac14$Age>=70, ] 

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = ageo, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=ageo, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allopm<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DFS'
  Group<-'Age>70'
  size = nrow(ageo)
})

sagec14dfs<-rbind(AUC_allypm, AUC_allopm)
resave(sagec14dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

################## sex ##########################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")

# female
female<-strac14[strac14$Gender=='Female', ] ### 
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = female, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=female, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'Female'
  size = nrow(female)
  
})


male<-strac14[strac14$Gender=='Male', ]

strac14[strac14$Gender=='Male', ]

### 382
####### OS
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = male, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=male, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'Male'
  size = nrow(male)
  
})

# bind results
sexc14os<-rbind(AUC_ally, AUC_allo)
resave(sexc14os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

########## PFS ##########
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")
strac14<-strac14[!is.na(strac14$recurr_cp), ]

## for female
female<-strac14[strac14$Gender=='Female', ] 

covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = female, y=TRUE, x = TRUE)})

overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=female, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allypm<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DFS'
  Group<-'Female'
  size = nrow(female)
  
})


# male
male<-strac14[strac14$Gender=='Male', ] 

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = male, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=male, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allopm<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DFS'
  Group<-'Male'
  size = nrow(male)
})

sexc14dfs<-rbind(AUC_allypm, AUC_allopm)
resave(sexc14dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

####### TNM stage ################### 
############# 3. stage I/II vs stage III/IV ###################################### 
#### stage I II vs, III vs IV #
### OS
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")
summary(strac14$Stage_at_diagnosis)

sg12<-subset(strac14, 
             Stage_at_diagnosis=='I'|Stage_at_diagnosis=='II') ### 

covariates<-c("WangY_2020_C14")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg12, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=sg12, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'I & II'
  size = nrow(sg12)
  
})


## III 
sg3<-subset(strac14, Stage_at_diagnosis=='III')
####### OS
covariates<-c("WangY_2020_C14")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg3, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=sg3, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'III'
  size = nrow(sg3)
  
  
})


## IV
sg4<-subset(strac14, Stage_at_diagnosis=='IV')
####### OS
covariates<-c("WangY_2020_C14")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg4, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=sg4, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_all4<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'IV'
  size = nrow(sg4)
  
})

stagec14os<-rbind(AUC_ally,AUC_allo, AUC_all4)
resave(stagec14os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")


## DFS
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")
strac14<-strac14[!is.na(strac14$recurr_cp), ]

sg12<-subset(strac14,  Stage_at_diagnosis=='I'|Stage_at_diagnosis=='II') ### 
## PFS
covariates<-c("WangY_2020_C14")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg12, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=sg12, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'I & II'
  size = nrow(sg12)
  
})


sg3<-subset(strac14, Stage_at_diagnosis=='III') ### 
## PFS
covariates<-c("WangY_2020_C14")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg3, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=sg3, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'III'
  size = nrow(sg3)
  
  
})

################### stage IV
sg4<-subset(strac14, Stage_at_diagnosis=='IV')

covariates<-c("WangY_2020_C14")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg4, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=sg4, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_all4<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'IV'
  size = nrow(sg4)
})

stagec14dfs<-rbind(AUC_ally,AUC_allo, AUC_all4)
resave(stagec14dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

####### 5. no treatment vs with treatment #####
### OS ##
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")

summary(strac14$chemradther)
treat<-subset(strac14, chemradther=='Yes') ### 
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = treat, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=treat, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'treat'
  size = nrow(treat)
  
})


### notreat
notreat<-subset(strac14, chemradther=='No')### 1769
####### OS
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = notreat, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=notreat, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'notreat'
  size = nrow(notreat)
  
  
})

# bind results
treatc14os<-rbind(AUC_ally, AUC_allo)
resave(treatc14os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")


########### DFS #####
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")
strac14<-strac14[!is.na(strac14$recurr_cp), ]
###  treat
treat<-subset(strac14, chemradther=='Yes') ### 
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = treat, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=treat, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'treat'
  size = nrow(treat)
  
})

### notreat
notreat<-subset(strac14, chemradther=='No') ### 
## PFS
covariates<-c("WangY_2020_C14")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = notreat, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=notreat, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'notreat'
  size = nrow(notreat)
  
  
})

# bind results
treatc14dfs<-rbind(AUC_ally, AUC_allo)
resave(treatc14dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

################## 6. MSI vs MSS ##########
### OS ##
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")
summary(strac14$MSI_gent)
# MSI
MSI<-subset(strac14, MSI_gent=='MSI') ### 223

covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = MSI, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=MSI, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'MSI'
  size = nrow(MSI)
  
})


## MSS
MSS<-subset(strac14, MSI_gent=='MSS')### 1769
####### OS
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = MSS, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=MSS, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'MSS'
  size = nrow(MSS)
  
  
})

# bind results
msic14os<-rbind(AUC_ally, AUC_allo)
resave(msic14os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")


## DSS
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")
strac14<-strac14[!is.na(strac14$recurr_cp), ]
# MSI
MSI<-subset(strac14, MSI_gent=='MSI') ### 

covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = MSI, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=MSI, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'MSI'
  size = nrow(MSI)
  
})

## MSS
MSS<-subset(strac14, MSI_gent=='MSS') ### 
## PFS
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = MSS, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=MSS, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'MSS'
  size = nrow(MSS)
  
})

msi13dfs<-rbind(AUC_ally, AUC_allo)
resave(msi13dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

############# 7. BRAF mutation ##########
### OS ##
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")
summary(strac14$brafmut)
## mut
mut<-subset(strac14, brafmut=='BRAF mut') ### 
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mut, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=mut, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'mut'
  size = nrow(mut)
  
})

## wt
wt<-subset(strac14, brafmut=='BRAF wt')### 1769
####### OS
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = wt, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=wt, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'wt'
  size = nrow(wt)
  
})

brafc14os<-rbind(AUC_ally, AUC_allo)
resave(brafc14os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")


## DFS
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")
strac14<-strac14[!is.na(strac14$recurr_cp), ]
## mut
mut<-subset(strac14, brafmut=='BRAF mut') ### 
## PFS
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mut, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=mut, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'mut'
  size = nrow(mut)
  
})

## wt
wt<-subset(strac14, brafmut=='BRAF wt') ### 
## PFS
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = wt, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=wt, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'wt'
  size = nrow(wt)
  
})

brafc14dfs<-rbind(AUC_ally, AUC_allo)
resave(brafc14dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")


########## 8. KRAS mutation ##########
### OS ##
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")
summary(strac14$krasmut)
## mut
mut<-subset(strac14, krasmut=='KRAS mut') ### 
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mut, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=mut, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'mut'
  size = nrow(mut)
  
})

## wt
wt<-subset(strac14, krasmut=='KRAS wt')### 1769
####### OS
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = wt, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=wt, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'wt'
  size = nrow(wt)
  
  
})

krasc14os<-rbind(AUC_ally, AUC_allo)
resave(krasc14os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

## DFS
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")
strac14<-strac14[!is.na(strac14$recurr_cp), ]
## mut
mut<-subset(strac14, krasmut=='KRAS mut') ### 
## PFS
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mut, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=mut, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'mut'
  size = nrow(mut)
  
})


## wt
wt<-subset(strac14, krasmut=='KRAS wt') ### 
## PFS
covariates<-c("WangY_2020_C14")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = wt, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=wt, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'wt'
  size = nrow(wt)
  
  
})

krasc14dfs<-rbind(AUC_ally, AUC_allo)
resave(krasc14dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

######### 9. CIMP mutation ##################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")
summary(strac14$cimphi)
high<-subset(strac14, cimphi=='CIMP-high') ### 
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = high, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=high, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'high'
  size = nrow(high)
  
})


######## low
low<-subset(strac14, cimphi=='CIMP-low')### 1769
####### OS
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = low, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=low, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'OS'
  Group<-'low'
  size = nrow(low)
  
  
})

cimpc14os<-rbind(AUC_ally, AUC_allo)
resave(cimpc14os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

## DFS
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/strac14.RData")
strac14<-strac14[!is.na(strac14$recurr_cp), ]
high<-subset(strac14, cimphi=='CIMP-high') ### 

covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = high, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=high, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_ally<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'high'
  size = nrow(high)
})

######## low
low<-subset(strac14, cimphi=='CIMP-low') ### 
####### PFS
covariates<-c("WangY_2020_C14")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = low, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=low, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_allo<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'low'
  size = nrow(low)
  
})

cimpc14dfs<-rbind(AUC_ally, AUC_allo)
resave(cimpc14dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

## bind all results together
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/stratify.RData")

os<-rbind(sagec14os,sexc14os, stagec14os,  treatc14os, msic14os, brafc14os, krasc14os, cimpc14os )
dfs<-rbind(sagec14dfs,sexc14dfs, stagec14dfs,  treatc14dfs, msi13dfs, brafc14dfs, krasc14dfs, cimpc14dfs )

c14<-cbind(os, dfs)

write.csv(c14, 
          "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/tablesfigures/stra_c14.csv")

######################### c14 ############################




