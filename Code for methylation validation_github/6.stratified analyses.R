## prepare data #####
# overall AUC and Brier score, the five models, unimputed 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/covariable_all.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/unicox_all.RData")

summary(covar)
summary(unicox)
dput(names(unicox))
score<-subset(unicox, select = c( "id", "TYang_2019_OS", "TGong_2020_OS", "TWangX_2020_OS_PFS", 
                                  "TWangY_2020_OS", "TXiang_2020_OS", "TChen_2021_OS", "THuang_2021_OS", 
                                  "THuang_2021_PFS", "TLi_2021_OS", "TGundert_2019_OS_DSS", "Gundert_2019_OS_DSS", 
                                  "Li_2021_OS", "Huang_2021_PFS", "Huang_2021_OS", "Chen_2021_OS", 
                                  "Xiang_2020_OS", "WangY_2020_OS", "WangX_2020_OS_PFS", "Gong_2020_OS", 
                                  "Yang_2019_OS"))
stradata<-merge(covar, score, by='id')
summary(stradata)
save(stradata, file = ("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData"))

# make the data for PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/PMunicox.RData")

stradatapm<-stradata[(stradata$id %in% unicoxpm$id), ]

save(stradatapm, file = ("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData"))

###################  1. Age<70 vs age>=70 ###################################### 

# OS 
# exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")

# Age <70
agey<-stradata[stradata$Age<70, ] ### 1928
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


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
ageo<-stradata[stradata$Age>=70, ] ### 382
####### OS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

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



## for PM classifier 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")

# Age <70
agey<-stradatapm[stradatapm$Age<70, ] ### 1447
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = agey, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=agey, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'Age<70'
  size = nrow(agey)
  
})



## for age > 79
ageo<-stradatapm[stradatapm$Age>=70, ] ### 291
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = ageo, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=ageo, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'Age>=70'
  size = nrow(ageo)
  
  
})


## bind all results together 
Sageaucos<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm)
Sageaucos<-Sageaucos[order(Sageaucos$model, Sageaucos$Group), ]


write_xlsx(Sageaucos, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Sageaucos.xlsx", col_names = T)


### DFS
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
stradata<-stradata[!is.na(stradata$recurr_cp), ]
# Age <70
agey<-stradata[stradata$Age<70, ] ### 1909

covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = agey, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=agey, metrics = c("auc", "brier"))

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
  Group<-'Age<70'
  size = nrow(agey)
  
})



## for age > 79
ageo<-stradata[stradata$Age>=70, ] ### 

## PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = ageo, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(recurr_timey, recurr_cp==1)~1,data=ageo, metrics = c("auc", "brier"))

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
  Outcome<-'DFS'
  Group<-'Age>=70'
  size = nrow(ageo)
  
  
})


## for PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")
stradatapm<-stradatapm[!is.na(stradatapm$recurr_cp), ]
# Age <70
agey<-stradatapm[stradatapm$Age<70, ] ### 1437
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = agey, y=TRUE, x = TRUE)})

# overall AUC
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


## for age > 79
ageo<-stradatapm[stradatapm$Age>=70, ] ### 290
## PFS
covariates<-c("Gundert_2019_OS_DSS")

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
  Group<-'Age>=70'
  size = nrow(ageo)
  
  
})


## bind all results together 
SageaucDSS<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm)
SageaucDSS<-SageaucDSS[order(SageaucDSS$model, SageaucDSS$Group), ]

write_xlsx(SageaucDSS, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/SageaucDFS.xlsx", col_names = T)


### bind all results together ###
Sageaucos <- read_excel("Validation_tablefigure/Sageaucos.xlsx")
Sageaucos<-Sageaucos[, c(1,3,2, 5)]
colnames(Sageaucos)[3:4]<-c('size OS', 'AUC (95%CI) OS')

SageaucDFS <- read_excel("Validation_tablefigure/SageaucDFS.xlsx")
dput(names(SageaucDFS))
SageaucDFS<-subset(SageaucDFS, select = c("size", "AUCCI"))
colnames(SageaucDFS)<-c('size DFS', 'AUC (95%CI) DFS')

Sage<-bind_cols(Sageaucos, SageaucDFS)

write_xlsx(Sage, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Strage70.xlsx", col_names = T) 


 ###################  2. Male vs Female ###################################### 
### OS ##
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
summary(stradata$Gender)
 #Male 
male<-stradata[stradata$Gender=='Male', ] ### 1346
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = male, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=male, metrics = c("auc", "brier"))

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
  Group<-'Male'
  size = nrow(male)
  
})



# female
female<-stradata[stradata$Gender=='Female', ] ### 964
## OS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")
)

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = female, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=female, metrics = c("auc", "brier"))

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
  Group<-'Female'
  size = nrow(female)
  
})




## for PM classifier 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")

male<-stradatapm[stradatapm$Gender =='Male', ] ### 1039
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = male, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=male, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'Male'
  size = nrow(male)
  
})



## female 
female<-stradatapm[stradatapm$Gender=='Female', ] ### 699
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = female, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=female, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'Female'
  size = nrow(female)
  
  
})



## bind all results together 
Ssexaucos<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm)
Ssexaucos<-Ssexaucos[order(Ssexaucos$model, Ssexaucos$Group), ]


write_xlsx(Ssexaucos, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Ssexaucos.xlsx", col_names = T)


## DFS
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
stradata<-stradata[!is.na(stradata$recurr_cp), ]
 # male
male<-stradata[stradata$Gender=='Male', ] ### 1332
## PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = male, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=male, metrics = c("auc", "brier"))

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
  Group<-'Male'
  size = nrow(male)
  
})



### female  
female<-stradata[stradata$Gender == 'Female', ] ### 955
## PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = female, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=female, metrics = c("auc", "brier"))

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
  Group<-'Female'
  size = nrow(female)
  
  
})


## for PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")
stradatapm<-stradatapm[!is.na(stradatapm$recurr_cp), ]
 ## male
male<-stradatapm[stradatapm$Gender =='Male', ] ### 1034
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = male, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=male, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'Male'
  size = nrow(male)
  
})



## female
female<-stradatapm[stradatapm$Gender =='Female', ] ### 693
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = female, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=female, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'Female'
  size = nrow(female)
  
  
})



## bind all results together 
SsexaucDSS<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm)
SsexaucDSS<-SsexaucDSS[order(SsexaucDSS$model, SsexaucDSS$Group), ]
write_xlsx(SsexaucDSS, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/SsexaucDFS.xlsx", col_names = T)


Ssexaucos <- read_excel("Validation_tablefigure/Ssexaucos.xlsx")
Ssexaucos<-Ssexaucos[, c(1,3,2, 5)]
colnames(Ssexaucos)[3:4]<-c('size OS', 'AUC (95%CI) OS')

SsexaucDFS <- read_excel("Validation_tablefigure/SsexaucDFS.xlsx")
dput(names(SageaucDFS))
SsexaucDFS<-subset(SsexaucDFS, select = c("size", "AUCCI"))
colnames(SsexaucDFS)<-c('size DFS', 'AUC (95%CI) DFS')

Ssex<-bind_cols(Ssexaucos, SsexaucDFS)
Ssex<-Ssex[order(Ssex$model, Ssex$Group), ]

write_xlsx(Ssex, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Ssex.xlsx", col_names = T)


############# 3. stage I/II vs stage III/IV ###################################### 
#### stage I II vs, III vs IV #
### OS
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
summary(stradata$Stage_at_diagnosis)

sg12<-subset(stradata, 
             Stage_at_diagnosis=='I'|Stage_at_diagnosis=='II') ### 

covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


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
sg3<-subset(stradata, Stage_at_diagnosis=='III')
####### OS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


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
sg4<-subset(stradata, Stage_at_diagnosis=='IV')
####### OS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


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


## for PM classifier 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")

sg12<-subset(stradatapm, Stage_at_diagnosis=='I'|Stage_at_diagnosis=='II') # 153

####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg12, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=sg12, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'I & II'
  size = nrow(sg12)
  
})



sg3<-subset(stradatapm, Stage_at_diagnosis=='III') 
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg3, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=sg3, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'III'
  size = nrow(sg3)
  
})


# stage IV 
sg4<-subset(stradatapm, Stage_at_diagnosis=='IV') 
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg4, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=sg4, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_all4pm<-within.data.frame(AUC_all, {
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


## bind all results together 
Stageaucos<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm, AUC_all4, AUC_all4pm)
Stageaucos<-Stageaucos[order(Stageaucos$model, Stageaucos$Group), ]


write_xlsx(Stageaucos, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Stageaucos.xlsx", col_names = T)



## DFS
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
stradata<-stradata[!is.na(stradata$recurr_cp), ]

sg12<-subset(stradata,  Stage_at_diagnosis=='I'|Stage_at_diagnosis=='II') ### 
## PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


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


sg3<-subset(stradata, Stage_at_diagnosis=='III') ### 
## PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


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
sg4<-subset(stradata, Stage_at_diagnosis=='IV')

covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


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



## for PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")
stradatapm<-stradatapm[!is.na(stradatapm$recurr_cp), ]

sg12<-subset(stradatapm, Stage_at_diagnosis=='I'|Stage_at_diagnosis=='II')### 
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg12, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=sg12, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'I & II'
  size = nrow(sg12)
  
})



sg3<-subset(stradatapm, Stage_at_diagnosis=='III') #
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg3, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=sg3, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'III'
  size = nrow(sg3)
  
})


### stage IV 
sg4<-subset(stradatapm, Stage_at_diagnosis=='IV') 
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg4, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=sg4, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_all4pm<-within.data.frame(AUC_all, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  AUC<-NULL
  lower<-NULL
  upper<-NULL
  times<-NULL
  Outcome<-'DSS'
  Group<-'IV'
  size = nrow(sg4)
  
})



## bind all results together 
Stageaucdds<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm, AUC_all4, AUC_all4pm)
Stageaucdds<-Stageaucdds[order(Stageaucdds$model, Stageaucdds$Group), ]

write_xlsx(Stageaucdds, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/StageaucDFS.xlsx", col_names = T)


### bind all results together ###
### bind all results together ###
Stageaucos <- read_excel("Validation_tablefigure/Stageaucos.xlsx")
Stageaucos<-Stageaucos[, c(1,3,2, 5)]
colnames(Stageaucos)[3:4]<-c('size OS', 'AUC (95%CI) OS')

StageaucDFS <- read_excel("Validation_tablefigure/StageaucDFS.xlsx")
dput(names(StageaucDFS))
StageaucDFS<-subset(StageaucDFS, select = c("size", "AUCCI"))
colnames(StageaucDFS)<-c('size DFS', 'AUC (95%CI) DFS')

Stage<-bind_cols(Stageaucos, StageaucDFS)


write_xlsx(Stage, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Stage.xlsx", col_names = T)


####### 4. location proximal vs distal vs rectum ################### 
### OS ##
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
summary(stradata$Location)

sg12<-subset(stradata, Location=='Distal colon') 

covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


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
  Group<-'Distal colon'
  size = nrow(sg12)
  
})



sg3<-subset(stradata, Location=='Proximal colon')
####### OS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


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
  Group<-'Proximal colon'
  size = nrow(sg3)
  
  
})


sg4<-subset(stradata, Location=='Rectum')
####### OS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

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


## for PM classifier 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")

sg12<-subset(stradatapm, Location=='Distal colon') # 153

####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg12, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=sg12, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'Distal colon'
  size = nrow(sg12)
  
})



sg3<-subset(stradatapm, Location=='Proximal colon') 
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg3, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=sg3, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'Proximal colon'
  size = nrow(sg3)
  
})


## location IV
sg4<-subset(stradata, Location=='Rectum')
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg4, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=sg4, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_all4pm<-within.data.frame(AUC_all, {
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



## bind all results together 
locationaucos<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm, AUC_all4, AUC_all4pm)
locationaucos<-locationaucos[order(locationaucos$model, locationaucos$Group), ]

write_xlsx(locationaucos, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/locationaucos.xlsx", col_names = T)

## DFS
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
stradata<-stradata[!is.na(stradata$recurr_cp), ]

sg12<-subset(stradata,  Location=='Distal colon') ### 
## PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


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
  Group<-'Distal colon'
  size = nrow(sg12)
  
})

sg3<-subset(stradata, Location=='Proximal colon') ### 

## PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

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
  Group<-'Proximal colon'
  size = nrow(sg3)
  
  
})

################### location IV
sg4<-subset(stradata, Location=='Rectum')

covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

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
  Outcome<-'DSS'
  Group<-'Rectum'
  size = nrow(sg4)
  
  
})


## for PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")
stradatapm<-stradatapm[!is.na(stradatapm$recurr_cp), ]

sg12<-subset(stradatapm, Location=='Distal colon')### 
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg12, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=sg12, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'Distal colon'
  size = nrow(sg12)
  
})

sg3<-subset(stradatapm, Location=='Proximal colon') #
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg3, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=sg3, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'Proximal colon'
  size = nrow(sg3)
  
})


## location IV
sg4<-subset(stradata, Location=='Rectum')
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = sg4, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=sg4, metrics = c("auc", "brier"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)

AUC_all4pm<-within.data.frame(AUC_all, {
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


## bind all results together 
locationaucdds<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm, AUC_all4, AUC_all4pm)
locationaucdds<-locationaucdds[order(locationaucdds$model, locationaucdds$Group), ]


write_xlsx(locationaucdds, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/locationaucdfs.xlsx", col_names = T)


### bind all results together ###
locationaucos <- read_excel("Validation_tablefigure/locationaucos.xlsx")
locationaucos<-locationaucos[, c(1,3,2, 5)]
colnames(locationaucos)[3:4]<-c('size OS', 'AUC (95%CI) OS')

locationaucdfs<- read_excel("Validation_tablefigure/locationaucdfs.xlsx")
dput(names(locationaucdfs))
locationaucdfs<-subset(locationaucdfs, select = c("size", "AUCCI"))
colnames(locationaucdfs)<-c('size DFS', 'AUC (95%CI) DFS')

Location3<-bind_cols(locationaucos, locationaucdfs)


write_xlsx(Location3, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Location3.xlsx", col_names = T)


####### 5. no treatment vs with treatment #####
### OS ##
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
summary(stradata$chemradther)
treat<-subset(stradata, chemradther=='Yes') ### 
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


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
notreat<-subset(stradata, chemradther=='No')### 1769
####### OS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

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


## for PM classifier 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")

treat<-subset(stradatapm, chemradther =='Yes') # 153

####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = treat, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=treat, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'treat'
  size = nrow(treat)
  
})



### notreat
notreat<-subset(stradatapm, chemradther=='No') 
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = notreat, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=notreat, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'notreat'
  size = nrow(notreat)
  
})


## bind all results together 
Streataucos<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm)
Streataucos<-Streataucos[order(Streataucos$model, Streataucos$Group), ]

write_xlsx(Streataucos, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Streataucos.xlsx", col_names = T)


## DFS
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
stradata<-stradata[!is.na(stradata$recurr_cp), ]
###  treat
treat<-subset(stradata, chemradther=='Yes') ### 
## PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = treat, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=treat, metrics = c("auc", "brier"))

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
notreat<-subset(stradata, chemradther=='No') ### 
## PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = notreat, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=notreat, metrics = c("auc", "brier"))

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


## for PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")
stradatapm<-stradatapm[!is.na(stradatapm$recurr_cp), ]

treat<-subset(stradatapm, chemradther=='Yes')### 
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = treat, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=treat, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'treat'
  size = nrow(treat)
  
})


### notreat
notreat<-subset(stradatapm, chemradther=='No') #
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = notreat, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=notreat, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'notreat'
  size = nrow(notreat)
  
})



## bind all results together 
StreataucDSS<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm)
StreataucDSS<-StreataucDSS[order(StreataucDSS$model, StreataucDSS$Group), ]

write_xlsx(StreataucDSS, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/StreataucDFS.xlsx", col_names = T)

### bind all results together ###
Streataucos <- read_excel("Validation_tablefigure/Streataucos.xlsx")
Streataucos<-Streataucos[, c(1,3,2, 5)]
colnames(Streataucos)[3:4]<-c('size OS', 'AUC (95%CI) OS')

StreataucDFS<- read_excel("Validation_tablefigure/StreataucDFS.xlsx")
dput(names(StreataucDFS))
StreataucDFS<-subset(StreataucDFS, select = c("size", "AUCCI"))
colnames(StreataucDFS)<-c('size DFS', 'AUC (95%CI) DFS')

Strtreat<-bind_cols(Streataucos, StreataucDFS)

write_xlsx(Strtreat, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Strtreat.xlsx", col_names = T)

################## 6. MSI vs MSS ##########
### OS ##
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
summary(stradata$MSI_gent)
# MSI
MSI<-subset(stradata, MSI_gent=='MSI') ### 223

covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

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
MSS<-subset(stradata, MSI_gent=='MSS')### 1769
####### OS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

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


## for PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")

MSI<-subset(stradatapm, MSI_gent =='MSI') # 153

####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = MSI, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=MSI, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'MSI'
  size = nrow(MSI)
  
})

## MSS
MSS<-subset(stradatapm, MSI_gent=='MSS') 
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = MSS, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=MSS, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'MSS'
  size = nrow(MSS)
  
})


## bind all results together 
Smsiaucos<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm)
Smsiaucos<-Smsiaucos[order(Smsiaucos$model, Smsiaucos$Group), ]

write_xlsx(Smsiaucos, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Smsiaucos.xlsx", col_names = T)

## DSS
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
stradata<-stradata[!is.na(stradata$recurr_cp), ]
# MSI
MSI<-subset(stradata, MSI_gent=='MSI') ### 

covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = MSI, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=MSI, metrics = c("auc", "brier"))

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
MSS<-subset(stradata, MSI_gent=='MSS') ### 
## PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = MSS, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=MSS, metrics = c("auc", "brier"))

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


## for PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")
stradatapm<-stradatapm[!is.na(stradatapm$recurr_cp), ]
# MSI
MSI<-subset(stradatapm, MSI_gent=='MSI')### 
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = MSI, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=MSI, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'MSI'
  size = nrow(MSI)
  
})



## MSS
MSS<-subset(stradatapm, MSI_gent=='MSS') #
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = MSS, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=MSS, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'MSS'
  size = nrow(MSS)
  
})

## bind all results together 
SmsiaucDSS<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm)
SmsiaucDSS<-SmsiaucDSS[order(SmsiaucDSS$model, SmsiaucDSS$Group), ]


write_xlsx(SmsiaucDSS, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/SmsiaucDFS.xlsx", col_names = T)


### bind all results together ###
### bind all results together ###
MSIaucos <- read_excel("Validation_tablefigure/Smsiaucos.xlsx")
MSIaucos<-MSIaucos[, c(1,3,2, 5)]
colnames(MSIaucos)[3:4]<-c('size OS', 'AUC (95%CI) OS')

MSIaucdfs<- read_excel("Validation_tablefigure/SmsiaucDFS.xlsx")
dput(names(MSIaucdfs))
MSIaucdfs<-subset(MSIaucdfs, select = c("size", "AUCCI"))
colnames(MSIaucdfs)<-c('size DFS', 'AUC (95%CI) DFS')

MSI<-bind_cols(MSIaucos, MSIaucdfs)


write_xlsx(MSI, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/strMSI.xlsx", col_names = T)

############# 7. BRAF mutation ##########
### OS ##
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
summary(stradata$brafmut)
## mut
mut<-subset(stradata, brafmut=='BRAF mut') ### 
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

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
wt<-subset(stradata, brafmut=='BRAF wt')### 1769
####### OS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

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


## for PM classifier 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")

mut<-subset(stradatapm, brafmut=='BRAF mut') # 153

####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mut, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=mut, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'mut'
  size = nrow(mut)
  
})

## wt
wt<-subset(stradatapm, brafmut=='BRAF wt') 
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = wt, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=wt, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'wt'
  size = nrow(wt)
  
})


## bind all results together 
Smutaucos<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm)
Smutaucos<-Smutaucos[order(Smutaucos$model, Smutaucos$Group), ]

write_xlsx(Smutaucos, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Smutaucos.xlsx", col_names = T)

## DFS
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
stradata<-stradata[!is.na(stradata$recurr_cp), ]
## mut
mut<-subset(stradata, brafmut=='BRAF mut') ### 
## PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mut, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=mut, metrics = c("auc", "brier"))

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
wt<-subset(stradata, brafmut=='BRAF wt') ### 
## PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = wt, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=wt, metrics = c("auc", "brier"))

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

## for PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")
stradatapm<-stradatapm[!is.na(stradatapm$recurr_cp), ]
## mut
mut<-subset(stradatapm, brafmut=='BRAF mut')### 
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mut, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=mut, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'mut'
  size = nrow(mut)
  
})

## wt
wt<-subset(stradatapm, brafmut=='BRAF wt') #
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = wt, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=wt, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'wt'
  size = nrow(wt)
  
})


## bind all results together 
SmutaucDSS<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm)
SmutaucDSS<-SmutaucDSS[order(SmutaucDSS$model, SmutaucDSS$Group), ]

write_xlsx(SmutaucDSS, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/SmutaucDFS.xlsx", col_names = T)

### bind all results together ###

brafaucos <- read_excel("Validation_tablefigure/Smutaucos.xlsx")
brafaucos<-brafaucos[, c(1,3,2, 5)]
colnames(brafaucos)[3:4]<-c('size OS', 'AUC (95%CI) OS')

brafaucdfs<- read_excel("Validation_tablefigure/SmutaucDFS.xlsx")
dput(names(brafaucdfs))
brafaucdfs<-subset(brafaucdfs, select = c("size", "AUCCI"))
colnames(brafaucdfs)<-c('size DFS', 'AUC (95%CI) DFS')

braf<-bind_cols(brafaucos, brafaucdfs)


write_xlsx(braf, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/strbraf.xlsx", col_names = T)

########## 8. KRAS mutation ##########
### OS ##
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
summary(stradata$krasmut)
## mut
mut<-subset(stradata, krasmut=='KRAS mut') ### 
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

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
wt<-subset(stradata, krasmut=='KRAS wt')### 1769
####### OS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

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


## for PM classifier 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")

mut<-subset(stradatapm, krasmut=='KRAS mut') # 153

####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mut, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=mut, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'mut'
  size = nrow(mut)
  
})

## wt
wt<-subset(stradatapm, krasmut=='KRAS wt') 
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = wt, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=wt, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'wt'
  size = nrow(wt)
  
})


## bind all results together 
Smutaucos<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm)
Smutaucos<-Smutaucos[order(Smutaucos$model, Smutaucos$Group), ]

write_xlsx(Smutaucos, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Smutaucos.xlsx", col_names = T)

## DFS
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
stradata<-stradata[!is.na(stradata$recurr_cp), ]
## mut
mut<-subset(stradata, krasmut=='KRAS mut') ### 
## PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mut, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=mut, metrics = c("auc", "brier"))

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
wt<-subset(stradata, krasmut=='KRAS wt') ### 
## PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = wt, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=wt, metrics = c("auc", "brier"))

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


## for PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")
stradatapm<-stradatapm[!is.na(stradatapm$recurr_cp), ]
## mut
mut<-subset(stradatapm, krasmut=='KRAS mut')### 
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = mut, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=mut, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'mut'
  size = nrow(mut)
  
})

## wt
wt<-subset(stradatapm, krasmut=='KRAS wt') #
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = wt, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=wt, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'wt'
  size = nrow(wt)
  
})


## bind all results together 
SmutaucDSS<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm)
SmutaucDSS<-SmutaucDSS[order(SmutaucDSS$model, SmutaucDSS$Group), ]


write_xlsx(SmutaucDSS, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/SmutaucDFS.xlsx", col_names = T)

### bind all results together ###
krasaucos <- read_excel("Validation_tablefigure/Smutaucos.xlsx")
krasaucos<-krasaucos[, c(1,3,2, 5)]
colnames(krasaucos)[3:4]<-c('size OS', 'AUC (95%CI) OS')

krasaucdfs<- read_excel("Validation_tablefigure/SmutaucDFS.xlsx")
dput(names(krasaucdfs))
krasaucdfs<-subset(krasaucdfs, select = c("size", "AUCCI"))
colnames(krasaucdfs)<-c('size DFS', 'AUC (95%CI) DFS')

kras<-bind_cols(krasaucos, krasaucdfs)


write_xlsx(kras, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/strkras.xlsx", col_names = T)


write_xlsx(Strkrasmut, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Strkrasmut.xlsx", col_names = T)

######### 9. CIMP mutation ##################
 OS
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
summary(stradata$cimphi)
 high 
high<-subset(stradata, cimphi=='CIMP-high') ### 
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

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
low<-subset(stradata, cimphi=='CIMP-low')### 1769
####### OS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

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

#################### for PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")

high<-subset(stradatapm, cimphi=='CIMP-high') # 153

####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = high, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=high, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'high'
  size = nrow(high)
  
})


######## low 
low<-subset(stradatapm, cimphi=='CIMP-low') 
####### OS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = low, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=low, metrics = c("auc", "brier"))

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
  Outcome<-'OS'
  Group<-'low'
  size = nrow(low)
  
})

#### bind all results together 
Shighaucos<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm)
Shighaucos<-Shighaucos[order(Shighaucos$model, Shighaucos$Group), ]

write_xlsx(Shighaucos, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Shighaucos.xlsx", col_names = T)

## DFS
#### exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydata.RData")
stradata<-stradata[!is.na(stradata$recurr_cp), ]
 high 
high<-subset(stradata, cimphi=='CIMP-high') ### 

covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = high, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=high, metrics = c("auc", "brier"))

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
low<-subset(stradata, cimphi=='CIMP-low') ### 
####### PFS
covariates<-c("Gong_2020_OS", "Huang_2021_PFS", "Huang_2021_OS",  "WangY_2020_OS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = low, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=low, metrics = c("auc", "brier"))

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

#################### for PM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/stratifydatapm.RData")
stradatapm<-stradatapm[!is.na(stradatapm$recurr_cp), ]
 high
high<-subset(stradatapm, cimphi=='CIMP-high')### 
####### PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = high, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=high, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'high'
  size = nrow(high)
  
})

######## low 
low<-subset(stradatapm, cimphi=='CIMP-low') #
## PFS
covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = low, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, recurr_cp==1)~1,data=low, metrics = c("auc", "brier"))

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
  Outcome<-'DSS'
  Group<-'low'
  size = nrow(low)
  
})


#### bind all results together 
ShighaucDSS<-bind_rows(AUC_allo, AUC_allopm, AUC_ally, AUC_allypm)
ShighaucDSS<-ShighaucDSS[order(ShighaucDSS$model, ShighaucDSS$Group), ]

write_xlsx(ShighaucDSS, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/ShighaucDFS.xlsx", col_names = T)

### bind all results together ###
Shighaucos <- read_excel("Validation_tablefigure/Shighaucos.xlsx")
Shighbrieros <- read_excel("Validation_tablefigure/Shighbrieros.xlsx")

OS<-bind_cols(Shighaucos, Shighbrieros[,5])

colnames(OS)[5:6]<-c("AUC (95%CI) OS", "Brier Score (95%CI) OS")

ShighaucPFS <- read_excel("Validation_tablefigure/ShighaucPFS.xlsx")
ShighbrierPFS <- read_excel("Validation_tablefigure/ShighbrierPFS.xlsx")

PFS<-bind_cols(ShighaucPFS, ShighbrierPFS[,5])

colnames(PFS)[5:6]<-c("AUC (95%CI) PFS", "Brier Score (95%CI) PFS")

Strhigh<-bind_cols(OS, PFS[, c(4, 2, 5, 6)])

ShighaucDSS <- read_excel("Validation_tablefigure/ShighaucDSS.xlsx")
ShighbrierDSS <- read_excel("Validation_tablefigure/ShighbrierDSS.xlsx")

DSS<-bind_cols(ShighaucDSS, ShighbrierDSS[,5])

colnames(DSS)[5:6]<-c("AUC (95%CI) DSS", "Brier Score (95%CI) DSS")

Strhigh<-bind_cols(Strhigh, DSS[, c(4, 2, 5, 6)])

Strcimp<-Strhigh

write_xlsx(Strcimp, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Strcimp.xlsx", col_names = T)










