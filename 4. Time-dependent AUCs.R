library(riskRegression)
library(survival)
library(dplyr)
library(plyr)
library(writexl)
library(tableone)
library(rlist)
library(cgwtools)
############################ 1. time-dependent AUCs #############################
######################### 1.1 OS ##############################################
# Colon, stage1-4 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c14.RData")
# continous 
covariates_con<-c("WangY_2020_C14","Yang_2019_C14","Xiang_2020_C14")
# distribution
summary(unicox_c14)


univ_formulas <- sapply(covariates_con,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_c14, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(timey, death_all)~1,data=unicox_c14, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_td<-td_aucbrier_alldeath$AUC$score
AUC_td<-AUC_td %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_td<-within.data.frame(AUC_td, {
  AUC_OS_C<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})

AUC_td$AUC<-NULL
AUC_td$lower<-NULL
AUC_td$upper<-NULL

# quartile
covariates_qu<-c("TWangY_2020_C14","TYang_2019_C14","TXiang_2020_C14")

univ_formulas <- sapply(covariates_qu,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_c14, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(timey, death_all)~1,data=unicox_c14, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_tdq<-td_aucbrier_alldeath$AUC$score
AUC_tdq<-AUC_tdq %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_tdq<-within.data.frame(AUC_tdq, {
  AUC_OS_Q<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})


AUC_td$AUC_OS_Q<-AUC_tdq$AUC_OS_Q

AUC_c14<-AUC_td

save(AUC_c14, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_OS.RData")

# Colon and rectum, stage1-2 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr12.RData")
dput(names(unicox_cr12))
# continous 
covariates_con<-c("Yu_2023_DFS_CR12")

univ_formulas <- sapply(covariates_con,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr12, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(timey, death_all)~1,data=unicox_cr12, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_td<-td_aucbrier_alldeath$AUC$score
AUC_td<-AUC_td %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_td<-within.data.frame(AUC_td, {
  AUC_OS_C<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})

AUC_td$AUC<-NULL
AUC_td$lower<-NULL
AUC_td$upper<-NULL

# quartile
covariates_qu<-c("TYu_2023_DFS_CR12")

univ_formulas <- sapply(covariates_qu,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr12, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(timey, death_all)~1,data=unicox_cr12, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_tdq<-td_aucbrier_alldeath$AUC$score
AUC_tdq<-AUC_tdq %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_tdq<-within.data.frame(AUC_tdq, {
  AUC_OS_Q<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})


AUC_td$AUC_OS_Q<-AUC_tdq$AUC_OS_Q

AUC_cr12<-AUC_td

resave(AUC_cr12, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_OS.RData")

# Colon and rectum, stage1-3

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr13.RData")
dput(names(unicox_cr13))
summary(unicox_cr13$Gundert_2019_CR13)

#continous 
covariates_con<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates_con,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr13, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(timey, death_all)~1,data=unicox_cr13, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_td<-td_aucbrier_alldeath$AUC$score
AUC_td<-AUC_td %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_td<-within.data.frame(AUC_td, {
  AUC_OS_C<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})

AUC_td$AUC<-NULL
AUC_td$lower<-NULL
AUC_td$upper<-NULL

# quartile
covariates_qu<-c("TGundert_2019_CR13")

univ_formulas <- sapply(covariates_qu,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr13, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(timey, death_all)~1,data=unicox_cr13, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_tdq<-td_aucbrier_alldeath$AUC$score
AUC_tdq<-AUC_tdq %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_tdq<-within.data.frame(AUC_tdq, {
  AUC_OS_Q<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})


AUC_td$AUC_OS_Q<-AUC_tdq$AUC_OS_Q

AUC_cr13<-AUC_td

resave(AUC_cr13, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_OS.RData")

# Colon and rectum, stage2-3 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr23.RData")
dput(names(unicox_cr23))

#continous 
covariates_con<-c("ChenF_2021_CR23")

univ_formulas <- sapply(covariates_con,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr23, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(timey, death_all)~1,data=unicox_cr23, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_td<-td_aucbrier_alldeath$AUC$score
AUC_td<-AUC_td %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_td<-within.data.frame(AUC_td, {
  AUC_OS_C<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})

AUC_td$AUC<-NULL
AUC_td$lower<-NULL
AUC_td$upper<-NULL

# quartile
covariates_qu<-c("TChenF_2021_CR23")

univ_formulas <- sapply(covariates_qu,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr23, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(timey, death_all)~1,data=unicox_cr23, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_tdq<-td_aucbrier_alldeath$AUC$score
AUC_tdq<-AUC_tdq %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_tdq<-within.data.frame(AUC_tdq, {
  AUC_OS_Q<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})


AUC_td$AUC_OS_Q<-AUC_tdq$AUC_OS_Q

AUC_cr23<-AUC_td

resave(AUC_cr23, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_OS.RData")

# Colon and rectum, stage1-4 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")
dput(names(unicox_cr14))

#continous 
covariates_con<-c("Gong_2020_CR14", "WangX_2020_CR14", "Huang_2021_OS_CR14", "Huang_2021_DFS_CR14", "Li_2021_CR14")

univ_formulas <- sapply(covariates_con,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr14, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(timey, death_all)~1,data=unicox_cr14, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_td<-td_aucbrier_alldeath$AUC$score
AUC_td<-AUC_td %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_td<-within.data.frame(AUC_td, {
  AUC_OS_C<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})

AUC_td$AUC<-NULL
AUC_td$lower<-NULL
AUC_td$upper<-NULL

# quartile
covariates_qu<-c("TGong_2020_CR14", "TWangX_2020_CR14", "THuang_2021_OS_CR14", "THuang_2021_DFS_CR14", "TLi_2021_CR14")

univ_formulas <- sapply(covariates_qu,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr14, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(timey, death_all)~1,data=unicox_cr14, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_tdq<-td_aucbrier_alldeath$AUC$score
AUC_tdq<-AUC_tdq %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_tdq<-within.data.frame(AUC_tdq, {
  AUC_OS_Q<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})


AUC_td$AUC_OS_Q<-AUC_tdq$AUC_OS_Q

AUC_cr14<-AUC_td

resave(AUC_cr14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_OS.RData")

# Colon and rectum, stage1-4, 2 studies from the DACHS 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_jiaso.RData")
dput(names(unicox_jiaso))

#continous 
covariates_con<-c("Jia_2019_CR14")

univ_formulas <- sapply(covariates_con,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_jiaso, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(timey, death_all)~1,data=unicox_jiaso, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_td<-td_aucbrier_alldeath$AUC$score
AUC_td<-AUC_td %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_td<-within.data.frame(AUC_td, {
  AUC_OS_C<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})

AUC_td$AUC<-NULL
AUC_td$lower<-NULL
AUC_td$upper<-NULL

# quartile
covariates_qu<-c("TJia_2019_CR14")

univ_formulas <- sapply(covariates_qu,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_jiaso, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(timey, death_all)~1,data=unicox_jiaso, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_tdq<-td_aucbrier_alldeath$AUC$score
AUC_tdq<-AUC_tdq %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_tdq<-within.data.frame(AUC_tdq, {
  AUC_OS_Q<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})


AUC_td$AUC_OS_Q<-AUC_tdq$AUC_OS_Q

AUC_jiaso<-AUC_td

resave(AUC_jiaso, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_OS.RData")

# bind all results together
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_OS.RData")
AUC_OS<-rbind(AUC_c14, AUC_cr12, AUC_cr13, AUC_cr14, AUC_cr23, AUC_jiaso)
save(AUC_OS, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_OS_comb.RData")

##################################### 1.2 DFS #############################################
# Colon Stage 1-4 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c14.RData")
unicox_c14<-unicox_c14[!is.na(unicox_c14$recurr_cp), ]
dput(names(unicox_c14))
# continous 
covariates_con<-c("WangY_2020_C14","Yang_2019_C14","Xiang_2020_C14")

univ_formulas <- sapply(covariates_con,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_c14, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_c14, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_td<-td_aucbrier_alldeath$AUC$score
AUC_td<-AUC_td %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_td<-within.data.frame(AUC_td, {
  AUC_DFS_C<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})

AUC_td$AUC<-NULL
AUC_td$lower<-NULL
AUC_td$upper<-NULL

# quartile
covariates_qu<-c("TWangY_2020_C14","TYang_2019_C14","TXiang_2020_C14")

univ_formulas <- sapply(covariates_qu,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_c14, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_c14, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_tdq<-td_aucbrier_alldeath$AUC$score
AUC_tdq<-AUC_tdq %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_tdq<-within.data.frame(AUC_tdq, {
  AUC_DFS_Q<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})


AUC_td$AUC_DFS_Q<-AUC_tdq$AUC_DFS_Q

AUC_c14<-AUC_td

save(AUC_c14, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_DFS.RData")

# Colon and rectum Stage 1-2

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr12.RData")
dput(names(unicox_cr12))
unicox_cr12<-unicox_cr12[!is.na(unicox_cr12$recurr_cp), ]
# continous 
covariates_con<-c("Yu_2023_DFS_CR12")

univ_formulas <- sapply(covariates_con,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr12, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_cr12, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_td<-td_aucbrier_alldeath$AUC$score
AUC_td<-AUC_td %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_td<-within.data.frame(AUC_td, {
  AUC_DFS_C<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})

AUC_td$AUC<-NULL
AUC_td$lower<-NULL
AUC_td$upper<-NULL

# quartile
covariates_qu<-c("TYu_2023_DFS_CR12")

univ_formulas <- sapply(covariates_qu,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr12, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_cr12, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_tdq<-td_aucbrier_alldeath$AUC$score
AUC_tdq<-AUC_tdq %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_tdq<-within.data.frame(AUC_tdq, {
  AUC_DFS_Q<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})


AUC_td$AUC_DFS_Q<-AUC_tdq$AUC_DFS_Q

AUC_cr12<-AUC_td

resave(AUC_cr12, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_DFS.RData")

# Colon and rectum Stage 1-3

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr13.RData")
dput(names(unicox_cr13))
unicox_cr13<-unicox_cr13[!is.na(unicox_cr13$recurr_cp), ]
#continous 
covariates_con<-c("Gundert_2019_CR13")

univ_formulas <- sapply(covariates_con,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr13, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_cr13, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_td<-td_aucbrier_alldeath$AUC$score
AUC_td<-AUC_td %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_td<-within.data.frame(AUC_td, {
  AUC_DFS_C<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})

AUC_td$AUC<-NULL
AUC_td$lower<-NULL
AUC_td$upper<-NULL

# quartile
covariates_qu<-c("TGundert_2019_CR13")

univ_formulas <- sapply(covariates_qu,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr13, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_cr13, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_tdq<-td_aucbrier_alldeath$AUC$score
AUC_tdq<-AUC_tdq %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_tdq<-within.data.frame(AUC_tdq, {
  AUC_DFS_Q<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})


AUC_td$AUC_DFS_Q<-AUC_tdq$AUC_DFS_Q

AUC_cr13<-AUC_td

resave(AUC_cr13, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_DFS.RData")

# Colon and rectum Stage 2-3

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr23.RData")
dput(names(unicox_cr23))
unicox_cr23<-unicox_cr23[!is.na(unicox_cr23$recurr_cp), ]

#continous 
covariates_con<-c("ChenF_2021_CR23")

univ_formulas <- sapply(covariates_con,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr23, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_cr23, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_td<-td_aucbrier_alldeath$AUC$score
AUC_td<-AUC_td %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_td<-within.data.frame(AUC_td, {
  AUC_DFS_C<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})

AUC_td$AUC<-NULL
AUC_td$lower<-NULL
AUC_td$upper<-NULL

# quartile
covariates_qu<-c("TChenF_2021_CR23")

univ_formulas <- sapply(covariates_qu,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr23, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_cr23, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_tdq<-td_aucbrier_alldeath$AUC$score
AUC_tdq<-AUC_tdq %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_tdq<-within.data.frame(AUC_tdq, {
  AUC_DFS_Q<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})


AUC_td$AUC_DFS_Q<-AUC_tdq$AUC_DFS_Q

AUC_cr23<-AUC_td

resave(AUC_cr23, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_DFS.RData")

# Colon and rectum Stage 1-4

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")
dput(names(unicox_cr14))
unicox_cr14<-unicox_cr14[!is.na(unicox_cr14$recurr_cp), ]

#continous 
covariates_con<-c("Gong_2020_CR14", "WangX_2020_CR14", "Huang_2021_DFS_CR14", "Huang_2021_OS_CR14", "Li_2021_CR14")

univ_formulas <- sapply(covariates_con,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr14, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_cr14, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_td<-td_aucbrier_alldeath$AUC$score
AUC_td<-AUC_td %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_td<-within.data.frame(AUC_td, {
  AUC_DFS_C<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})

AUC_td$AUC<-NULL
AUC_td$lower<-NULL
AUC_td$upper<-NULL

# quartile
covariates_qu<-c("TGong_2020_CR14", "TWangX_2020_CR14", "THuang_2021_DFS_CR14", "THuang_2021_OS_CR14", "TLi_2021_CR14")

univ_formulas <- sapply(covariates_qu,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr14, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_cr14, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_tdq<-td_aucbrier_alldeath$AUC$score
AUC_tdq<-AUC_tdq %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_tdq<-within.data.frame(AUC_tdq, {
  AUC_DFS_Q<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})


AUC_td$AUC_DFS_Q<-AUC_tdq$AUC_DFS_Q

AUC_cr14<-AUC_td

resave(AUC_cr14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_DFS.RData")

# Colon and rectum Stage 1-4, two studies from the DACHS
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_jiaso.RData")
dput(names(unicox_jiaso))
unicox_jiaso<-unicox_jiaso[!is.na(unicox_jiaso$recurr_cp), ]


#continous 
covariates_con<-c("Jia_2019_CR14")

univ_formulas <- sapply(covariates_con,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_jiaso, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_jiaso, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_td<-td_aucbrier_alldeath$AUC$score
AUC_td<-AUC_td %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_td<-within.data.frame(AUC_td, {
  AUC_DFS_C<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})

AUC_td$AUC<-NULL
AUC_td$lower<-NULL
AUC_td$upper<-NULL

# quartile
covariates_qu<-c("TJia_2019_CR14")

univ_formulas <- sapply(covariates_qu,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_jiaso, y=TRUE, x = TRUE)})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_jiaso, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_tdq<-td_aucbrier_alldeath$AUC$score
AUC_tdq<-AUC_tdq %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_tdq<-within.data.frame(AUC_tdq, {
  AUC_DFS_Q<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  times<-paste0(times,' ', 'year')
  
})

AUC_td$AUC_DFS_Q<-AUC_tdq$AUC_DFS_Q

AUC_jiaso<-AUC_td

resave(AUC_jiaso, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_DFS.RData")

# bind all results together
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_DFS.RData")
AUC_DFS<-rbind(AUC_c14, AUC_cr12, AUC_cr13, AUC_cr14, AUC_cr23, AUC_jiaso)
save(AUC_DFS, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_DFS_comb.RData")

# combine all results together 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_DFS_comb.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUC_OS_comb.RData")

AUC<-merge(AUC_OS, AUC_DFS, by = c('model', 'times'))

write.csv(AUC, 
          file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/tablesfigures/AUC.csv")

## distribution
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_jiaso.RData")
summary(unicox_jiaso$Jia_2019_CR14)

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr12.RData")
summary(unicox_cr12$Yu_2023_DFS_CR12)



















