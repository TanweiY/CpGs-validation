####################### 1. unicox  #####################################################
## 1.1 OS ####
# exclude ProM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/unicox_all.RData")
covariates<-c("TYang_2019_OS", 
              "TGong_2020_OS", "TWangX_2020_OS_PFS", "TWangY_2020_OS", "TXiang_2020_OS", 
              "TChen_2021_OS", "THuang_2021_OS", "THuang_2021_PFS", "TLi_2021_OS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timeD, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

#bind results
OS<-ldply(univ_results, rbind)
OS$.id<-NULL
colnames(OS)<-c('Score', 'HR_OS', 'p_OS')

# add ProM classifier

covariates<-c("TGundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timeD, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicoxpm)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

#将一个list中所有dataframe合并

OSpm<-ldply(univ_results, rbind)
OSpm$.id<-NULL
colnames(OSpm)<-c('Score', 'HR_OS', 'p_OS')
unicoxscore_os<-bind_rows(OS, OSpm)

## 1.2 DFS ####
# exclude ProM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/unicox_all.RData")
unicox_crc<-unicox[!is.na(unicox$recurr_cp), ]
covariates<-c("TYang_2019_OS",
              "TGong_2020_OS", "TWangX_2020_OS_PFS", "TWangY_2020_OS", "TXiang_2020_OS",
              "TChen_2021_OS", "THuang_2021_OS", "THuang_2021_PFS", "TLi_2021_OS")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_crc)})
uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)
DFS<-ldply(univ_results, rbind)
DFS$.id<-NULL
colnames(DFS)<-c('Score', 'HR_DFS', 'p_DFS')

### validation ProM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/PMunicox.RData")
unicoxpm_crc<-unicoxpm[!is.na(unicoxpm$recurr_cp), ]
covariates<-c("TGundert_2019_OS_DSS")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicoxpm_crc)})
uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)
DFSpm<-ldply(univ_results, rbind)
DFSpm$.id<-NULL
colnames(DFSpm)<-c('Score', 'HR_DFS', 'p_DFS')
unicoxscore_dfs<-bind_rows(DFS, DFSpm)

######################################## 2. multicox_main adjustment  ##################################################################################
## 2.1 OS ####
## exclude  PM classifier ##
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicox_impu20list.RData")
dput(names(multi_impu[[1]]))

scoret<-c("TYang_2019_OS", 
          "TGong_2020_OS", "TWangX_2020_OS_PFS", "TWangY_2020_OS", "TXiang_2020_OS", 
          "TChen_2021_OS", "THuang_2021_OS", "THuang_2021_PFS", "TLi_2021_OS")


multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(scoret)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  Score<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(Score, HR, LL, UL, p)
  
}

OS<-ldply(result, rbind)
OS$Score<-str_sub(OS$Score, start = 2, end = 30)
multicoxscore_os<-bind_rows(OS, OSpm)

## 2.2 DFS ####
# exclude CpGs in PM classifier ##
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicox_impu20list.RData")
dput(names(multi_impu[[1]]))
scoret<-c("TYang_2019_OS",
          "TGong_2020_OS", "TWangX_2020_OS_PFS", "TWangY_2020_OS", "TXiang_2020_OS",
          "TChen_2021_OS", "THuang_2021_OS", "THuang_2021_PFS", "TLi_2021_OS")
multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))
n<-length(scoret)
result<-vector(n, mode="list")
for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu[[j]])
  }
  p<-summary(pool(mice.fit))
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  Score<-rownames(HRCI)[ c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  result[[i]]<-data.frame(Score, HR, LL, UL, p)
}
DFS<-ldply(result, rbind)
DFS$Score<-str_sub(DFS$Score, start = 2, end = 30)
# keep twp decimals
DFS[,c(2:4)]<-round(DFS[,c(2:4)], digits = 2)
DFS$p_round<-round(DFS$p, digits = 3)
DFS$HRCI<-paste0(DFS$HR, ' ', '(', DFS$LL, ',', ' ', DFS$UL, ')' )

# validate PM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicoxPM_impu20list.RData")
dput(names(multipm_impu[[1]]))
scoret<-c("TGundert_2019_OS_DSS")
multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))
n<-length(scoret)
result<-vector(n, mode="list")
for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multipm_impu[[j]])
  }
  p<-summary(pool(mice.fit))
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  Score<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  result[[i]]<-data.frame(Score, HR, LL, UL, p)
}

DFSpmt<-ldply(result, rbind)
DFSpmt$Score<-str_sub(DFSpmt$Score, start = 2, end = 30)
DFSpm<-DFSpmt
# keep twp decimals
DFSpm[,c(2:4)]<-round(DFSpm[,c(2:4)], digits = 2)
DFSpm$p_round<-round(DFSpm$p, digits = 3)
DFSpm$HRCI<-paste0(DFSpm$HR, ' ', '(', DFSpm$LL, ',', ' ', DFSpm$UL, ')' )

multicoxscore_dfs<-bind_rows(DFS, DFSpm)

######################################## 3. multicox supplement adjustment add MSI  ##################################################################################
## 3.1 OS ####
# exclude CpGs in PM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicox_impu20list.RData")
dput(names(multi_impu[[1]]))

scoret<-c("TYang_2019_OS", 
          "TGong_2020_OS", "TWangX_2020_OS_PFS", "TWangY_2020_OS", "TXiang_2020_OS", 
          "TChen_2021_OS", "THuang_2021_OS", "THuang_2021_PFS", "TLi_2021_OS")


multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+MSI_gent+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(scoret)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  Score<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(Score, HR, LL, UL, p)
  
}

OS<-ldply(result, rbind)
OS$Score<-str_sub(OS$Score, start = 2, end = 30)

### validate  ProM classifier

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicoxPM_impu20list.RData")

dput(names(multipm_impu[[1]]))

scoret<-c("TGundert_2019_OS_DSS")


multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+MSI_gent+chemradther+Stage_at_diagnosis+Location+', x)))


n<-length(scoret)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multipm_impu[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  Score<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(Score, HR, LL, UL, p)
  
}

OSpm<-ldply(result, rbind)
OSpm$Score<-str_sub(OSpm$Score, start = 2, end = 30)

multicox2score_os<-bind_rows(OS, OSpm)

## 3.2 DFS ####
# exclude PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicox_impu20list.RData")
scoret<-c("TYang_2019_OS",
          "TGong_2020_OS", "TWangX_2020_OS_PFS", "TWangY_2020_OS", "TXiang_2020_OS",
          "TChen_2021_OS", "THuang_2021_OS", "THuang_2021_PFS", "TLi_2021_OS")
multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+MSI_gent+chemradther+Stage_at_diagnosis+Location+', x)))
n<-length(scoret)
result<-vector(n, mode="list")
for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu[[j]])
  }
  p<-summary(pool(mice.fit))
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  Score<-rownames(HRCI)[ c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  result[[i]]<-data.frame(Score, HR, LL, UL, p)
}
DFS<-ldply(result, rbind)
DFS$Score<-str_sub(DFS$Score, start = 2, end = 30)

# validate PM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicoxPM_impu20list.RData")
dput(names(multipm_impu[[1]]))
scoret<-c("TGundert_2019_OS_DSS")
multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+MSI_gent+chemradther+Stage_at_diagnosis+Location+', x)))
n<-length(scoret)
result<-vector(n, mode="list")
for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multipm_impu[[j]])
  }
  p<-summary(pool(mice.fit))
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  Score<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  result[[i]]<-data.frame(Score, HR, LL, UL, p)
}
DFSpm<-ldply(result, rbind)
DFSpm$Score<-str_sub(DFSpm$Score, start = 2, end = 30)

multicox2score_os<-bind_rows(DFS, DFSpm)

######################################## 4. overall and time-dependent AUC  ##################################################################################
## 4.1 OS ####
# exclude  PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/unicox_all.RData")
covariates<-c( "Li_2021_OS", "Huang_2021_PFS", "Huang_2021_OS", "Chen_2021_OS", 
               "Xiang_2020_OS", "WangY_2020_OS", "WangX_2020_OS_PFS", "Gong_2020_OS", 
               "Yang_2019_OS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))



univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox, y=TRUE, x = TRUE)})


# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=unicox, metrics = c("auc"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_all<-within.data.frame(AUC_all, {
  times ='Overall'
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  Outcome<-'OS'
  
})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(timey, death_all)~1,data=unicox, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_td<-td_aucbrier_alldeath$AUC$score
AUC_td<-AUC_td %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_td<-within.data.frame(AUC_td, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  Outcome<-'OS'
  times<-paste0(times,' ', 'year')
  
})

AUCOS<-bind_rows(AUC_all, AUC_td)
AUCOS<-AUCOS[order(AUCOS$model),]

# same procedure for PM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/PMunicox.RData")

covariates<-c("Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(timey, death_all)~', x)))


univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicoxpm, y=TRUE, x = TRUE)})


# overall AUC
overall_aucbrier_alldeath<-Score(univ_models,
                                 formula=Surv(timey, death_all)~1,data=unicoxpm, metrics = c("auc"))

AUC_all<-overall_aucbrier_alldeath$AUC$score
AUC_all<-AUC_all %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_all<-within.data.frame(AUC_all, {
  times ='Overall'
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  Outcome<-'OS'
  
})

# time-dependent AUC

td_aucbrier_alldeath<-Score(univ_models,
                            formula=Surv(timey, death_all)~1,data=unicoxpm, metrics = c("auc", "brier"), 
                            times=c(1, 3, 5, 8))


AUC_td<-td_aucbrier_alldeath$AUC$score
AUC_td<-AUC_td %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_td<-within.data.frame(AUC_td, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  Outcome<-'OS'
  times<-paste0(times,' ', 'year')
  
})

AUCOSpm<-bind_rows(AUC_all, AUC_td)
AUCOS<-bind_rows(AUCOSpm, AUCOS)
colnames(AUCOS)[7]<-'AUC (95%CI)'

write_xlsx(AUCOS, 
           path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/AUC_OS.xlsx", col_names = T)

## 4.2 DFS ####
# exclude PM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/unicox_all.RData")

dput(names(unicox))
unicox_DFS<-unicox[!is.na(unicox$recurr_cp),]

covariates<-c( "Li_2021_OS", "Huang_2021_PFS", "Huang_2021_OS", "Chen_2021_OS", 
               "Xiang_2020_OS", "WangY_2020_OS", "WangX_2020_OS_PFS", "Gong_2020_OS", 
               "Yang_2019_OS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_DFS, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_DFS<-Score(univ_models,
                            formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_DFS, metrics = c("auc"))

AUC_DFS<-overall_aucbrier_DFS$AUC$score
AUC_DFS<-AUC_DFS %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_DFS<-within.data.frame(AUC_DFS, {
  times ='Overall'
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  Outcome<-'DFS'
  
})

# time-dependent AUC

td_aucbrier_DFS<-Score(univ_models,
                       formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_DFS, metrics = c("auc"), 
                       times=c(1, 3, 5, 8))


AUC_tdDFS<-td_aucbrier_DFS$AUC$score
AUC_tdDFS<-AUC_tdDFS %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_tdDFS<-within.data.frame(AUC_tdDFS, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  Outcome<-'DFS'
  times<-paste0(times,' ', 'year')
  
})

AUCDFS<-bind_rows(AUC_DFS, AUC_tdDFS)
AUCDFS<-AUCDFS[order(AUCDFS$model),]

# same procedure for PM classifer 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/PMunicox.RData")

unicox_DFS<-unicoxpm[!is.na(unicoxpm$recurr_cp),]
summary(unicox_DFS$recurr_cp==1)
covariates<-c( "Gundert_2019_OS_DSS")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_DFS, y=TRUE, x = TRUE)})

# overall AUC
overall_aucbrier_DFS<-Score(univ_models,
                            formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_DFS, metrics = c("auc"))

AUC_DFS<-overall_aucbrier_DFS$AUC$score
AUC_DFS<-AUC_DFS %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_DFS<-within.data.frame(AUC_DFS, {
  times ='Overall'
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  Outcome<-'DFS'
  
})

# time-dependent AUC
td_aucbrier_DFS<-Score(univ_models,
                       formula=Surv(recurr_timey, recurr_cp==1)~1,data=unicox_DFS, metrics = c("auc"), 
                       times=c(1, 3, 5, 8))


AUC_tdDFS<-td_aucbrier_DFS$AUC$score
AUC_tdDFS<-AUC_tdDFS %>% 
  mutate_if(is.numeric, round, digits=2)
AUC_tdDFS<-within.data.frame(AUC_tdDFS, {
  AUCCI<-paste0(AUC, ' ', '(', lower, ',', ' ', upper, ')')
  se<-NULL
  Outcome<-'DFS'
  times<-paste0(times,' ', 'year')
  
})

AUCDFSpm<-bind_rows(AUC_DFS, AUC_tdDFS)
AUCDFS<-bind_rows(AUCDFSpm, AUCDFS)























