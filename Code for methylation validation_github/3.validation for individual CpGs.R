####################### 1. Unicox  #####################################################
## 1.1 OS ####

# CpGs excluding CpGs in ProM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/unicox_all.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/cpgcanvalidate.RData")

cpgd <- read_excel("CpG_Data extraction.xlsx",
                   sheet = "all CpGs", col_types = c("text",
                                                     "skip", "text", "skip", "skip", "skip",
                                                     "skip", "skip", "skip", "skip"))

cpgd<-cpgd$CpGs[cpgd$Study=='Gündert et al, 2019']


cpg<-cpg[!(cpg %in% cpgd)]
rm(cpgd)

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

# bind results

OS<-ldply(univ_results, rbind)
OS$rownames.b.<-NULL
colnames(OS)<-c('cpg', 'HR_OS', 'p_OS')

# validate CpGs for ProM classifier 
cpg <- read_excel("CpG_Data extraction.xlsx",
                  sheet = "all CpGs", col_types = c("text",
                                                    "skip", "text", "skip", "skip", "skip",
                                                    "skip", "skip", "skip", "skip"))

cpg<-cpg$CpGs[cpg$Study=='Gündert et al, 2019']

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/PMunicox.RData")

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicoxpm)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

#bind results

OSpm<-ldply(univ_results, rbind)
OSom$rownames.b.<-NULL
colnames(OSpm)<-c('cpg', 'HR_OS', 'p_OS')

# bind results for all CpGs
unicoxos<-bind_rows(OS, OSpm)

## 1.2 DFS ####
# exclude ProM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/unicox_all.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/cpgcanvalidate.RData")

library(readxl)
cpgd <- read_excel("CpG_Data extraction.xlsx", 
                   sheet = "all CpGs", col_types = c("text", 
                                                     "skip", "skip", "text", "skip", "skip", 
                                                     "skip", "skip", "skip", "skip", "skip"))

cpgd<-cpgd$CpGs[cpgd$Study=='Gündert et al, 2019']


cpg<-cpg[!(cpg %in% cpgd)]
rm(cpgd)

unicox_crc<-unicox[!is.na(unicox$recurr_cp), ]
univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_crc)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

DFS<-ldply(univ_results, rbind)
DFS$rownames.b.<-NULL
colnames(DFS)<-c('cpg', 'HR_DFS', 'p_DFS')

### validation ProM classifier
cpg <-  read_excel("CpG_Data extraction.xlsx", 
                   sheet = "all CpGs", col_types = c("text", 
                                                     "skip", "skip", "text", "skip", "skip", 
                                                     "skip", "skip", "skip", "skip", "skip"))

cpg<-cpg$CpGs[cpg$Study=='Gündert et al, 2019']

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/PMunicox.RData")

unicoxpm_crc<-unicoxpm[!is.na(unicoxpm$recurr_cp), ]
univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicoxpm_crc)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

DFSpm<-ldply(univ_results, rbind)
DFSpm$rownames.b.<-NULL
colnames(DFSpm)<-c('cpg', 'HR_DFS', 'p_DFS')

unicoxdfs<-bind_rows(DFS, DFSpm)

######################################## 2. Multicox_main adjustment  ##################################################################################
## 2.1 OS ####
## exclude cpgs in PM classifier ##
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicox_impu20list.RData")

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/cpgcanvalidate.RData")

cpgd <- read_excel("CpG_Data extraction.xlsx",
                   sheet = "all CpGs", col_types = c("text",
                                                     "skip", "text", "skip", "skip", "skip",
                                                     "skip", "skip", "skip", "skip"))

cpgd<-cpgd$CpGs[cpgd$Study=='Gündert et al, 2019']

cpg<-cpg[!(cpg %in% cpgd)]
rm(cpgd)

dput(names(multi_impu[[1]]))
multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))


n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[nrow(p), ncol(p)]
  
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[nrow(HRCI), 1])
  LL<-exp(HRCI[nrow(HRCI), 3])
  UL<-exp(HRCI[nrow(HRCI), 4])
  cpgname<-rownames(HRCI)[nrow(HRCI)]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

OS<-ldply(result, rbind)
# keep twp decimals
OS[,c(2:4)]<-round(OS[,c(2:4)], digits = 2)
OS$p_round<-round(OS$p, digits = 3)
OS$HRCI<-paste0(OS$HR, ' ', '(', OS$LL, ',', ' ', OS$UL, ')' )

# add results of PM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicoxPM_impu20list.RData")

cpg <- read_excel("CpG_Data extraction.xlsx",
                  sheet = "all CpGs", col_types = c("text",
                                                    "skip", "text", "skip", "skip", "skip",
                                                    "skip", "skip", "skip", "skip"))

cpgpm<-cpg$CpGs[cpg$Study== 'Gündert et al, 2019']
rm(cpg)
multi_formulas <- sapply(cpgpm,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))


n<-length(cpgpm)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multipm_impu[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[nrow(p), ncol(p)]
  
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[nrow(HRCI), 1])
  LL<-exp(HRCI[nrow(HRCI), 3])
  UL<-exp(HRCI[nrow(HRCI), 4])
  cpgname<-rownames(HRCI)[nrow(HRCI)]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

OSpm<-ldply(result, rbind)
# keep twp decimals
OSpm[,c(2:4)]<-round(OSpm[,c(2:4)], digits = 2)
OSpm$p_round<-round(OSpm$p, digits = 3)
OSpm$HRCI<-paste0(OSpm$HR, ' ', '(', OSpm$LL, ',', ' ', OSpm$UL, ')' )

## merge results together
multicoxos<-bind_rows(OS, OSpm)

## 2.2 DFS ####
# exclude CpGs in PM classifier ##
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/cpgcanvalidate.RData")

cpgd <-read_excel("CpG_Data extraction.xlsx", 
                  sheet = "all CpGs", col_types = c("text", 
                                                    "skip", "skip", "text", "skip", "skip", 
                                                    "skip", "skip", "skip", "skip", "skip"))

cpgd<-cpgd$CpGs[cpgd$Study=='Gündert et al, 2019']

cpg<-cpg[!(cpg %in% cpgd)]
rm(cpgd)

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))


n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[nrow(p), ncol(p)]
  
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[nrow(HRCI), 1])
  LL<-exp(HRCI[nrow(HRCI), 3])
  UL<-exp(HRCI[nrow(HRCI), 4])
  cpgname<-rownames(HRCI)[nrow(HRCI)]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

DFS<-ldply(result, rbind)
# keep twp decimals
DFS[,c(2:4)]<-round(DFS[,c(2:4)], digits = 2)
DFS$p_round<-round(DFS$p, digits = 3)
DFS$HRCI<-paste0(DFS$HR, ' ', '(', DFS$LL, ',', ' ', DFS$UL, ')' )

# add results of PM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicoxPM_impu20list.RData")

cpg <- read_excel("CpG_Data extraction.xlsx", 
                  sheet = "all CpGs", col_types = c("text", 
                                                    "skip", "skip", "text", "skip", "skip", 
                                                    "skip", "skip", "skip", "skip", "skip"))

cpgpm<-cpg$CpGs[cpg$Study=='Gündert et al, 2019']


multi_formulas <- sapply(cpgpm,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))


n<-length(cpgpm)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multipm_impu[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[nrow(p), ncol(p)]
  
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[nrow(HRCI), 1])
  LL<-exp(HRCI[nrow(HRCI), 3])
  UL<-exp(HRCI[nrow(HRCI), 4])
  cpgname<-rownames(HRCI)[nrow(HRCI)]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

DFSpm<-ldply(result, rbind)
# keep twp decimals
DFSpm[,c(2:4)]<-round(DFSpm[,c(2:4)], digits = 2)
DFSpm$p_round<-round(DFSpm$p, digits = 3)
DFSpm$HRCI<-paste0(DFSpm$HR, ' ', '(', DFSpm$LL, ',', ' ', DFSpm$UL, ')' )

## merge results together

multicoxdfs<-bind_rows(DFS, DFSpm)

##################################### 3. Multicox supplement adjustment add MSI  ##################################################################################
## 3.1 OS ####
## exclude cpgs in PM classifier ##
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicox_impu20list.RData")

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/cpgcanvalidate.RData")

cpgd <- read_excel("CpG_Data extraction.xlsx",
                   sheet = "all CpGs", col_types = c("text",
                                                     "skip", "text", "skip", "skip", "skip",
                                                     "skip", "skip", "skip", "skip"))

cpgd<-cpgd$CpGs[cpgd$Study=='Gündert et al, 2019']

cpg<-cpg[!(cpg %in% cpgd)]
rm(cpgd)

dput(names(multi_impu[[1]]))
multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+MSI_gent+chemradther+Stage_at_diagnosis+Location+', x)))


n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[nrow(p), ncol(p)]
  
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[nrow(HRCI), 1])
  LL<-exp(HRCI[nrow(HRCI), 3])
  UL<-exp(HRCI[nrow(HRCI), 4])
  cpgname<-rownames(HRCI)[nrow(HRCI)]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

OS<-ldply(result, rbind)
# keep twp decimals
OS[,c(2:4)]<-round(OS[,c(2:4)], digits = 2)
OS$p_round<-round(OS$p, digits = 3)
OS$HRCI<-paste0(OS$HR, ' ', '(', OS$LL, ',', ' ', OS$UL, ')' )

### add results of PM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicoxPM_impu20list.RData")

cpg <- read_excel("CpG_Data extraction.xlsx",
                  sheet = "all CpGs", col_types = c("text",
                                                    "skip", "text", "skip", "skip", "skip",
                                                    "skip", "skip", "skip", "skip"))

cpgpm<-cpg$CpGs[cpg$Study== 'Gündert et al, 2019']
rm(cpg)
multi_formulas <- sapply(cpgpm,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+MSI_gent+chemradther+Stage_at_diagnosis+Location+', x)))


n<-length(cpgpm)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multipm_impu[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[nrow(p), ncol(p)]
  
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[nrow(HRCI), 1])
  LL<-exp(HRCI[nrow(HRCI), 3])
  UL<-exp(HRCI[nrow(HRCI), 4])
  cpgname<-rownames(HRCI)[nrow(HRCI)]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

OSpm<-ldply(result, rbind)
# keep twp decimals
OSpm[,c(2:4)]<-round(OSpm[,c(2:4)], digits = 2)
OSpm$p_round<-round(OSpm$p, digits = 3)
OSpm$HRCI<-paste0(OSpm$HR, ' ', '(', OSpm$LL, ',', ' ', OSpm$UL, ')' )

## merge results together

multicox2os<-bind_rows(OS, OSpm)

## 3.2 DFS ####
# exclude CpGs in PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/cpgcanvalidate.RData")

cpgd <- read_excel("CpG_Data extraction.xlsx",
                   sheet = "all CpGs", col_types = c("text",
                                                     "skip", "text", "skip", "skip", "skip",
                                                     "skip", "skip", "skip", "skip"))

cpgd<-cpgd$CpGs[cpgd$Study=='Gündert et al, 2019']

cpg<-cpg[!(cpg %in% cpgd)]
rm(cpgd)

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age+Gender+MSI_gent+chemradther+Stage_at_diagnosis+Location+', x)))


n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[nrow(p), ncol(p)]
  
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[nrow(HRCI), 1])
  LL<-exp(HRCI[nrow(HRCI), 3])
  UL<-exp(HRCI[nrow(HRCI), 4])
  cpgname<-rownames(HRCI)[nrow(HRCI)]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

DFS<-ldply(result, rbind)
# keep twp decimals
DFS[,c(2:4)]<-round(DFS[,c(2:4)], digits = 2)
DFS$p_round<-round(DFS$p, digits = 3)
DFS$HRCI<-paste0(DFS$HR, ' ', '(', DFS$LL, ',', ' ', DFS$UL, ')' )

### add results of PM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicoxPM_impu20list.RData")

cpg <- read_excel("CpG_Data extraction.xlsx",
                  sheet = "all CpGs", col_types = c("text",
                                                    "skip", "text", "skip", "skip", "skip",
                                                    "skip", "skip", "skip", "skip"))

cpgpm<-cpg$CpGs[cpg$Study=='Gündert et al, 2019']


multi_formulas <- sapply(cpgpm,
                         function(x) as.formula(paste('Surv(timeD, death_crccp==1)~Age+Gender++MSI_gent+chemradther+Stage_at_diagnosis+Location+', x)))


n<-length(cpgpm)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multipm_impu[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[nrow(p), ncol(p)]
  
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[nrow(HRCI), 1])
  LL<-exp(HRCI[nrow(HRCI), 3])
  UL<-exp(HRCI[nrow(HRCI), 4])
  cpgname<-rownames(HRCI)[nrow(HRCI)]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

DFSpm<-ldply(result, rbind)
# keep twp decimals
DFSpm[,c(2:4)]<-round(DFSpm[,c(2:4)], digits = 2)
DFSpm$p_round<-round(DFSpm$p, digits = 3)
DFSpm$HRCI<-paste0(DFSpm$HR, ' ', '(', DFSpm$LL, ',', ' ', DFSpm$UL, ')' )

## merge results together

multicox2dfs<-bind_rows(DFS, DFSpm)
































