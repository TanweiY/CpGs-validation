library(survival)
library(dplyr)
library(plyr)
library(writexl)
library(tableone)
library(rlist)
library(cgwtools)
library(mice)
library(mitools)
library(survival)
library(dplyr)
library(plyr)
library(writexl)
library(tableone)
library(rlist)
library(cgwtools)
library(mice)
library(mitools)
####################### 1. prepare the data  #####################################################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/covariable_all.RData")

dput(names(covar))
summary(covar)
outcome<-subset(covar, select = c("id", "Stage_at_diagnosis", "Location", "timeD", "timey", "timeD_recurr", "recurr_timey", "death_all", 
                                  "death_crc", "recurr", "death_crccp", "recurr_cp", "PFS", "timeD_PFS", 
                                  "timey_PFS"))

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/beta_new288.RData")

beta_new288$uid<-NULL

# merge clinical variables with CpGs 
unicox<-merge(outcome, beta_new288, by='id')
summary(unicox)

######################### 1.1 calculate the score for each prognostic model ##################################

unicox<-within.data.frame(unicox, {
  
  Yang_2019_C14<-(0.12*cg02196655+1.35*cg03763616+0.73*cg03944089+0.73*cg06117855+0.76*cg07173760-3.96*cg07293947-0.76*cg07509155+0.58*cg09244244+0.4*cg10451565+0.28*cg12582008+1.99*cg13796218+3.6*cg20247048+1.34*cg21481775+0.42*cg23829949-0.28*cg23964386+0.96*cg24127989-0.45*cg24674703+0.84*cg24938727)
  Gong_2020_CR14<-(0.09*cg14660573-0.04*cg09353563 + 2.06*cg00110724)
  WangX_2020_CR14<-(0.919*cg03091331 + 0.963*cg06884352 + 0.703*cg07707546 + 0.721*cg08081805 + 0.587*cg21347353 + 0.528*cg25164589)
  WangY_2020_C14<-(38.52*cg00177496-4.13*cg01963906+ 2.574*cg05165940-79.32*cg12921795+2.31*cg19414598+6.061*cg25783173)
  Xiang_2020_C14<-(3.02*cg03017653) + (22.84*cg03977782) + (3.00*cg05417950) + (36.54*cg06250108) + (2.18*cg09893305) + (3.12*cg10414946) + (3.06*cg15170424) + (5.77*cg15639045) + (-9.51*cg15786837) + (22.02*cg17329249) + (8.12*cg21212956) + (-23.31*cg24206256) 
  ChenF_2021_CR23<- (1.87*cg11621464 + 1.11*cg13565656 -1.74*cg18976437 + 2.40*cg20505223 -1.97*cg20528583)
  Huang_2021_OS_CR14<-(0.754*cg21614638 + 1.230*cg21770617 + 0.701*cg12751565)
  Huang_2021_DFS_CR14<-(0.573*cg21614638 + 0.888*cg21770617 + 0.477*cg12751565)
  Li_2021_CR14<-(-2.1*cg01408654) + (-2.2*cg04035209) + (-3*cg10196720) + (-5*cg10379890) + (2.3*cg11097433) + (4.86*cg14675211) + (2.93*cg15428578) + (-4.7*cg19343464) + (2.7*cg21384402) + (3.7*cg27404023)
  Gundert_2019_CR13<-(0.18*cg23750514) + (0.15*cg01131395) + (0.26*cg12510999) + (0.19*cg18195165) + (0.3*cg05646575) + (0.24*cg11056055) + (0.21*cg10758824) + (0.31*cg16336556) + (0.22*cg19184885) + (0.12*cg19340296) + (0.22*cg22522598) + (0.29*cg17431888) + (0.15*cg18736676) + (0.23*cg08804626) + (0.19*cg08617020) + (0.31*cg14270346) + (0.07*cg16399624) + (0.22*cg14983135) + (0.24*cg00832644)
  Jia_2019_CR14<-(-1.24524*cg16935707) + (-1.32717*cg05481217) + (-1.23953*cg08044454) + (-1.15191*cg01552551) + (-0.62350*cg24311416) + (-0.53487*cg02425108) + (-0.17845*cg15659052)
  Yu_2023_DFS_CR12<-(1.635*cg00561674) + (1.525*cg06887407) + (0.634*cg13598109)  # cg16747321 is not in the array
  
  TGundert_2019_CR13<-as.factor(ntile(Gundert_2019_CR13, 4))
  TLi_2021_CR14<-as.factor(ntile(Li_2021_CR14, 4)) 
  THuang_2021_DFS_CR14<-as.factor(ntile(Huang_2021_DFS_CR14, 4))
  THuang_2021_OS_CR14<-as.factor(ntile(Huang_2021_OS_CR14, 4))
  TChenF_2021_CR23<-as.factor(ntile(ChenF_2021_CR23, 4))
  TXiang_2020_C14<-as.factor(ntile(Xiang_2020_C14, 4))
  TWangY_2020_C14<-as.factor(ntile(WangY_2020_C14, 4))
  TWangX_2020_CR14<-as.factor(ntile(WangX_2020_CR14, 4))
  TGong_2020_CR14<-as.factor(ntile(Gong_2020_CR14, 4))
  TYang_2019_C14<-as.factor(ntile(Yang_2019_C14, 4))
  TJia_2019_CR14<-as.factor(ntile(Jia_2019_CR14, 4))
  TYu_2023_DFS_CR12<-as.factor(ntile(Yu_2023_DFS_CR12, 4))
  
})

# get the name of the 288 cpgs could be validated
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/cpg_validate.RData")
cpg<-cpg_validate

#  normalization Z score for each CpGs for Cox models
unicox[cpg]<-sapply(unicox[cpg], function(data) (data-mean(data))/sd(data))

save(unicox, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_all.RData")


######################### 1.2 prepare the validation cohort for each study ##################################
#1. colon, stage 1-4
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_all.RData")

#subset validation
summary(unicox$Location)
unicox_c14<-subset(unicox, Location=='Distal colon'|Location=='Proximal colon')

summary(unicox$Stage_at_diagnosis)

#subset cpg and score

cpg <- read_excel("data_extraction_newest.xlsx", sheet = "all CpGs")

cpgs<-cpg$CpGs[cpg$Location=='Colon' & cpg$Stage== 'I-IV']

unicox_c14<-subset(unicox_c14, select=c('id', cpgs, 'Yang_2019_C14', 'TYang_2019_C14', 
                                        'WangY_2020_C14', 'TWangY_2020_C14', 'Xiang_2020_C14', 'TXiang_2020_C14',
                                        "timeD", "timey", "timeD_recurr", "recurr_timey", "death_all", 
                                        "death_crc", "recurr", "death_crccp", "recurr_cp", "PFS", "timeD_PFS", "timey_PFS"))

summary(unicox_c14)

save(unicox_c14, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c14.RData")

#2. colon, stage2
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_all.RData")

#subset validation
summary(unicox$Location)
summary(unicox$Stage_at_diagnosis)
unicox_c2<-subset(unicox, Location=='Distal colon'|Location=='Proximal colon')
unicox_c2<-subset(unicox_c2, Stage_at_diagnosis=='II')

#subset cpg and score
unicox_c2<-subset(unicox_c2, select=c('id', 'cg14353573', "timeD", "timey", "timeD_recurr", "recurr_timey", "death_all", "death_crc", "recurr", 
                                      "death_crccp", "recurr_cp", "PFS", "timeD_PFS", "timey_PFS"))

save(unicox_c2, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c2.RData")

#3. colorectal stage1-4
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_all.RData")

cpg <- read_excel("data_extraction_newest.xlsx", sheet = "all CpGs")

cpgs<-subset(cpg, Location=='Colorect' & Stage== 'I-IV')

# exclude jia and neumeyer
cpgs<-cpgs$CpGs[!(cpgs$Study == 'Jia et al, 2019'|cpgs$Study == 'Neumeyer et al, 2019')]

unicox_cr14<-subset(unicox, select=c('id', cpgs, 'Gong_2020_CR14', 'WangX_2020_CR14', 
                                     'Huang_2021_OS_CR14', 'Huang_2021_DFS_CR14', 'Li_2021_CR14', 
                                     'TGong_2020_CR14', 'TWangX_2020_CR14', 
                                     'THuang_2021_OS_CR14', 'THuang_2021_DFS_CR14', 'TLi_2021_CR14',
                                     "timeD", "timey", "timeD_recurr", "recurr_timey", "death_all", "death_crc", "recurr", 
                                     "death_crccp", "recurr_cp", "PFS", "timeD_PFS", "timey_PFS"))



save(unicox_cr14, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")

#4. colorectal stageI-II
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_all.RData")

#subset validation
summary(unicox$Stage_at_diagnosis)
unicox_cr12<-subset(unicox, Stage_at_diagnosis=='I'|Stage_at_diagnosis=='II')

#subset cpg and score

unicox_cr12<-subset(unicox_cr12, select=c('id', 'cg00561674', 'cg06887407', 'cg13598109',
                                          'Yu_2023_DFS_CR12', 'TYu_2023_DFS_CR12',
                                          "timeD", "timey", "timeD_recurr", "recurr_timey", "death_all", "death_crc", "recurr", 
                                          "death_crccp", "recurr_cp", "PFS", "timeD_PFS", "timey_PFS"))

summary(unicox_cr12)

save(unicox_cr12, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr12.RData")

#colorectal stageII-III
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_all.RData")

#subset validation
summary(unicox$Stage_at_diagnosis)
unicox_cr23<-subset(unicox, Stage_at_diagnosis=='II'|Stage_at_diagnosis=='III')

#subset cpg and score
cpg <- read_excel("data_extraction_newest.xlsx", sheet = "all CpGs")

cpgs<-cpg$CpGs[cpg$Location=='Colorect' & cpg$Stage== 'II/III']

unicox_cr23<-subset(unicox_cr23, select=c('id', cpgs, 'ChenF_2021_CR23', 'TChenF_2021_CR23',
                                          "timeD", "timey", "timeD_recurr", "recurr_timey", "death_all", "death_crc", "recurr", 
                                          "death_crccp", "recurr_cp", "PFS", "timeD_PFS", "timey_PFS"))

save(unicox_cr23, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr23.RData")

#5. colorectal stageI-III

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_all.RData")

## exclude duplicate patient ID 
library(readr)
proid <- read_csv("dachs_cohort/ProMprobe.csv")
proid[nrow(proid)+1, ]<-colnames(proid)
colnames(proid)<-'id'

inter<-Reduce(intersect, list(proid$id, unicox$id)) ## 572 overlapping patients

unicox_cr13<-unicox[!(unicox$id %in% inter ), ] # 1738

#subset validation
summary(unicox_cr13$Stage_at_diagnosis)

unicox_cr13<-subset(unicox_cr13, Stage_at_diagnosis=='I'|Stage_at_diagnosis=='III'|Stage_at_diagnosis=='II')

cpg <- read_excel("data_extraction_newest.xlsx", sheet = "all CpGs")

cpgs<-cpg$CpGs[cpg$Location=='Colorect' & cpg$Stage== 'I-III']

unicox_cr13<-subset(unicox_cr13, select=c('id', cpgs, 
                                          'Gundert_2019_CR13', 'TGundert_2019_CR13',
                                          "timeD", "timey", "timeD_recurr", "recurr_timey", "death_all", "death_crc", "recurr", 
                                          "death_crccp", "recurr_cp", "PFS", "timeD_PFS", "timey_PFS"))

save(unicox_cr13, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr13.RData")

## prepare the data for another two studies originating from the DACHS (Jia et al, 2019 and Neumeyer et al, 2019)

# 
library(haven)
jia<-read_sas('/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/clinical_molecular_data.sas7bdat')

summary(as.factor(jia$set))

jia<-subset(jia, select = c('dachs_id', 'set'))
# Min's analysis were focusing on set 1  

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_all.RData")

colnames(jia)[1]<-'id'

jia<-subset(jia, set ==1) 

unicox_jia<-unicox[!(unicox$id %in% jia$id), ] # 1216

cpg <- read_excel("data_extraction_newest.xlsx", sheet = "all CpGs")

cpgs<-cpg$CpGs[cpg$Study == 'Jia et al, 2019'|cpg$Study == 'Neumeyer et al, 2019']

unicox_jiaso<-subset(unicox_jia, select=c('id',  cpgs, 'Jia_2019_CR14', 'TJia_2019_CR14',
                                          "timeD", "timey", "timeD_recurr", "recurr_timey", "death_all", 
                                          "death_crc", "recurr", "death_crccp", "recurr_cp", "PFS", "timeD_PFS", 
                                          "timey_PFS"))
save(unicox_jiaso, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_jiaso.RData")


################################## 2. Univariable Cox for each CpGs and prognostic scores ###################
################################## 2.1 OS ###################
# 1. Colon, stage 2 OS 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c2.RData")

cpg<-colnames(unicox_c2)[2]

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_c2)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)}
univ_results <- lapply(univ_models, uni_coxresult)



OS_c2<-ldply(univ_results, rbind)
OS_c2$rownames.b.<-NULL
colnames(OS_c2)<-c('cpg', 'HR_OS', 'p_OS')
OS_c2$group<-'colon_2'

save(OS_c2, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox.RData")

#2. Colon, stage 1-4 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c14.RData")

cpg<-colnames(unicox_c14)[2: (length(unicox_c14)-12)]

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_c14)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)



OS_c14<-ldply(univ_results, rbind)
OS_c14$.id<-NULL
colnames(OS_c14)<-c('cpg', 'HR_OS', 'p_OS')
OS_c14$group<-'colon_14'

resave(OS_c14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox.RData")

#3. Colon and rectum, stage 1-2 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr12.RData")

cpg<-colnames(unicox_cr12)[2: (length(unicox_cr12)-12)]

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr12)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)



OS_cr12<-ldply(univ_results, rbind)
OS_cr12$.id<-NULL
colnames(OS_cr12)<-c('cpg', 'HR_OS', 'p_OS')
OS_cr12$group<-'cr12'

resave(OS_cr12, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox.RData")

# 4. Colon and rectum, stage 1-3 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr13.RData")

cpg<-colnames(unicox_cr13)[2: (length(unicox_cr13)-12)]

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr13)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

OS_cr13<-ldply(univ_results, rbind)
OS_cr13$.id<-NULL
colnames(OS_cr13)<-c('cpg', 'HR_OS', 'p_OS')
OS_cr13$group<-'cr13'

resave(OS_cr13, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox.RData")

# Colon and rectum, stage 1-4 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")

cpg<-colnames(unicox_cr14)[2: (length(unicox_cr14)-12)]

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr14)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)



OS_cr14<-ldply(univ_results, rbind)
OS_cr14$.id<-NULL
colnames(OS_cr14)<-c('cpg', 'HR_OS', 'p_OS')
OS_cr14$group<-'cr14'

resave(OS_cr14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox.RData")

# Colon and rectum + stage 2-3 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr23.RData")

cpg<-colnames(unicox_cr23)[2: (length(unicox_cr23)-12)]

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr23)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)



OS_cr23<-ldply(univ_results, rbind)
OS_cr23$.id<-NULL
colnames(OS_cr23)<-c('cpg', 'HR_OS', 'p_OS')
OS_cr23$group<-'cr23'

resave(OS_cr23, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox.RData")

# Colon and rectum, stage 2-3 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_jiaso.RData")

cpg<-colnames(unicox_jiaso)[2: (length(unicox_jiaso)-12)]

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_jiaso)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

OS_jiaso<-ldply(univ_results, rbind)
OS_jiaso$.id<-NULL
colnames(OS_jiaso)<-c('cpg', 'HR_OS', 'p_OS')
OS_jiaso$group<-'jiaso'

resave(OS_jiaso, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox.RData")

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox.RData")
ls()
unicox_comb<-rbind(OS_c14,
                   OS_c2, OS_cr12, OS_cr13,OS_cr14,OS_cr23, OS_jiaso)

unicox_comb<-unicox_comb[order(unicox_comb$cpg), ]

save(unicox_comb,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox_OS_combine.RData")

########################### 2.2 DFS #######################################
#1. Colon, stage2 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c2.RData")
unicox_c2<-unicox_c2[!is.na(unicox_c2$recurr_cp), ]

cpg<-colnames(unicox_c2)[2]

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_c2)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)



DFS_c2<-ldply(univ_results, rbind)
DFS_c2$.id<-NULL
colnames(DFS_c2)<-c('cpg', 'HR_DFS', 'p_DFS')
DFS_c2$group<-'colon_2'

save(DFS_c2, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox_DFS.RData")

#2. Colon, stage1-4  
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c14.RData")
unicox_c14<-unicox_c14[!is.na(unicox_c14$recurr_cp), ]

cpg<-colnames(unicox_c14)[2: (length(unicox_c14)-12)]

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_c14)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)



DFS_c14<-ldply(univ_results, rbind)
DFS_c14$.id<-NULL
colnames(DFS_c14)<-c('cpg', 'HR_DFS', 'p_DFS')
DFS_c14$group<-'c14'

resave(DFS_c14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox_DFS.RData")

# 3. Colon and rectum, stage 1-2 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr12.RData")
unicox_cr12<-unicox_cr12[!is.na(unicox_cr12$recurr_cp), ]

cpg<-colnames(unicox_cr12)[2: (length(unicox_cr12)-12)]

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr12)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)



DFS_cr12<-ldply(univ_results, rbind)
DFS_cr12$.id<-NULL
colnames(DFS_cr12)<-c('cpg', 'HR_DFS', 'p_DFS')
DFS_cr12$group<-'cr_14'

resave(DFS_cr12, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox_DFS.RData")

#4. Colon and rectum, stage 1-3 ########
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr13.RData")
unicox_cr13<-unicox_cr13[!is.na(unicox_cr13$recurr_cp), ]

cpg<-colnames(unicox_cr13)[2: (length(unicox_cr13)-12)]

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr13)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)



DFS_cr13<-ldply(univ_results, rbind)
DFS_cr13$.id<-NULL
colnames(DFS_cr13)<-c('cpg', 'HR_DFS', 'p_DFS')
DFS_cr13$group<-'cr13'

resave(DFS_cr13, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox_DFS.RData")

#5. Colon and rectum, stage 1-4 ########
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")
unicox_cr14<-unicox_cr14[!is.na(unicox_cr14$recurr_cp), ]

cpg<-colnames(unicox_cr14)[2: (length(unicox_cr14)-12)]

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr14)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)



DFS_cr14<-ldply(univ_results, rbind)
DFS_cr14$.id<-NULL
colnames(DFS_cr14)<-c('cpg', 'HR_DFS', 'p_DFS')
DFS_cr14$group<-'cr14'

resave(DFS_cr14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox_DFS.RData")

#6. Colon and rectum, stage 2-4 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr23.RData")
unicox_cr23<-unicox_cr23[!is.na(unicox_cr23$recurr_cp), ]

cpg<-colnames(unicox_cr23)[2: (length(unicox_cr23)-12)]

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_cr23)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

DFS_cr23<-ldply(univ_results, rbind)
DFS_cr23$.id<-NULL
colnames(DFS_cr23)<-c('cpg', 'HR_DFS', 'p_DFS')
DFS_cr23$group<-'cr23'

resave(DFS_cr23, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox_DFS.RData")

# 7. Colon and rectum, two studies from DACHS  

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_jiaso.RData")
unicox_jiaso<-unicox_jiaso[!is.na(unicox_jiaso$recurr_cp), ]

cpg<-colnames(unicox_jiaso)[2: (length(unicox_jiaso)-12)]

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = unicox_jiaso)})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)


DFS_jiaso<-ldply(univ_results, rbind)
DFS_jiaso$.id<-NULL
colnames(DFS_jiaso)<-c('cpg', 'HR_DFS', 'p_DFS')
DFS_jiaso$group<-'jiaso'

resave(DFS_jiaso, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox_DFS.RData")

# combine all the results 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox_DFS.RData")
ls()
unicox_comb<-rbind(DFS_c14,
                   DFS_c2, DFS_cr12, DFS_cr13, DFS_cr14, DFS_cr23, DFS_jiaso)

unicox_comb<-unicox_comb[order(unicox_comb$cpg), ]

save(unicox_comb,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/unicox_DFS_combine.RData")

################################## 3. Multivariable Cox for each CpGs and prognostic scores ###################
################################## 3.1 Primary multivariable model ######################################
################################## 3.1.1 OS ######################################
# 1. Colon, stage2 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c2.RData")

multi_impu_c2<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_c2[[i]]<-merge(unicox_c2, multi_impu[[i]], by='id')
}

dput(names(multi_impu_c2[[1]]))

cpg<-colnames(multi_impu_c2[[1]])[2: (length(multi_impu_c2[[1]])-22)]

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_c2[[j]])
  }
  
  p<-summary(pool(mice.fit))
  
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

OS_c2<-ldply(result, rbind)

OS_c2[,c(2:4)]<-round(OS_c2[,c(2:4)], digits = 2)
OS_c2$p_round<-round(OS_c2$p, digits = 3)
OS_c2$HRCI<-paste0(OS_c2$HR, ' ', '(', OS_c2$LL, ',', ' ', OS_c2$UL, ')' )

colnames(OS_c2)[7]<-'HRCI_OS'

save(OS_c2, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_OS.RData")

#2. Colon, stage 1-4 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c14.RData")

multi_impu_c14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_c14[[i]]<-merge(unicox_c14, multi_impu[[i]], by='id')
}

dput(names(multi_impu_c14[[1]]))

cpg<-colnames(multi_impu_c14[[1]])[2: (length(multi_impu_c14[[1]])-22)]

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_c14[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

OS_c14<-ldply(result, rbind)

OS_c14[,c(2:4)]<-round(OS_c14[,c(2:4)], digits = 2)
OS_c14$p_round<-round(OS_c14$p, digits = 3)
OS_c14$HRCI<-paste0(OS_c14$HR, ' ', '(', OS_c14$LL, ',', ' ', OS_c14$UL, ')' )

colnames(OS_c14)[7]<-'HRCI_OS'

resave(OS_c14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_OS.RData")

#Colon and rectum stage 1-2 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr12.RData")

multi_impu_cr12<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr12[[i]]<-merge(unicox_cr12, multi_impu[[i]], by='id')
}

dput(names(multi_impu_cr12[[1]]))

cpg<-colnames(multi_impu_cr12[[1]])[2: (length(multi_impu_cr12[[1]])-22)]

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_cr12[[j]])
  }
  
  p<-summary(pool(mice.fit))
  
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

OS_cr12<-ldply(result, rbind)

OS_cr12[,c(2:4)]<-round(OS_cr12[,c(2:4)], digits = 2)
OS_cr12$p_round<-round(OS_cr12$p, digits = 3)
OS_cr12$HRCI<-paste0(OS_cr12$HR, ' ', '(', OS_cr12$LL, ',', ' ', OS_cr12$UL, ')' )

colnames(OS_cr12)[7]<-'HRCI_OS'

resave(OS_cr12, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_OS.RData")

# Colon and rectum, stage 1-3 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr13.RData")

multi_impu_cr13<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr13[[i]]<-merge(unicox_cr13, multi_impu[[i]], by='id')
}

dput(names(multi_impu_cr13[[1]]))

cpg<-colnames(multi_impu_cr13[[1]])[2: (length(multi_impu_cr13[[1]])-22)]

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_cr13[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

OS_cr13<-ldply(result, rbind)

OS_cr13[,c(2:4)]<-round(OS_cr13[,c(2:4)], digits = 2)
OS_cr13$p_round<-round(OS_cr13$p, digits = 3)
OS_cr13$HRCI<-paste0(OS_cr13$HR, ' ', '(', OS_cr13$LL, ',', ' ', OS_cr13$UL, ')' )

colnames(OS_cr13)[7]<-'HRCI_OS'

resave(OS_cr13, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_OS.RData")

# Colon and rectum, stage 1-4
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")

multi_impu_cr14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr14[[i]]<-merge(unicox_cr14, multi_impu[[i]], by='id')
}

dput(names(multi_impu_cr14[[1]]))

cpg<-colnames(multi_impu_cr14[[1]])[2: (length(multi_impu_cr14[[1]])-22)]

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_cr14[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

OS_cr14<-ldply(result, rbind)

OS_cr14[,c(2:4)]<-round(OS_cr14[,c(2:4)], digits = 2)
OS_cr14$p_round<-round(OS_cr14$p, digits = 3)
OS_cr14$HRCI<-paste0(OS_cr14$HR, ' ', '(', OS_cr14$LL, ',', ' ', OS_cr14$UL, ')' )

colnames(OS_cr14)[7]<-'HRCI_OS'

resave(OS_cr14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_OS.RData")

## Colon and rectum, stage 2-3 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr23.RData")

multi_impu_cr23<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr23[[i]]<-merge(unicox_cr23, multi_impu[[i]], by='id')
}

dput(names(multi_impu_cr23[[1]]))

cpg<-colnames(multi_impu_cr23[[1]])[2: (length(multi_impu_cr23[[1]])-22)]

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_cr23[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

OS_cr23<-ldply(result, rbind)

OS_cr23[,c(2:4)]<-round(OS_cr23[,c(2:4)], digits = 2)
OS_cr23$p_round<-round(OS_cr23$p, digits = 3)
OS_cr23$HRCI<-paste0(OS_cr23$HR, ' ', '(', OS_cr23$LL, ',', ' ', OS_cr23$UL, ')' )

colnames(OS_cr23)[7]<-'HRCI_OS'

resave(OS_cr23, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_OS.RData")

# two studies from DACHS
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_jiaso.RData")

multi_impu_jiaso<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_jiaso[[i]]<-merge(unicox_jiaso, multi_impu[[i]], by='id')
}

dput(names(multi_impu_jiaso[[1]]))

cpg<-colnames(multi_impu_jiaso[[1]])[2: (length(multi_impu_jiaso[[1]])-22)]

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_jiaso[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

OS_jiaso<-ldply(result, rbind)

OS_jiaso[,c(2:4)]<-round(OS_jiaso[,c(2:4)], digits = 2)
OS_jiaso$p_round<-round(OS_jiaso$p, digits = 3)
OS_jiaso$HRCI<-paste0(OS_jiaso$HR, ' ', '(', OS_jiaso$LL, ',', ' ', OS_jiaso$UL, ')' )

colnames(OS_jiaso)[7]<-'HRCI_OS'

resave(OS_jiaso, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_OS.RData")

# combine all the results 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_OS.RData")
ls()

multicox_comb<-rbind(OS_c14,
                     OS_c2, OS_cr12, OS_cr13, OS_cr14, OS_cr23, OS_jiaso)

## remove needless stuff
multicox_comb <- multicox_comb[multicox_comb[, "cpgname"] != 	'LocationProximal colon', ]
multicox_comb <- multicox_comb[multicox_comb[, "cpgname"] != 	'LocationRectum', ]

save(multicox_comb,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_OS_combine.RData")

######################################### 3.1.2 DFS #########################################
# Colon, stage2 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c2.RData")

multi_impu_c2<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_c2[[i]]<-merge(unicox_c2, multi_impu[[i]], by='id')
  multi_impu_c2[[i]]<-multi_impu_c2[[i]][!is.na(multi_impu_c2[[i]]$recurr_cp), ]
  
}

dput(names(multi_impu_c2[[1]]))

cpg<-colnames(multi_impu_c2[[1]])[2: (length(multi_impu_c2[[1]])-22)]

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_c2[[j]])
  }
  
  p<-summary(pool(mice.fit))
  
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

DFS_c2<-ldply(result, rbind)

DFS_c2[,c(2:4)]<-round(DFS_c2[,c(2:4)], digits = 2)
DFS_c2$p_round<-round(DFS_c2$p, digits = 3)
DFS_c2$HRCI<-paste0(DFS_c2$HR, ' ', '(', DFS_c2$LL, ',', ' ', DFS_c2$UL, ')' )

colnames(DFS_c2)[7]<-'HRCI_DFS'

save(DFS_c2, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_DFS.RData")

# Colon, stage 1-4 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c14.RData")

multi_impu_c14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_c14[[i]]<-merge(unicox_c14, multi_impu[[i]], by='id')
  multi_impu_c14[[i]]<-multi_impu_c14[[i]][!is.na(multi_impu_c14[[i]]$recurr_cp), ]
}

dput(names(multi_impu_c14[[1]]))

cpg<-colnames(multi_impu_c14[[1]])[2: (length(multi_impu_c14[[1]])-22)]

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_c14[[j]])
  }
  
  p<-summary(pool(mice.fit))
  
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

DFS_c14<-ldply(result, rbind)

DFS_c14[,c(2:4)]<-round(DFS_c14[,c(2:4)], digits = 2)
DFS_c14$p_round<-round(DFS_c14$p, digits = 3)
DFS_c14$HRCI<-paste0(DFS_c14$HR, ' ', '(', DFS_c14$LL, ',', ' ', DFS_c14$UL, ')' )

colnames(DFS_c14)[7]<-'HRCI_DFS'

resave(DFS_c14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_DFS.RData")

# Colon and rectum, stage 1-2 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr12.RData")

multi_impu_cr12<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr12[[i]]<-merge(unicox_cr12, multi_impu[[i]], by='id')
  multi_impu_cr12[[i]]<-multi_impu_cr12[[i]][!is.na(multi_impu_cr12[[i]]$recurr_cp), ]
  
}

dput(names(multi_impu_cr12[[1]]))

cpg<-colnames(multi_impu_cr12[[1]])[2: (length(multi_impu_cr12[[1]])-22)]

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_cr12[[j]])
  }
  
  p<-summary(pool(mice.fit))
  
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

DFS_cr12<-ldply(result, rbind)

DFS_cr12[,c(2:4)]<-round(DFS_cr12[,c(2:4)], digits = 2)
DFS_cr12$p_round<-round(DFS_cr12$p, digits = 3)
DFS_cr12$HRCI<-paste0(DFS_cr12$HR, ' ', '(', DFS_cr12$LL, ',', ' ', DFS_cr12$UL, ')' )

colnames(DFS_cr12)[7]<-'HRCI_DFS'

resave(DFS_cr12, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_DFS.RData")

# Colon and rectum, stage 1-3 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr13.RData")
# 1369

multi_impu_cr13<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr13[[i]]<-merge(unicox_cr13, multi_impu[[i]], by='id', all.x = TRUE)
  multi_impu_cr13[[i]]<-multi_impu_cr13[[i]][!is.na(multi_impu_cr13[[i]]$recurr_cp), ]
  
}

dput(names(multi_impu_cr13[[1]]))

cpg<-colnames(multi_impu_cr13[[1]])[2: (length(multi_impu_cr13[[1]])-22)]

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_cr13[[j]])
  }
  
  p<-summary(pool(mice.fit))
  
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

DFS_cr13<-ldply(result, rbind)

DFS_cr13[,c(2:4)]<-round(DFS_cr13[,c(2:4)], digits = 2)
DFS_cr13$p_round<-round(DFS_cr13$p, digits = 3)
DFS_cr13$HRCI<-paste0(DFS_cr13$HR, ' ', '(', DFS_cr13$LL, ',', ' ', DFS_cr13$UL, ')' )

colnames(DFS_cr13)[7]<-'HRCI_DFS'

resave(DFS_cr13, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_DFS.RData")

# Colon and rectum, stage 1-4 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")

multi_impu_cr14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr14[[i]]<-merge(unicox_cr14, multi_impu[[i]], by='id')
  multi_impu_cr14[[i]]<-multi_impu_cr14[[i]][!is.na(multi_impu_cr14[[i]]$recurr_cp), ]
  
}


summary(multi_impu_cr14[[1]]$timeD_recurr)

cpg<-colnames(multi_impu_cr14[[1]])[2: (length(multi_impu_cr14[[1]])-22)]

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_cr14[[j]])
  }
  
  p<-summary(pool(mice.fit))
  
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

DFS_cr14<-ldply(result, rbind)

DFS_cr14[,c(2:4)]<-round(DFS_cr14[,c(2:4)], digits = 2)
DFS_cr14$p_round<-round(DFS_cr14$p, digits = 3)
DFS_cr14$HRCI<-paste0(DFS_cr14$HR, ' ', '(', DFS_cr14$LL, ',', ' ', DFS_cr14$UL, ')' )

colnames(DFS_cr14)[7]<-'HRCI_DFS'

resave(DFS_cr14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_DFS.RData")


# Colon and rectum, stage 2-3 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr23.RData")

multi_impu_cr23<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr23[[i]]<-merge(unicox_cr23, multi_impu[[i]], by='id')
  multi_impu_cr23[[i]]<-multi_impu_cr23[[i]][!is.na(multi_impu_cr23[[i]]$recurr_cp), ]
  
}

summary(multi_impu_cr23[[1]]$timeD_recurr)

cpg<-colnames(multi_impu_cr23[[1]])[2: (length(multi_impu_cr23[[1]])-22)]

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_cr23[[j]])
  }
  
  p<-summary(pool(mice.fit))
  
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

DFS_cr23<-ldply(result, rbind)

DFS_cr23[,c(2:4)]<-round(DFS_cr23[,c(2:4)], digits = 2)
DFS_cr23$p_round<-round(DFS_cr23$p, digits = 3)
DFS_cr23$HRCI<-paste0(DFS_cr23$HR, ' ', '(', DFS_cr23$LL, ',', ' ', DFS_cr23$UL, ')' )

colnames(DFS_cr23)[7]<-'HRCI_DFS'

resave(DFS_cr23, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_DFS.RData")

# Colon and rectum, stage 1-4, two studies from DACHS 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_jiaso.RData")

multi_impu_jiaso<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_jiaso[[i]]<-merge(unicox_jiaso, multi_impu[[i]], by='id')
  multi_impu_jiaso[[i]]<-multi_impu_jiaso[[i]][!is.na(multi_impu_jiaso[[i]]$recurr_cp), ]
  
}

summary(multi_impu_jiaso[[1]]$timeD_recurr)

cpg<-colnames(multi_impu_jiaso[[1]])[2: (length(multi_impu_jiaso[[1]])-22)]

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_jiaso[[j]])
  }
  
  p<-summary(pool(mice.fit))
  
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

DFS_jiaso<-ldply(result, rbind)

DFS_jiaso[,c(2:4)]<-round(DFS_jiaso[,c(2:4)], digits = 2)
DFS_jiaso$p_round<-round(DFS_jiaso$p, digits = 3)
DFS_jiaso$HRCI<-paste0(DFS_jiaso$HR, ' ', '(', DFS_jiaso$LL, ',', ' ', DFS_jiaso$UL, ')' )

colnames(DFS_jiaso)[7]<-'HRCI_DFS'

resave(DFS_jiaso, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_DFS.RData")


# combine all the results 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_DFS.RData")
ls()

multicox_comb<-rbind(DFS_c14,
                     DFS_c2, DFS_cr12, DFS_cr13, DFS_cr14, DFS_cr23, DFS_jiaso)

## remove needless stuff
multicox_comb <- multicox_comb[multicox_comb[, "cpgname"] != 	'LocationProximal colon', ]
multicox_comb <- multicox_comb[multicox_comb[, "cpgname"] != 	'LocationRectum', ]

save(multicox_comb,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/multicox_DFS_combine.RData")


################################## 3.2 Sensitivity multivariable model, add MSI ######################################
# repeat the process of 3.1, aside from adding MSI onto the multivariable Cox Model


load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")

################################## 3.3 validation location-specific cpgs in Xin et al ######################################

## left side CRC --> include both distal and rectum #####

multi_impu_cr14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr14[[i]]<-merge(unicox_cr14, multi_impu[[i]], by='id')
  multi_impu_cr14[[i]]<-subset(multi_impu_cr14[[i]], Location == 'Distal colon'|Location == 'Rectum')
}

summary(multi_impu_cr14[[1]]$Location)

cpga <- read_excel("data_extraction_newest.xlsx",sheet = "all CpGs")
###left

cpgl<-subset(cpga, Note =='Left')
cpg<-cpgl$CpGs

# OS 
#unicox
univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD, death_all)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = multi_impu_cr14[[1]])})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)

OS_unil<-ldply(univ_results, rbind)
OS_unil$.id<-NULL
colnames(OS_unil)<-c('cpg', 'HR_OS_uni', 'p_OS_uni')
save(OS_unil, 
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/xin_left.RData')

# primary multicox
multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_cr14[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

OS_cr14<-ldply(result, rbind)

OS_cr14[,c(2:4)]<-round(OS_cr14[,c(2:4)], digits = 2)
OS_cr14$p_round<-round(OS_cr14$p, digits = 3)
OS_cr14$HRCI<-paste0(OS_cr14$HR, ' ', '(', OS_cr14$LL, ',', ' ', OS_cr14$UL, ')' )

colnames(OS_cr14)[7]<-'HRCI_OS'
OS_cr14 <- OS_cr14[OS_cr14[, "cpgname"] != 	'LocationProximal colon', ]
OS_cr14 <- OS_cr14[OS_cr14[, "cpgname"] != 	'LocationRectum', ]
OS_cr14<-subset(OS_cr14, select = c('cpgname', 'HRCI_OS', 'p_round'))
colnames(OS_cr14)<-c('cpg', 'HR_OS_multi1', 'p_OS_multi1')
OS_multi1l<-OS_cr14

resave(OS_multi1l, 
       file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/xin_left.RData')

# sensitivity multicox, add MSI 
multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD, death_all)~Age+Gender+MSI_gent+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_cr14[[j]])
  }
  
  p<-summary(pool(mice.fit))
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

OS_cr14<-ldply(result, rbind)

OS_cr14[,c(2:4)]<-round(OS_cr14[,c(2:4)], digits = 2)
OS_cr14$p_round<-round(OS_cr14$p, digits = 3)
OS_cr14$HRCI<-paste0(OS_cr14$HR, ' ', '(', OS_cr14$LL, ',', ' ', OS_cr14$UL, ')' )

colnames(OS_cr14)[7]<-'HRCI_OS'
OS_cr14 <- OS_cr14[OS_cr14[, "cpgname"] != 	'LocationProximal colon', ]
OS_cr14 <- OS_cr14[OS_cr14[, "cpgname"] != 	'LocationRectum', ]
OS_cr14<-subset(OS_cr14, select = c('cpgname', 'HRCI_OS', 'p_round'))
colnames(OS_cr14)<-c('cpg', 'HR_OS_multi1', 'p_OS_multi1')
OS_multi2l<-OS_cr14

resave(OS_multi2l, 
       file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/xin_left.RData')

# DFS 
for (i in 1:20) {
  multi_impu_cr14[[i]]<-multi_impu_cr14[[i]][!is.na(multi_impu_cr14[[i]]$recurr_cp), ]
}

univ_formulas <- sapply(cpg,
                        function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = multi_impu_cr14[[1]])})

uni_coxresult<-function(x){
  b<-as.data.frame(ShowRegTable(x, exp = TRUE, digits = 2, pDigits = 3,
                                printToggle = TRUE, ciFun = confint), optional=T, fix.empty.names=T);b<-data.frame(rownames(b),b)
                                return(b)
}
univ_results <- lapply(univ_models, uni_coxresult)


DFS_unil<-ldply(univ_results, rbind)
DFS_unil$.id<-NULL
colnames(DFS_unil)<-c('cpg', 'HR_DFS_uni', 'p_DFS_uni')

resave(DFS_unil, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/xin_left.RData")

# primary multicox
multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_cr14[[j]])
  }
  
  p<-summary(pool(mice.fit))
  
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

DFS_cr14<-ldply(result, rbind)

DFS_cr14[,c(2:4)]<-round(DFS_cr14[,c(2:4)], digits = 2)
DFS_cr14$p_round<-round(DFS_cr14$p, digits = 3)
DFS_cr14$HRCI<-paste0(DFS_cr14$HR, ' ', '(', DFS_cr14$LL, ',', ' ', DFS_cr14$UL, ')' )

colnames(DFS_cr14)[7]<-'HRCI_DFS'
DFS_cr14 <- DFS_cr14[DFS_cr14[, "cpgname"] != 	'LocationProximal colon', ]
DFS_cr14 <- DFS_cr14[DFS_cr14[, "cpgname"] != 	'LocationRectum', ]
DFS_cr14<-subset(DFS_cr14, select = c('cpgname', 'HRCI_DFS', 'p_round'))
colnames(DFS_cr14)<-c('cpg', 'HR_DFS_multi1', 'p_DFS_multi1')
DFS_multi1l<-DFS_cr14

resave(DFS_multi1l, 
       file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/xin_left.RData')

# sensitivity multicox, add MSI

multi_formulas <- sapply(cpg,
                         function(x) as.formula(paste('Surv(timeD_recurr, recurr_cp==1)~Age+Gender+MSI_gent+chemradther+Stage_at_diagnosis+Location+', x)))

n<-length(cpg)

result<-vector(n, mode="list")

for (i in 1:n) {
  n_impu<-20
  mice.fit<-vector(n_impu,mode="list")
  
  for (j in 1:n_impu) {
    mice.fit[[j]]<-coxph(multi_formulas[[i]], data = multi_impu_cr14[[j]])
  }
  
  p<-summary(pool(mice.fit))
  
  p<-p[c(nrow(p)-2, nrow(p)-1, nrow(p)), ncol(p)]
  HRCI<-summary(MIcombine(mice.fit))
  HR<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 1])
  LL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 3])
  UL<-exp(HRCI[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI)), 4])
  cpgname<-rownames(HRCI)[c(nrow(HRCI)-2, nrow(HRCI)-1, nrow(HRCI))]
  
  result[[i]]<-data.frame(cpgname, HR, LL, UL, p)
  
}

DFS_cr14<-ldply(result, rbind)

DFS_cr14[,c(2:4)]<-round(DFS_cr14[,c(2:4)], digits = 2)
DFS_cr14$p_round<-round(DFS_cr14$p, digits = 3)
DFS_cr14$HRCI<-paste0(DFS_cr14$HR, ' ', '(', DFS_cr14$LL, ',', ' ', DFS_cr14$UL, ')' )

colnames(DFS_cr14)[7]<-'HRCI_DFS'
DFS_cr14 <- DFS_cr14[DFS_cr14[, "cpgname"] != 	'LocationProximal colon', ]
DFS_cr14 <- DFS_cr14[DFS_cr14[, "cpgname"] != 	'LocationRectum', ]
DFS_cr14<-subset(DFS_cr14, select = c('cpgname', 'HRCI_DFS', 'p_round'))
colnames(DFS_cr14)<-c('cpg', 'HR_DFS_multi1', 'p_DFS_multi1')
DFS_multi2l<-DFS_cr14

resave(DFS_multi2l,
       file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/xin_left.RData')


## bind all results 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/xin_left.RData")

os<-cbind(OS_unil, OS_multi1l, OS_multi2l)
dfs<-cbind(DFS_unil, DFS_multi1l, DFS_multi2l)
xinleft<-cbind(os, dfs)
write.csv(xinleft,
          '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/tablesfigures/xin_left.csv')




