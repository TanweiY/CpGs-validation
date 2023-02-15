######################################## 1. DNA methylation preprocessing  ##################################################################################
library(ChAMP)
gc()

## 1.1 load the methylation array ####
myLoad <- champ.import(directory = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Data_methylation450")

# Set two filters: 1) CpG with p>0.01 set to NA; 2) CpGs with less thab 3 beads

myfilter2<-champ.filter(beta=myLoad$beta,
                        pd=myLoad$pd,
                        filterDetP=TRUE,
                        detPcut=0.01,  ## The detection p-value threshhold. Probes about this cutoff will be filtered out. 
                        detP=myLoad$detP,
                        ProbeCutoff=1,
                        SampleCutoff=1,
                        
                        beadcount=myLoad$beadcount, 
                        filterBeads=TRUE, ## probes with less then 3 beads would be set NA
                        beadCutoff = 1,
                        
                        autoimpute=FALSE,
                        filterNoCG = FALSE,
                        filterSNPs = FALSE,
                        filterMultiHit = FALSE,
                        filterXY = FALSE,
                        fixOutlier = FALSE,
                        arraytype = "450K")

save(myfilter2, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/myfilter2.RData')


pd<-myfilter2$pd
pd$Sample_Group<-substring(pd$Sample_ID, first = 11, last = 12)
summary(as.factor(pd$Sample_Group))
pd$Sample_Group[!(pd$Sample_Group=='NG' | pd$Sample_Group=='TU')]<-'TU'
save(pd, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/pd.RData")

## 1.2 Normalization ####
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/myfilter2.RData")

myLoad<-myfilter2

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=100)

## 1.3 check and correct batch effect ####

# Scree and SVD plot to check components needed to be considered in correcting batch effect (80%)

champ.SVD(beta=myNorm, pd=myLoad$pd)

pd<-myLoad$pd
dput(names(pd))

myCombat <- champ.runCombat(beta=myNorm, pd=myLoad$pd, 
                            batchname=c("Array", "Sample_Plate" ))


# check the corrected result

champ.SVD(beta=myCombat,pd=myLoad$pd)

## 1.4 Select the targeted CpGs for further analysis ####
library(readxl)
CpG_Data_extraction <- read_excel("CpG_Data extraction.xlsx", 
                                  sheet = "all CpGs")


cpg_orig<-CpG_Data_extraction$CpGs
cpg_unique<-unique(cpg_orig) ## 290 unique

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/myCombat.RData")
beta<-myCombat

beta<-as.data.frame(beta)
methy<-rownames(beta)

intersect<-Reduce(intersect, list(cpg_unique, methy)) # only 288 CpGs were contained

# two cpgs cannot be validated:
setdiff(cpg_unique, intersect)

"cg10022744" "cg16747321"

# select CpGs
beta_vali<-beta[intersect, ]
rm(beta)
beta_vali$cpg<-rownames(beta_vali)
beta_vali[nrow(beta_vali)+1, ]<-colnames(beta_vali)
library(data.table)
beta_vali<-transpose(beta_vali)
colnames(beta_vali)<-beta_vali[2681, ]

names(beta_vali)[289]<-'id'

library(dplyr)
beta_vali <- beta_vali %>%
  select(id, everything())
beta_vali<-beta_vali[-2681, ]

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/pd.RData")
sum(beta_vali$id == pd$Sample_Name) ## 2680 
beta_vali$id<-pd$Sample_ID
beta_vali$uid<-pd$Sample_Name
beta_vali <- beta_vali %>%
  select(uid, everything())

# slice and create id, remove normal tissue sample
beta_vali$bat <-strtrim(beta_vali$id, 2)
summary(as.factor(beta_vali$bat))

beta_vali$case<-substring(beta_vali$id, first = 11, last = 12)

summary(as.factor(beta_vali$case))

beta_vali$id[!(beta_vali$case=='NG' | beta_vali$case=='TU')] # all tumor

beta_vali$case[!(beta_vali$case=='NG' | beta_vali$case=='TU')]<-'TU'

summary(as.factor(beta_vali$case))

beta_vali<-subset(beta_vali, case=='TU') ## 

save(beta_vali, 
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/allmethbeforeremovedup.RData')


# process hipo first
hipo<-subset(beta_vali, bat=='H0')
rm(beta_vali)

# remove duplicate
dupid<-hipo[duplicated(hipo$id), ]$id 
dupid<-unique(dupid)
dp<-hipo[hipo$id %in% dupid, ]
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/pd.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/myCombat.RData")
beta<-myCombat[, pd$Sample_Name %in% dp$uid]
pd<-pd[pd$Sample_Name %in% dp$uid, ]

QC.GUI(beta = beta, pheno=pd$Sample_Group)

## select 3998568005_R01C01, and exclude others
hipo[duplicated(hipo$id), ]$uid

exclud<-dp[-1,1]
save(exclud, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/duplicatestoexclude.RData")

hipo<-subset(hipo, !(hipo$uid %in% exclud))

dput(names(hipo))
hipo$bat<-NULL
hipo$case<-NULL

hipo[,-c(1,2)]<-lapply(hipo[,-c(1,2)], as.numeric)
summary(hipo)

# matach names with DACHS id
library(readxl)
hipmatch <- read_excel("validation_processeddata/hipmatch.xlsx")
View(hipmatch)

hipo_final<-merge(hipmatch, hipo, by='id')
hipo_final[duplicated(hipo_final$id), ]$id 

hipo_final[duplicated(hipo_final$tn_id), ]$tn_id
hipo_final$id<-NULL
colnames(hipo_final)[1]<-'id'

save(hipo_final, 
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/hipo_processed.RData')

# then process the other 2 arrays
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/allmethbeforeremovedup.RData")
summary(as.factor(beta_vali$bat))
beta_vali<-subset(beta_vali, bat=='DA')

beta_vali$tn_id <-strtrim(beta_vali$id, 9)

beta_vali <- beta_vali %>%
  select(tn_id, everything())

# check and exclude duplicate # tn_id, rather than id
# we plot the distribution plot using QC.GUI(), and exclude the ID with abnormal distribution
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/duplicatestoexclude2array.RData")

beta_vali_final<-subset(beta_vali, !(beta_vali$uid %in% exclude))

sum(duplicated(beta_vali_final$tn_id)) ## 0

dput(names(beta_vali_final))
beta_vali_final$bat<-NULL
beta_vali_final$case<-NULL
beta_vali_final$id<-NULL
colnames(beta_vali_final)[1]<-'id'

beta_vali_final[,-c(1,2)]<-lapply(beta_vali_final[,-c(1,2)], as.numeric)

save(beta_vali_final, 
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/array2_processed.RData')

# save all the excluded duplicate ID
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/duplicatestoexclude.RData")
exclude<-append(exclude, exclud)

save(exclude, 
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/excludedduplicatesall.RData')

# stack three batachs all together 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/hipo_processed.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/array2_processed.RData")

beta_all<-bind_rows(beta_vali_final, hipo_final)
# further exclude duplicates
beta_all[duplicated(beta_all$id), ]$id 

beta_all<-subset(beta_all, !(beta_all$uid %in% excl)) ### 

save(beta_all, file="/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/allarray_processed.RData")

######################################## 2. Clinical variable description ##################################################################################
## 2.1 Data preparation ####
library(haven)
follow<-read_sas('/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/follow.sas7bdat')

save(follow, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/follow_raw.RData')

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/follow.RData")
# select variables
library(readr)
follow_variable <- read_csv("dachs_cohort/follow_variable.csv")
var<-follow_variable$var

df<-subset(follow, select = var)
summary(df)
df$id[duplicated(df$tn_id)] # no duplicate

# encode variables

df<-within.data.frame(df,{
  P01SEX<- factor(P01SEX,levels=c(1, 2), labels=c("Female", "Male"))
  
  chemradther<- as.factor(ifelse(chemradther=='Ja', 'Yes', 
                                 ifelse(chemradther == 'Nein', 'No', NA)))
  
  stagediag<- factor(stagediag,
                     levels=c(1,2,3,4),
                     labels=c("I", "II", "III", "IV"))
  
  Location = as.factor(ifelse(crc2sites=='rectum', 'Rectum', ifelse(crcprox == 1, 'Proximal colon', 'Distal colon')))
  recurr_what<-as.factor(recurr_what)
  
  crc2sites<-NULL
  crcprox<-NULL
  crcdist<-NULL
  
})

colnames(df)<-follow_variable$Names
colnames(df)[1]<-'id'

# drop the one without following up N = 1
df<-df[!is.na(df$timeD), ]
summary(df)

save(df, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/follow_orig.RData')

# molecular characteristics
library(haven)
mole<-read_sas('/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/tumchars_set12_20160802.sas7bdat')
save(mole, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/mole.RData')
write.csv(mole, '/home/t532n/molecular.csv')

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/mole_orig.RData")

mole<-subset(mole, select = c('tn_id', 'MSI_gent' ,'brafmut', 'krasmut', 'cimphi'))

# ecode varaibles

mole<-within.data.frame(mole,{
  MSI_gent<- factor(MSI_gent,levels=c(1, 0), labels=c("Yes", "No"))
  brafmut<- factor(brafmut,levels=c(1, 0), labels=c("Yes", "No"))
  krasmut<- factor(krasmut,levels=c(1, 0), labels=c("Yes", "No"))
  cimphi<- factor(cimphi,levels=c(1, 0), labels=c("Yes", "No"))
  
})

colnames(mole)[1]<-'id'

colnames(mole)[2:5]<-c('')

save(mole, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/mole_orig.RData')

# matach ID 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/follow_orig.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/mole_orig.RData")

library(dplyr)

intersect<-Reduce(intersect, list(beta_all$id, df$id))

covar<-df[df$id %in% intersect, ]
beta_all<-beta_all[beta_all$id %in% intersect, ]

save(beta_all,
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/allarray_matched.RData')

covar<-merge(covar, mole, by='id', all.x = T)
summary(covar)

library(naniar)
miss<-as.data.frame(miss_var_summary(covar))
summary(covar)
library(dplyr)
covar$MSI_gent<-recode(covar$MSI_gent, Yes = 'MSI', No = 'MSS')
covar$brafmut<-recode(covar$brafmut, Yes = 'BRAF mut', No = 'BRAF wt')
covar$krasmut<-recode(covar$krasmut, Yes = 'KRAS mut', No = 'KRAS wt')
covar$cimphi<-recode(covar$cimphi, Yes = 'CIMP-high', No = 'CIMP-low')
covar$timey<-covar$timeD/365.25
covar$recurr_timey<-covar$timeD_recurr/365.25

# consider competing risk
library(cmprsk)
summary(covar)
attach(covar)
covar$death_crccp<-ifelse(death_all==1&death_crc==1, 1, 
                          ifelse(death_all==1&death_crc==0, 2, 
                                 ifelse(is.na(death_crc), NA, 0)))

summary(covar$death_crccp)

covar$recurr_cp<-ifelse(death_all==1&recurr==1, 1, 
                        ifelse(death_all==1&recurr==0, 2, 
                               ifelse(is.na(recurr), NA, 0)))


# create PFS and corresponding time
covar$PFS<-ifelse(death_all==1 | recurr==1, 1,
                  ifelse(death_all==0&is.na(recurr), NA, 0))


# recurrence could happen before death from all causes 
covar$timeD_PFS<-ifelse(recurr==1 & !is.na(recurr), timeD_recurr, 
                        ifelse(death_all==1, timeD,
                               ifelse(death_all==0 & is.na(recurr), NA, timeD)))


covar$timey_PFS<-covar$timeD_PFS/365.25
save(covar, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/covariable_all.RData')

## 2.2 make table one ####
libray(tableone)
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/covariable_all.RData")

dput(names(covar))
cols<-c("Diagnosis_year", "Age", "Gender", "Stage_at_diagnosis", "chemradther", 
        "MSI_gent", "brafmut", "krasmut", "cimphi", "timey")

vars<-c("Diagnosis_year", "Age", "Gender", "Stage_at_diagnosis", "Location", "chemradther", 
        "MSI_gent", "brafmut", "krasmut", "cimphi", "timey")

nonNormalVars<-c("Diagnosis_year", "Age", "timey")

table1<- CreateTableOne(vars = vars,  data = covar, includeNA =F, test = T)

b<-print(table1, nonnormal = nonNormalVars,  catDigits=0,  contDigits=0, showAllLevels=T, missing=T, quote = TRUE, noSpaces=T )
b<-as.data.frame(b)
total<-data.frame(rownames(b), b)
total$rownames.b.<-substring(total$rownames.b., first = 3)

library(writexl)

write_xlsx(total, path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Table1_updated.xlsx", col_names = T)

## 2.3 outcome description ####

### 2.3.1 all cause mortality #####

# median follow-up 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/covariable_all.RData")

summary(covar)
library(prodlim)
quantile(prodlim(Hist(timey, death_all)~1, data = covar, reverse = T))

# mortality per 1000 person-years (95% CI)
library(survival)
summary(baselineend$time)

death<-pyears(Surv(timeD, death_all)~Gender, data=covar)

sum<-summary(death, header = T, n=T, event=T, pyears = T, rate=T, ci.r = T, totals = T, rr=T, ci.rr=T,
             scale = 1000, digits=4)

# cumulative mortality at 3, 5, and 10 years 

fit_development <- survfit(Surv(timey, death_all) ~ 1, data=covar)
survival_development <- 0
survival_development$time <- fit_development$time
survival_development$surv <- fit_development$surv

S_3 <- min(survival_development$surv[survival_development$time <= 3])
S_5 <- min(survival_development$surv[survival_development$time <= 5])
S_10 <- min(survival_development$surv[survival_development$time <= 10])


(1-S_3)*100  

(1-S_5)*100 

(1-S_10)*100 

### 2.3.2  disease free survival/Relapse free survival (DFS) #####
# median follow-up 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/covariable_all.RData")

quantile(prodlim(Hist(recurr_timey, recurr)~1, data = covar, reverse = T))

# mortality per 1000 person-years (95% CI)

death<-pyears(Surv(timeD_recurr_recurr, recurr)~Gender, data=covar)

sum<-summary(death, header = T, n=T, event=T, pyears = T, rate=T, ci.r = T, totals = T, rr=T, ci.rr=T,
             scale = 1000, digits=4)

# cumulative mortality at 5, 10, and 15 years 
recur<-covar[!is.na(covar$recurr), ]

fit_development <- cuminc(ftime = recur$recurr_timey, fstatus =  recur$recurr_cp, cencode = 0)

survival_development <- 0
survival_development$time <- fit_development$`1 1`$time
survival_development$surv <- fit_development$`1 1`$est

S_3 <- max(survival_development$surv[survival_development$time <= 3])
S_3*100
S_5 <- max(survival_development$surv[survival_development$time <= 5])
S_5*100
S_10 <- max(survival_development$surv[survival_development$time <= 10])
S_10*100


## 2.4. make complete case dataset, and compare the outcome with overall sample ####
### 2.4.1  all cause mortality #####
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/covariable_all.RData")
dput(names(covar))
OS_complete <- subset(covar, select = c("id", "Diagnosis_year", "Age", "Gender", "chemradther",
                                        "Stage_at_diagnosis", "timeD", "timeM", "death_all", 
                                        "Location", "MSI_gent", "brafmut", "krasmut", "cimphi", "timey"))

# drop rows with any missing

OS_complete <- OS_complete[complete.cases(OS_complete),]

save(OS_complete, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/OS_complete.RData')


# describe outcome ##

# median follow-up

quantile(prodlim(Hist(timey, death_all)~1, data = OS_complete, reverse = T))

# mortality per 1000 person-years (95% CI) 

death<-pyears(Surv(timeD, death_all)~Gender, data=OS_complete)

sum<-summary(death, header = T, n=T, event=T, pyears = T, rate=T, ci.r = T, totals = T, rr=T, ci.rr=T,
             scale = 1000, digits=4)

# cumulative mortality at 3, 5, 10, and 15 years 

fit_development <- survfit(Surv(timey, death_all) ~ 1, data=OS_complete)
survival_development <- 0
survival_development$time <- fit_development$time
survival_development$surv <- fit_development$surv

S_3 <- min(survival_development$surv[survival_development$time <= 3])
(1-S_3)*100 


S_5 <- min(survival_development$surv[survival_development$time <= 5])
S_10 <- min(survival_development$surv[survival_development$time <= 10])
S_15 <- min(survival_development$surv[survival_development$time <= 15])

(1-S_5)*100 

(1-S_10)*100  

(1-S_15)*100  

# compare K=M curves
# complete case
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/OS_complete.RData")
dput(names(OS_complete))
osc<-subset(OS_complete, select = c("id",  "death_all",  "timey"))
osc$Group<-'Complete cases (N = 1734)'
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/covariable_all.RData")
oso<-subset(covar, select = c("id",  "death_all",  "timey"))
oso$Group<-'Overall sample (N = 2310)'
library(dplyr)
oskmc<-bind_rows(osc, oso)

oskmc$Group<-as.factor(oskmc$Group)
oskmc$Group<-factor(oskmc$Group, levels = rev(levels(oskmc$Group)))

summary(oskmc)

# plot figure
library (survminer)
deathc<-survfit(Surv(timey, death_all)~Group, data=oskmc)

ggsurvplot(deathc, data = oskmc, fun="event",
           censor.size=0.5,  conf.int=T, palette ="jco",
           risk.table = T,pval.size=4, pval.method.size=4,
           xlab="Time (years)", ylab="Cumulative mortality (%)",
           y.offset= FALSE, x.offset= FALSE,xlim=c(0, 15), break.x.by=1,
           legend = c(0.35, 0.9),  pval.coord=c(11, 0),
           legend.title="Overall survival",  font.legend=c(15), 
           legend.labs= c("Overall sample (N = 2310)", "Complete cases (N = 1734)"),
           font.x=c(10),font.xlab=c(10),ylim=c(0, 0.8),
           risk.table.col = "strata",
           font.y=c(10), font.ylab=c(10), surv.scale=c('percent'))

### 2.4.2  DFS #####
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/covariable_all.RData")

dput(names(covar))

recur_complete <- subset(covar, select = c("id", "Diagnosis_year", "Age", "Gender", "chemradther",
                                           "Stage_at_diagnosis","recurr",
                                           "Location", "MSI_gent", "brafmut", "krasmut", "cimphi", 
                                           "timeD_recurr_recurr", "recurr_timey", "recurr_cp"))

# drop rows with any missing

recur_complete <- recur_complete[complete.cases(recur_complete),]

save(recur_complete, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/recur_complete.RData')

# describe outcome 

# median follow-up

quantile(prodlim(Hist(recurr_timey, recurr)~1, data = recur_complete, reverse = T))

# mortality per 1000 person-years (95% CI) 

death<-pyears(Surv(timeD_recurr_recurr, recurr)~Gender, data=recur_complete)

sum<-summary(death, header = T, n=T, event=T, pyears = T, rate=T, ci.r = T, totals = T, rr=T, ci.rr=T,
             scale = 1000, digits=4)

## cumulative mortality at 5, 10, and 15 years 

fit_development <- cuminc(ftime = recur_complete$recurr_timey,
                          fstatus =  recur_complete$recurr_cp, cencode = 0)


survival_development <- 0
survival_development$time <- fit_development$`1 1`$time
survival_development$surv <- fit_development$`1 1`$est

S_3 <- max(survival_development$surv[survival_development$time <= 3])
S_3*100
S_5 <- max(survival_development$surv[survival_development$time <= 5])
S_5*100
S_10 <- max(survival_development$surv[survival_development$time <= 10])
S_10*100

# compare cumulative incidence curves
library (survminer)
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/recur_complete.RData")
dput(names(recur_complete))
rsc<-subset(recur_complete, select = c("id",  "recurr_cp",  "recurr_timey" ))
rsc$Group<-'Complete'
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/covariable_all.RData")
rso<-subset(covar, select = c("id",  "recurr_cp",  "recurr_timey"))
rso<-rso[complete.cases(rso), ]
rso$Group<-'Overall'
library(dplyr)
rskmc<-bind_rows(rsc, rso)

rskmc$Group<-factor(rskmc$Group, levels = c("Overall", 
                                            "Complete"))

save(rskmc, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/recur_KMcompare.RData')

summary(rskmc)

# plot figure
recur<-cuminc(ftime = rskmc$recurr_timey, 
              fstatus = rskmc$recurr_cp, cencode = 0, 
              group=rskmc$Group)

ggcompetingrisks(recur,  xlab = "Follow up time (years)",  
                 ylab = "Cumulative incidence", 
                 title = "",
                 multiple_panels=F)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), limits = c(0,0.3))+
  scale_x_continuous(limits = c(0,10), breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9,10))+
  scale_linetype_discrete(name  ="",
                          breaks=c("Overall", "Complete" ),
                          labels=c("Overall sample (N = 2290)", "Complete cases (N = 1718)"))+
  scale_colour_discrete(name  ="",
                        breaks=c("1", "2"),
                        labels=c("Recurrence", 
                                 "Death without recurrence"))+
  theme(legend.position = c(0.82, 0.28), legend.box = "vertical")

######################################## 3. Prepare data for unicox ##################################################################################
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/cpg_validate.RData")
cpg<-cpg_unique
rm(cpg_orig)
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/covariable_all.RData")

dput(names(covar))

outcome<-subset(covar, select = c("id", "timeD", "timey", "timeD_recurr", "recurr_timey", "death_all", 
                                  "death_crc", "recurr", "death_crccp", "recurr_cp", "PFS", "timeD_PFS", 
                                  "timey_PFS"))

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/allarray_matched.RData")

unicox<-merge(outcome, beta_all, by='id')
summary(unicox)
unicox$uid<-NULL


## calculate score and quartiles

unicox<-within.data.frame(unicox, {
  Yang_2019_OS<-(0.12*cg02196655+1.35*cg03763616+0.73*cg03944089+0.73*cg06117855+0.76*cg07173760-3.96*cg07293947-0.76*cg07509155+0.58*cg09244244+0.4*cg10451565+0.28*cg12582008+1.99*cg13796218+3.6*cg20247048+1.34*cg21481775+0.42*cg23829949-0.28*cg23964386+0.96*cg24127989-0.45*cg24674703+0.84*cg24938727)
  Gong_2020_OS<-(0.09*cg14660573-0.04*cg09353563 + 2.06*cg00110724)
  WangX_2020_OS_PFS<-(0.919*cg03091331 + 0.963*cg06884352 + 0.703*cg07707546 + 0.721*cg08081805 + 0.587*cg21347353 + 0.528*cg25164589)
  WangY_2020_OS<-(38.52*cg00177496-4.13*cg01963906+ 2.574*cg05165940-79.32*cg12921795+2.31*cg19414598+6.061*cg25783173)
  Xiang_2020_OS<-(3.02*cg03017653) + (22.84*cg03977782) + (3.00*cg05417950) + (36.54*cg06250108) + (2.18*cg09893305) + (3.12*cg10414946) + (3.06*cg15170424) + (5.77*cg15639045) + (-9.51*cg15786837) + (22.02*cg17329249) + (8.12*cg21212956) + (-23.31*cg24206256) 
  Chen_2021_OS<- (1.87*cg11621464 + 1.11*cg13565656 -1.74*cg18976437 + 2.40*cg20505223 -1.97*cg20528583)
  Huang_2021_OS<-(0.754*cg21614638 + 1.230*cg21770617 + 0.701*cg12751565)
  Huang_2021_PFS<-(0.573*cg21614638 + 0.888*cg21770617 + 0.477*cg12751565)
  Li_2021_OS<-(-2.1*cg01408654) + (-2.2*cg04035209) + (-3*cg10196720) + (-5*cg10379890) + (2.3*cg11097433) + (4.86*cg14675211) + (2.93*cg15428578) + (-4.7*cg19343464) + (2.7*cg21384402) + (3.7*cg27404023)
  #Li_2021_PFS<- (-1.6363*cg10541864 -1.05*cg02774439 + 0.3498*cg06223767 -1.6098*cg19505136) 
  Gundert_2019_OS_DSS<-(0.18*cg23750514) + (0.15*cg01131395) + (0.26*cg12510999) + (0.19*cg18195165) + (0.3*cg05646575) + (0.24*cg11056055) + (0.21*cg10758824) + (0.31*cg16336556) + (0.22*cg19184885) + (0.12*cg19340296) + (0.22*cg22522598) + (0.29*cg17431888) + (0.15*cg18736676) + (0.23*cg08804626) + (0.19*cg08617020) + (0.31*cg14270346) + (0.07*cg16399624) + (0.22*cg14983135) + (0.24*cg00832644)
  
  TGundert_2019_OS_DSS<-as.factor(ntile(Gundert_2019_OS_DSS, 4))
  TLi_2021_OS<-as.factor(ntile(Li_2021_OS, 4)) 
  THuang_2021_PFS<-as.factor(ntile(Huang_2021_PFS, 4))
  THuang_2021_OS<-as.factor(ntile(Huang_2021_OS, 4))
  TChen_2021_OS<-as.factor(ntile(Chen_2021_OS, 4))
  TXiang_2020_OS<-as.factor(ntile(Xiang_2020_OS, 4))
  TWangY_2020_OS<-as.factor(ntile(WangY_2020_OS, 4))
  TWangX_2020_OS_PFS<-as.factor(ntile(WangX_2020_OS_PFS, 4))
  TGong_2020_OS<-as.factor(ntile(Gong_2020_OS, 4))
  TYang_2019_OS<-as.factor(ntile(Yang_2019_OS, 4))
  
})

#  normalization beta value to Z score

unicox[cpg]<-sapply(unicox[cpg], function(data) (data-mean(data))/sd(data))

save(unicox, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/unicox_all.RData")

# prepare the dataset without ProM classifier
library(readxl)
cpg <- read_excel("CpG_Data extraction.xlsx", 
                  sheet = "all CpGs", col_types = c("text", 
                                                    "skip", "skip", "text", "skip", "skip", 
                                                    "skip", "skip", "skip", "skip", "skip"))
cpgd<-cpg$CpGs[cpg$Study=='Gündert et al, 2019']

## subset needed columns
dput(names(unicox))
unicoxpm<-subset(unicox, select=c("id", "timeD", "timey", "timeD_recurr", "recurr_timey", "death_all", 
                                  "death_crc", "recurr", "death_crccp", "recurr_cp", "PFS", "timeD_PFS", 
                                  "timey_PFS", cpgd, "Gundert_2019_OS_DSS",  "TGundert_2019_OS_DSS"))


summary(unicoxpm)

# exclude the 572 patients used to develop ProMClassifier
library(readr)

proid <- read_csv("dachs_cohort/ProMprobe.csv")
proid[nrow(proid)+1, ]<-colnames(proid)
colnames(proid)<-'id'

inter<-Reduce(intersect, list(proid$id, unicox$id)) ## 572

unicoxpm<-unicoxpm[!(unicoxpm$id %in% proid$id ), ] # 1738

summary(unicoxpm)

save(unicoxpm,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/PMunicox.RData")

# describe distribution of risk score
# exclude PM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/unicox_all.RData")
cols<-c("Li_2021_OS", "Huang_2021_PFS", "Huang_2021_OS", "Chen_2021_OS", 
        "Xiang_2020_OS", "WangY_2020_OS", "WangX_2020_OS_PFS", "Gong_2020_OS", 
        "Yang_2019_OS")

vars<-c("Li_2021_OS", "Huang_2021_PFS", "Huang_2021_OS", "Chen_2021_OS", 
        "Xiang_2020_OS", "WangY_2020_OS", "WangX_2020_OS_PFS", "Gong_2020_OS", 
        "Yang_2019_OS")

nonNormalVars<-c("Li_2021_OS", "Huang_2021_PFS", "Huang_2021_OS", "Chen_2021_OS", 
                 "Xiang_2020_OS", "WangY_2020_OS", "WangX_2020_OS_PFS", "Gong_2020_OS", 
                 "Yang_2019_OS")

table1<- CreateTableOne(vars = vars,  data = unicox, includeNA =T, test = T)

b<-print(table1, nonnormal = nonNormalVars,  catDigits=1,  contDigits=1, showAllLevels=T, missing=T, quote = TRUE, noSpaces=T )
b<-as.data.frame(b)
total<-data.frame(rownames(b), b)
total$X.level.<-NULL
total<-total[-1, ]
total$X.Missing.<-NULL
colnames(total)<-c('Score', 'Median [IQR]')
library(stringr)
total$Score<-str_sub(total$Score, start = 2, end = -17)
total$Score[2]
rownames(total)<-(1:9)

# add the distribution of PM classifier #####
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/PMunicox.RData")
summary(unicoxpm)
cols<-c("Gundert_2019_OS_DSS")

vars<-c("Gundert_2019_OS_DSS")

nonNormalVars<-c("Gundert_2019_OS_DSS")

table1<- CreateTableOne(vars = vars,  data = unicoxpm, includeNA =T, test = T)

b<-print(table1, nonnormal = nonNormalVars,  catDigits=1,  contDigits=1, showAllLevels=T, missing=T, quote = TRUE, noSpaces=T )
b<-as.data.frame(b)
pm<-data.frame(rownames(b), b)
pm$X.level.<-NULL
pm<-pm[-1, ]
pm$X.Missing.<-NULL
colnames(pm)<-c('Score', 'Median [IQR]')
library(stringr)
pm$Score<-str_sub(pm$Score, start = 2, end = -17)
pm$Score[2]
rownames(pm)<-1

distr<-bind_rows(total, pm)

write_xlsx(distr, path = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_tablefigure/Scoredistribution.xlsx", col_names = T)





