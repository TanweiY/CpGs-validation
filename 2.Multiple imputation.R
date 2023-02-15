######################################## 3. Multiple imputation ############################################################
# install mice 3.1.3
packlocfit<-file.choose()
install.packages(packlocfit, repos = NULL)

library(mice)
library(survival)
library(mitools)


## 3.1 imputation for clincial variables ####

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/covariable_all.RData")
dput(names(covar))
# calculate cumulative hazard###########change the outcome
clinimput<-subset(covar, select = c("id", "Diagnosis_year", "Age", "Gender", "chemradther", 
                                    "Stage_at_diagnosis","Location", "death_all", "timey"))
summary(clinimput)
HT1 <- summary(survival::survfit(Surv(timey, death_all)~1,data=clinimput))
clinimput$haz_os <-  approx(c(0, HT1$time), -log(c(1,HT1$surv)),xout=clinimput$time,method="constant",f=0,rule=2)$y

summary(clinimput)

# set outcome as factor
clinimput$death_all<- factor(clinimput$death_all)

save(clinimput, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/clinimpute_orig.RData')

# see all the default settings for imputation
impu_default <- mice(clinimput, maxit = 0)
summary(impu_default)

# see the predictor structure
pred <- quickpred(clinimput, exclude = c("id","timey"))
pred

# imputation method
meth <- impu_default$meth
meth
meth[6]<-'polr'    # stage is ordered

# multiple imputation for 20 times

clin_imputation_20 <- mice(clinimput, maxit = 10, m = 20, seed = 1234, pred = pred, meth = meth, print = TRUE)

save(clin_imputation_20, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/clinimputed.RData')

K <- 20

clin_imputated <- vector(K,mode="list")

for (i in 1:K) {
  clin_imputated[[i]] <- mice::complete(clin_imputation_20, i)
}

summary(clin_imputated[[1]])

save(clin_imputated, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/clinimputed20list.RData')

## 3.2 impute molecular characteristics ####

# selec the 1852 CpGs related with CIMP genes and use it for imputation
cpgall1852 <- read_csv("imputationcpg/cpgall1852.csv")
cpgimput<-cpgall1852$VAR1
save(cpgimput, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/imputationcpg/cpgimput.RData')

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/imputationcpg/cpgimput.RData")

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/myCombat.RData")
beta<-myCombat
rm(myCombat)
beta<-as.data.frame(beta)
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/imputationcpg/cpgimput.RData")

# select and clean
beta<-beta[cpgimput, ]

beta$cpg<-rownames(beta)
beta[nrow(beta)+1, ]<-colnames(beta)

library(data.table)
beta<-transpose(beta)
colnames(beta)<-beta[2681, ]

colnames(beta)[1853]<-'id'

library(dplyr)
beta <- beta %>%
  select(id, everything())
beta<-beta[-2681, ]

### select the id needed
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/allarray_matched.RData")
intersect<-Reduce(intersect, list(beta_all$uid, beta$id))

beta<-beta[beta$id %in% intersect, ]


save(beta, 
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/imputationcpg/beta_imputprocessed.RData')

library(mice)
library(survival)

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/covariable_all.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/imputationcpg/beta_imputprocessed.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation_processeddata/allarray_matched.RData")
colnames(beta)[1]<-'uid'
ID<-beta_all[,c(1,2)]
beta<-merge(beta, ID, by='uid')
beta <- beta %>%
  select(id, everything())
beta$uid<-NULL

beta[,-1]<-lapply(beta[,-1], as.numeric)

# calculate cumulative hazard
clinimput<-subset(covar, select = c("id", "Age", "Gender","Stage_at_diagnosis", "Location", "death_all", "timey",
                                    "MSI_gent", "brafmut", "krasmut", "cimphi"))

clinimput<-merge(clinimput, beta, by='id')

summary(clinimput)
HT1 <- summary(survival::survfit(Surv(timey, death_all)~1,data=clinimput))
clinimput$haz_os <-  approx(c(0, HT1$time), -log(c(1,HT1$surv)),xout=clinimput$time,method="constant",f=0,rule=2)$y

summary(clinimput)

## set outcome as factor
clinimput$death_all<- factor(clinimput$death_all)
moleimput<-clinimput
save(moleimput, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/moleimpute_data.RData')

# see all the default settings for imputation
impu_default <- mice(moleimput, maxit = 0)
summary(impu_default)

# see the predictor structure
pred <- quickpred(moleimput, exclude = c("id","timey"))
pred

# imputation method
meth <- impu_default$meth
meth
meth[4]<-'polr'    # stage is ordered
meth[8:11]<-'rf'

# 2 imputation

# multiple imputation for 20 times

mole_imputation_20 <- mice(moleimput, maxit = 10, m = 20, seed = 1234, pred = pred, meth = meth, print = TRUE)

save(mole_imputation_20, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/molemputed.RData')

K <- 20

mole_imputated <- vector(K,mode="list")

for (i in 1:K) {
  mole_imputated[[i]] <- mice::complete(mole_imputation_20, i)
}

summary(mole_imputated[[1]])

save(mole_imputated, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/moleimputed20list.RData')

# prepare the imputed data for downward analyses

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/clinimputed.RData")
nrow(clin_imputation_20[[1]])
dput(names(clin_imputation_20[[1]]))

clin_imputated <- vector(20,mode="list")

for (i in 1:20) {
  clin_imputated[[i]] <- mice::complete(clin_imputation_20, i)
}

save(clin_imputated, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/clinimputed20list.RData')

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/clinimputed20list.RData")
nrow(clin_imputated[[1]])

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/moleimputed20list.RData")
summary(mole_imputated[[1]])
nrow(mole_imputated[[1]])
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/unicox_all.RData")
summary(unicox)

n_impu <- 20
multi_impu <- vector(n_impu,mode="list")
for (i in 1:n_impu) {
  clin_imputated[[i]]$death_all<-NULL
  clin_imputated[[i]]$timey <-NULL
  clin_imputated[[i]]$haz_os <-NULL
  
  multi_impu[[i]]<-merge(clin_imputated[[i]], mole_imputated[[i]], by='id')
  multi_impu[[i]]<-merge(multi_impu[[i]], unicox, by = 'id')
  
}

save(multi_impu, file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicox_impu20list.RData')





