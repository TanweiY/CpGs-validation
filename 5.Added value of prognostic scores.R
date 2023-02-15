####################### 1. Likelihood ratio test  #####################################################
####################### 1.1 OS  #####################################################
# Colon and rectum, stage 1-4
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")

multi_impu_cr14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr14[[i]]<-merge(unicox_cr14, multi_impu[[i]], by='id')
}

dput(names(multi_impu_cr14[[1]]))

scoret<-c("Huang_2021_OS_CR14", "Huang_2021_DFS_CR14", "Gong_2020_CR14")


multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+', x)))


P<-vector(20, mode="list")

for (i in 1:20){
  
  coxref<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu_cr14[[i]], y=TRUE, x = TRUE)
  
  
  cox1<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+Huang_2021_OS_CR14,
              data=multi_impu_cr14[[i]], y=TRUE, x = TRUE)
  
  cox2<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+Huang_2021_DFS_CR14,
              data=multi_impu_cr14[[i]], y=TRUE, x = TRUE)
  
  cox3<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+Gong_2020_CR14,
              data=multi_impu_cr14[[i]], y=TRUE, x = TRUE)
  
  p1<-anova(coxref, cox1)$`P(>|Chi|)`[2]
  p2<-anova(coxref, cox2)$`P(>|Chi|)`[2]
  p3<-anova(coxref, cox3)$`P(>|Chi|)`[2]
  
  
  P[[i]]<-data.frame(p1, p2, p3)
  
}

Pos<-ldply(P, rbind)
colnames(Pos)<-scoret

poscr14<-Pos

save(poscr14, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/likelihoodp.RData")

# Colon and rectum, stage 1-3

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr13.RData")

multi_impu_cr13<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr13[[i]]<-merge(unicox_cr13, multi_impu[[i]], by='id')
}

dput(names(multi_impu_cr13[[1]]))

scoret<-c("Gundert_2019_CR13")


multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+', x)))


P<-vector(20, mode="list")

for (i in 1:20){
  
  coxref<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu_cr13[[i]], y=TRUE, x = TRUE)
  
  
  cox1<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+Gundert_2019_CR13,
              data=multi_impu_cr13[[i]], y=TRUE, x = TRUE)
  
  p1<-anova(coxref, cox1)$`P(>|Chi|)`[2]
  
  P[[i]]<-data.frame(p1)
  
}

Pos<-ldply(P, rbind)
colnames(Pos)<-scoret

poscr13<-Pos

resave(poscr13, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/likelihoodp.RData")

# Colon, stage 1-4

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c14.RData")

multi_impu_c14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_c14[[i]]<-merge(unicox_c14, multi_impu[[i]], by='id')
}

dput(names(multi_impu_c14[[1]]))

scoret<-c("WangY_2020_C14")


multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+', x)))


P<-vector(20, mode="list")

for (i in 1:20){
  
  coxref<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu_c14[[i]], y=TRUE, x = TRUE)
  
  
  cox1<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+WangY_2020_C14,
              data=multi_impu_c14[[i]], y=TRUE, x = TRUE)
  
  p1<-anova(coxref, cox1)$`P(>|Chi|)`[2]
  
  P[[i]]<-data.frame(p1)
  
}

Pos<-ldply(P, rbind)
colnames(Pos)<-scoret

posc14<-Pos

resave(posc14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/likelihoodp.RData")

## bind results together 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/likelihoodp.RData")

pos<-cbind(posc14, poscr13, poscr14)

resave(pos, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/likelihoodp.RData")

####################### 1.2 DFS  #####################################################
# Colon and rectum, stage 1-4
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")

multi_impu_cr14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr14[[i]]<-merge(unicox_cr14, multi_impu[[i]], by='id')
  multi_impu_cr14[[i]]<-multi_impu_cr14[[i]][!is.na(multi_impu_cr14[[i]]$recurr_cp), ]
  
}

dput(names(multi_impu_cr14[[1]]))

scoret<-c("Huang_2021_OS_CR14", "Huang_2021_DFS_CR14", "Gong_2020_CR14")


multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+', x)))


P<-vector(20, mode="list")

for (i in 1:20){
  
  coxref<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu_cr14[[i]], y=TRUE, x = TRUE)
  
  
  cox1<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+Huang_2021_OS_CR14,
              data=multi_impu_cr14[[i]], y=TRUE, x = TRUE)
  
  cox2<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+Huang_2021_DFS_CR14,
              data=multi_impu_cr14[[i]], y=TRUE, x = TRUE)
  
  cox3<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+Gong_2020_CR14,
              data=multi_impu_cr14[[i]], y=TRUE, x = TRUE)
  
  p1<-anova(coxref, cox1)$`P(>|Chi|)`[2]
  p2<-anova(coxref, cox2)$`P(>|Chi|)`[2]
  p3<-anova(coxref, cox3)$`P(>|Chi|)`[2]
  
  
  P[[i]]<-data.frame(p1, p2, p3)
  
}

pdfs<-ldply(P, rbind)
colnames(pdfs)<-scoret

pdfscr14<-pdfs

resave(pdfscr14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/likelihoodp.RData")

# Colon and rectum, stage 1-3
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr13.RData")

multi_impu_cr13<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr13[[i]]<-merge(unicox_cr13, multi_impu[[i]], by='id')
  multi_impu_cr13[[i]]<-multi_impu_cr13[[i]][!is.na(multi_impu_cr13[[i]]$recurr_cp), ]
  
}

dput(names(multi_impu_cr13[[1]]))

scoret<-c("Gundert_2019_CR13")


multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+', x)))


P<-vector(20, mode="list")

for (i in 1:20){
  
  coxref<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu_cr13[[i]], y=TRUE, x = TRUE)
  
  
  cox1<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+Gundert_2019_CR13,
              data=multi_impu_cr13[[i]], y=TRUE, x = TRUE)
  
  p1<-anova(coxref, cox1)$`P(>|Chi|)`[2]
  
  P[[i]]<-data.frame(p1)
  
}

pdfs<-ldply(P, rbind)
colnames(pdfs)<-scoret

pdfscr13<-pdfs

resave(pdfscr13, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/likelihoodp.RData")

# Colon, stage 1-4

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c14.RData")

multi_impu_c14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_c14[[i]]<-merge(unicox_c14, multi_impu[[i]], by='id')
  multi_impu_c14[[i]]<-multi_impu_c14[[i]][!is.na(multi_impu_c14[[i]]$recurr_cp), ]
  
}

dput(names(multi_impu_c14[[1]]))

scoret<-c("WangY_2020_C14")


multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+', x)))


P<-vector(20, mode="list")

for (i in 1:20){
  
  coxref<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu_c14[[i]], y=TRUE, x = TRUE)
  
  
  cox1<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+WangY_2020_C14,
              data=multi_impu_c14[[i]], y=TRUE, x = TRUE)
  
  p1<-anova(coxref, cox1)$`P(>|Chi|)`[2]
  
  P[[i]]<-data.frame(p1)
  
}

pdfs<-ldply(P, rbind)
colnames(pdfs)<-scoret

pdfsc14<-pdfs

resave(pdfsc14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/likelihoodp.RData")

# bind os results together 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/likelihoodp.RData")

pdfs<-cbind(pdfsc14, pdfscr13, pdfscr14)

resave(pdfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/likelihoodp.RData")

# bind two outcomes together 
pos$outcome<-'OS'

likelip<-cbind(pos, pdfs)


write.csv(likelip, 
          "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/tablesfigures/likelip.csv")

####################### 2. AUC difference  #####################################################
####################### 2.1 OS  ################################

# Colon and rectum, stage 1-4 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")

multi_impu_cr14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr14[[i]]<-merge(unicox_cr14, multi_impu[[i]], by='id')
}

dput(names(multi_impu_cr14[[1]]))

scoret<-c("Huang_2021_OS_CR14", "Huang_2021_DFS_CR14", "Gong_2020_CR14")


multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+', x)))


# add reference
diff_AUC<-vector(20, mode="list")
diff_AUC_se<-vector(20, mode="list")

for (i in 1:20){
  
  multi_models <- lapply(multi_formulas, function(x){coxph(x, data = multi_impu_cr14[[i]], y=TRUE, x = TRUE)})
  
  coxref<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu_cr14[[i]], y=TRUE, x = TRUE)
  
  all_models<-c(list(coxref), multi_models)
  
  
  aucscore<-Score(all_models,
                  formula=Surv(timey, death_all)~1, data= multi_impu_cr14[[i]], conf.int=TRUE, metrics = 'AUC',
                  time=c(1, 3, 5, 8))
  
  model<-aucscore$AUC$score[, 1]
  time<-aucscore$AUC$score[, 2]
  
  contr<-subset(aucscore$AUC$contrasts, reference == 'coxph')
  modelc<-contr$model
  timec<-contr$times
  ref<-contr$reference
  diff_AUC[[i]]<-contr$delta.AUC
  diff_AUC_se[[i]]<-contr$se
  
}

diff_AUC<-bind_cols(diff_AUC)
diff_AUC<-rowMeans(diff_AUC)

diff_AUC_se<-bind_cols(diff_AUC_se)
diff_AUC_se<-rowMeans(diff_AUC_se)

AUC_diff<-data.frame(timec, modelc, ref, diff_AUC, diff_AUC_se)

AUC_diff$ll<-AUC_diff$diff_AUC-qnorm(0.975)*AUC_diff$diff_AUC_se
AUC_diff$ul<-AUC_diff$diff_AUC+qnorm(0.975)*AUC_diff$diff_AUC_se


AUC_diff$diffAUC_CI<-paste0(round(AUC_diff$diff_AUC, digits = 3), ' ', '(',
                            round(AUC_diff$ll, digits = 3), ',', ' ', round(AUC_diff$ul, digits = 3), ')' )


AUCdiff_os_cr14<-AUC_diff
save(AUCdiff_os_cr14, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUCdiff.RData")

# Colon and rectum, stage 1-3

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr13.RData")

multi_impu_cr13<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr13[[i]]<-merge(unicox_cr13, multi_impu[[i]], by='id')
}

dput(names(multi_impu_cr13[[1]])) # Gundert_2019_CR13
summary(multi_impu_cr13[[1]]$timey)


diff_AUCpm<-vector(20, mode="list")
diff_AUCpm_se<-vector(20, mode="list")

for (i in 1:20){
  
  # otherwise error
  multi_impu_cr13[[i]]$Stage_at_diagnosis<-as.factor(as.character(multi_impu_cr13[[i]]$Stage_at_diagnosis))
  
  coxref<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu_cr13[[i]], y=TRUE, x = TRUE)
  
  coxpm<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+Gundert_2019_CR13,
               data=multi_impu_cr13[[i]], y=TRUE, x = TRUE)
  
  AUCpmscore<-Score(list("ref"=coxref,"Gundert_2019_CR13"=coxpm),
                    formula=Surv(timey, death_all)~1, data= multi_impu_cr13[[i]], conf.int=TRUE,
                    time=c(1, 3, 5, 8),
                    metrics = 'AUC')
  
  score<-subset(AUCpmscore$AUC$score, model == 'Gundert_2019_CR13')
  
  model<-score$model
  times<-score$times
  
  contr<-AUCpmscore$AUC$contrasts
  modelc<-contr$model
  timec<-contr$times
  ref<-contr$reference
  diff_AUCpm[[i]]<-contr$delta.AUC
  diff_AUCpm_se[[i]]<-contr$se
  
}

# individual AUCpm
diff_AUC<-bind_cols(diff_AUCpm)
diff_AUC<-rowMeans(diff_AUC)

diff_AUC_se<-bind_cols(diff_AUCpm_se)
diff_AUC_se<-rowMeans(diff_AUC_se)

AUC_diff<-data.frame(timec, modelc, ref, diff_AUC, diff_AUC_se)

AUC_diff$ll<-AUC_diff$diff_AUC-qnorm(0.975)*AUC_diff$diff_AUC_se
AUC_diff$ul<-AUC_diff$diff_AUC+qnorm(0.975)*AUC_diff$diff_AUC_se


AUC_diff$diffAUC_CI<-paste0(round(AUC_diff$diff_AUC, digits = 3), ' ', '(',
                            round(AUC_diff$ll, digits = 3), ',', ' ', round(AUC_diff$ul, digits = 3), ')' )


AUCdiff_os_cr13<-AUC_diff
resave(AUCdiff_os_cr13, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUCdiff.RData")

# Colon, stage 1-4 

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c14.RData")

multi_impu_c14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_c14[[i]]<-merge(unicox_c14, multi_impu[[i]], by='id')
}

dput(names(multi_impu_c14[[1]])) 
summary(multi_impu_c14[[1]]$Stage_at_diagnosis)


diff_AUCpm<-vector(20, mode="list")
diff_AUCpm_se<-vector(20, mode="list")

for (i in 1:20){
  
  # otherwise error
  multi_impu_c14[[i]]$Stage_at_diagnosis<-as.factor(as.character(multi_impu_c14[[i]]$Stage_at_diagnosis))
  
  coxref<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu_c14[[i]], y=TRUE, x = TRUE)
  
  coxpm<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+WangY_2020_C14,
               data=multi_impu_c14[[i]], y=TRUE, x = TRUE)
  
  AUCpmscore<-Score(list("ref"=coxref,"WangY_2020_C14"=coxpm),
                    formula=Surv(timey, death_all)~1, data= multi_impu_c14[[i]], conf.int=TRUE,
                    time=c(1, 3, 5, 8),
                    metrics = 'AUC')
  
  score<-subset(AUCpmscore$AUC$score, model == 'WangY_2020_C14')
  
  model<-score$model
  times<-score$times
  
  contr<-AUCpmscore$AUC$contrasts
  modelc<-contr$model
  timec<-contr$times
  ref<-contr$reference
  diff_AUCpm[[i]]<-contr$delta.AUC
  diff_AUCpm_se[[i]]<-contr$se
  
}

### individual AUCpm
diff_AUC<-bind_cols(diff_AUCpm)
diff_AUC<-rowMeans(diff_AUC)

diff_AUC_se<-bind_cols(diff_AUCpm_se)
diff_AUC_se<-rowMeans(diff_AUC_se)

AUC_diff<-data.frame(timec, modelc, ref, diff_AUC, diff_AUC_se)

AUC_diff$ll<-AUC_diff$diff_AUC-qnorm(0.975)*AUC_diff$diff_AUC_se
AUC_diff$ul<-AUC_diff$diff_AUC+qnorm(0.975)*AUC_diff$diff_AUC_se


AUC_diff$diffAUC_CI<-paste0(round(AUC_diff$diff_AUC, digits = 3), ' ', '(',
                            round(AUC_diff$ll, digits = 3), ',', ' ', round(AUC_diff$ul, digits = 3), ')' )


AUCdiff_os_c14<-AUC_diff
resave(AUCdiff_os_c14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUCdiff.RData")

# bind all results together 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUCdiff.RData")
AUC_diff_os<-rbind(AUCdiff_os_c14, AUCdiff_os_cr13, AUCdiff_os_cr14)

AUC_diff_os<-AUC_diff_os[order(AUC_diff_os$modelc), ]

resave(AUC_diff_os, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUCdiff.RData")

############################## 2.2 DFS ##################
# Colon and rectum, stage 1-4
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")

multi_impu_cr14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr14[[i]]<-merge(unicox_cr14, multi_impu[[i]], by='id')
  multi_impu_cr14[[i]]<-multi_impu_cr14[[i]][!is.na(multi_impu_cr14[[i]]$recurr_cp), ]
  
}

dput(names(multi_impu_cr14[[1]]))

scoret<-c("Huang_2021_OS_CR14", "Huang_2021_DFS_CR14", "Gong_2020_CR14")


multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+', x)))


# add reference
diff_AUC<-vector(20, mode="list")
diff_AUC_se<-vector(20, mode="list")

for (i in 1:20){
  
  multi_models <- lapply(multi_formulas, function(x){coxph(x, data = multi_impu_cr14[[i]], y=TRUE, x = TRUE)})
  
  coxref<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu_cr14[[i]], y=TRUE, x = TRUE)
  
  all_models<-c(list(coxref), multi_models)
  
  
  aucscore<-Score(all_models,
                  formula=Surv(recurr_timey, recurr_cp==1)~1, data= multi_impu_cr14[[i]], conf.int=TRUE, metrics = 'AUC',
                  time=c(1, 3, 5, 8))
  
  model<-aucscore$AUC$score[, 1]
  time<-aucscore$AUC$score[, 2]
  
  contr<-subset(aucscore$AUC$contrasts, reference == 'coxph')
  modelc<-contr$model
  timec<-contr$times
  ref<-contr$reference
  diff_AUC[[i]]<-contr$delta.AUC
  diff_AUC_se[[i]]<-contr$se
  
}

diff_AUC<-bind_cols(diff_AUC)
diff_AUC<-rowMeans(diff_AUC)

diff_AUC_se<-bind_cols(diff_AUC_se)
diff_AUC_se<-rowMeans(diff_AUC_se)

AUC_diff<-data.frame(timec, modelc, ref, diff_AUC, diff_AUC_se)

AUC_diff$ll<-AUC_diff$diff_AUC-qnorm(0.975)*AUC_diff$diff_AUC_se
AUC_diff$ul<-AUC_diff$diff_AUC+qnorm(0.975)*AUC_diff$diff_AUC_se


AUC_diff$diffAUC_CI<-paste0(round(AUC_diff$diff_AUC, digits = 3), ' ', '(',
                            round(AUC_diff$ll, digits = 3), ',', ' ', round(AUC_diff$ul, digits = 3), ')' )


AUCdiff_dfs_cr14<-AUC_diff
resave(AUCdiff_dfs_cr14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUCdiff.RData")

# Colon and rectum, stage 1-3
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr13.RData")

multi_impu_cr13<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr13[[i]]<-merge(unicox_cr13, multi_impu[[i]], by='id')
  multi_impu_cr13[[i]]<-multi_impu_cr13[[i]][!is.na(multi_impu_cr13[[i]]$recurr_cp), ]
  
}

dput(names(multi_impu_cr13[[1]])) # Gundert_2019_CR13
summary(multi_impu_cr13[[1]]$timey)


diff_AUCpm<-vector(20, mode="list")
diff_AUCpm_se<-vector(20, mode="list")

for (i in 1:20){
  
  # otherwise error
  multi_impu_cr13[[i]]$Stage_at_diagnosis<-as.factor(as.character(multi_impu_cr13[[i]]$Stage_at_diagnosis))
  
  coxref<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu_cr13[[i]], y=TRUE, x = TRUE)
  
  coxpm<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+Gundert_2019_CR13,
               data=multi_impu_cr13[[i]], y=TRUE, x = TRUE)
  
  AUCpmscore<-Score(list("ref"=coxref,"Gundert_2019_CR13"=coxpm),
                    formula=Surv(recurr_timey, recurr_cp==1)~1, data= multi_impu_cr13[[i]], conf.int=TRUE,
                    time=c(1, 3, 5, 8),
                    metrics = 'AUC')
  
  score<-subset(AUCpmscore$AUC$score, model == 'Gundert_2019_CR13')
  
  model<-score$model
  times<-score$times
  
  contr<-AUCpmscore$AUC$contrasts
  modelc<-contr$model
  timec<-contr$times
  ref<-contr$reference
  diff_AUCpm[[i]]<-contr$delta.AUC
  diff_AUCpm_se[[i]]<-contr$se
  
}

### individual AUCpm
diff_AUC<-bind_cols(diff_AUCpm)
diff_AUC<-rowMeans(diff_AUC)

diff_AUC_se<-bind_cols(diff_AUCpm_se)
diff_AUC_se<-rowMeans(diff_AUC_se)

AUC_diff<-data.frame(timec, modelc, ref, diff_AUC, diff_AUC_se)

AUC_diff$ll<-AUC_diff$diff_AUC-qnorm(0.975)*AUC_diff$diff_AUC_se
AUC_diff$ul<-AUC_diff$diff_AUC+qnorm(0.975)*AUC_diff$diff_AUC_se


AUC_diff$diffAUC_CI<-paste0(round(AUC_diff$diff_AUC, digits = 3), ' ', '(',
                            round(AUC_diff$ll, digits = 3), ',', ' ', round(AUC_diff$ul, digits = 3), ')' )


AUCdiff_dfs_cr13<-AUC_diff
resave(AUCdiff_dfs_cr13, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUCdiff.RData")

# Colon, stage 1-4 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_c14.RData")

multi_impu_c14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_c14[[i]]<-merge(unicox_c14, multi_impu[[i]], by='id')
  multi_impu_c14[[i]]<-multi_impu_c14[[i]][!is.na(multi_impu_c14[[i]]$recurr_cp), ]
  
}

dput(names(multi_impu_c14[[1]])) 
summary(multi_impu_c14[[1]]$Stage_at_diagnosis)


diff_AUCpm<-vector(20, mode="list")
diff_AUCpm_se<-vector(20, mode="list")

for (i in 1:20){
  
  # otherwise error
  #multi_impu_c14[[i]]$Stage_at_diagnosis<-as.factor(as.character(multi_impu_c14[[i]]$Stage_at_diagnosis))
  
  coxref<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu_c14[[i]], y=TRUE, x = TRUE)
  
  coxpm<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+WangY_2020_C14,
               data=multi_impu_c14[[i]], y=TRUE, x = TRUE)
  
  AUCpmscore<-Score(list("ref"=coxref,"WangY_2020_C14"=coxpm),
                    formula=Surv(recurr_timey, recurr_cp==1)~1, data= multi_impu_c14[[i]], conf.int=TRUE,
                    time=c(1, 3, 5, 8),
                    metrics = 'AUC')
  
  score<-subset(AUCpmscore$AUC$score, model == 'WangY_2020_C14')
  
  model<-score$model
  times<-score$times
  
  contr<-AUCpmscore$AUC$contrasts
  modelc<-contr$model
  timec<-contr$times
  ref<-contr$reference
  diff_AUCpm[[i]]<-contr$delta.AUC
  diff_AUCpm_se[[i]]<-contr$se
  
}

### individual AUCpm
diff_AUC<-bind_cols(diff_AUCpm)
diff_AUC<-rowMeans(diff_AUC)

diff_AUC_se<-bind_cols(diff_AUCpm_se)
diff_AUC_se<-rowMeans(diff_AUC_se)

AUC_diff<-data.frame(timec, modelc, ref, diff_AUC, diff_AUC_se)

AUC_diff$ll<-AUC_diff$diff_AUC-qnorm(0.975)*AUC_diff$diff_AUC_se
AUC_diff$ul<-AUC_diff$diff_AUC+qnorm(0.975)*AUC_diff$diff_AUC_se


AUC_diff$diffAUC_CI<-paste0(round(AUC_diff$diff_AUC, digits = 3), ' ', '(',
                            round(AUC_diff$ll, digits = 3), ',', ' ', round(AUC_diff$ul, digits = 3), ')' )


AUCdiff_dfs_c14<-AUC_diff
resave(AUCdiff_dfs_c14, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUCdiff.RData")

## bind all results together 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUCdiff.RData")
AUC_diff_dfs<-rbind(AUCdiff_dfs_c14, AUCdiff_dfs_cr13, AUCdiff_dfs_cr14)

AUC_diff_dfs<-AUC_diff_dfs[order(AUC_diff_dfs$modelc), ]

resave(AUC_diff_dfs, 
       file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUCdiff.RData")

# make the tables 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/AUCdiff.RData")
dput(names(AUC_diff_os))
os<-select(AUC_diff_os, select = c("modelc", "timec", "diffAUC_CI"))
colnames(os)<-c('Study', 'Times (year)', 'OS')

dfs<-select(AUC_diff_dfs, select = c("modelc", "timec", "diffAUC_CI"))
colnames(dfs)<-c('Study', 'Times (year)', 'DFS')

AUCdiff<-merge(os, dfs, by =c('Study', 'Times (year)'))

write.csv(AUCdiff, 
          '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/tablesfigures/AUCdiff.csv')












