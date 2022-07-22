####################### 1. Likelihood ratio test  #####################################################
## 1.1 OS ####
# exclude PM classifer

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicox_impu20list.RData")

scoret<-c("Huang_2021_OS", "Huang_2021_PFS", "WangY_2020_OS", "Gong_2020_OS")


multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+', x)))


P<-vector(20, mode="list")

for (i in 1:20){
  
  coxref<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu[[i]], y=TRUE, x = TRUE)
  
  
  cox1<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+Huang_2021_OS,
              data=multi_impu[[i]], y=TRUE, x = TRUE)
  
  cox2<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+Huang_2021_PFS,
              data=multi_impu[[i]], y=TRUE, x = TRUE)
  
  cox3<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+WangY_2020_OS,
              data=multi_impu[[i]], y=TRUE, x = TRUE)
  
  cox4<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+Gong_2020_OS,
              data=multi_impu[[i]], y=TRUE, x = TRUE)
  
  p1<-anova(coxref, cox1)$`P(>|Chi|)`[2]
  p2<-anova(coxref, cox2)$`P(>|Chi|)`[2]
  p3<-anova(coxref, cox3)$`P(>|Chi|)`[2]
  p4<-anova(coxref, cox4)$`P(>|Chi|)`[2]
  
  P[[i]]<-data.frame(p1, p2, p3, p4)
  
}


Pos<-ldply(P, rbind)
colnames(Pos)<-scoret

# add PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicoxPM_impu20list.RData")

P<-vector(20, mode="list")

for (i in 1:20){
  
  coxref<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis,
                data=multipm_impu[[i]], y=TRUE, x = TRUE)
  
  coxpm<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+Gundert_2019_OS_DSS,
               data=multipm_impu[[i]], y=TRUE, x = TRUE)
  
  Gundert_2019_OS_DSS<-anova(coxref, coxpm)$`P(>|Chi|)`[2]
  
  P[[i]]<-data.frame(Gundert_2019_OS_DSS)
  
}

ppm<-ldply(P, rbind)

Pos<-bind_cols(Pos, ppm)

## 1.2 DFS ####
# exclude PM classifer
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicox_impu20list.RData")
for (i in 1:20){
  multi_impu[[i]]<- subset( multi_impu[[i]], !is.na(recurr_cp))
}

scoret<-c("Huang_2021_OS", "Huang_2021_PFS", "WangY_2020_OS", "Gong_2020_OS")

multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+', x)))

P<-vector(20, mode="list")

for (i in 1:20){
  
  coxref<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu[[i]], y=TRUE, x = TRUE)
  
  cox1<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+Huang_2021_OS,
              data=multi_impu[[i]], y=TRUE, x = TRUE)
  
  cox2<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+Huang_2021_PFS,
              data=multi_impu[[i]], y=TRUE, x = TRUE)
  
  cox3<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+WangY_2020_OS,
              data=multi_impu[[i]], y=TRUE, x = TRUE)
  
  cox4<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+Gong_2020_OS,
              data=multi_impu[[i]], y=TRUE, x = TRUE)
  
  p1<-anova(coxref, cox1)$`P(>|Chi|)`[2]
  p2<-anova(coxref, cox2)$`P(>|Chi|)`[2]
  p3<-anova(coxref, cox3)$`P(>|Chi|)`[2]
  p4<-anova(coxref, cox4)$`P(>|Chi|)`[2]
  
  P[[i]]<-data.frame(p1, p2, p3, p4)
  
}


Pfs<-ldply(P, rbind)
colnames(Pdfs)<-scoret

### add PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicoxPM_impu20list.RData")

for (i in 1:20){
  multipm_impu[[i]]<- subset( multipm_impu[[i]], !is.na(recurr_cp))
}

P<-vector(20, mode="list")

for (i in 1:20){
  
  coxref<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis,
                data=multipm_impu[[i]], y=TRUE, x = TRUE)
  
  coxpm<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+Gundert_2019_OS_DSS,
               data=multipm_impu[[i]], y=TRUE, x = TRUE)
  
  Gundert_2019_OS_DSS<-anova(coxref, coxpm)$`P(>|Chi|)`[2]
  
  P[[i]]<-data.frame(Gundert_2019_OS_DSS)
  
}

ppm<-ldply(P, rbind)

Pdfs<-bind_cols(Pdfs, ppm)

####################### 2. AUC difference #####################################################
## 2.1 OS ####
### exclude PM classifier

## time-dependent AUC difference
scoret<-c( "Li_2021_OS", "Huang_2021_PFS", "Huang_2021_OS", "Chen_2021_OS", 
           "Xiang_2020_OS", "WangY_2020_OS", "WangX_2020_OS_PFS", "Gong_2020_OS", 
           "Yang_2019_OS")

multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+', x)))

# add reference
diff_AUC<-vector(20, mode="list")
diff_AUC_se<-vector(20, mode="list")

for (i in 1:20){
  
  
  multi_models <- lapply(multi_formulas, function(x){coxph(x, data = multi_impu[[i]], y=TRUE, x = TRUE)})
  
  coxref<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu[[i]], y=TRUE, x = TRUE)
  
  all_models<-c(list(coxref), multi_models)
  
  
  aucscore<-Score(all_models,
                  formula=Surv(timey, death_all)~1, data= multi_impu[[i]], conf.int=TRUE, metrics = 'AUC',
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

## overall AUC difference
diff_AUC<-vector(20, mode="list")
diff_AUC_se<-vector(20, mode="list")


for (i in 1:20){
  
  
  multi_models <- lapply(multi_formulas, function(x){coxph(x, data = multi_impu[[i]], y=TRUE, x = TRUE)})
  
  coxref<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu[[i]], y=TRUE, x = TRUE)
  
  all_models<-c(list(coxref), multi_models)
  
  
  aucscore<-Score(all_models,
                  formula=Surv(timey, death_all)~1, data= multi_impu[[i]], conf.int=TRUE, metrics = 'AUC'
  )
  
  model<-aucscore$AUC$score[, 1]
  times<-'Overall'
  
  contr<-subset(aucscore$AUC$contrasts, reference == 'coxph')
  modelc<-contr$model
  timec<-'Overall'
  ref<-contr$reference
  diff_AUC[[i]]<-contr$delta.AUC
  diff_AUC_se[[i]]<-contr$se
  
}


### differences in AUC
diff_AUC<-bind_cols(diff_AUC)
diff_AUC<-rowMeans(diff_AUC)

diff_AUC_se<-bind_cols(diff_AUC_se)
diff_AUC_se<-rowMeans(diff_AUC_se)

AUC_diff_all<-data.frame(timec, modelc, ref, diff_AUC, diff_AUC_se)

### add PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicoxPM_impu20list.RData")

diff_AUCpm<-vector(20, mode="list")
diff_AUCpm_se<-vector(20, mode="list")


for (i in 1:20){
  
  coxref<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis,
                data=multipm_impu[[1]], y=TRUE, x = TRUE)
  
  coxpm<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+Gundert_2019_OS_DSS,
               data=multipm_impu[[1]], y=TRUE, x = TRUE)
  
  
  AUCpmscore<-Score(list("ref"=coxref,"Gundert_2019_OS_DSS"=coxpm),
                    formula=Surv(timey, death_all)~1, data= multipm_impu[[1]], conf.int=TRUE, metrics = 'AUC',
                    time=c(1, 3, 5, 8))
  
  score<-subset(AUCpmscore$AUC$score, model == 'Gundert_2019_OS_DSS')
  
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
AUC<-bind_cols(AUCpm)
AUC<-rowMeans(AUC)
AUC_se<-bind_cols(AUCpm_se)
AUC_se<-rowMeans(AUC_se)
ll<-AUC-qnorm(0.975)*AUC_se
ul<-AUC+qnorm(0.975)*AUC_se

AUCpm_model<-data.frame(times, model, AUC, AUC_se, ll, ul)

diff_AUC<-bind_cols(diff_AUCpm)
diff_AUC<-rowMeans(diff_AUC)

diff_AUC_se<-bind_cols(diff_AUCpm_se)
diff_AUC_se<-rowMeans(diff_AUC_se)

AUCpm_diff<-data.frame(timec, modelc, ref, diff_AUC, diff_AUC_se)

# overall AUC 
diff_AUCpm<-vector(20, mode="list")
diff_AUCpm_se<-vector(20, mode="list")

for (i in 1:20){
  multipm_impu[[i]]<- subset( multipm_impu[[i]], !is.na(recurr_cp==1))
}

for (i in 1:20){
  
  coxref<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis,
                data=multipm_impu[[i]], y=TRUE, x = TRUE)
  
  coxpm<-coxph(Surv(timey, death_all)~Age+Gender+Stage_at_diagnosis+Gundert_2019_OS_DSS,
               data=multipm_impu[[i]], y=TRUE, x = TRUE)
  
  
  AUCpmscore<-Score(list("ref"=coxref,"Gundert_2019_OS_DSS"=coxpm),
                    formula=Surv(timey, death_all)~1, data= multipm_impu[[i]], conf.int=TRUE, metrics = 'AUC')
  
  score<-subset(AUCpmscore$AUC$score, model == 'Gundert_2019_OS_DSS')
  
  model<-score$model
  times<-'Overall'
  
  contr<-AUCpmscore$AUC$contrasts
  modelc<-contr$model
  timec<-'Overall'
  ref<-contr$reference
  diff_AUCpm[[i]]<-contr$delta.AUC
  diff_AUCpm_se[[i]]<-contr$se
  
}

### differences in AUCpm
diff_AUC<-bind_cols(diff_AUCpm)
diff_AUC<-rowMeans(diff_AUC)

diff_AUC_se<-bind_cols(diff_AUCpm_se)
diff_AUC_se<-rowMeans(diff_AUC_se)

AUCpm_diffall<-data.frame(timec, modelc, ref, diff_AUC, diff_AUC_se)

# bind results of AUC 

AUC_diff$timec<-as.character(AUC_diff$timec)
AUCpm_diff$timec<-as.character(AUCpm_diff$timec)

AUC_diff<-bind_rows(AUC_diff, AUC_diff_all, AUCpm_diff, AUCpm_diffall)

AUC_diff<-AUC_diff[order(AUC_diff$modelc), ]

AUC_diff$ll<-AUC_diff$diff_AUC-qnorm(0.975)*AUC_diff$diff_AUC_se
AUC_diff$ul<-AUC_diff$diff_AUC+qnorm(0.975)*AUC_diff$diff_AUC_se


AUC_diff$diffAUC_CI<-paste0(round(AUC_diff$diff_AUC, digits = 3), ' ', '(',
                            round(AUC_diff$ll, digits = 3), ',', ' ', round(AUC_diff$ul, digits = 3), ')' )

## 2.2 DFS ####
# exclude PM classifier ###
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicox_impu20list.RData")
options(scipen = 10000)
### AUC ###

### exclude PM classifier ###

scoret<-c( "Li_2021_OS", "Huang_2021_PFS", "Huang_2021_OS", "Chen_2021_OS", 
           "Xiang_2020_OS", "WangY_2020_OS", "WangX_2020_OS_PFS", "Gong_2020_OS", 
           "Yang_2019_OS")

multi_formulas <- sapply(scoret,
                         function(x) as.formula(paste('Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+', x)))

# add reference
diff_AUC<-vector(20, mode="list")
diff_AUC_se<-vector(20, mode="list")

for (i in 1:20){
  multi_impu[[i]]<- subset( multi_impu[[i]], !is.na(recurr_cp==1))
}

summary(multi_impu[[1]]$recurr_cp==1)

for (i in 1:20){
  
  
  multi_models <- lapply(multi_formulas, function(x){coxph(x, data = multi_impu[[i]], y=TRUE, x = TRUE)})
  
  coxref<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis,
                data=multi_impu[[i]], y=TRUE, x = TRUE)
  
  all_models<-c(list(coxref), multi_models)
  
  
  aucscore<-Score(all_models,
                  formula=Surv(recurr_timey, recurr_cp==1)~1, data= multi_impu[[i]], conf.int=TRUE, metrics = 'AUC',
                  time=c(1, 3, 5, 8),
                  cause = 1)
  
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

AUC_diff_all<-data.frame(timec, modelc, ref, diff_AUC, diff_AUC_se)

# add PM classifier 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicoxPM_impu20list.RData")

# time dependent AUC difference
diff_AUCpm<-vector(20, mode="list")
diff_AUCpm_se<-vector(20, mode="list")

for (i in 1:20){
  multipm_impu[[i]]<- subset( multipm_impu[[i]], !is.na(recurr_cp==1))
}

for (i in 1:20){
  
  coxref<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis,
                data=multipm_impu[[i]], y=TRUE, x = TRUE)
  
  coxpm<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+Gundert_2019_OS_DSS,
               data=multipm_impu[[i]], y=TRUE, x = TRUE)
  
  
  AUCpmscore<-Score(list("ref"=coxref,"Gundert_2019_OS_DSS"=coxpm),
                    formula=Surv(recurr_timey, recurr_cp==1)~1, data= multipm_impu[[i]], conf.int=TRUE, metrics = 'AUC',
                    time=c(1, 3, 5, 8), cause = 1)
  
  score<-subset(AUCpmscore$AUC$score, model == 'Gundert_2019_OS_DSS')
  
  model<-score$model
  times<-score$times

  contr<-AUCpmscore$AUC$contrasts
  modelc<-contr$model
  timec<-contr$times
  ref<-contr$reference
  diff_AUCpm[[i]]<-contr$delta.AUC
  diff_AUCpm_se[[i]]<-contr$se
  
}

### differences in AUCpm
diff_AUC<-bind_cols(diff_AUCpm)
diff_AUC<-rowMeans(diff_AUC)

diff_AUC_se<-bind_cols(diff_AUCpm_se)
diff_AUC_se<-rowMeans(diff_AUC_se)

AUCpm_diff<-data.frame(timec, modelc, ref, diff_AUC, diff_AUC_se)

# Overall AUC difference
diff_AUCpm<-vector(20, mode="list")
diff_AUCpm_se<-vector(20, mode="list")

for (i in 1:20){
  multipm_impu[[i]]<- subset( multipm_impu[[i]], !is.na(recurr_cp==1))
}

for (i in 1:20){
  
  coxref<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis,
                data=multipm_impu[[i]], y=TRUE, x = TRUE)
  
  coxpm<-coxph(Surv(recurr_timey, recurr_cp==1)~Age+Gender+Stage_at_diagnosis+Gundert_2019_OS_DSS,
               data=multipm_impu[[i]], y=TRUE, x = TRUE)
  
  
  AUCpmscore<-Score(list("ref"=coxref,"Gundert_2019_OS_DSS"=coxpm),
                    formula=Surv(recurr_timey, recurr_cp==1)~1, data= multipm_impu[[i]], conf.int=TRUE, metrics = 'AUC',
                    cause = 1)
  
  score<-subset(AUCpmscore$AUC$score, model == 'Gundert_2019_OS_DSS')
  
  model<-score$model
  times<-'Overall'

  contr<-AUCpmscore$AUC$contrasts
  modelc<-contr$model
  timec<-'Overall'
  ref<-contr$reference
  diff_AUCpm[[i]]<-contr$delta.AUC
  diff_AUCpm_se[[i]]<-contr$se
  
}

diff_AUC<-bind_cols(diff_AUCpm)
diff_AUC<-rowMeans(diff_AUC)

diff_AUC_se<-bind_cols(diff_AUCpm_se)
diff_AUC_se<-rowMeans(diff_AUC_se)

AUCpm_diffall<-data.frame(timec, modelc, ref, diff_AUC, diff_AUC_se)

###### bind results of AUC 
AUC_diff$timec<-as.character(AUC_diff$timec)
AUCpm_diff$timec<-as.character(AUCpm_diff$timec)

AUC_diff<-bind_rows(AUC_diff, AUC_diff_all, AUCpm_diff, AUCpm_diffall)

AUC_diff<-AUC_diff[order(AUC_diff$modelc), ]

AUC_diff$ll<-AUC_diff$diff_AUC-qnorm(0.975)*AUC_diff$diff_AUC_se
AUC_diff$ul<-AUC_diff$diff_AUC+qnorm(0.975)*AUC_diff$diff_AUC_se


AUC_diff$diffAUC_CI<-paste0(round(AUC_diff$diff_AUC, digits = 3), ' ', '(',
                            round(AUC_diff$ll, digits = 3), ',', ' ', round(AUC_diff$ul, digits = 3), ')' )

















