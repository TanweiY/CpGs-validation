######################################## 1. Dose response curves  ##################################################################################
## 1.1 prepare data for dose response curves ####
#  model1
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/unicox_all.RData")
# determine the best number of nodes
dd <- datadist(unicox)
options(datadist="dd")
for (i in 3:7) {
  fit <- cph(Surv(timeD, death_all)~rcs(Huang_2021_OS)+Age+Gender+chemradther+Stage_at_diagnosis+Location,
             data=multi_impu[[1]], x=TRUE)
  tmp <- extractAIC(fit)
  if(i == 3) {AIC = tmp[2]; nk = 3}
  if(tmp[2] < AIC) {AIC = tmp[2]; nk = i} 
}

print(nk) ## 3


### make the dataset for multicox
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/unicox_all.RData")

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicox_impu20list.RData")

n<-20

x<-vector(n, mode="list")
y<-vector(n, mode="list")
lc<-vector(n, mode="list")
uc<-vector(n, mode="list")

for (i in 1: n){
  
  dd <- datadist(multi_impu[[i]])
  options(datadist="dd")
  
  fit <- cph(Surv(timeD, death_all)~ rcs(Huang_2021_OS, 3)+Age+Gender+chemradther+Stage_at_diagnosis+Location,
             data=multi_impu[[i]], x=TRUE) #
  
  surv2<-Predict(fit, name='Huang_2021_OS', ref.zero=TRUE, fun = exp)
  
  x[[i]] = surv2$Huang_2021_OS
  y[[i]] = surv2$yhat
  lc[[i]] = surv2$lower
  uc[[i]] = surv2$uppe
}

x<-colMeans(ldply(x, rbind))
y<-colMeans(ldply(y, rbind))
lc<-colMeans(ldply(lc, rbind))
uc<-colMeans(ldply(uc, rbind))

multidata1<-data.frame(x, y, lc, uc)

# model2
# determine the best number of nodes

for (i in 3:7) {
  fit <- cph(Surv(timeD, death_all)~rcs(Huang_2021_PFS)+Age+Gender+chemradther+Stage_at_diagnosis+Location,
             data=multi_impu[[1]], x=TRUE)
  tmp <- extractAIC(fit)
  if(i == 3) {AIC = tmp[2]; nk = 3}
  if(tmp[2] < AIC) {AIC = tmp[2]; nk = i} 
}

print(nk) ## 3

### make the dataset for multicox
n<-20

x<-vector(n, mode="list")
y<-vector(n, mode="list")
lc<-vector(n, mode="list")
uc<-vector(n, mode="list")

for (i in 1: n){
  
  dd <- datadist(multi_impu[[i]])
  options(datadist="dd")
  
  fit <- cph(Surv(timeD, death_all)~ rcs(Huang_2021_PFS, 3)+Age+Gender+chemradther+Stage_at_diagnosis+Location,
             data=multi_impu[[i]], x=TRUE) #
  
  surv2<-Predict(fit, name='Huang_2021_PFS', ref.zero=TRUE, fun = exp)
  
  x[[i]] = surv2$Huang_2021_PFS
  y[[i]] = surv2$yhat
  lc[[i]] = surv2$lower
  uc[[i]] = surv2$uppe
}

x<-colMeans(ldply(x, rbind))
y<-colMeans(ldply(y, rbind))
lc<-colMeans(ldply(lc, rbind))
uc<-colMeans(ldply(uc, rbind))

multidata2<-data.frame(x, y, lc, uc)

# model3
# determine the best number of nodes

for (i in 3:7) {
  fit <- cph(Surv(timeD, death_all)~rcs(WangY_2020_OS)+Age+Gender+chemradther+Stage_at_diagnosis+Location,
             data=multi_impu[[1]], x=TRUE)
  tmp <- extractAIC(fit)
  if(i == 3) {AIC = tmp[2]; nk = 3}
  if(tmp[2] < AIC) {AIC = tmp[2]; nk = i} 
}

print(nk) ## 3

### make the dataset for multicox
n<-20

x<-vector(n, mode="list")
y<-vector(n, mode="list")
lc<-vector(n, mode="list")
uc<-vector(n, mode="list")

for (i in 1: n){
  
  dd <- datadist(multi_impu[[i]])
  options(datadist="dd")
  
  fit <- cph(Surv(timeD, death_all)~ rcs(WangY_2020_OS, 3)+Age+Gender+chemradther+Stage_at_diagnosis+Location,
             data=multi_impu[[i]], x=TRUE) #
  
  surv2<-Predict(fit, name='WangY_2020_OS', ref.zero=TRUE, fun = exp)
  
  x[[i]] = surv2$WangY_2020_OS
  y[[i]] = surv2$yhat
  lc[[i]] = surv2$lower
  uc[[i]] = surv2$uppe
}

x<-colMeans(ldply(x, rbind))
y<-colMeans(ldply(y, rbind))
lc<-colMeans(ldply(lc, rbind))
uc<-colMeans(ldply(uc, rbind))

multidata3<-data.frame(x, y, lc, uc)

save(multidata3,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/dp_WangY_2020_OS.RData")
# model4
# determine the best number of nodes

for (i in 3:7) {
  fit <- cph(Surv(timeD, death_all)~rcs(Gong_2020_OS)+Age+Gender+chemradther+Stage_at_diagnosis+Location,
             data=multi_impu[[1]], x=TRUE)
  tmp <- extractAIC(fit)
  if(i == 3) {AIC = tmp[2]; nk = 3}
  if(tmp[2] < AIC) {AIC = tmp[2]; nk = i} 
}

print(nk) ## 3

### make the dataset for multicox
n<-20

x<-vector(n, mode="list")
y<-vector(n, mode="list")
lc<-vector(n, mode="list")
uc<-vector(n, mode="list")

for (i in 1: n){
  
  dd <- datadist(multi_impu[[i]])
  options(datadist="dd")
  
  fit <- cph(Surv(timeD, death_all)~ rcs(Gong_2020_OS, 3)+Age+Gender+chemradther+Stage_at_diagnosis+Location,
             data=multi_impu[[i]], x=TRUE) #
  
  surv2<-Predict(fit, name='Gong_2020_OS', ref.zero=TRUE, fun = exp)
  
  x[[i]] = surv2$Gong_2020_OS
  y[[i]] = surv2$yhat
  lc[[i]] = surv2$lower
  uc[[i]] = surv2$uppe
}

x<-colMeans(ldply(x, rbind))
y<-colMeans(ldply(y, rbind))
lc<-colMeans(ldply(lc, rbind))
uc<-colMeans(ldply(uc, rbind))

multidata4<-data.frame(x, y, lc, uc)

# model5 - PM classifier
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/PMunicox.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/multicoxPM_impu20list.RData")

# determine the best number of nodes

for (i in 3:7) {
  fit <- cph(Surv(timeD, death_all)~rcs(Gundert_2019_OS_DSS)+Age+Gender+chemradther+Stage_at_diagnosis+Location,
             data=multipm_impu[[1]], x=TRUE)
  tmp <- extractAIC(fit)
  if(i == 3) {AIC = tmp[2]; nk = 3}
  if(tmp[2] < AIC) {AIC = tmp[2]; nk = i} 
}

print(nk) ## 3

### make the dataset for multicox
n<-20

x<-vector(n, mode="list")
y<-vector(n, mode="list")
lc<-vector(n, mode="list")
uc<-vector(n, mode="list")

for (i in 1: n){
  
  dd <- datadist(multipm_impu[[i]])
  options(datadist="dd")
  
  fit <- cph(Surv(timeD, death_all)~ rcs(Gundert_2019_OS_DSS, 3)+Age+Gender+chemradther+Stage_at_diagnosis+Location,
             data=multipm_impu[[i]], x=TRUE) #
  
  surv2<-Predict(fit, name='Gundert_2019_OS_DSS', ref.zero=TRUE, fun = exp)
  
  x[[i]] = surv2$Gundert_2019_OS_DSS
  y[[i]] = surv2$yhat
  lc[[i]] = surv2$lower
  uc[[i]] = surv2$uppe
}

x<-colMeans(ldply(x, rbind))
y<-colMeans(ldply(y, rbind))
lc<-colMeans(ldply(lc, rbind))
uc<-colMeans(ldply(uc, rbind))

multidata5<-data.frame(x, y, lc, uc)

save(multidata5,
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/dp_Gundert_2019_OS_DSS.RData")

## 1.2 plot dose response curves ####
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/dachs_cohort/unicox_all.RData")

load("~/Desktop/Projects/PhD projects/Biomarkers External validation /Data analysis/biomarker validation/dp_Gong_2020_OS.RData")
t<-quantile(unicox$Gong_2020_OS, probs = seq(0, 1, 1/4))[c(2, 3, 4)]
t
## plot figure
os1<-ggplot(multidata4, aes(x = x, y = y))+
  geom_ribbon(aes(ymin=lc, ymax=uc), alpha=0.2, fill = "steelblue3")+
  geom_line(colour ="steelblue3")+
  labs(x="Score", y = "Hazard ratio", title='Gong et al, 2020 (OS)')+
  scale_fill_manual(values=c("steelblue3", "darkred"))+
  scale_color_manual(values=c("steelblue3", "darkred"))+
  ylim(0,2.5)+
  theme(axis.text = element_text(size = 12),
        axis.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(color="black", size=15, face="bold", hjust = 0.5),
        # control legend
        legend.title=element_blank(),
        legend.position="top",
        legend.spacing.x = unit(1.2, 'cm'),
        legend.text=element_text(size=11),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm") ))

load("~/Desktop/Projects/PhD projects/Biomarkers External validation /Data analysis/biomarker validation/dp_Gundert_2019_OS_DSS.RData")
t<-quantile(unicoxpm$Gundert_2019_OS_DSS, probs = seq(0, 1, 1/4))[c(2, 3, 4)]
t

os2<-ggplot(multidata5, aes(x = x, y = y))+
  geom_ribbon(aes(ymin=lc, ymax=uc), alpha=0.2, fill = "steelblue3")+
  geom_line(colour ="steelblue3")+
  labs(x="Score", y = "Hazard ratio", title='Gündert et al, 2019 (OS)')+
  scale_fill_manual(values=c("steelblue3", "darkred"))+
  scale_color_manual(values=c("steelblue3", "darkred"))+
  ylim(0,2.5)+
  xlim(-7, 5)+
  theme(axis.text = element_text(size = 12),
        axis.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(color="black", size=15, face="bold", hjust = 0.5),
        # control legend
        legend.title=element_blank(),
        legend.position="top",
        legend.spacing.x = unit(1.2, 'cm'),
        legend.text=element_text(size=11),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm") ))


load("~/Desktop/Projects/PhD projects/Biomarkers External validation /Data analysis/biomarker validation/dp_Huang_2021_OS.RData")
t<-quantile(unicox$Huang_2021_OS, probs = seq(0, 1, 1/4))[c(2, 3, 4)]
t
os3<-ggplot(multidata1, aes(x = x, y = y))+
  geom_ribbon(aes(ymin=lc, ymax=uc), alpha=0.2, fill = "steelblue3")+
  geom_line(colour ="steelblue3")+
  labs(x="Score", y = "Hazard ratio", title='Huang et al, 2021 (OS)')+
  ylim(0,2.5)+
  theme(axis.text = element_text(size = 12),
        axis.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(color="black", size=15, face="bold", hjust = 0.5),
        # control legend
        legend.title=element_blank(),
        legend.position="top",
        legend.spacing.x = unit(1.2, 'cm'),
        legend.text=element_text(size=11),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm") ))


load("~/Desktop/Projects/PhD projects/Biomarkers External validation /Data analysis/biomarker validation/dp_Huang_2021_DFS.RData")

t<-quantile(unicox$Huang_2021_PFS, probs = seq(0, 1, 1/4))[c(2, 3, 4)]
t
## plot figure
os4<-ggplot(multidata2, aes(x = x, y = y))+
  geom_ribbon(aes(ymin=lc, ymax=uc), alpha=0.2, fill = "steelblue3")+
  geom_line(colour ="steelblue3")+
  labs(x="Score", y = "Hazard ratio", title='Huang et al, 2021 (DFS)')+
  scale_fill_manual(values=c("steelblue3", "darkred"))+
  scale_color_manual(values=c("steelblue3", "darkred"))+
  ylim(0,2.5)+
  theme(axis.text = element_text(size = 12),
        axis.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(color="black", size=15, face="bold", hjust = 0.5),
        # control legend
        legend.title=element_blank(),
        legend.position="top",
        legend.spacing.x = unit(1.2, 'cm'),
        legend.text=element_text(size=11),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm") ))


## fg5
load("~/Desktop/Projects/PhD projects/Biomarkers External validation /Data analysis/biomarker validation/dp_WangY_2020_OS.RData")
t<-quantile(unicox$WangY_2020_OS, probs = seq(0, 1, 1/4))[c(2, 3, 4)]
t

## plot figure
os5<-ggplot(multidata3, aes(x = x, y = y))+
  geom_ribbon(aes(ymin=lc, ymax=uc), alpha=0.2, fill = "steelblue3")+
  geom_line(colour ="steelblue3")+
  labs(x="Score", y = "Hazard ratio", title='Wang Y et al, 2020 (OS)')+
  scale_fill_manual(values=c("steelblue3", "darkred"))+
  scale_color_manual(values=c("steelblue3", "darkred"))+
  ylim(0,2.5)+
  theme(axis.text = element_text(size = 12),
        axis.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(color="black", size=15, face="bold", hjust = 0.5),
        # control legend
        legend.title=element_blank(),
        legend.position="top",
        legend.spacing.x = unit(1.2, 'cm'),
        legend.text=element_text(size=11),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm") ))

os5

ggarrange(os1, os2, os3, os4, os5, ncol = 1, nrow = 5)

######################################## 2. forest plot for multicox  ##################################################################################
gong <- read_excel("coxforest.xlsx",
                   sheet = "all")
# text
gong_table<-cbind(
  c('Score (outcome)', gong$`Score (outcome)`),
  c('Group', gong$Quartiles),
  c('aHR (95%CI)', gong$HR),
  c('p value', gong$p))

# data
dput(names(gong))
gong_data<-subset(gong, select=c("coef", "low", "high"))
gong_data<-gong_data %>% add_row(coef = NA, low = NA, high = NA, .before = 1)

forestplot(gong_table,
           gong_data$coef, gong_data$low, gong_data$high,
           align='l', xlab='Hazard ratio', graphwidth= unit(10, 'cm') ,
           is.summary=c(TRUE, rep(FALSE,28)),
           line.margin = unit(0.8, 'cm'),
           zero=1,
           col=fpColors(box="#D73027",
                        line="#4575B4"),
           txt_gp = fpTxtGp(label=gpar(fontsize=15, cex =0.8, fontfamily='Arial'),
                            ticks=gpar(fontsize=15, cex=0.8, fontfamily='Arial'),
                            xlab=gpar(fontsize=15,cex = 0.8, fontfamily='Arial')),
           xticks=c(0.5, 1, 1.5, 2),
           lwd.ci=1.3, lwd.zero=1,
           ci.vertices = TRUE, ci.vertices.height=0.17,
           boxsize=0.3,
           colgap = unit(1, 'cm'))

######################################## 3. forest plot for AUC  ##################################################################################
AUC_OS <- read_excel("AUC_OS.xlsx")
# only plot significant
AUC_OS<-subset(AUC_OS, significant=='yes')
dput(names(AUC_OS))
AUC_OS$model<-factor(AUC_OS$model, levels = c( 'Gong et al, 2020',
                                               'Gündert et al, 2019','Huang et al, 2021 (OS)', 'Huang et al, 2021 (DFS)',
                                               'Wang Y et al, 2020'))
summary(AUC_OS$model)

## handpick color
display.brewer.pal(7, 'Set1')
brewer.pal(n = 7, name = "Set1")

compare3Cgzplot<-ggplot(AUC_OS, aes(times, AUC)) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper, color = model),
    position = position_dodge(0.3), width = 0.2, size=0.6
  )+
  geom_point(aes(color = model), position = position_dodge(0.3), size=2) +
  scale_color_manual(values = c("#A65628" , "#E31A1C", "#1F78B4", "#A6CEE3", "#33A02C"))+
  labs(x="Follow-up time points", y = "AUC")+
  geom_hline(yintercept=0.5, linetype="dashed", color = "gray")+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        legend.position='top',
        legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill = "white"),
        legend.key.size = unit(1.5, "cm"),
        legend.key.width = unit(0.6,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=13),
        legend.spacing.x = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.5, 'cm'),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm") ))+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))+
  scale_y_continuous(limits=c(0,1), breaks = seq(0,1,0.1))

compare3Cgzplot

# add a table under figure
# make table
data_table <- ggplot(AUC_OS, aes(x = times, y = model,
                                 label = format(`AUC (95%CI)`, nsmall = 1), colour = model)) +
  geom_text(size = 5) + theme_bw() +
  scale_y_discrete(limits = rev(c('Gong et al, 2020',
                                  'Gündert et al, 2019','Huang et al, 2021 (OS)', 'Huang et al, 2021 (DFS)',
                                  'Wang Y et al, 2020') ),
                   labels = rev(c('Gong',
                                  'Günd9','Hu)', 'Huan)','Wan0')))+
  scale_color_manual(values = c("black" , "black", "black", "black", "black"))+
  theme(axis.text = element_text(size = 12), axis.title=element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.title.x = element_blank(),
        #axis.text.y  = element_text(colour=rev(c("#E41A1C" , "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
        axis.title.y = element_blank()
  )


# combine table and figure
Layout <- grid.layout(nrow = 2, ncol = 1, heights = unit(c(2,0.5), c("null", "null")))
grid.show.layout(Layout)
vplayout <- function(...) {
  grid.newpage()
  pushViewport(viewport(layout = Layout))
}
subplot <- function(x, y) viewport(layout.pos.row = x,
                                   layout.pos.col = y)

mmplot <- function(a, b) {
  vplayout()
  print(a, vp = subplot(1, 1))
  print(b, vp = subplot(2, 1))
}

mmplot(compare3Cgzplot, data_table)

### Difference in AUC ###
library(readxl)
AUC_OS <- read_excel("DifAUC_imput_OS.xlsx")
# only plot significant
AUC_OS<-subset(AUC_OS, significant=='yes')
dput(names(AUC_OS))
AUC_OS$model<-factor(AUC_OS$model, levels = c( 'Gong et al, 2020',
                                               'Gündert et al, 2019','Huang et al, 2021 (OS)', 'Huang et al, 2021 (DFS)',
                                               'Wang Y et al, 2020'))

summary(AUC_OS$model)

compare3Cgzplot<-ggplot(AUC_OS, aes(times, diff_AUC)) +
  geom_errorbar(
    aes(ymin = ll, ymax = ul, color = model),
    position = position_dodge(0.3), width = 0.2, size=0.6
  )+
  geom_point(aes(color = model), position = position_dodge(0.3), size=2.5) +
  scale_color_manual(values = c("#A65628" , "#E31A1C", "#1F78B4", "#A6CEE3", "#33A02C"))+
  labs(x="Follow-up time points", y = "AUC difference")+
  geom_hline(yintercept=0, linetype="dashed", color = "gray")+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm") ))+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))+
  scale_y_continuous(limits=c(-0.05,0.05), breaks = seq(-0.05,0.05,0.01))

# add a table under figure
# make table
data_table <- ggplot(AUC_OS, aes(x = times, y = model,
                                 label = format(diffAUC_CI, nsmall = 1), colour = model)) +
  geom_text(size = 5) + theme_bw() +
  scale_y_discrete(limits = rev(c('Gong et al, 2020',
                                  'Gündert et al, 2019','Huang et al, 2021 (OS)', 'Huang et al, 2021 (DFS)',
                                  'Wang Y et al, 2020') ),
                   labels = c('Gong',
                              'Günd9','Hu)', 'Huan)','Wan0'))+
  scale_color_manual(values = c("black" , "black", "black", "black", "black"))+
  theme(axis.text = element_text(size = 12), axis.title=element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


# combine table and figure
Layout <- grid.layout(nrow = 2, ncol = 1, heights = unit(c(2,0.5), c("null", "null")))
grid.show.layout(Layout)
vplayout <- function(...) {
  grid.newpage()
  pushViewport(viewport(layout = Layout))
}
subplot <- function(x, y) viewport(layout.pos.row = x,
                                   layout.pos.col = y)

mmplot <- function(a, b) {
  vplayout()
  print(a, vp = subplot(1, 1))
  print(b, vp = subplot(2, 1))
}

mmplot(compare3Cgzplot, data_table)













