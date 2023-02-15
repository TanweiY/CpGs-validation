library(rms)
library(rlist)
library(ggpubr)
## OS ###
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")

multi_impu_cr14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr14[[i]]<-merge(unicox_cr14, multi_impu[[i]], by='id')
}

covariates_con<-c("Gong_2020_CR14", "WangX_2020_CR14", "Huang_2021_DFS_CR14", "Huang_2021_OS_CR14", "Li_2021_CR14")

# model1
baseriskav<-vector(20, mode="list")
basehazardtime<-vector(20,mode="list")

for (i in 1:20) {
  final_model <- cph(Surv(timey, death_all)~Gong_2020_CR14, data=multi_impu_cr14[[i]],init=1, iter.max=0,x=TRUE,y=TRUE)
  basehazard<-basehaz(final_model, centered=FALSE)
  baseriskav[i]<-as.data.frame(basehazard$hazard)
  basehazardtime[i]<-as.data.frame(basehazard$time)
}

baserisk<-rowMeans(list.cbind(baseriskav))
basehazardtime<-rowMeans(list.cbind(basehazardtime))
baseriskvali<-data.frame(baserisk, basehazardtime)
baseriskvali$basesurv<-exp(-baseriskvali$baserisk)

# Calibration plot
n_impu <- 20
sample_size<-2310

S2_nomogram <- seq(1:sample_size)
S5_nomogram <- seq(1:sample_size)
S10_nomogram <- seq(1:sample_size)

survival_predicted_2_combine <- seq(1:20)
survival_predicted_5_combine <- seq(1:20)
survival_predicted_10_combine <- seq(1:20)

survival_observed_2_combine <- seq(1:20)
survival_observed_5_combine <- seq(1:20)
survival_observed_10_combine <- seq(1:20)

survival_observed_2_var_combine <- seq(1:20)
survival_observed_5_var_combine <- seq(1:20)
survival_observed_10_var_combine <- seq(1:20)

for (i in 1:n_impu) {
  # calculate predicted survival probability
  multi_impu_cr14[[i]]$S2_nomogram <- 0.9341356^exp(multi_impu_cr14[[i]]$Gong_2020_CR14)
  multi_impu_cr14[[i]]$S5_nomogram <- 0.8385566^exp(multi_impu_cr14[[i]]$Gong_2020_CR14)
  multi_impu_cr14[[i]]$S10_nomogram <- 0.7244421^exp(multi_impu_cr14[[i]]$Gong_2020_CR14)
  
  # evenly divide patients into 10 groups
  
  multi_impu_cr14[[i]]$group10<-cut(multi_impu_cr14[[i]]$Gong_2020_CR14, quantile(multi_impu_cr14[[i]]$Gong_2020_CR14, seq(0,1,0.1)), right=FALSE, labels=c(1:10))
  multi_impu_cr14[[i]]$group10[multi_impu_cr14[[i]]$Gong_2020_CR14==max(multi_impu_cr14[[i]]$Gong_2020_CR14)] <- 10
  
  
  survival_predicted_2 <- 0
  survival_predicted_5 <- 0
  survival_predicted_10 <- 0
  
  
  # 2-year predicted survival
  survival_predicted_2 <- aggregate(multi_impu_cr14[[i]]$S2_nomogram, list(multi_impu_cr14[[i]]$group10), mean)
  survival_predicted_2_combine <- data.frame(survival_predicted_2_combine,survival_predicted_2$x)
  
  # 5-year predicted survival
  survival_predicted_5 <- aggregate(multi_impu_cr14[[i]]$S5_nomogram, list(multi_impu_cr14[[i]]$group10), mean)
  survival_predicted_5_combine <- data.frame(survival_predicted_5_combine,survival_predicted_5$x)
  
  # 10-year predicted survival
  survival_predicted_10 <- aggregate(multi_impu_cr14[[i]]$S10_nomogram, list(multi_impu_cr14[[i]]$group10), mean)
  survival_predicted_10_combine <- data.frame(survival_predicted_10_combine,survival_predicted_10$x)
  
  # observed survival
  survival_observed_2 <- 0
  survival_observed_2_var <- 0
  
  survival_observed_5 <- 0
  survival_observed_5_var <- 0
  
  survival_observed_10 <- 0
  survival_observed_10_var <- 0
  
  for (j in 1:10) {
    
    data_temp <- subset(multi_impu_cr14[[i]],multi_impu_cr14[[i]]$group10==j)
    
    fit_calibration <- survfit(Surv(timey, death_all) ~ 1, data=data_temp)
    
    survival_observed_2[j] <- min(fit_calibration$surv[fit_calibration$time <= 2])
    survival_observed_2_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 2],1))^2
    
    survival_observed_5[j] <- min(fit_calibration$surv[fit_calibration$time <= 5])
    survival_observed_5_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 5],1))^2
    
    survival_observed_10[j] <- min(fit_calibration$surv[fit_calibration$time <= 10])
    survival_observed_10_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 10],1))^2
  }
  
  
  survival_observed_2_combine <- data.frame(survival_observed_2_combine,survival_observed_2)
  survival_observed_5_combine <- data.frame(survival_observed_5_combine,survival_observed_5)
  survival_observed_10_combine <- data.frame(survival_observed_10_combine,survival_observed_10)
  
  survival_observed_2_var_combine <- data.frame(survival_observed_2_var_combine,survival_observed_2_var)
  survival_observed_5_var_combine <- data.frame(survival_observed_5_var_combine,survival_observed_5_var)
  survival_observed_10_var_combine <- data.frame(survival_observed_10_var_combine,survival_observed_10_var)
  
}


# plot calibration plot

# 2-year

survival_predicted_2_final <- exp(rowMeans(log(survival_predicted_2_combine[,-1])))
survival_observed_2_final <- exp(rowMeans(log(survival_observed_2_combine[,-1])))
survival_observed_2_var_final <- rowMeans(survival_observed_2_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_2_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower2_final<-exp(log(survival_observed_2_final) - qnorm(0.975)*survival_observed_2_var_final^0.5)
survival_upper2_final<-exp(log(survival_observed_2_final) + qnorm(0.975)*survival_observed_2_var_final^0.5)
survival_comparison2 <- data.frame(survival_predicted_2_final, survival_observed_2_final,
                                   survival_lower2_final,survival_upper2_final)

#survival_comparison2$survival_upper2_final<-ifelse(survival_comparison2$survival_upper2_final>1, 1,survival_comparison2$survival_upper1_final)
#survival_comparison2$underestimate<-(survival_comparison2$survival_observed_2_final-survival_comparison2$survival_predicted_2_final)/survival_comparison2$survival_observed_2_final

c1<-ggplot(data=survival_comparison2, aes(x=survival_predicted_2_final, y=survival_observed_2_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison2, mapping=aes(x=survival_predicted_2_final, ymin=survival_lower2_final,
                                                       ymax=survival_upper2_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="A. Calibration curve at 2 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme( axis.text = element_text(size = 15), axis.title=element_text(size=15),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=15),
         axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
         axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.line = element_line(color='black'),
  )


# 5-year
survival_predicted_5_final <- exp(rowMeans(log(survival_predicted_5_combine[,-1])))
survival_observed_5_final <- exp(rowMeans(log(survival_observed_5_combine[,-1])))
survival_observed_5_var_final <- rowMeans(survival_observed_5_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_5_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower5_final<-exp(log(survival_observed_5_final) - qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_upper5_final<-exp(log(survival_observed_5_final) + qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_comparison5 <- data.frame(survival_predicted_5_final, survival_observed_5_final,
                                   survival_lower5_final,survival_upper5_final)

survival_comparison5$survival_upper5_final<-ifelse(survival_comparison5$survival_upper5_final>1, 1,survival_comparison5$survival_upper5_final)
survival_comparison5$underestimate<-(survival_comparison5$survival_observed_5_final-survival_comparison5$survival_predicted_5_final)/survival_comparison5$survival_observed_5_final

c2<-ggplot(data=survival_comparison5, aes(x=survival_predicted_5_final, y=survival_observed_5_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison5, mapping=aes(x=survival_predicted_5_final, ymin=survival_lower5_final,
                                                       ymax=survival_upper5_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="B. Calibration curve at 5 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme( axis.text = element_text(size = 15), axis.title=element_text(size=15),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=15),
         axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
         axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.line = element_line(color='black'),
  )


# 10-year
survival_predicted_10_final <- exp(rowMeans(log(survival_predicted_10_combine[,-1])))
survival_observed_10_final <- exp(rowMeans(log(survival_observed_10_combine[,-1])))
survival_observed_10_var_final <- rowMeans(survival_observed_10_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_10_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower10_final<-exp(log(survival_observed_10_final) - qnorm(0.975)*survival_observed_10_var_final^0.5)
survival_upper10_final<-exp(log(survival_observed_10_final) + qnorm(0.975)*survival_observed_10_var_final^0.5)
survival_comparison10 <- data.frame(survival_predicted_10_final, survival_observed_10_final,
                                    survival_lower10_final,survival_upper10_final)

survival_comparison10$survival_upper10_final<-ifelse(survival_comparison10$survival_upper10_final>1, 1,survival_comparison10$survival_upper10_final)
survival_comparison10$underestimate<-(survival_comparison10$survival_observed_10_final-survival_comparison10$survival_predicted_10_final)/survival_comparison10$survival_observed_10_final

c3<-ggplot(data=survival_comparison10, aes(x=survival_predicted_10_final, y=survival_observed_10_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison10, mapping=aes(x=survival_predicted_10_final, ymin=survival_lower10_final,
                                                        ymax=survival_upper10_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="C. Calibration curve at 10 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.line = element_line(color='black'),
  )

Gong_2020_OS<-ggarrange(c1, c2, c3, ncol = 3, nrow = 1)
save(Gong_2020_OS, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/calibraGong.RData")

## Huang OS 
## OS ###
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")

multi_impu_cr14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr14[[i]]<-merge(unicox_cr14, multi_impu[[i]], by='id')
}


# model1
baseriskav<-vector(20, mode="list")
basehazardtime<-vector(20,mode="list")

for (i in 1:20) {
  final_model <- cph(Surv(timey, death_all)~Huang_2021_OS_CR14, data=multi_impu_cr14[[i]],init=1, iter.max=0,x=TRUE,y=TRUE)
  basehazard<-basehaz(final_model, centered=FALSE)
  baseriskav[i]<-as.data.frame(basehazard$hazard)
  basehazardtime[i]<-as.data.frame(basehazard$time)
}

baserisk<-rowMeans(list.cbind(baseriskav))
basehazardtime<-rowMeans(list.cbind(basehazardtime))
baseriskvali<-data.frame(baserisk, basehazardtime)
baseriskvali$basesurv<-exp(-baseriskvali$baserisk)

# Calibration plot
n_impu <- 20
sample_size<-2310

S2_nomogram <- seq(1:sample_size)
S5_nomogram <- seq(1:sample_size)
S10_nomogram <- seq(1:sample_size)

survival_predicted_2_combine <- seq(1:20)
survival_predicted_5_combine <- seq(1:20)
survival_predicted_10_combine <- seq(1:20)

survival_observed_2_combine <- seq(1:20)
survival_observed_5_combine <- seq(1:20)
survival_observed_10_combine <- seq(1:20)

survival_observed_2_var_combine <- seq(1:20)
survival_observed_5_var_combine <- seq(1:20)
survival_observed_10_var_combine <- seq(1:20)

for (i in 1:n_impu) {
  # calculate predicted survival probability
  multi_impu_cr14[[i]]$S2_nomogram <- 0.9738218^exp(multi_impu_cr14[[i]]$Huang_2021_OS_CR14)
  multi_impu_cr14[[i]]$S5_nomogram <- 0.9326499^exp(multi_impu_cr14[[i]]$Huang_2021_OS_CR14)
  multi_impu_cr14[[i]]$S10_nomogram <- 0.8771890^exp(multi_impu_cr14[[i]]$Huang_2021_OS_CR14)
  
  # evenly divide patients into 10 groups
  
  multi_impu_cr14[[i]]$group10<-cut(multi_impu_cr14[[i]]$Huang_2021_OS_CR14, quantile(multi_impu_cr14[[i]]$Huang_2021_OS_CR14, seq(0,1,0.1)), right=FALSE, labels=c(1:10))
  multi_impu_cr14[[i]]$group10[multi_impu_cr14[[i]]$Huang_2021_OS_CR14==max(multi_impu_cr14[[i]]$Huang_2021_OS_CR14)] <- 10
  
  
  survival_predicted_2 <- 0
  survival_predicted_5 <- 0
  survival_predicted_10 <- 0
  
  
  # 2-year predicted survival
  survival_predicted_2 <- aggregate(multi_impu_cr14[[i]]$S2_nomogram, list(multi_impu_cr14[[i]]$group10), mean)
  survival_predicted_2_combine <- data.frame(survival_predicted_2_combine,survival_predicted_2$x)
  
  # 5-year predicted survival
  survival_predicted_5 <- aggregate(multi_impu_cr14[[i]]$S5_nomogram, list(multi_impu_cr14[[i]]$group10), mean)
  survival_predicted_5_combine <- data.frame(survival_predicted_5_combine,survival_predicted_5$x)
  
  # 10-year predicted survival
  survival_predicted_10 <- aggregate(multi_impu_cr14[[i]]$S10_nomogram, list(multi_impu_cr14[[i]]$group10), mean)
  survival_predicted_10_combine <- data.frame(survival_predicted_10_combine,survival_predicted_10$x)
  
  # observed survival
  survival_observed_2 <- 0
  survival_observed_2_var <- 0
  
  survival_observed_5 <- 0
  survival_observed_5_var <- 0
  
  survival_observed_10 <- 0
  survival_observed_10_var <- 0
  
  for (j in 1:10) {
    
    data_temp <- subset(multi_impu_cr14[[i]],multi_impu_cr14[[i]]$group10==j)
    
    fit_calibration <- survfit(Surv(timey, death_all) ~ 1, data=data_temp)
    
    survival_observed_2[j] <- min(fit_calibration$surv[fit_calibration$time <= 2])
    survival_observed_2_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 2],1))^2
    
    survival_observed_5[j] <- min(fit_calibration$surv[fit_calibration$time <= 5])
    survival_observed_5_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 5],1))^2
    
    survival_observed_10[j] <- min(fit_calibration$surv[fit_calibration$time <= 10])
    survival_observed_10_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 10],1))^2
  }
  
  
  survival_observed_2_combine <- data.frame(survival_observed_2_combine,survival_observed_2)
  survival_observed_5_combine <- data.frame(survival_observed_5_combine,survival_observed_5)
  survival_observed_10_combine <- data.frame(survival_observed_10_combine,survival_observed_10)
  
  survival_observed_2_var_combine <- data.frame(survival_observed_2_var_combine,survival_observed_2_var)
  survival_observed_5_var_combine <- data.frame(survival_observed_5_var_combine,survival_observed_5_var)
  survival_observed_10_var_combine <- data.frame(survival_observed_10_var_combine,survival_observed_10_var)
  
}

# plot calibration plot

# 2-year

survival_predicted_2_final <- exp(rowMeans(log(survival_predicted_2_combine[,-1])))
survival_observed_2_final <- exp(rowMeans(log(survival_observed_2_combine[,-1])))
survival_observed_2_var_final <- rowMeans(survival_observed_2_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_2_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower2_final<-exp(log(survival_observed_2_final) - qnorm(0.975)*survival_observed_2_var_final^0.5)
survival_upper2_final<-exp(log(survival_observed_2_final) + qnorm(0.975)*survival_observed_2_var_final^0.5)
survival_comparison2 <- data.frame(survival_predicted_2_final, survival_observed_2_final,
                                   survival_lower2_final,survival_upper2_final)

survival_comparison2$survival_upper2_final<-ifelse(survival_comparison2$survival_upper2_final>1, 1,survival_comparison2$survival_upper1_final)
#survival_comparison2$underestimate<-(survival_comparison2$survival_observed_2_final-survival_comparison2$survival_predicted_2_final)/survival_comparison2$survival_observed_2_final

c1<-ggplot(data=survival_comparison2, aes(x=survival_predicted_2_final, y=survival_observed_2_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison2, mapping=aes(x=survival_predicted_2_final, ymin=survival_lower2_final,
                                                       ymax=survival_upper2_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="A. Calibration curve at 2 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme( axis.text = element_text(size = 15), axis.title=element_text(size=15),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=15),
         axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
         axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.line = element_line(color='black'),
  )

c1
# 5-year

survival_predicted_5_final <- exp(rowMeans(log(survival_predicted_5_combine[,-1])))
survival_observed_5_final <- exp(rowMeans(log(survival_observed_5_combine[,-1])))
survival_observed_5_var_final <- rowMeans(survival_observed_5_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_5_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower5_final<-exp(log(survival_observed_5_final) - qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_upper5_final<-exp(log(survival_observed_5_final) + qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_comparison5 <- data.frame(survival_predicted_5_final, survival_observed_5_final,
                                   survival_lower5_final,survival_upper5_final)

survival_comparison5$survival_upper5_final<-ifelse(survival_comparison5$survival_upper5_final>1, 1,survival_comparison5$survival_upper5_final)
survival_comparison5$underestimate<-(survival_comparison5$survival_observed_5_final-survival_comparison5$survival_predicted_5_final)/survival_comparison5$survival_observed_5_final

c2<-ggplot(data=survival_comparison5, aes(x=survival_predicted_5_final, y=survival_observed_5_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison5, mapping=aes(x=survival_predicted_5_final, ymin=survival_lower5_final,
                                                       ymax=survival_upper5_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="B. Calibration curve at 5 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme( axis.text = element_text(size = 15), axis.title=element_text(size=15),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=15),
         axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
         axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.line = element_line(color='black'),
  )
c2

# 10-year
survival_predicted_10_final <- exp(rowMeans(log(survival_predicted_10_combine[,-1])))
survival_observed_10_final <- exp(rowMeans(log(survival_observed_10_combine[,-1])))
survival_observed_10_var_final <- rowMeans(survival_observed_10_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_10_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower10_final<-exp(log(survival_observed_10_final) - qnorm(0.975)*survival_observed_10_var_final^0.5)
survival_upper10_final<-exp(log(survival_observed_10_final) + qnorm(0.975)*survival_observed_10_var_final^0.5)
survival_comparison10 <- data.frame(survival_predicted_10_final, survival_observed_10_final,
                                    survival_lower10_final,survival_upper10_final)

survival_comparison10$survival_upper10_final<-ifelse(survival_comparison10$survival_upper10_final>1, 1,survival_comparison10$survival_upper10_final)
survival_comparison10$underestimate<-(survival_comparison10$survival_observed_10_final-survival_comparison10$survival_predicted_10_final)/survival_comparison10$survival_observed_10_final

c3<-ggplot(data=survival_comparison10, aes(x=survival_predicted_10_final, y=survival_observed_10_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison10, mapping=aes(x=survival_predicted_10_final, ymin=survival_lower10_final,
                                                        ymax=survival_upper10_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="C. Calibration curve at 10 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.line = element_line(color='black'),
  )

c3

Huang_2021_OS<-ggarrange(c1, c2, c3, ncol = 3, nrow = 1)
save(Huang_2021_OS, 
     file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/calibrHuangOS.RData")

## Huang DFS ##
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr14.RData")

multi_impu_cr14<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr14[[i]]<-merge(unicox_cr14, multi_impu[[i]], by='id')
}

# model1
baseriskav<-vector(20, mode="list")
basehazardtime<-vector(20,mode="list")

for (i in 1:20) {
  final_model <- cph(Surv(timey, death_all)~Huang_2021_DFS_CR14, data=multi_impu_cr14[[i]],init=1, iter.max=0,x=TRUE,y=TRUE)
  basehazard<-basehaz(final_model, centered=FALSE)
  baseriskav[i]<-as.data.frame(basehazard$hazard)
  basehazardtime[i]<-as.data.frame(basehazard$time)
}

baserisk<-rowMeans(list.cbind(baseriskav))
basehazardtime<-rowMeans(list.cbind(basehazardtime))
baseriskvali<-data.frame(baserisk, basehazardtime)
baseriskvali$basesurv<-exp(-baseriskvali$baserisk)

# Calibration plot
n_impu <- 20
sample_size<-2310

S2_nomogram <- seq(1:sample_size)
S5_nomogram <- seq(1:sample_size)
S10_nomogram <- seq(1:sample_size)

survival_predicted_2_combine <- seq(1:20)
survival_predicted_5_combine <- seq(1:20)
survival_predicted_10_combine <- seq(1:20)

survival_observed_2_combine <- seq(1:20)
survival_observed_5_combine <- seq(1:20)
survival_observed_10_combine <- seq(1:20)

survival_observed_2_var_combine <- seq(1:20)
survival_observed_5_var_combine <- seq(1:20)
survival_observed_10_var_combine <- seq(1:20)

for (i in 1:n_impu) {
  # calculate predicted survival probability
  multi_impu_cr14[[i]]$S2_nomogram <- 0.9738218^exp(multi_impu_cr14[[i]]$Huang_2021_DFS_CR14)
  multi_impu_cr14[[i]]$S5_nomogram <- 0.9326499^exp(multi_impu_cr14[[i]]$Huang_2021_DFS_CR14)
  multi_impu_cr14[[i]]$S10_nomogram <- 0.8771890^exp(multi_impu_cr14[[i]]$Huang_2021_DFS_CR14)
  
  # evenly divide patients into 10 groups
  
  multi_impu_cr14[[i]]$group10<-cut(multi_impu_cr14[[i]]$Huang_2021_DFS_CR14, quantile(multi_impu_cr14[[i]]$Huang_2021_DFS_CR14, seq(0,1,0.1)), right=FALSE, labels=c(1:10))
  multi_impu_cr14[[i]]$group10[multi_impu_cr14[[i]]$Huang_2021_DFS_CR14==max(multi_impu_cr14[[i]]$Huang_2021_DFS_CR14)] <- 10
  
  
  survival_predicted_2 <- 0
  survival_predicted_5 <- 0
  survival_predicted_10 <- 0
  
  
  # 2-year predicted survival
  survival_predicted_2 <- aggregate(multi_impu_cr14[[i]]$S2_nomogram, list(multi_impu_cr14[[i]]$group10), mean)
  survival_predicted_2_combine <- data.frame(survival_predicted_2_combine,survival_predicted_2$x)
  
  # 5-year predicted survival
  survival_predicted_5 <- aggregate(multi_impu_cr14[[i]]$S5_nomogram, list(multi_impu_cr14[[i]]$group10), mean)
  survival_predicted_5_combine <- data.frame(survival_predicted_5_combine,survival_predicted_5$x)
  
  # 10-year predicted survival
  survival_predicted_10 <- aggregate(multi_impu_cr14[[i]]$S10_nomogram, list(multi_impu_cr14[[i]]$group10), mean)
  survival_predicted_10_combine <- data.frame(survival_predicted_10_combine,survival_predicted_10$x)
  
  # observed survival
  survival_observed_2 <- 0
  survival_observed_2_var <- 0
  
  survival_observed_5 <- 0
  survival_observed_5_var <- 0
  
  survival_observed_10 <- 0
  survival_observed_10_var <- 0
  
  for (j in 1:10) {
    
    data_temp <- subset(multi_impu_cr14[[i]],multi_impu_cr14[[i]]$group10==j)
    
    fit_calibration <- survfit(Surv(timey, death_all) ~ 1, data=data_temp)
    
    survival_observed_2[j] <- min(fit_calibration$surv[fit_calibration$time <= 2])
    survival_observed_2_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 2],1))^2
    
    survival_observed_5[j] <- min(fit_calibration$surv[fit_calibration$time <= 5])
    survival_observed_5_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 5],1))^2
    
    survival_observed_10[j] <- min(fit_calibration$surv[fit_calibration$time <= 10])
    survival_observed_10_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 10],1))^2
  }
  
  
  survival_observed_2_combine <- data.frame(survival_observed_2_combine,survival_observed_2)
  survival_observed_5_combine <- data.frame(survival_observed_5_combine,survival_observed_5)
  survival_observed_10_combine <- data.frame(survival_observed_10_combine,survival_observed_10)
  
  survival_observed_2_var_combine <- data.frame(survival_observed_2_var_combine,survival_observed_2_var)
  survival_observed_5_var_combine <- data.frame(survival_observed_5_var_combine,survival_observed_5_var)
  survival_observed_10_var_combine <- data.frame(survival_observed_10_var_combine,survival_observed_10_var)
  
}

# plot calibration plot

# 2-year

survival_predicted_2_final <- exp(rowMeans(log(survival_predicted_2_combine[,-1])))
survival_observed_2_final <- exp(rowMeans(log(survival_observed_2_combine[,-1])))
survival_observed_2_var_final <- rowMeans(survival_observed_2_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_2_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower2_final<-exp(log(survival_observed_2_final) - qnorm(0.975)*survival_observed_2_var_final^0.5)
survival_upper2_final<-exp(log(survival_observed_2_final) + qnorm(0.975)*survival_observed_2_var_final^0.5)
survival_comparison2 <- data.frame(survival_predicted_2_final, survival_observed_2_final,
                                   survival_lower2_final,survival_upper2_final)

#survival_comparison2$survival_upper2_final<-ifelse(survival_comparison2$survival_upper2_final>1, 1,survival_comparison2$survival_upper1_final)
survival_comparison2$underestimate<-(survival_comparison2$survival_observed_2_final-survival_comparison2$survival_predicted_2_final)/survival_comparison2$survival_observed_2_final

c1<-ggplot(data=survival_comparison2, aes(x=survival_predicted_2_final, y=survival_observed_2_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison2, mapping=aes(x=survival_predicted_2_final, ymin=survival_lower2_final,
                                                       ymax=survival_upper2_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="A. Calibration curve at 2 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme( axis.text = element_text(size = 15), axis.title=element_text(size=15),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=15),
         axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
         axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.line = element_line(color='black'),
  )

c1
# 5-year

survival_predicted_5_final <- exp(rowMeans(log(survival_predicted_5_combine[,-1])))
survival_observed_5_final <- exp(rowMeans(log(survival_observed_5_combine[,-1])))
survival_observed_5_var_final <- rowMeans(survival_observed_5_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_5_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower5_final<-exp(log(survival_observed_5_final) - qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_upper5_final<-exp(log(survival_observed_5_final) + qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_comparison5 <- data.frame(survival_predicted_5_final, survival_observed_5_final,
                                   survival_lower5_final,survival_upper5_final)

survival_comparison5$survival_upper5_final<-ifelse(survival_comparison5$survival_upper5_final>1, 1,survival_comparison5$survival_upper5_final)
survival_comparison5$underestimate<-(survival_comparison5$survival_observed_5_final-survival_comparison5$survival_predicted_5_final)/survival_comparison5$survival_observed_5_final

c2<-ggplot(data=survival_comparison5, aes(x=survival_predicted_5_final, y=survival_observed_5_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison5, mapping=aes(x=survival_predicted_5_final, ymin=survival_lower5_final,
                                                       ymax=survival_upper5_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="B. Calibration curve at 5 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme( axis.text = element_text(size = 15), axis.title=element_text(size=15),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=15),
         axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
         axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.line = element_line(color='black'),
  )
c2

# 10-year
survival_predicted_10_final <- exp(rowMeans(log(survival_predicted_10_combine[,-1])))
survival_observed_10_final <- exp(rowMeans(log(survival_observed_10_combine[,-1])))
survival_observed_10_var_final <- rowMeans(survival_observed_10_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_10_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower10_final<-exp(log(survival_observed_10_final) - qnorm(0.975)*survival_observed_10_var_final^0.5)
survival_upper10_final<-exp(log(survival_observed_10_final) + qnorm(0.975)*survival_observed_10_var_final^0.5)
survival_comparison10 <- data.frame(survival_predicted_10_final, survival_observed_10_final,
                                    survival_lower10_final,survival_upper10_final)

#survival_comparison10$survival_upper10_final<-ifelse(survival_comparison10$survival_upper10_final>1, 1,survival_comparison10$survival_upper10_final)
survival_comparison10$underestimate<-(survival_comparison10$survival_observed_10_final-survival_comparison10$survival_predicted_10_final)/survival_comparison10$survival_observed_10_final

c3<-ggplot(data=survival_comparison10, aes(x=survival_predicted_10_final, y=survival_observed_10_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison10, mapping=aes(x=survival_predicted_10_final, ymin=survival_lower10_final,
                                                        ymax=survival_upper10_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="C. Calibration curve at 10 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.line = element_line(color='black'),
  )

c3

Huang_2021_PFS<-ggarrange(c1, c2, c3, ncol = 3, nrow = 1)
save(Huang_2021_PFS, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/calibrHuangDFS.RData")

## Grundent

load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/multicox_impu20list.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/unicox_cr13.RData")

multi_impu_cr13<-vector(20, mode="list")

for (i in 1:20) {
  multi_impu_cr13[[i]]<-merge(unicox_cr13, multi_impu[[i]], by='id')
  multi_impu_cr13[[i]]$Gundert_2019_CR1<-(-multi_impu_cr13[[i]]$Gundert_2019_CR1)
}


summary(multi_impu_cr13[[1]]$Gundert_2019_CR1)

baseriskav<-vector(20, mode="list")
basehazardtime<-vector(20,mode="list")
library(rms)

for (i in 1:20) {
  final_model <- cph(Surv(timey, death_all)~Gundert_2019_CR1, data=multi_impu_cr13[[i]],init=1, iter.max=0,x=TRUE,y=TRUE)
  basehazard<-basehaz(final_model, centered=FALSE)
  baseriskav[i]<-as.data.frame(basehazard$hazard)
  basehazardtime[i]<-as.data.frame(basehazard$time)
}

baserisk<-rowMeans(list.cbind(baseriskav))
basehazardtime<-rowMeans(list.cbind(basehazardtime))
baseriskvali<-data.frame(baserisk, basehazardtime)
baseriskvali$basesurv<-exp(-baseriskvali$baserisk)

# Calibration plot
n_impu <- 20
sample_size<-1738

S2_nomogram <- seq(1:sample_size)
S5_nomogram <- seq(1:sample_size)
S10_nomogram <- seq(1:sample_size)

survival_predicted_2_combine <- seq(1:20)
survival_predicted_5_combine <- seq(1:20)
survival_predicted_10_combine <- seq(1:20)

survival_observed_2_combine <- seq(1:20)
survival_observed_5_combine <- seq(1:20)
survival_observed_10_combine <- seq(1:20)

survival_observed_2_var_combine <- seq(1:20)
survival_observed_5_var_combine <- seq(1:20)
survival_observed_10_var_combine <- seq(1:20)

for (i in 1:n_impu) {
  # calculate predicted survival probability
  multi_impu_cr13[[i]]$S2_nomogram <- 0.9996193^exp(multi_impu_cr13[[i]]$Gundert_2019_CR1)
  multi_impu_cr13[[i]]$S5_nomogram <- 0.9983834^exp(multi_impu_cr13[[i]]$Gundert_2019_CR1)
  multi_impu_cr13[[i]]$S10_nomogram <- 0.9954706^exp(multi_impu_cr13[[i]]$Gundert_2019_CR1)
  
  # evenly divide patients into 10 groups
  
  multi_impu_cr13[[i]]$group10<-cut(multi_impu_cr13[[i]]$Gundert_2019_CR1, quantile(multi_impu_cr13[[i]]$Gundert_2019_CR1, seq(0,1,0.1)), right=FALSE, labels=c(1:10))
  multi_impu_cr13[[i]]$group10[multi_impu_cr13[[i]]$Gundert_2019_CR1==max(multi_impu_cr13[[i]]$Gundert_2019_CR1)] <- 10
  
  
  survival_predicted_2 <- 0
  survival_predicted_5 <- 0
  survival_predicted_10 <- 0
  
  
  # 2-year predicted survival
  survival_predicted_2 <- aggregate(multi_impu_cr13[[i]]$S2_nomogram, list(multi_impu_cr13[[i]]$group10), mean)
  survival_predicted_2_combine <- data.frame(survival_predicted_2_combine,survival_predicted_2$x)
  
  # 5-year predicted survival
  survival_predicted_5 <- aggregate(multi_impu_cr13[[i]]$S5_nomogram, list(multi_impu_cr13[[i]]$group10), mean)
  survival_predicted_5_combine <- data.frame(survival_predicted_5_combine,survival_predicted_5$x)
  
  # 10-year predicted survival
  survival_predicted_10 <- aggregate(multi_impu_cr13[[i]]$S10_nomogram, list(multi_impu_cr13[[i]]$group10), mean)
  survival_predicted_10_combine <- data.frame(survival_predicted_10_combine,survival_predicted_10$x)
  
  # observed survival
  survival_observed_2 <- 0
  survival_observed_2_var <- 0
  
  survival_observed_5 <- 0
  survival_observed_5_var <- 0
  
  survival_observed_10 <- 0
  survival_observed_10_var <- 0
  
  for (j in 1:10) {
    
    data_temp <- subset(multi_impu_cr13[[i]],multi_impu_cr13[[i]]$group10==j)
    
    fit_calibration <- survfit(Surv(timey, death_all) ~ 1, data=data_temp)
    
    survival_observed_2[j] <- min(fit_calibration$surv[fit_calibration$time <= 2])
    survival_observed_2_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 2],1))^2
    
    survival_observed_5[j] <- min(fit_calibration$surv[fit_calibration$time <= 5])
    survival_observed_5_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 5],1))^2
    
    survival_observed_10[j] <- min(fit_calibration$surv[fit_calibration$time <= 10])
    survival_observed_10_var[j] <- (tail(fit_calibration$std.err[fit_calibration$time <= 10],1))^2
  }
  
  
  survival_observed_2_combine <- data.frame(survival_observed_2_combine,survival_observed_2)
  survival_observed_5_combine <- data.frame(survival_observed_5_combine,survival_observed_5)
  survival_observed_10_combine <- data.frame(survival_observed_10_combine,survival_observed_10)
  
  survival_observed_2_var_combine <- data.frame(survival_observed_2_var_combine,survival_observed_2_var)
  survival_observed_5_var_combine <- data.frame(survival_observed_5_var_combine,survival_observed_5_var)
  survival_observed_10_var_combine <- data.frame(survival_observed_10_var_combine,survival_observed_10_var)
  
}

# plot calibration plot


# 2-year

survival_predicted_2_final <- exp(rowMeans(log(survival_predicted_2_combine[,-1])))
survival_observed_2_final <- exp(rowMeans(log(survival_observed_2_combine[,-1])))
survival_observed_2_var_final <- rowMeans(survival_observed_2_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_2_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower2_final<-exp(log(survival_observed_2_final) - qnorm(0.975)*survival_observed_2_var_final^0.5)
survival_upper2_final<-exp(log(survival_observed_2_final) + qnorm(0.975)*survival_observed_2_var_final^0.5)
survival_comparison2 <- data.frame(survival_predicted_2_final, survival_observed_2_final,
                                   survival_lower2_final,survival_upper2_final)

#survival_comparison2$survival_upper2_final<-ifelse(survival_comparison2$survival_upper2_final>1, 1,survival_comparison2$survival_upper1_final)
survival_comparison2$underestimate<-(survival_comparison2$survival_observed_2_final-survival_comparison2$survival_predicted_2_final)/survival_comparison2$survival_observed_2_final

c1<-ggplot(data=survival_comparison2, aes(x=survival_predicted_2_final, y=survival_observed_2_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison2, mapping=aes(x=survival_predicted_2_final, ymin=survival_lower2_final,
                                                       ymax=survival_upper2_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="A. Calibration curve at 2 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme( axis.text = element_text(size = 15), axis.title=element_text(size=15),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=15),
         axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
         axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.line = element_line(color='black'),
  )

c1
# 5-year

survival_predicted_5_final <- exp(rowMeans(log(survival_predicted_5_combine[,-1])))
survival_observed_5_final <- exp(rowMeans(log(survival_observed_5_combine[,-1])))
survival_observed_5_var_final <- rowMeans(survival_observed_5_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_5_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower5_final<-exp(log(survival_observed_5_final) - qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_upper5_final<-exp(log(survival_observed_5_final) + qnorm(0.975)*survival_observed_5_var_final^0.5)
survival_comparison5 <- data.frame(survival_predicted_5_final, survival_observed_5_final,
                                   survival_lower5_final,survival_upper5_final)

survival_comparison5$survival_upper5_final<-ifelse(survival_comparison5$survival_upper5_final>1, 1,survival_comparison5$survival_upper5_final)
survival_comparison5$underestimate<-(survival_comparison5$survival_observed_5_final-survival_comparison5$survival_predicted_5_final)/survival_comparison5$survival_observed_5_final

c2<-ggplot(data=survival_comparison5, aes(x=survival_predicted_5_final, y=survival_observed_5_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison5, mapping=aes(x=survival_predicted_5_final, ymin=survival_lower5_final,
                                                       ymax=survival_upper5_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="B. Calibration curve at 5 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme( axis.text = element_text(size = 15), axis.title=element_text(size=15),
         panel.background = element_rect(fill = "white"),
         plot.title = element_text(size=15),
         axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
         axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
         plot.margin = margin(1, 1, 1, 1, "cm"),
         axis.line = element_line(color='black'),
  )
c2

# 10-year
survival_predicted_10_final <- exp(rowMeans(log(survival_predicted_10_combine[,-1])))
survival_observed_10_final <- exp(rowMeans(log(survival_observed_10_combine[,-1])))
survival_observed_10_var_final <- rowMeans(survival_observed_10_var_combine[,-1]) + (1+1/n_impu)*apply(survival_observed_10_var_combine[,-1], MARGIN=1, FUN=var, na.rm=TRUE)

survival_lower10_final<-exp(log(survival_observed_10_final) - qnorm(0.975)*survival_observed_10_var_final^0.5)
survival_upper10_final<-exp(log(survival_observed_10_final) + qnorm(0.975)*survival_observed_10_var_final^0.5)
survival_comparison10 <- data.frame(survival_predicted_10_final, survival_observed_10_final,
                                    survival_lower10_final,survival_upper10_final)

survival_comparison10$survival_upper10_final<-ifelse(survival_comparison10$survival_upper10_final>1, 1,survival_comparison10$survival_upper10_final)
survival_comparison10$underestimate<-(survival_comparison10$survival_observed_10_final-survival_comparison10$survival_predicted_10_final)/survival_comparison10$survival_observed_10_final

c3<-ggplot(data=survival_comparison10, aes(x=survival_predicted_10_final, y=survival_observed_10_final)) +
  geom_line(size=1, colour="dodgerblue3")+
  geom_errorbar(data=survival_comparison10, mapping=aes(x=survival_predicted_10_final, ymin=survival_lower10_final,
                                                        ymax=survival_upper10_final),
                colour="dodgerblue3", size=1,alpha=0.5, linetype=1)+
  geom_point(size=3, colour="dodgerblue3")+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1,lty=2)+
  labs(title="C. Calibration curve at 10 year",x="Predicted Survival Probability", y = "Observed Survival Probability")+
  theme(axis.text = element_text(size = 15), axis.title=element_text(size=15),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size=15),
        axis.title.y = element_text(margin = margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        axis.title.x= element_text(margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm")),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.line = element_line(color='black'),
  )

c3
library(ggpubr)
Gundert_2019_OS_DSS<-ggarrange(c1, c2, c3, ncol = 3, nrow = 1)
save(Gundert_2019_OS_DSS, file = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/calibrGundert.RData")

## arrange multiple figures together 
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/calibrGundert.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/calibrHuangDFS.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/calibrHuangOS.RData")
load("/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/validation2023_processresults/results/calibraGong.RData")

ggarrange(Gong_2020_OS, Gundert_2019_OS_DSS,
          Huang_2021_OS, Huang_2021_PFS, ncol = 1, nrow = 4)

# 14* 15
install.packages('car', dependencies = TRUE)

options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/centos7/latest"))






































