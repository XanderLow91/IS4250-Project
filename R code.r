survey <- read.csv("H1N1_surveys.csv")

# Get wave starting and ending date: "MM dd to MM dd"
library(plyr)
mat <- ddply(survey,.(survey$wave,survey$date),nrow)
mat$V1 <- NULL
colnames(mat) <- c("wave","date")
mat <- ddply(mat, .(wave), mutate, start = min(date),end=max(date))
mat$date <- NULL
mat <- unique(mat) # unique start and ending
rownames(mat) <- mat$wave
mat$start <- as.Date(as.character(mat$start), format='%Y%m%d')
mat$end <- as.Date(as.character(mat$end), format='%Y%m%d')
mat$text <- paste(format(mat$start, format='%b %d')," to ",format(mat$end, format='%b %d'))
# Calculate anxiety score
survey$ax1_1.1 <- 5-survey$ax1_1; survey$ax1_2.1 <- 5-survey$ax1_2; survey$ax1_3.1 <- 5-survey$ax1_3
survey$ax1_4.1 <- 5-survey$ax1_4; survey$ax1_5.1 <- 5-survey$ax1_5
survey$ascore <- rowMeans(survey[c("ax1_1.1","ax1_2.1","ax1_3.1","ax1_4.1","ax1_5.1","ax1_6","ax1_7","ax1_8","ax1_9","ax1_10")],na.rm=T)

# Replication1:Plot the anxiety score graph
library(ggplot2)
dat <- aggregate(ascore ~ wave, survey, mean)
dat$ref <- c(1)
dat <- merge(dat,mat,by="wave")
ggplot(data=dat, aes(x=factor(wave), y=ascore,group=ref)) + 
  # geom_line(size=1) +     # Set linetype by sex
  geom_point(size=3) +geom_line()+    # Use larger points
  expand_limits(y=4) +    # Set y range to include 4
  expand_limits(x=14) +    # Set x range to include 14
  xlab("Wave") + ylab("STAI score") + # Set axis labels
  ggtitle("Anxiety Level")  +   # Set title
  geom_text(angle = 45,hjust=-0.05,vjust = -1,aes(label=text))

# Replication 2:worry, absolute & relative susceptibility, severity compared to SARS
dat_bf1 <- aggregate(bf1 ~ wave, survey, mean)#"Absolute susceptibility"
dat_bf2 <- aggregate(bf2 ~ wave, survey, mean)#"Relative susceptibility"
dat_bf4 <- aggregate(bf4 ~ wave, survey, mean)#"Perceived severity compared to SARS"
dat_bf5 <- aggregate(bf5 ~ wave, survey, mean)#"Worry if developed ILI"
dat_bf1$bf1  <- dat_bf1$bf1/max(na.omit(survey$"bf1"))
dat_bf2$bf2  <- dat_bf2$bf2/max(na.omit(survey$"bf2"))
dat_bf4$bf4  <- dat_bf4$bf4/max(na.omit(survey$"bf4"))
dat_bf5$bf5  <- dat_bf5$bf5/max(na.omit(survey$"bf5"))
dat2 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(dat_bf1,dat_bf2,dat_bf4,dat_bf5,mat))
colnames(dat2) <- c("wave","Absolute susceptibility", "Relative susceptibility", "Perceived severity compared to SARS","Worry if developed ILI","start","end","text" )
library(reshape2)
data <- melt(dat2,id=c("wave","text","start","end"))
ggplot(data=data, aes(x=factor(wave), y=value,group=variable,colour=variable)) + 
  geom_line(size=1) +     # Set linetype by sex
  geom_point(size=3) +geom_line()+    # Use larger points
  expand_limits(y=0.8) +    # Set y range to include 4
  xlab("Wave") + ylab("Proportion") + # Set axis labels
  ggtitle("Worry, absolute & relative susceptibility, severity compared to SARS")     # Set title
  
  
  
setwd("C:/Users/User/Dropbox/1. healthcare analytics/Project R")

# Data Preparation
survey <- read.csv("H1N1_surveys.csv") #12965 observations
subset <- survey[survey$wave>2,] # Survey waves in mitigation phase, 10940 observations
# Suvery waves with >=1000 success respondents, 10436 observations
for (waveIndex in 2:13){
  if(sum(survey$wave==waveIndex)<1000)
    subset <- subset[subset$wave!=waveIndex,] 
}
subset$sick <- 1*((subset$sm1_1==1&subset$sm1_5==1)|(subset$sm1_1==1&subset$sm1_9==1)) # New attribute: sick=fever+cough OR fever+sore throat
subset <- subset[subset$sick==0&!is.na(subset$sick),] # Exclude sick/unknown status, 10334 observations
# Data validity check: Convert unauthorized data to NA
set12 <- c(1:2) #1-Yes; 2-No
subset$bf3_1[!is.element(subset$bf3_1,set12)] <- NA # Only Response 1 and 2 are valid for qn bf3_1
subset$bf3_3[!is.element(subset$bf3_3,set12)] <- NA
subset$bf3_4[!is.element(subset$bf3_4,set12)] <- NA
subset$bf3_5[!is.element(subset$bf3_5,set12)] <- NA
subset$bf3_6[!is.element(subset$bf3_6,set12)] <- NA
set14<-c(1:4) #1-Always; 2-Usually; 3-Sometimes; 4-Never
subset$pm2[!is.element(subset$pm2,set14)] <- NA 
subset$pm3[!is.element(subset$pm3,set14)] <- NA
subset$pm3a[!is.element(subset$pm3a,set14)] <- NA
subset$pm4[!is.element(subset$pm4,set14)] <- NA
subset$pm5[!is.element(subset$pm5,set14)] <- NA
subset$pm7[!is.element(subset$pm7,set14)] <- NA
subset$pm7b[!is.element(subset$pm7b,set14)] <- NA
set13<-c(1:3) #1-Yes, due to swine flu; 2-Yes, but not due to swine flu; 3-No
subset$pm10_1[!is.element(subset$pm10_1,set13)] <- NA
subset$pm10_2[!is.element(subset$pm10_2,set13)] <- NA
subset$pm10_3[!is.element(subset$pm10_3,set13)] <- NA
subset$pm10_5[!is.element(subset$pm10_5,set13)] <- NA
subset$pm10_6[!is.element(subset$pm10_6,set13)] <- NA
rm(set12,set13,set14)
# For invalid data, cary out Multiple Imputation using Additive Regression, Bootstrapping, and Predictive Mean Matching
library(Hmisc)
set.seed(12345) # same as the Original R script
subset.i <- aregImpute( ~ I(pm3)+I(pm5)+sex+I(age_gp)+I(edu)+I(ax1_1)+I(ax1_2)+I(ax1_3)+I(ax1_4)+I(ax1_5)+I(ax1_6)+I(ax1_7)+I(ax1_8)+I(ax1_9)+I(ax1_10)
                        +I(ph1)+I(bf1)+I(bf2)+I(pm2)+I(pm3a)+I(pm4)+I(pm7)+I(pm7b)+I(pm10_1)+I(pm10_2)+I(pm10_3)+I(pm10_5)+I(pm10_6)
                        +I(bf5)+I(bf4)+I(bf3_1)+I(bf3_3)+I(bf3_4)+wave,#AsIs for each category
                        data=subset, 
                        n.impute=10) #number of multiple imputations per missing value
subset.nomiss <- list(subset, subset, subset, subset, subset, subset, subset, subset, subset, subset) # List of df
for(i in 1:10){
  subset.nomiss[[i]]$pm3[is.na(subset.nomiss[[i]]$pm3)] <- subset.i$imputed$pm3[,i]
  subset.nomiss[[i]]$pm3a[is.na(subset.nomiss[[i]]$pm3a)] <- subset.i$imputed$pm3a[,i]
  subset.nomiss[[i]]$pm4[is.na(subset.nomiss[[i]]$pm4)] <- subset.i$imputed$pm4[,i]
  subset.nomiss[[i]]$pm5[is.na(subset.nomiss[[i]]$pm5)] <- subset.i$imputed$pm5[,i]
  subset.nomiss[[i]]$pm7[is.na(subset.nomiss[[i]]$pm7)] <- subset.i$imputed$pm7[,i]
  subset.nomiss[[i]]$pm7b[is.na(subset.nomiss[[i]]$pm7b)] <- subset.i$imputed$pm7b[,i]
  subset.nomiss[[i]]$sex[is.na(subset.nomiss[[i]]$sex)] <- subset.i$imputed$sex[,i]
  subset.nomiss[[i]]$age_gp[is.na(subset.nomiss[[i]]$age_gp)] <- subset.i$imputed$age_gp[,i]
  subset.nomiss[[i]]$edu[is.na(subset.nomiss[[i]]$edu)] <- subset.i$imputed$edu[,i]
  subset.nomiss[[i]]$ax1_1[is.na(subset.nomiss[[i]]$ax1_1)] <- subset.i$imputed$ax1_1[,i]
  subset.nomiss[[i]]$ax1_2[is.na(subset.nomiss[[i]]$ax1_2)] <- subset.i$imputed$ax1_2[,i]
  subset.nomiss[[i]]$ax1_3[is.na(subset.nomiss[[i]]$ax1_3)] <- subset.i$imputed$ax1_3[,i]
  subset.nomiss[[i]]$ax1_4[is.na(subset.nomiss[[i]]$ax1_4)] <- subset.i$imputed$ax1_4[,i]
  subset.nomiss[[i]]$ax1_5[is.na(subset.nomiss[[i]]$ax1_5)] <- subset.i$imputed$ax1_5[,i]
  subset.nomiss[[i]]$ax1_6[is.na(subset.nomiss[[i]]$ax1_6)] <- subset.i$imputed$ax1_6[,i]
  subset.nomiss[[i]]$ax1_7[is.na(subset.nomiss[[i]]$ax1_7)] <- subset.i$imputed$ax1_7[,i]
  subset.nomiss[[i]]$ax1_8[is.na(subset.nomiss[[i]]$ax1_8)] <- subset.i$imputed$ax1_8[,i]
  subset.nomiss[[i]]$ax1_9[is.na(subset.nomiss[[i]]$ax1_9)] <- subset.i$imputed$ax1_9[,i]
  subset.nomiss[[i]]$ax1_10[is.na(subset.nomiss[[i]]$ax1_10)] <- subset.i$imputed$ax1_10[,i]
  subset.nomiss[[i]]$ph1[is.na(subset.nomiss[[i]]$ph1)] <- subset.i$imputed$ph1[,i]
  subset.nomiss[[i]]$bf1[is.na(subset.nomiss[[i]]$bf1)] <- subset.i$imputed$bf1[,i]
  subset.nomiss[[i]]$bf2[is.na(subset.nomiss[[i]]$bf2)] <- subset.i$imputed$bf2[,i]
  subset.nomiss[[i]]$pm10_1[is.na(subset.nomiss[[i]]$pm10_1)] <- subset.i$imputed$pm10_1[,i]
  subset.nomiss[[i]]$pm10_2[is.na(subset.nomiss[[i]]$pm10_2)] <- subset.i$imputed$pm10_2[,i]
  subset.nomiss[[i]]$pm10_3[is.na(subset.nomiss[[i]]$pm10_3)] <- subset.i$imputed$pm10_3[,i]
  subset.nomiss[[i]]$pm10_5[is.na(subset.nomiss[[i]]$pm10_5)] <- subset.i$imputed$pm10_5[,i]
  subset.nomiss[[i]]$pm10_6[is.na(subset.nomiss[[i]]$pm10_6)] <- subset.i$imputed$pm10_6[,i]
  subset.nomiss[[i]]$bf5[is.na(subset.nomiss[[i]]$bf5)] <- subset.i$imputed$bf5[,i]
  subset.nomiss[[i]]$bf4[is.na(subset.nomiss[[i]]$bf4)] <- subset.i$imputed$bf4[,i]
  subset.nomiss[[i]]$bf3_1[is.na(subset.nomiss[[i]]$bf3_1)] <- subset.i$imputed$bf3_1[,i]
  subset.nomiss[[i]]$bf3_3[is.na(subset.nomiss[[i]]$bf3_3)] <- subset.i$imputed$bf3_3[,i]
  subset.nomiss[[i]]$bf3_4[is.na(subset.nomiss[[i]]$bf3_4)] <- subset.i$imputed$bf3_4[,i]  
}

# set ref X and combine categories for easier interpretation
for (i in 1:10){
  # anxiety, with new calculated attribute ax1_score
  subset.nomiss[[i]]$ax1_1 <- 5-subset.nomiss[[i]]$ax1_1
  subset.nomiss[[i]]$ax1_2 <- 5-subset.nomiss[[i]]$ax1_2
  subset.nomiss[[i]]$ax1_3 <- 5-subset.nomiss[[i]]$ax1_3
  subset.nomiss[[i]]$ax1_4 <- 5-subset.nomiss[[i]]$ax1_4
  subset.nomiss[[i]]$ax1_5 <- 5-subset.nomiss[[i]]$ax1_5  
  subset.nomiss[[i]]$ax1_score <- (subset.nomiss[[i]]$ax1_1+subset.nomiss[[i]]$ax1_2+subset.nomiss[[i]]$ax1_3+subset.nomiss[[i]]$ax1_4+subset.nomiss[[i]]$ax1_5+subset.nomiss[[i]]$ax1_6+subset.nomiss[[i]]$ax1_7+subset.nomiss[[i]]$ax1_8+subset.nomiss[[i]]$ax1_9+subset.nomiss[[i]]$ax1_10)/10  
  subset.nomiss[[i]]$ax1_score[subset.nomiss[[i]]$ax1_score<2] <- 1
  subset.nomiss[[i]]$ax1_score[subset.nomiss[[i]]$ax1_score>=2&subset.nomiss[[i]]$ax1_score<2.5] <- 2
  subset.nomiss[[i]]$ax1_score[subset.nomiss[[i]]$ax1_score>=2.5&subset.nomiss[[i]]$ax1_score<=4] <- 3  
  # sex ref:M
  subset.nomiss[[i]]$sex[subset.nomiss[[i]]$sex==1] <- 0
  subset.nomiss[[i]]$sex[subset.nomiss[[i]]$sex==2] <- 1
  # age ref: 35-44
  subset.nomiss[[i]]$age_gp[subset.nomiss[[i]]$age_gp==3] <- 0  
  # education ref:primary, grouping 1,2,3;4,5,6;7
  subset.nomiss[[i]]$edu[subset.nomiss[[i]]$edu<=3] <- 1
  subset.nomiss[[i]]$edu[subset.nomiss[[i]]$edu<=6&subset.nomiss[[i]]$edu>=4] <- 2
  subset.nomiss[[i]]$edu[subset.nomiss[[i]]$edu==7] <- 3  
  # perceived health
  subset.nomiss[[i]]$ph1[subset.nomiss[[i]]$ph1==3] <- 0  
  # absolute susceptibility, grouping 6,7
  subset.nomiss[[i]]$bf1[subset.nomiss[[i]]$bf1==4] <- 0
  subset.nomiss[[i]]$bf1[subset.nomiss[[i]]$bf1==7] <- 6  
  # relative susceptibility, grouping 6,7
  subset.nomiss[[i]]$bf2[subset.nomiss[[i]]$bf2==4] <- 0
  subset.nomiss[[i]]$bf2[subset.nomiss[[i]]$bf2==7] <- 6  
  # severity vs SARS
  subset.nomiss[[i]]$bf5[subset.nomiss[[i]]$bf5==3] <- 0  
  # how worry if developed ILI tomorrow
  subset.nomiss[[i]]$bf4[subset.nomiss[[i]]$bf4==4] <- 0  
  # knowledge
  subset.nomiss[[i]]$bf3_1[subset.nomiss[[i]]$bf3_1==2] <- 0
  subset.nomiss[[i]]$bf3_3[subset.nomiss[[i]]$bf3_3==2] <- 0
  subset.nomiss[[i]]$bf3_4[subset.nomiss[[i]]$bf3_4==2] <- 0    
  # handwashing after sneezing
  subset.nomiss[[i]]$pm3[subset.nomiss[[i]]$pm3<=2] <- 1
  subset.nomiss[[i]]$pm3[subset.nomiss[[i]]$pm3>=3&subset.nomiss[[i]]$pm3<=4] <- 0  
  # handwashing-use liquid soup
  subset.nomiss[[i]]$pm4[subset.nomiss[[i]]$pm4<=2] <- 1
  subset.nomiss[[i]]$pm4[subset.nomiss[[i]]$pm4>=3&subset.nomiss[[i]]$pm4<=4] <- 0  
  # handwashing after home
  subset.nomiss[[i]]$pm3a[subset.nomiss[[i]]$pm3a<=2] <- 1
  subset.nomiss[[i]]$pm3a[subset.nomiss[[i]]$pm3a>=3&subset.nomiss[[i]]$pm3a<=4] <- 0 
  # adopt any preventive measures when touching common objects (pm7)
  subset.nomiss[[i]]$pm7[subset.nomiss[[i]]$pm7<=2] <- 1
  subset.nomiss[[i]]$pm7[subset.nomiss[[i]]$pm7>=3&subset.nomiss[[i]]$pm7<=4] <- 0  
  # wash hands after toucning common objects (pm7b)
  subset.nomiss[[i]]$pm7b[subset.nomiss[[i]]$pm7b<=2] <- 1
  subset.nomiss[[i]]$pm7b[subset.nomiss[[i]]$pm7b>=3&subset.nomiss[[i]]$pm7b<=4] <- 0
  
}

# OR and IC calculation for the multiple logistic regression
combine.mi <- function(model, n.impute){
  betas <- matrix(c(model[[1]][[4]]$fixed, model[[2]][[4]]$fixed, model[[3]][[4]]$fixed,
                    model[[4]][[4]]$fixed, model[[5]][[4]]$fixed,model[[6]][[4]]$fixed, model[[7]][[4]]$fixed,
                    model[[8]][[4]]$fixed, model[[9]][[4]]$fixed, model[[10]][[4]]$fixed),
                  byrow=FALSE, ncol=n.impute)  # coefficients
  vars <- matrix(c(diag(model[[1]][[5]]), diag(model[[2]][[5]]), diag(model[[3]][[5]]),
                   diag(model[[4]][[5]]), diag(model[[5]][[5]]),diag(model[[6]][[5]]), diag(model[[7]][[5]]),
                   diag(model[[8]][[5]]), diag(model[[9]][[5]]), diag(model[[10]][[5]])),
                 byrow=FALSE, ncol=n.impute) # variance (diagonal of the variance-covariance matrix of the fixed effects)
  coef.names <- names(model[[1]][[4]]$fixed)
  mean.coefs <- rowMeans(betas) # mean of coefficients for all 10 imputed data without NA
  Ubar <- rowMeans(vars) #  mean of variance 
  B <- rowSums((betas-mean.coefs)^2 /(n.impute-1)) #sigma square
  T <- (1 + 1/n.impute)*B+Ubar # s square
  degf <- (n.impute-1)*(1+Ubar/((1+1/n.impute)*B))*(1+Ubar/((1+1/n.impute)*B)) # degree of freedom
  data.frame(OR = exp(mean.coefs), # odds ratio for a change in X (Xi vs ref )
             lowerCI = exp(mean.coefs-qt(0.975, df=degf)*sqrt(T)) #CI, alpha=0.05
             ,upperCI = exp(mean.coefs + qt(0.975, df=degf)*sqrt(T))
             ,row.names=coef.names)
}
# Multiple logistic regression for 5 interested attributes in the Table Replicated
# handwashing after sneezing, coughing or touching nose
fit.m <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
library(MASS)
for (i in 1:10){ #Fit Generalized Linear Mixed Models via PQL
  fit.m[[i]] <- glmmPQL(pm3~sex+factor(age_gp)+factor(edu)+factor(ax1_score)+
                          factor(ph1)+factor(bf1)+factor(bf2)+factor(bf5)+factor(bf4)+bf3_1+bf3_3+bf3_4,
                        random = ~ 1 | wave, family = binomial,data=subset.nomiss[[i]])
}
round(combine.mi(fit.m,10),2) # return OR and CI
# use liquid soup when washing hands
fit.m <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
for (i in 1:10){
  fit.m[[i]] <- glmmPQL(pm4~sex+factor(age_gp)+factor(edu)+factor(ax1_score)+
                          factor(ph1)+factor(bf1)+factor(bf2)+factor(bf5)+factor(bf4)+bf3_1+bf3_3+bf3_4,
                        random = ~ 1 | wave, #Random effects
                        family = binomial,data=subset.nomiss[[i]])
}
round(combine.mi(fit.m,10),2)
# wash hands after returning home
fit.m <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
for (i in 1:10){
  fit.m[[i]] <- glmmPQL(pm3a~sex+factor(age_gp)+factor(edu)+factor(ax1_score)+
                          factor(ph1)+factor(bf1)+factor(bf2)+factor(bf5)+factor(bf4)+bf3_1+bf3_3+bf3_4,
                        random = ~ 1 | wave, family = binomial,data=subset.nomiss[[i]])
}
round(combine.mi(fit.m,10),2)
# wash hands after toucning common objects (pm7b)
fit.m <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
for (i in 1:10){
  fit.m[[i]] <- glmmPQL(pm7b~sex+factor(age_gp)+factor(edu)+factor(ax1_score)+
                          factor(ph1)+factor(bf1)+factor(bf2)+factor(bf5)+factor(bf4)+bf3_1+bf3_3+bf3_4,
                        random = ~ 1 | wave, family = binomial,data=subset.nomiss[[i]])
}
round(combine.mi(fit.m,10),2)
# clean or disinfect house more often
fit.m <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
for (i in 1:10){
  fit.m[[i]] <- glmmPQL(pm10_6~sex+factor(age_gp)+factor(edu)+factor(ax1_score)+
                          factor(ph1)+factor(bf1)+factor(bf2)+factor(bf5)+factor(bf4)+bf3_1+bf3_3+bf3_4,
                        random = ~ 1 | wave, family = binomial,data=subset.nomiss[[i]])
}
round(combine.mi(fit.m,10),2)

#Get the number of successful respondents in each category in the actual survey
library(plyr)
# Gender
ddply(subset,.(subset$sex),nrow)
#age group
ddply(subset,.(subset$age),nrow)
#educational level
ddply(subset,.(subset$edu),nrow)
(sum(284,332,1185)) #primary and below
(sum(1487,2864,966)) #Secondary
# ax1_score, with calculation
subset2 <- subset
subset2$ax1_score <- (25-subset2$ax1_1-subset2$ax1_2-subset2$ax1_3-subset2$ax1_4-subset2$ax1_5+subset2$ax1_6+subset2$ax1_7+subset2$ax1_8+subset2$ax1_9+subset2$ax1_10)/10
subset2$ax1_score[subset2$ax1_score<2] <- 1
subset2$ax1_score[subset2$ax1_score>=2&subset2$ax1_score<2.5] <- 2
subset2$ax1_score[subset2$ax1_score>=2.5&subset2$ax1_score<=4] <- 3
# perceived health
#1-Excellent; 2-Very good; 3-Good; 4-Fair; 5-Poor
ddply(subset,.(subset$ph1),nrow)
# absolute susceptibility
#1-Never; 2-Very unlikely; 3-Unlikely; 4-Evens; 5-Likely; 6-Very likely; 7-Certain.
ddply(subset,.(subset$bf1),nrow)
sum(64,15) #combine 6 and 7 together
# relative susceptibility
#1-Not at all; 2-Much less; 3- Less; 4-Evens; 5-More; 6-Much more; 7-Certain.
ddply(subset,.(subset$bf2),nrow)
sum(44,44) #combine 6 and 7 together
# severity vs SARS
#1-Much less; 2-Less; 3-About the same; 4-More; 5-Much more.
ddply(subset,.(subset$bf5),nrow)
# how worry if developed ILI tomorrow
#1-Not at all worried; 2-Much less worried than normal; 3-Worried less than normal; 4-About same; 5-Worried more than normal; 6-Worried much more than normal; 7-Extremely worried.
ddply(subset,.(subset$bf4),nrow)
# knowledge
ddply(subset,.(subset$bf3_1),nrow) #droplets
ddply(subset,.(subset$bf3_3),nrow) #indirect hand contact
ddply(subset,.(subset$bf3_4),nrow) #oral-faecal