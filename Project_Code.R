# importing required libraries
library(formula.tools)
library(nlme)
library(ggplot2)
library(clubSandwich)

# reading data
lead <- read.table("lead.full.txt", header = F)
colnames(lead) = c("id", "ind.age", "sex", "week", "blood", "trt")
head(lead)

# defining variables
id <- lead$id
Y <- lead$blood
t <- lead$week
sex <-lead$sex
ind.age <-lead$ind.age
trt <-lead$trt
Trt1 = as.numeric(lead$trt==1)
Trt2 = as.numeric(lead$trt==2)
Trt3 = as.numeric(lead$trt==3)
tfact <- as.numeric( factor(t, labels = 1:5))
lead_new <- data.frame(id,t,sex,ind.age,Trt1,Trt2,Trt3,tfact,Y)

# fitting different lme models
meanform <- Y ~ -1 + Trt1 + Trt1:sex + Trt1:ind.age + Trt1:t + Trt1:sex:ind.age + Trt1:sex:t + Trt1:ind.age:t + Trt1:sex:ind.age:t +
  Trt2 + Trt2:sex + Trt2:ind.age + Trt2:t + Trt2:sex:ind.age + Trt2:sex:t + Trt2:ind.age:t + Trt2:sex:ind.age:t +
  Trt3 + Trt3:sex + Trt3:ind.age + Trt3:t + Trt3:sex:ind.age + Trt3:sex:t + Trt3:ind.age:t + Trt3:sex:ind.age:t


fit1 = lme(fixed = meanform
          ,random =  ~ t|id
          ,method = "ML"
          ,control = lmeControl(opt='optim')
          ,data=lead_new)
AIC1 = AIC(fit1)
BIC1 = BIC(fit1)

fit2 = lme(fixed = meanform
           ,random =  ~ t|id
           ,method = "ML"
           ,control = lmeControl(opt='optim')
           ,weights = varIdent(form = ~ 1 | t)
           ,data=lead_new) 
AIC2 = AIC(fit2)
BIC2 = BIC(fit2)

fit3 = lme(fixed = meanform
           ,random =  ~ t|id
           ,method = "ML"
           ,control = lmeControl(opt='optim')
           ,correlation = corAR1(form = ~ t | id)
           ,data=lead_new)
AIC3 = AIC(fit3)
BIC3 = BIC(fit3)

fit4 = lme(fixed = meanform
           ,random =  ~ t|id
           ,method = "ML"
           ,control = lmeControl(opt='optim')
           ,correlation = corAR1(form = ~ t | id)
           ,weights = varIdent(form = ~ 1 | t)
           ,data=lead_new) 
AIC4 = AIC(fit4)
BIC4 = BIC(fit4)

fit5 = lme(fixed = meanform
           ,random =  ~ t|id
           ,method = "ML"
           ,control = lmeControl(opt='optim')
           ,correlation = corSymm(form = ~ tfact | id)
           ,data=lead_new)
AIC5 = AIC(fit5)
BIC5 = BIC(fit5)

fit6 = lme(fixed = meanform
           ,random =  ~ t|id
           ,method = "ML"
           ,control = lmeControl(opt='optim')
           ,correlation = corSymm(form = ~ tfact | id)
           ,weights = varIdent(form = ~ 1 | tfact)
           ,data=lead_new) 
AIC6 = AIC(fit6)
BIC6 = BIC(fit6)

table = data.frame(cbind(c(AIC1,AIC2,AIC3,AIC4,AIC5,AIC6)
                         ,c(BIC1,BIC2,BIC3,BIC4,BIC5,BIC6)))
colnames(table) = c("AIC","BIC")
rownames(table) = c("       Independent w/ equal error variances (fit1)"
                    ," Independent w/ diff. error var. over weeks (fit2)"
                    ,"             AR(1) w/ equal error variances (fit3)"
                    ,"       AR(1) w/ diff. error var. over weeks (fit4)"
                    ,"      Unstructured w/ equal error variances (fit5)"
                    ,"Unstructured w/ diff. error var. over weeks (fit6)")
table

# choosing bestfit from AIC/BIC table
bestfit = fit1
anova.lme(bestfit)

# testing sex main effect and interactions
t.test.reg(bestfit,get.id.contain(bestfit,"sex")) # t-test for all 12 sex coefficients
L = create_L(bestfit,get.id.contain(bestfit,"sex")) # create L vector for age terms
hypo_test(L, bestfit, "Wald")   # combined Wald test

# testing age main effect and interactions
t.test.reg(bestfit,get.id.contain(bestfit,"age")) # t-test for all 12 age coefficients
L = create_L(bestfit,get.id.contain(bestfit,"age")) # create L vector for age terms
hypo_test(L, bestfit, "Wald")   # combined Wald test

# removing sex main effect and interactions
meanform2 <- Y ~ -1 + Trt1 + Trt1:ind.age + Trt1:t + Trt1:ind.age:t +
  Trt2 + Trt2:ind.age + Trt2:t + Trt2:ind.age:t +
  Trt3 + Trt3:ind.age + Trt3:t + Trt3:ind.age:t 
bestfit_reduced = lme(fixed = meanform2
           ,random =  ~ t|id
           ,method = "ML"
           ,control = lmeControl(opt='optim')
           ,data=lead_new)
anova.lme(bestfit_reduced)

# testing for ind.age:t and derivative interactions
t.test.reg(bestfit_reduced,get.id.contain(bestfit_reduced,"ind.age:t"))
L = create_L(bestfit_reduced,get.id.contain(bestfit_reduced,"ind.age:t"))
hypo_test(L, bestfit_reduced, "F-test") 
hypo_test(L, bestfit_reduced, "Wald")

# removing sex main effect and interactions
# and ind.age:t and derivative interactions
meanform3 <- Y ~ -1 + Trt1 + Trt1:ind.age + Trt1:t +
  Trt2 + Trt2:ind.age + Trt2:t +
  Trt3 + Trt3:ind.age + Trt3:t  

# reduced model
meanform3 <- Y ~ ind.age + Trt1:t + Trt2:t + Trt3:t  
bestfit_reduced2 = lme(fixed = meanform3
                      ,random =  ~ t|id
                      ,method = "ML"
                      ,control = lmeControl(opt='optim')
                      ,data=lead_new)

# checked if reduce model is sufficient
anova.lme(bestfit,bestfit_reduced2)


# checking whether mean trends for treatments are same in reduced model
L1 <- c(1,0,-1,0,0,0,0,0,0)
L2 <- c(0,0,0,1,0,0,0,-1,0)
L3 <- c(0,0,0,0,1,0,0,0,-1)
anova.lme(bestfit_reduced2, L=L1,adjustSigma = TRUE)
anova.lme(bestfit_reduced2, L=L2,adjustSigma = TRUE)
anova.lme(bestfit_reduced2, L=L3,adjustSigma = TRUE)


L <- rep(0,length(fixed.effects(bestfit)))
L[c(5,8,11)] <- c(1,0,-1)
anova.lme(bestfit,L=L)

############################ Plotting Code #########################
### plotting actual vs. predicted
actual <- lead_new$Y
predicted <- predict(bestfit_reduced2,level=1)
png("predvsact.png",width=5,height=5,units="in",res=1200)
plot(predicted,actual,pch=19,xlim=c(0,max(actual,predicted))
     ,ylim=c(0,max(actual,predicted)),col='#666666'
     ,xlab="Predicted Values",ylab="Actual Response",main="Actual Response v/s Predicted Values")
plot_fit <- lm(predicted~actual)
abline(0,1,lwd=2.5,lty=2,col='blue')
dev.off()


### plotting residuals 
res <- resid(bestfit_reduced2, level = 1, type = "normalized")

index <- 1:nrow(lead_new)
index_plac <- index[lead_new$Trt1==1]
index_sucl <- index[lead_new$Trt2==1]
index_such <- index[lead_new$Trt3==1]

rand_intercept <- random.effects(bestfit_reduced2)

png('randomeffects.png',width=6,height=5.5,units="in",res=1200)
par(mfrow=c(2,2))
plot(rand_intercept[1:40,1],rand_intercept[1:40,2],pch=19
     ,xlab='Intercept',ylab='week',main='Treatment: Placebo'
     ,xlim=c(-2,2),ylim=c(-2,2))
plot(rand_intercept[41:80,1],rand_intercept[41:80,2],pch=19
     ,xlab='Intercept',ylab='week',main='Treatment: Succimer Low'
     ,xlim=c(-2,2),ylim=c(-2,2))
plot(rand_intercept[81:120,1],rand_intercept[81:120,2],pch=19
     ,xlab='Intercept',ylab='week',main='Treatment: Succimer High'
     ,xlim=c(-2,2),ylim=c(-2,2))
dev.off()

png('randomqq.png',width=9,height=6,units="in",res=800)
par(mfrow=c(2,3))
qqnorm(rand_intercept[1:40,1],pch=19,main='Q-Q plot, Intercept random effect, Placebo')
qqline(rand_intercept[1:40,1],pch=19)
qqnorm(rand_intercept[1:40,2],pch=19,main='Q-Q plot, week random effect, Placebo')
qqline(rand_intercept[1:40,2],pch=19)
qqnorm(rand_intercept[41:80,1],pch=19,main='Q-Q plot, Intercept random effect, Succ. Low')
qqline(rand_intercept[41:80,1],pch=19)
qqnorm(rand_intercept[41:80,2],pch=19,main='Q-Q plot, week random effect, Succ. Low')
qqline(rand_intercept[41:80,2],pch=19)
qqnorm(rand_intercept[81:120,1],pch=19,main='Q-Q plot, Intercept random effect, Succ. High')
qqline(rand_intercept[81:120,1],pch=19)
qqnorm(rand_intercept[81:120,2],pch=19,main='Q-Q plot, week random effect, Succ. High')
qqline(rand_intercept[81:120,2],pch=19)
dev.off()

age_levels <- c("Age <= 24 months","Age > 24 months")
treat_levels <- c("Placebo","Low Succimer","High Succimer")
sex_levels <- c("Female","Male")
l=1
for (i in 1:length(age_levels)){
  for (j in 1:length(sex_levels)){
    for (k in 1:length(treat_levels)){
      name = paste0("resid",l,".png")
      l=l+1
      if (k==1){
        index_new <- index[lead_new$Trt1==1 & lead_new$ind.age==i-1 & lead_new$sex==j-1]
      }
      if (k==2){
        index_new <- index[lead_new$Trt2==1 & lead_new$ind.age==i-1 & lead_new$sex==j-1]
      }
      if (k==3){
        index_new <- index[lead_new$Trt3==1 & lead_new$ind.age==i-1 & lead_new$sex==j-1]
      }
      title = paste0(age_levels[i],", ",treat_levels[k],", ",sex_levels[j])
      save_resid_plot(res[index_new],name,title,0.1)
      
    }
  }
}

### plotting qq plots
l=1
for (i in 1:length(age_levels)){
  for (j in 1:length(sex_levels)){
    for (k in 1:length(treat_levels)){
      name = paste0("qqplot",l,".png")
      l=l+1
      if (k==1){
        index_new <- index[lead_new$Trt1==1 & lead_new$ind.age==i-1 & lead_new$sex==j-1]
      }
      if (k==2){
        index_new <- index[lead_new$Trt2==1 & lead_new$ind.age==i-1 & lead_new$sex==j-1]
      }
      if (k==3){
        index_new <- index[lead_new$Trt3==1 & lead_new$ind.age==i-1 & lead_new$sex==j-1]
      }
      title = paste0(age_levels[i],", ",treat_levels[k],", ",sex_levels[j])
      save_qq_plot(res[index_new],name,title)
    }
  }
}

############### Plotting for comparing treatments for every combo###################################
index_plac <- index[lead_new$Trt1==1 & lead_new$ind.age==1 & lead_new$sex==0]
index_sucl <- index[lead_new$Trt2==1 & lead_new$ind.age==1 & lead_new$sex==0]
index_such <- index[lead_new$Trt3==1 & lead_new$ind.age==1 & lead_new$sex==0]

x <- c(0,2,4,6,8)
y_plac <- treat_means(lead_new,index_plac)
y_sucl <- treat_means(lead_new,index_sucl)
y_such <- treat_means(lead_new,index_such)

png("treatment_plot4.png",width=6,height=5.5,units="in",res=1200)

plot(lead_new$t[index_plac],lead_new$Y[index_plac],pch=16
     ,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=1.5,ylab = 'Blood Level'
     ,xlab = "Week", main = "Mean plot : Age > 24 months, Female")
points(lead_new$t[index_sucl],lead_new$Y[index_sucl],pch=16
       ,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=1.5)
points(lead_new$t[index_such],lead_new$Y[index_such],pch=16
       ,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=1.5)
lines(x=x,y=y_plac,col="blue",type='b',pch=16,lty=1,lwd=3,cex=1.5)
lines(x=x,y=y_sucl,col="red",type='b',pch=16,lty=1,lwd=3,cex=1.5)
lines(x=x,y=y_such,col="darkgreen",type='b',pch=16,lty=1,lwd=3,cex=1.5)
legend(x=5,y=41,legend=c("Placebo","Sucimmer-Low","Succimer-High"), 
       lwd=c(3,3,3), col=c("blue","red","darkgreen"), pch=c(19,19,19))

dev.off()
############### Plotting for comparing age/sex for every combo###################################
index_age0f <- index[lead_new$Trt3==1 & lead_new$ind.age==0 & lead_new$sex==0]
index_age1f <- index[lead_new$Trt3==1 & lead_new$ind.age==1 & lead_new$sex==0]
index_age0m <- index[lead_new$Trt3==1 & lead_new$ind.age==0 & lead_new$sex==1]
index_age1m <- index[lead_new$Trt3==1 & lead_new$ind.age==1 & lead_new$sex==1]

x <- c(0,2,4,6,8)
y_age0f <- treat_means(lead_new,index_age0f)
y_age1f <- treat_means(lead_new,index_age1f)
y_age0m <- treat_means(lead_new,index_age0m)
y_age1m <- treat_means(lead_new,index_age1m)

png("agesex_plot3.png",width=6,height=5.5,units="in",res=1200)

plot(lead_new$t[index_age0f],lead_new$Y[index_age0f],pch=16
     ,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=1.5,ylab = 'Blood Level'
     ,xlab = "Week", main = "Mean plot : Succimer High Dose", ylim=c(10,40))
points(lead_new$t[index_age1f],lead_new$Y[index_age1f],pch=16
       ,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=1.5)
points(lead_new$t[index_age0m],lead_new$Y[index_age0m],pch=16
       ,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=1.5)
points(lead_new$t[index_age1m],lead_new$Y[index_age1m],pch=16
       ,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=1.5)
lines(x=x,y=y_age0f,col="blue",type='b',pch=16,lty=1,lwd=3,cex=1.5)
lines(x=x,y=y_age1f,col="red",type='b',pch=16,lty=1,lwd=3,cex=1.5)
lines(x=x,y=y_age0m,col="blue",type='b',pch=16,lty=2,lwd=3,cex=1.5)
lines(x=x,y=y_age1m,col="red",type='b',pch=16,lty=2,lwd=3,cex=1.5)
legend(x=3.2,y=40,legend=c("Female, Age <= 24 months","Female, Age > 24 months"
                             ,"Male, Age <= 24 months","Male, Age > 24 months"), 
       lwd=c(3,3,3,3), col=c("blue","red","blue","red"), lty=c(1,1,2,2))

dev.off()
############### Plotting for comparing treatments###################################
index_plac <- index[lead_new$Trt1==1]
index_sucl <- index[lead_new$Trt2==1]
index_such <- index[lead_new$Trt3==1]

x <- c(0,2,4,6,8)
y_plac <- treat_means(lead_new,index_plac)
y_sucl <- treat_means(lead_new,index_sucl)
y_such <- treat_means(lead_new,index_such)

png("treatment_main.png",width=6,height=5.5,units="in",res=1200)

plot(lead_new$t[index_plac],lead_new$Y[index_plac],pch=16
     ,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=1.5,ylab = 'Blood Level'
     ,xlab = "Week", main = "Mean plot of Treatments", ylim=c(15,35))
points(lead_new$t[index_sucl],lead_new$Y[index_sucl],pch=16
       ,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=1.5)
points(lead_new$t[index_such],lead_new$Y[index_such],pch=16
       ,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=1.5)
lines(x=x,y=y_plac,col="blue",type='b',pch=16,lty=1,lwd=3,cex=1.5)
lines(x=x,y=y_sucl,col="red",type='b',pch=16,lty=1,lwd=3,cex=1.5)
lines(x=x,y=y_such,col="darkgreen",type='b',pch=16,lty=1,lwd=3,cex=1.5)
legend(x=5,y=35,legend=c("Placebo","Sucimmer-Low","Succimer-High"), 
       lwd=c(3,3,3), col=c("blue","red","darkgreen"), pch=c(19,19,19))

dev.off()


############### Plotting for comparing age ###################################
index_age0 <- index[lead_new$ind.age==0]
index_age1 <- index[lead_new$ind.age==1]
index_sex0 <- index[lead_new$sex==0]
index_sex1 <- index[lead_new$sex==1]

x <- c(0,2,4,6,8)
y_age0 <- treat_means(lead_new,index_age0)
y_age1 <- treat_means(lead_new,index_age1)
y_sex0 <- treat_means(lead_new,index_sex0)
y_sex1 <- treat_means(lead_new,index_sex1)

png("sex_main.png",width=6,height=5.5,units="in",res=1200)

plot(lead_new$t[index_sex0],lead_new$Y[index_sex0],pch=16
     ,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=1.5,ylab = 'Blood Level'
     ,xlab = "Week", main = "Mean plot of Sex groups", ylim=c(10,40))
points(lead_new$t[index_sex1],lead_new$Y[index_sex1],pch=16
       ,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=1.5)
#points(lead_new$t[index_sex0],lead_new$Y[index_sex0],pch=16
       #,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=1.5)
#points(lead_new$t[index_sex1],lead_new$Y[index_sex1],pch=16
       #,col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),cex=1.5)
lines(x=x,y=y_sex0,col="blue",type='b',pch=16,lty=1,lwd=3,cex=1.5)
lines(x=x,y=y_sex1,col="red",type='b',pch=16,lty=1,lwd=3,cex=1.5)
#lines(x=x,y=y_sex0,col="blue",type='b',pch=16,lty=2,lwd=3,cex=1.5)
#lines(x=x,y=y_sex1,col="red",type='b',pch=16,lty=2,lwd=3,cex=1.5)
legend(x=5.5,y=40,legend=c("Female","Male"), 
       lwd=c(3,3), col=c("blue","red"), lty=c(1,1))

dev.off()
############### Plotting for comparing treatments###################################
index_plac <- index[lead_new$Trt1==1]
index_sucl <- index[lead_new$Trt2==1]
index_such <- index[lead_new$Trt3==1]

x <- c(0,2,4,6,8)
plac_mat <- plot_mat(lead_new,index_plac)
sucl_mat <- plot_mat(lead_new,index_sucl)
such_mat <- plot_mat(lead_new,index_such)


pdf("succimerhigh_ind.pdf",width=6,height=5.5)

data_plot <- such_mat
plot(x=x,y=data_plot[1,],pch=16
     ,col=rgb(red = 1, green = 1, blue = 1, alpha = 0.1),cex=1.5,ylab = 'Blood Level'
     ,xlab = "Week", main = "Individual Trajectories - Succimer High", xlim=c(-0.2,8.2),ylim=c(0,50))
for (i in 1:nrow(data_plot)){
  temp =c()
  for (j in 1:ncol(data_plot)){
    if (is.na(data_plot[i,j])){
      temp <- append(temp,x[j])
    }
  }
  x_new <- x[!x %in% temp]
  y_new <- data_plot[i,][!data_plot[i,]%in% NA] 
  lines(x=x_new,y= y_new ,col=rgb(red=0,green=0.4,blue=0,alpha=0.3)
        ,type='b',pch=16,lty=1,lwd=1,cex=1.5)
}
y_mean <- colMeans(data_plot,na.rm = TRUE)
lines(x=x,y= y_mean ,col='black'
      ,type='b',pch=16,lty=1,lwd=3,cex=1.5)
legend(x=0,y=50,legend=c("Individual Trajectories","Mean Trajectory"), 
       lwd=c(1,3), col=c("darkgreen","black"), lty=c(1,1))
text(x=x,y= y_mean, labels=round(y_mean,2), cex=1.2, font=2, pos=1)

dev.off()
############### Plotting for comparing treatments###################################
index_plac <- index[lead_new$Trt1==1]
index_sucl <- index[lead_new$Trt2==1]
index_such <- index[lead_new$Trt3==1]

x <- c(0,2,4,6,8)
plac_mat <- plot_mat(lead_new,index_plac)
sucl_mat <- plot_mat(lead_new,index_sucl)
such_mat <- plot_mat(lead_new,index_such)

plac_sd <- apply(plac_mat, 2, sd, na.rm=TRUE)
sucl_sd <- apply(sucl_mat, 2, sd, na.rm=TRUE)
such_sd <- apply(such_mat, 2, sd, na.rm=TRUE)

png("sd_profiles.png",width=6,height=5.5,units="in",res=800)


plot(x=x,y=plac_sd,pch=16
     ,col=rgb(red = 1, green = 1, blue = 1, alpha = 0.1),cex=1.5,ylab = 'Std. Dev. of Blood Level'
     ,xlab = "Week", main = "Standard Deviation Profiles", xlim=c(-0.2,8.2),ylim=c(4,10))
lines(x=x,y= plac_sd ,col='blue'
      ,type='b',pch=16,lty=1,lwd=3,cex=1.5)
lines(x=x,y= sucl_sd ,col='red'
      ,type='b',pch=16,lty=1,lwd=3,cex=1.5)
lines(x=x,y= such_sd ,col='darkgreen'
      ,type='b',pch=16,lty=1,lwd=3,cex=1.5)
legend(x=0,y=10,legend=c("Placebo","Succimer Low","Succimer High"), 
       lwd=c(3,3,3), col=c("blue","red","darkgreen"), lty=c(1,1,1))
text(x=x,y= y_mean, labels=round(y_mean,2), cex=1.2, font=2, pos=1)

dev.off()

############### Plotting measurements frequency ##############
count0 <- nrow(lead_new[lead_new$t==0,])
count2 <- nrow(lead_new[lead_new$t==2,])
count4 <- nrow(lead_new[lead_new$t==4,])
count6 <- nrow(lead_new[lead_new$t==6,])
count8 <- nrow(lead_new[lead_new$t==8,])
count <- c(count0,count2,count4,count6,count8)

heat <- matrix(0, nrow=length(unique(lead_new$t)), ncol=length(unique(lead_new$id)))
for (k in 1:nrow(lead_new)){
  row <- lead_new[k,]
  heat[row$tfact,row$id]=1
}

rownames(heat) = c(0,2,4,6,8)
colnames(heat) = 1:ncol(heat)

png("missing_heatmap.png",width=6,height=5.5,units="in",res=800)
par(mfrow=c(1,1),oma = c(0.5, 0.5, 0.5, 0.5))
aspectHeatmap(heat,Rowv=NA, col=c('yellow','navyblue'),Colv=NA, 
        scale="none",xlab="Individual ID",ylab="Week"
        ,margins=c(3,3),hExp = 1, wExp = 2)
dev.off()

png("missing_freequency.png",width=6,height=5.5,units="in",res=800)
xx <-barplot(count, names.arg=c(0,2,4,6,8), xlab='Week',ylab='Frequency'
             ,main='Frequency of Measurements',border = NA)
text(x = xx, y = count, label = count, pos = 1, cex = 1.5, col = "black")
dev.off()

measure_count <- colSums(heat)
count2 <- c(sum(measure_count==5),
            sum(measure_count==4),
            sum(measure_count==3))
png("measurements_frequency.png",width=6,height=5.5,units="in",res=300)
xx2 <-barplot(count2, names.arg=c(5,4,3), xlab='No. of Measurements',ylab='Frequency'
             ,main='Frequency of No. of Measurements per Individual',border = NA)
text(x = xx2, y = count2, label = count2, pos = 1, cex = 1.5, col = "black")
dev.off()

### supporting functions######################################

get.id.contain <- function(fitobject,string){ 
  return (grep(string,names(fixed.effects(fitobject))))
}

t.test.reg <- function(fitobject,L){
  df <- fitobject$fixDF$X
  df_L <- df[L]
  betahat_L <- fixed.effects(fitobject)[L]
  SE_L <- sqrt( diag(fitobject$varFix)[L] )
  t.stat <- betahat_L/SE_L
  p.value <- round( 2*pt(q = abs(t.stat), df = df_L, lower.tail = FALSE), 4 )
  out_sex <- data.frame(betahat_L, SE_L, df_L,p.value)
  colnames(out_sex) <- c("Coefficients","SE","DF","P-value")
  return (round(out_sex, 3))
}

hypo_test <- function(contrast,lme_object,test) {
  betahat <- fixed.effects(lme_object)
  V.robust <- vcovCR(lme_object, type = "CR0")
  cc <- nrow(L)
  df <- length(eval(lhs(eval(lme_object$call$fixed)))) - length(betahat)
  # estimate and covariance matrix of L\beta
  est <- contrast %*% betahat
  SE <- contrast %*% V.robust %*% t(contrast)
  varmat <- L %*% V.robust %*% t(L)
  if (test == "Wald") {
    # Wald test
    Wald <- c( t(est) %*% solve(varmat) %*% (est) )
    p.value <- pchisq(q = Wald, df = cc, lower.tail=FALSE)
    return(data.frame(Wald, p.value))
  }
  if (test == "F-test") {
    # F-test
    Fstat <- c( t(est) %*% solve(varmat) %*% (est) ) / cc
    p.value <- pf(q = Fstat, df1 = cc, df2 = df, lower.tail=FALSE)
    return(data.frame(Fstat, p.value))
  }
}

create_L <- function(lme_object,vector){
  L = matrix(0L, nrow = length(vector), ncol = length(fixed.effects(lme_object)))
  for (i in 1:length(vector)){
    L[i,vector[i]] = 1
  }
  return(L)
}

save_resid_plot <- function(residual,name,title,margin){
  png(name,width=4.5,height=4.5,units="in",res=1200)
  range = max(residual)-min(residual)
  ul = ceiling(max(residual)+(margin*range))
  ll = floor(min(residual)-(margin*range))
  l = max(abs(ll),abs(ul))
  plot(residual,pch=19,ylim = c(-l,l),ylab="Standardized Residuals",xlab="Fitted Values",main=title)
  dev.off()
  print("Plot saved succesfully!")
}

save_qq_plot <- function(residual,name,title){
  png(name,width=4.5,height=4.5,units="in",res=1200)
  qqnorm(residual,pch=19,main=title)
  qqline(residual)
  dev.off()
  print("Plot saved succesfully!")
}

treat_means <-function(data,index){
  vec <- data$Y[index]
  time <- data$t[index]
  iter <- sort(unique(time))
  result <- c()
  for (i in 1:length(iter)){
    temp <- vec[time==iter[i]]
    result[i] <- mean(temp)
  } 
  return (result)
}

plot_mat <- function(data,index){
  data_iter <- data[index,]
  x <- sort(unique(data_iter$t))
  iter <- sort(unique(data_iter$id))
  matrix <- matrix(NA,nrow=length(iter),ncol=length(x))
  for (i in 1:length(iter)){
    for (j in 1:length(x)){
      replace <- data_iter[data_iter$id == iter[i] & data_iter$t == x[j],]$Y
      if (length(replace) == 1){
        matrix[i,j] <- replace
      }
    }
  }
  return (matrix)
}

################################################################



