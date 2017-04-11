---
title: "landscape and rice field structure"
author: "Dipl.-Biol. Jenny Schellenberg"
date: "22 März 2017"
output: html_document
---

Summary:

# 1: Species richness (S) negatively correlated with local land-use intensity
land-use intensitiy: which predictor to use for that? 
LUI preds

# 2: S is enhanced by increasing landscape heterogeneity (categories "complex"/"simple" are used)
tested in the separated analyses; yes, but strong regional patterns

# 3: S is enhanced by increasing landscape heterogeneity  (area of different land-use types (% in 300m radius) are used)
yes, SDiv is best explaining, but not better explaining than Region
yes, fruit and forest is influencing S but also not better than region

# 5:  The patterns found are independent from the study region
is to be rejected for all analyses yet. No overall patterns. Too strong influences of Region, predictors in Region too distinct.



#################    Preparation      #################################

setwd("G:/Arbeit/Oli Statistik")
setwd("/AAA/Dropbox/WORK/RData_OF/Paper3")

```{r packages and functions}
####### packages ########

library(lme4)
library(lattice)
library(MASS)
library(DHARMa)
library(sjPlot)
library(psych)
library(lmerTest)



##### function to test for overdispersion (from http://glmm.wikidot.com/faq) #######
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}



##### correlation panel function #####
# collinearity test panels
panel.cor.spearman<-function(x,y,digits=2, prefix="", cex.cor, ...)
{
  usr<-par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method="spearman")
  txt<-format(c(r, 0.123456789), digits=digits)[1]
  txt<-paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1/strwidth(txt)
  text(0.5, 0.5, txt, cex=cex.cor*r)
}
panel.cor.pearson<-function(x,y,digits=2, prefix="", cex.cor, ...)
{
  usr<-par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method="pearson"))
  txt<-format(c(r, 0.123456789), digits=digits)[1]
  txt<-paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1.25/strwidth(txt)
  text(0.5, 0.5, txt, cex=cex.cor*r)
}
panel.hist<-function(x,...)
{
  usr<-par("usr"); on.exit(par(usr))
  par(usr=c(usr[1:2], 0, 1.5))
  h<-hist(x, plot=FALSE)
  breaks=h$breaks;nB<-length(breaks)
  y<-h$counts; y<-y/max(y)
  rect(breaks[-nB], 0, breaks[-1],y, col="red", ...)
}

panel.smooth<-function (x, y, col = "black", bg = NA, pch = 18, 
                        cex = 0.8, col.smooth = "red", span = 3/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}



```

```{r data}
### env data "bhead"
bhead<-read.csv("bheadall.txt", header=T, sep="\t")
#Inspect the results
names(bhead)
str(bhead)
summary(bhead)

# select relevant columns
b<-bhead[,c(1:3, 5, 28:54)]    # only choosing landscape predictors

# combine all crop areas (excl. meadows)
b$LU_crops <- b$LU_cereals+b$LU_fruit+b$LU_rice+b$LU_vegetable
# combine all agrarian areas (incl. meadows)
b$LU_agrice <- b$LU_cereals+b$LU_fruit+b$LU_grass+b$LU_vegetable+b$LU_rice
# combine all agrarian areas (excl. rice fields)
b$LU_agrar <- b$LU_cereals+b$LU_fruit+b$LU_grass+b$LU_vegetable


summary(b)
```

far too many predictors... reducing....

```{r collinearity of predictors and setting subdataframes}
##### correlation panel function #####
# collinearity test panels
panel.cor.spearman<-function(x,y,digits=2, prefix="", cex.cor, ...)
{
  usr<-par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method="spearman")
  txt<-format(c(r, 0.123456789), digits=digits)[1]
  txt<-paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1/strwidth(txt)
  text(0.5, 0.5, txt, cex=cex.cor*r)
}
panel.cor.pearson<-function(x,y,digits=2, prefix="", cex.cor, ...)
{
  usr<-par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method="pearson"))
  txt<-format(c(r, 0.123456789), digits=digits)[1]
  txt<-paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1.25/strwidth(txt)
  text(0.5, 0.5, txt, cex=cex.cor*r)
}
panel.hist<-function(x,...)
{
  usr<-par("usr"); on.exit(par(usr))
  par(usr=c(usr[1:2], 0, 1.5))
  h<-hist(x, plot=FALSE)
  breaks=h$breaks;nB<-length(breaks)
  y<-h$counts; y<-y/max(y)
  rect(breaks[-nB], 0, breaks[-1],y, col="red", ...)
}

panel.smooth<-function (x, y, col = "black", bg = NA, pch = 18, 
                        cex = 0.8, col.smooth = "red", span = 3/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}


##### bls: land use composition and rice field structure predictors #####

names(b)
bls<-b[,c(1:16, 27:34)]
bls<-bls[complete.cases(bls),]
names(bls)
summary(bls)


pairs(~S+LU_cereals+LU_forest+LU_fruit+LU_grass+LU_noveg+LU_rice+LU_vegetable+LU_water+LU_crops+LU_agrice+LU_agrar+LS_Rich+LS_SDiv+LS_Even+LS_Dom+MPS_ricepatch+NP_Diss+ED_Diss, data=bls,upper.panel=panel.cor.pearson, diag.panel=panel.hist, lower.panel=panel.smooth, main="colinearity between predictors land use structures", na.action=na.omit)
```



#################### bls: landscape composition and rice field structure################################

```{r bls S in regions}
###### species richness in regions: histograms and density plots #####

frame<-bls
levels(frame$Region)

# create colors for regions: red/orange for PH, blue/bluegreen for VN, 
colour<-c("red", "orangered", "orangered3", "royalblue1", "lightseagreen", "midnightblue", "lightblue")

par(mfrow=c(2,2))
hist(frame$S, breaks=20)
plot(density(frame$S, na.rm=TRUE), main="density only PH", col="grey", lwd=3, xlab="species richness", ylab="density", ylim=c(0,0.15))
for (i in 1:3) { 
  lines(density(frame[which(frame$Region==levels(frame$Region)[i]),]$S, na.rm=TRUE), col=colour[i], lwd=3) }
legend("topright", lwd=3, col=c("grey", colour[1:3]), legend=c("all regions", levels(frame$Region)[1:3]))

plot(density(frame$S, na.rm=TRUE), main="density only VN", col="grey", lwd=3, xlab="species richness", ylab="density", ylim=c(0,0.15))
for (i in 4:7) { lines(density(frame[which(frame$Region==levels(frame$Region)[i]),]$S, na.rm=TRUE), col=colour[i], lwd=3)  }
legend("topright", lwd=3, col=c("grey", colour[4:7]), legend=c("all regions", levels(frame$Region)[4:7])) 

boxplot(frame$S~factor(frame$Region), whisklty=1, ylim=c(0,70), main= "boxplot with sample sizes")
for (i in 1:7) { 
  text( i ,68,labels=length(frame[frame$Region==levels(frame$Region)[i],]$Region)) 
       }
# failing normal distribution
# acceptable for PH_1, PH_2 
# PH_3, VN_1, VN_2, VN_4 borderline... risk of artefact in modeling
# hard risk of artefact in VN_3 - is it really necessary???

# rest ok; sample size low but even 
# variations in -/+ same range 
```

```{r bls predictor overview}
###### is there influence of LocStruc? #######
par(mfrow=c(2,4)) 
boxplot(bls$S~bls$LocStruct, main="all Regions", whisklty=1, pch=19, cex=0.8)    
for (i in 1:7) {
boxplot(bls[bls$Region==levels(bls$Region)[i],]$S~ bls[bls$Region==levels(bls$Region)[i],]$LocStruct, main=levels(bls$Region)[i], whisklty=1, pch=19, cex=0.8)
}
# yes, strong influence in some regions (PH_1, PH_2, PH_3, VN_4)
# other weaker but existing
# strong differences in patterns of influence between regions


###### landscape structures predicitors: histograms and distr. in areas #######
names(bls)
### only LU % area
par(mfrow=c(3,4)) 
for (i in c(5:12,22:24)) {
  hist(bls[,i], main=names(bls)[i]) }
# beware of cereals and vegetables: they are very unbalanced 

# preds in regions
par(mfrow=c(3,4)) 
for (i in c(5:12,22:24)) {
boxplot(bls[,i]~factor(bls$Region), whisklty=1, main=names(bls)[i])
}
# strange differences in area proportions between regions
# not to be analysed because there is risk of artefact:
# cereals

### LS 
par(mfrow=c(2,2)) 
for (i in 13:16) {
  hist(bls[,i], main=names(bls)[i]) }

par(mfrow=c(2,2)) 
for (i in 13:16) {
boxplot(bls[,i]~factor(bls$Region), whisklty=1, main=names(bls)[i])
}

# these 4 predictors are correlated
# taking SDiv only for analysis

```

```{r test of model technique}
# tested with LU_forest  and LS_Sdiv
model<-m1
plotConventionalResiduals(model) 
simR<-simulateResiduals(model)
plotSimulatedResiduals(simR)    

m1<-lmer(S~LU_forest + (1|Region/Landscape), data=bls, REML=F)
m2<-glmer(S~LU_forest + (1|Region/Landscape), data=bls, family=poisson)
overdisp_fun(m2)
m3<-glmer.nb(S~LU_forest + (1|Region/Landscape), data=bls)
# lmer is best

m4<-lmer(S~LS_SDiv + (1|Region/Landscape), data=bls, REML=F)
m5<-glmer(S~LS_SDiv + (1|Region/Landscape), data=bls, family=poisson)
overdisp_fun(m5)
m6<-glmer.nb(S~LS_SDiv + (1|Region/Landscape), data=bls)
# lmer is best

# taking lmer for modeling
```

```{r nonlinear relationships}
# function for diagnostic anova and plots for testing non-linear relationship of predictors to response
polytest<-function(pred,dataframe) {
a1<-lmer(S ~ pred + (1|Region/Landscape), data=dataframe, REML=F)
a2<-lmer(S ~ poly(pred,2) + (1|Region/Landscape), data=dataframe, REML=F)
a3<-lmer(S ~ poly(pred,3) + (1|Region/Landscape), data=dataframe, REML=F)
a4<-lmer(S ~ poly(pred,4) + (1|Region/Landscape), data=dataframe, REML=F) 


plot(pred,dataframe$S, pch=19,cex=0.8, col="grey")
lines(loess.smooth(pred,predict(a1)), col="black", lty=1, lwd=2)
lines(loess.smooth(pred,predict(a2)), col="blue", lty=1, lwd=2) 
lines(loess.smooth(pred,predict(a3)), col="orangered", lty=1, lwd=2) 
lines(loess.smooth(pred,predict(a4)), col="blueviolet", lty=1, lwd=2) 
legend("topright", legend=seq(1,4,1), title="polynomial grade", lty=1, lwd=2, col=c("black", "blue", "orangered", "blueviolet"))

print(anova(a1,a2,a3,a4))
}

par(mfrow=c(1,1))
polytest(bls$LU_cereals,bls)
polytest(bls$LU_forest,bls)
polytest(bls$LU_fruit,bls)
# 3rd grade polynomial, but not considered
polytest(bls$LU_grass,bls)
# quadratic, but not considered
polytest(bls$LU_noveg,bls)
polytest(bls$LU_rice,bls)
polytest(bls$LU_vegetable,bls)
polytest(bls$LU_water,bls)
# quadratic
polytest(bls$LU_crops,bls)
polytest(bls$LU_agrice,bls)
polytest(bls$LU_agrar,bls)
# 3rd grade but not really better
polytest(bls$LS_Rich,bls)
polytest(bls$LS_SDiv,bls)
polytest(bls$LS_Even,bls)
polytest(bls$LS_Dom,bls)
polytest(bls$N_ricepatch,bls)
polytest(bls$MPS_ricepatch,bls)
polytest(bls$NP_Diss,bls)
polytest(bls$ED_Diss,bls)
```

```{r bls forward selection}
# model selection from null model
m0<-lmer(S~1 + (1|Region/Landscape), data=bls, REML=F)

m7<-lmer(S~LU_fruit + (1|Region/Landscape), data=bls, REML=F)
m8<-lmer(S~LU_grass + (1|Region/Landscape), data=bls, REML=F)
m9<-lmer(S~LU_noveg + (1|Region/Landscape), data=bls, REML=F)
m10<-lmer(S~LU_rice + (1|Region/Landscape), data=bls, REML=F)
m11<-lmer(S~LU_vegetable + (1|Region/Landscape), data=bls, REML=F)
m12<-lmer(S~poly(LU_water,2) + (1|Region/Landscape), data=bls, REML=F)
m13<-lmer(S~LS_Rich + (1|Region/Landscape), data=bls, REML=F)
m14<-lmer(S~LS_Even + (1|Region/Landscape), data=bls, REML=F)
m15<-lmer(S~LS_Dom + (1|Region/Landscape), data=bls, REML=F)
m16<-lmer(S~LU_agrice + (1|Region/Landscape), data=bls, REML=F)
m17<-lmer(S~LU_crops + (1|Region/Landscape), data=bls, REML=F)
m18<-lmer(S~LU_agrar + (1|Region/Landscape), data=bls, REML=F)
m19<-lmer(S~N_ricepatch + (1|Region/Landscape), data=bls, REML=F)
m20<-lmer(S~MPS_ricepatch + (1|Region/Landscape), data=bls, REML=F)
m21<-lmer(S~NP_Diss + (1|Region/Landscape), data=bls, REML=F)
m22<-lmer(S~ED_Diss + (1|Region/Landscape), data=bls, REML=F)

anova(m0,m1,m7,m8,m9,m10,m11,m12,m13,m4,m14,m15,m16,m17,m18,m19,m20,m21,m22)
anova(m0,m10,m12)
anova(m7,m10)
# with rice it is significantly better than null model

# going on adding preds
# categories LU_crops, LU_agrice and LU_agrar not usable (products of LU_rice)
# LU_fruit, LU_crops, LU_agrice, LU_agrar, LS_SDiv, LS_Even, MPS_ricepatch not usable (are colinear with LU_rice)

# testing additive effects with remaining preds

m23<-lmer(S~scale(LU_rice) + scale(LU_grass) + (1|Region/Landscape), data=bls, REML=F)
m24<-lmer(S~scale(LU_rice) + scale(LU_noveg) + (1|Region/Landscape), data=bls, REML=F)
m25<-lmer(S~scale(LU_rice) + scale(LU_vegetable) + (1|Region/Landscape), data=bls, REML=F)
m26<-lmer(S~scale(LU_rice) + scale(poly(LU_water,2)) + (1|Region/Landscape), data=bls, REML=F)
m27<-lmer(S~scale(LU_rice) + scale(LS_Rich) + (1|Region/Landscape), data=bls, REML=F)
m28<-lmer(S~scale(LU_rice) + scale(N_ricepatch) + (1|Region/Landscape), data=bls, REML=F)
m29<-lmer(S~scale(LU_rice) + scale(NP_Diss) + (1|Region/Landscape), data=bls, REML=F)
m30<-lmer(S~scale(LU_rice) + scale(ED_Diss) + (1|Region/Landscape), data=bls, REML=F)

anova(m10,m23,m24,m25,m26,m27,m28,m29,m30)
anova(m10,m26)

m31<-lmer(S~scale(LU_rice) * scale(LU_grass) + (1|Region/Landscape), data=bls, REML=F)
m32<-lmer(S~scale(LU_rice) * scale(LU_noveg) + (1|Region/Landscape), data=bls, REML=F)
m33<-lmer(S~scale(LU_rice) * scale(LU_vegetable) + (1|Region/Landscape), data=bls, REML=F)
m34<-lmer(S~scale(LU_rice) * scale(poly(LU_water,2)) + (1|Region/Landscape), data=bls, REML=F)
m35<-lmer(S~scale(LU_rice) * scale(LS_Rich) + (1|Region/Landscape), data=bls, REML=F)
m36<-lmer(S~scale(LU_rice) * scale(N_ricepatch) + (1|Region/Landscape), data=bls, REML=F)
m37<-lmer(S~scale(LU_rice) * scale(NP_Diss) + (1|Region/Landscape), data=bls, REML=F)
m38<-lmer(S~scale(LU_rice) * scale(ED_Diss) + (1|Region/Landscape), data=bls, REML=F)

anova(m10,m26,m31,m32,m33,m34,m35,m36,m37,m38)
# no improvement with including interaction effects

# further testing additive effects

m39<-lmer(S~scale(LU_rice) + scale(poly(LU_water,2)) + scale(LU_grass) + (1|Region/Landscape), data=bls, REML=F)
m40<-lmer(S~scale(LU_rice) + scale(poly(LU_water,2)) + scale(LU_noveg) + (1|Region/Landscape), data=bls, REML=F)
m41<-lmer(S~scale(LU_rice) + scale(poly(LU_water,2)) + scale(LU_vegetable) + (1|Region/Landscape), data=bls, REML=F)
m42<-lmer(S~scale(LU_rice) + scale(poly(LU_water,2)) + scale(LS_Rich)  + (1|Region/Landscape), data=bls, REML=F)
m43<-lmer(S~scale(LU_rice) + scale(poly(LU_water,2)) + scale(N_ricepatch) + (1|Region/Landscape), data=bls, REML=F)
m44<-lmer(S~scale(LU_rice) + scale(poly(LU_water,2)) + scale(NP_Diss) + (1|Region/Landscape), data=bls, REML=F)
m45<-lmer(S~scale(LU_rice) + scale(poly(LU_water,2)) + scale(ED_Diss) + (1|Region/Landscape), data=bls, REML=F)

anova(m26,m39,m40,m41,m42,m43,m44,m45)
anova(m26,m45)
# nothing is better

m46<-lmer(S~(scale(LU_rice) + scale(poly(LU_water,2))) * scale(LU_grass) + (1|Region/Landscape), data=bls, REML=F)
m47<-lmer(S~(scale(LU_rice) + scale(poly(LU_water,2))) * scale(LU_noveg) + (1|Region/Landscape), data=bls, REML=F)
m48<-lmer(S~(scale(LU_rice) + scale(poly(LU_water,2))) * scale(LU_vegetable) + (1|Region/Landscape), data=bls, REML=F)
m49<-lmer(S~(scale(LU_rice) + scale(poly(LU_water,2))) * scale(LS_Rich)  + (1|Region/Landscape), data=bls, REML=F)
m50<-lmer(S~(scale(LU_rice) + scale(poly(LU_water,2))) * scale(N_ricepatch) + (1|Region/Landscape), data=bls, REML=F)
m51<-lmer(S~(scale(LU_rice) + scale(poly(LU_water,2))) * scale(NP_Diss) + (1|Region/Landscape), data=bls, REML=F)
m52<-lmer(S~(scale(LU_rice) + scale(poly(LU_water,2))) * scale(ED_Diss) + (1|Region/Landscape), data=bls, REML=F)

anova(m26,m46,m47,m48,m49,m50,m51,m52)
# nothing better

# final model is m26 with additive effects of LU_rice and LU_water
```

```{r bls model diagnostics}
model<-m26
plotConventionalResiduals(model) 
simR<-simulateResiduals(model)
plotSimulatedResiduals(simR)  
# ok

summary(model)

sjp.lmer(model, type="fe", showIntercept=FALSE, title="")
```

```{r bls generalized predictions}

# check min/max, variation... for setting up general gradients
sort(unique(bls$LU_rice)) 
sort(unique(bls$LU_water)) 

# gradients
rice<-seq(0.04,1,0.01)
water<-seq(0,0.2,0.01)

97*21*70    # schätzung der Gesamtdatensätze für alle Kombinationen

# Create new data.frame from original (fixed) combinations of Region/Landscape 
subset<-data.frame(Region=bls$Region, Landscape=bls$Landscape)
subset<-subset[sample(1:nrow(subset), 150000, replace=TRUE),]    # resample 

# Create data.frame with synthetic gradient
new<-expand.grid(LU_rice=rice, LU_water=water)     # create all combinations
new<-data.frame(new[sample(1:nrow(new), 150000, replace=TRUE),])
general<-data.frame(subset, new)   # zusammenfügen
# generalized model prediction
general$pred<-predict(model, newdata=general, type="response", allow.new.levels=T) 
hist(general$pred)
hist(bls$S)
# very good


###### only rice ########
# plot predicted species richness of the regions, depending on prop rice only
    
colour<-c("red", "orangered", "orangered3", "royalblue1", "lightseagreen", "midnightblue", "lightblue")
par(mfrow=c(1,1))
plot(general$LU_rice, general$pred, col=factor(general$Region), type="n")
for (i in 1:7) {
  lines(names(tapply(general[factor(general$Region)==levels(factor(general$Region))[i],]$pred, general[factor(general$Region)==levels(factor(general$Region))[i],]$LU_rice, mean)),
        tapply(general[factor(general$Region)==levels(factor(general$Region))[i],]$pred, general[factor(general$Region)==levels(factor(general$Region))[i],]$LU_rice,mean), col=colour[i]) 
}
legend("topleft", legend=levels(factor(general$Region)), lwd=1, col=colour)

# plot predicted species richness of the regions (with loess)
par(mfrow=c(1,1))
plot(general$LU_rice, general$pred, col=factor(general$Region), type="n")
for (i in 1:7) {
  lines(loess.smooth(general[factor(general$Region)==levels(factor(general$Region))[i],]$LU_rice, general[factor(general$Region)==levels(factor(general$Region))[i],]$pred), col=colour[i])
}
legend("topleft", legend=levels(factor(general$Region)), lwd=1, col=colour)

# regions in original ranges 
par(mfrow=c(1,1))
plot(general$LU_rice, general$pred, col=factor(general$Region), type="n")
for (i in 1:7) {
  lines(loess.smooth(general[factor(general$Region)==levels(factor(general$Region))[i] &
                             general$LU_rice <= max(bls[bls$Region==levels(factor(general$Region))[i],]$LU_rice) &
                             general$LU_rice >= min(bls[bls$Region==levels(factor(general$Region))[i],]$LU_rice),]$LU_rice, 
                     general[factor(general$Region)==levels(factor(general$Region))[i] &
                             general$LU_rice <= max(bls[bls$Region==levels(factor(general$Region))[i],]$LU_rice) &
                             general$LU_rice >= min(bls[bls$Region==levels(factor(general$Region))[i],]$LU_rice),]$pred), col=colour[i])
}
legend("topright", legend=levels(factor(general$Region)), lwd=1, col=colour)

###### only water ########
# plot predicted species richness of the regions, depending on prop water only
    
colour<-c("red", "orangered", "orangered3", "royalblue1", "lightseagreen", "midnightblue", "lightblue")
par(mfrow=c(1,1))
plot(general$LU_water, general$pred, col=factor(general$Region), type="n")
for (i in 1:7) {
  lines(names(tapply(general[factor(general$Region)==levels(factor(general$Region))[i],]$pred, general[factor(general$Region)==levels(factor(general$Region))[i],]$LU_water, mean)),
        tapply(general[factor(general$Region)==levels(factor(general$Region))[i],]$pred, general[factor(general$Region)==levels(factor(general$Region))[i],]$LU_water,mean), col=colour[i]) 
}
legend("topleft", legend=levels(factor(general$Region)), lwd=1, col=colour)

# plot predicted species richness of the regions (with loess)
par(mfrow=c(1,1))
plot(general$LU_water, general$pred, col=factor(general$Region), type="n")
for (i in 1:7) {
  lines(loess.smooth(general[factor(general$Region)==levels(factor(general$Region))[i] &
                             general$LU_water <= max(bls[bls$Region==levels(factor(general$Region))[i],]$LU_water) &
                             general$LU_water >= min(bls[bls$Region==levels(factor(general$Region))[i],]$LU_water),]$LU_water, 
                     general[factor(general$Region)==levels(factor(general$Region))[i] &
                             general$LU_water <= max(bls[bls$Region==levels(factor(general$Region))[i],]$LU_water) &
                             general$LU_water >= min(bls[bls$Region==levels(factor(general$Region))[i],]$LU_water),]$pred), col=colour[i])
}
legend("topright", legend=levels(factor(general$Region)), lwd=1, col=colour)

```

```{r bls combination plot for additive effects}
# plot for additive effects of rice and water proportion in each region
# grey=low water proportion
# blue=high water proportion
# only original region ranges
par(mfrow=c(2,4))
for (k in 1:7) {
region <- levels(factor(general$Region))[k]
colour<-colorRampPalette(c("grey", "darkblue")) (length(unique(general$LU_water)))
palette(colour)
plot(general[general$Region==region,]$LU_rice, general[general$Region==region,]$pred, type="n", main=region, xlab="% rice", ylab="species richness")
for (i in 1:length(unique(general[general$Region==region,]$LU_water))) {
  lines(general[general$Region==region & 
                general$LU_rice <= max(bls[bls$Region==region,]$LU_rice) &
                general$LU_rice >= min(bls[bls$Region==region,]$LU_rice) & 
                general$LU_water <= max(bls[bls$Region==region,]$LU_water) &
                general$LU_water >= min(bls[bls$Region==region,]$LU_water) &  
                general$LU_water==unique(general[general$Region==region,]$LU_water)[i],]$LU_rice,
                general[general$Region==region & 
                general$LU_rice <= max(bls[bls$Region==region,]$LU_rice) &
                general$LU_rice >= min(bls[bls$Region==region,]$LU_rice) & 
                general$LU_water <= max(bls[bls$Region==region,]$LU_water) &
                general$LU_water >= min(bls[bls$Region==region,]$LU_water) & 
                general$LU_water==unique(general[general$Region==region,]$LU_water)[i],]$pred,
      col=colour[unique(general$LU_water)==general[general$Region==region &
                             general$LU_rice <= max(bls[bls$Region==region,]$LU_rice) &
                             general$LU_rice >= min(bls[bls$Region==region,]$LU_rice) &
                             general$LU_water <= max(bls[bls$Region==region,]$LU_water) &
                             general$LU_water >= min(bls[bls$Region==region,]$LU_water) &   
                             general$LU_water==unique(general[general$Region==region,]$LU_water)[i],]$LU_water], lwd=2) 
    }
}
plot(general$LU_rice, general$pred, type="n", axes=F, xlab="", ylab="")
legend("left", legend=c(min(water)*100, rep("",length(water)-2), max(water)*100), lty=1, lwd=2, col=colour, bty="n", y.intersp=0.5)
text(0.6,20, labels="% water", cex=1.5)

# check
sort(bls[bls$Region==region[1],]$LU_water)


par(mfrow=c(2,4))
for (k in 1:7) {
region <- levels(factor(general$Region))[k]
colour<-colorRampPalette(c("grey", "red")) (length(unique(general[general$Region==region,]$LU_rice)))
plot(general[general$Region==region,]$LU_water, general[general$Region==region,]$pred, type="n", main=region, xlab="% water", ylab="species richness")
for (i in seq(1,90,10)) {
  points(general[general$Region==region &
                             general$LU_water <= max(bls[bls$Region==region,]$LU_water) &
                             general$LU_water >= min(bls[bls$Region==region,]$LU_water) &
                             general$LU_rice==unique(general[general$Region==region,]$LU_rice)[i],]$LU_water,
                     general[general$Region==region & 
                             general$LU_water <= max(bls[bls$Region==region,]$LU_water) &
                             general$LU_water >= min(bls[bls$Region==region,]$LU_water) &
                             general$LU_rice==unique(general[general$Region==region,]$LU_rice)[i],]$pred,
      col=colour[i], pch=19) 
    }
}
plot(general$LU_water, general$pred, type="n", axes=F, xlab="", ylab="")
legend("left", legend=c(min(rice)*100, rep("",length(rice)-2), max(rice)*100), pch=19, col=colour, bty="n", y.intersp=0.05)
text(0.1,20, labels="% rice", cex=1.5)

par(mfrow=c(1,1))
colour<-colorRampPalette(c("white", "red")) (length(unique(round(general$pred, digit=0))))
palette(colour)
plot(general$LU_rice, general$LU_water, col=round(general$pred, digit=0), pch=15, cex=2)
# no clear pattern

# general over all regions
par(mfrow=c(1,2))
colour<-colorRampPalette(c("grey", "blue")) (length(unique(general$LU_water)))
plot(general$LU_rice, general$pred, type="n", xlab="% rice", ylab="species richness", ylim=c(20,50))
for (i in 1:length(unique(general$LU_water))) {
  lines(loess.smooth(general[general$LU_water==unique(general$LU_water)[i],]$LU_rice,
                     general[general$LU_water==unique(general$LU_water)[i],]$pred),
      col=colour[i], lwd=2) 
    }
plot(general$LU_rice, general$pred, type="n", axes=F, xlab="", ylab="")
legend("left", legend=c(min(water), rep("",length(water)-2), max(water)), lty=1, col=colour, bty="n", y.intersp=0.5)
text(0.4,80, labels="% grass", cex=1.5)

# overall effect not clear

```

water effect is artefact of only one region (vn1) and therefore not valid

```{r bls model diagnostics without water}
model<-m10
plotConventionalResiduals(model) 
simR<-simulateResiduals(model)
plotSimulatedResiduals(simR)  
# ok

summary(model)

sjp.lmer(model, type="fe", showIntercept=FALSE, title="")
```

```{r bls generalized predictions}

# check min/max, variation... for setting up general gradients
sort(unique(bls$LU_rice)) 
sort(unique(bls$LU_water)) 

# gradients
rice<-seq(0.04,1,0.01)

97*70    # schätzung der Gesamtdatensätze für alle Kombinationen

# Create new data.frame from original (fixed) combinations of Region/Landscape 
subset<-data.frame(Region=bls$Region, Landscape=bls$Landscape)
subset<-subset[sample(1:nrow(subset), 7000, replace=TRUE),]    # resample 

# Create data.frame with synthetic gradient
new<-data.frame(LU_rice=sample(rice, 7000, replace=TRUE))
general<-data.frame(subset, new)   # zusammenfügen
# generalized model prediction
general$pred<-predict(model, newdata=general, type="response", allow.new.levels=T) 
hist(general$pred)
hist(bls$S)
# very good


###### only rice ########
# plot predicted species richness of the regions, depending on prop rice only
    
colour<-c("red", "orangered", "orangered3", "royalblue1", "lightseagreen", "midnightblue", "lightblue")
par(mfrow=c(1,1))
plot(general$LU_rice, general$pred, col=factor(general$Region), type="n")
for (i in 1:7) {
  lines(names(tapply(general[factor(general$Region)==levels(factor(general$Region))[i],]$pred, general[factor(general$Region)==levels(factor(general$Region))[i],]$LU_rice, mean)),
        tapply(general[factor(general$Region)==levels(factor(general$Region))[i],]$pred, general[factor(general$Region)==levels(factor(general$Region))[i],]$LU_rice,mean), col=colour[i]) 
}
legend("topleft", legend=levels(factor(general$Region)), lwd=1, col=colour)

# plot predicted species richness of the regions (with loess)
par(mfrow=c(1,1))
plot(general$LU_rice, general$pred, col=factor(general$Region), type="n")
for (i in 1:7) {
  lines(loess.smooth(general[factor(general$Region)==levels(factor(general$Region))[i],]$LU_rice, general[factor(general$Region)==levels(factor(general$Region))[i],]$pred), col=colour[i])
}
legend("topleft", legend=levels(factor(general$Region)), lwd=1, col=colour)

# regions in original ranges 
par(mfrow=c(1,1))
plot(general$LU_rice, general$pred, col=factor(general$Region), type="n")
for (i in 1:7) {
  lines(loess.smooth(general[factor(general$Region)==levels(factor(general$Region))[i] &
                             general$LU_rice <= max(bls[bls$Region==levels(factor(general$Region))[i],]$LU_rice) &
                             general$LU_rice >= min(bls[bls$Region==levels(factor(general$Region))[i],]$LU_rice),]$LU_rice, 
                     general[factor(general$Region)==levels(factor(general$Region))[i] &
                             general$LU_rice <= max(bls[bls$Region==levels(factor(general$Region))[i],]$LU_rice) &
                             general$LU_rice >= min(bls[bls$Region==levels(factor(general$Region))[i],]$LU_rice),]$pred), col=colour[i])
}
legend("topright", legend=levels(factor(general$Region)), lwd=1, col=colour)
```

model with only rice is best; but region effect still dominating


