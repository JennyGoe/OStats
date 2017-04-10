---
title: "management influence on species richness"
author: "Jenny Schellenberg"
date: "6 April 2017"
output: html_document
---


# 1: Species richness (S) negatively correlated with local land-use intensity
land-use intensitiy: which predictor to use for that? 
LUI preds

# management is influencing species richness patterns

# management is influencing life form composition of rice field vegetation

#################    Preparation      #################################

setwd("G:/Arbeit/Oli Statistik")
setwd("C:/Users/jschell/Dropbox/Uni/Oli")


setwd("/AAA/Dropbox/WORK/RData_OF/Paper3")

```{r packages and functions}
####### packages ########

library(lme4)
library(lattice)
library(MASS)
library(DHARMa)
library(sjPlot)
library(lmerTest)
library(gplots)



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

###### Create new Land use index 'LUI' ################################
# (incl. herbicide, cutting/handweeding, number of crops/year, ...)
# with scaled/normalized covariates before summarizing them

LUI_A1<- as.matrix((bhead$CropsPerYear/mean(bhead$CropsPerYear,na.rm=T))+(bhead$Herbicide/mean(bhead$Herbicide,na.rm=T))+(bhead$RemoveSoil/mean(bhead$RemoveSoil,na.rm=T))+(bhead$Cutting/mean(bhead$Cutting,na.rm=T))+(bhead$Handweeding/mean(bhead$Handweeding,na.rm=T))+(bhead$Grazing/mean(bhead$Grazing,na.rm=T)))

LUI_A2<- scale(bhead$CropsPerYear)+scale(bhead$Herbicide)+scale(bhead$RemoveSoil)+scale(bhead$Cutting)+scale(bhead$Handweeding)+scale(bhead$Grazing)

LUI_B1<- as.matrix((bhead$CropsPerYear/mean(bhead$CropsPerYear,na.rm=T))+(bhead$Herbicide/mean(bhead$Herbicide,na.rm=T))+(bhead$RemoveSoil/mean(bhead$RemoveSoil,na.rm=T))+(0.25*(bhead$Cutting/mean(bhead$Cutting,na.rm=T)))+(0.25*(bhead$Handweeding/mean(bhead$Handweeding,na.rm=T)))+(0.25*(bhead$Grazing/mean(bhead$Grazing,na.rm=T))))

LUI_B2<- scale(bhead$CropsPerYear)+scale(bhead$Herbicide)+scale(bhead$RemoveSoil)+scale(0.25*bhead$Cutting)+scale(0.25*bhead$Handweeding)+scale(0.25*bhead$Grazing)

# setting new working data frame
b <- data.frame(bhead[,c(1:3, 17:18,22:27,54)], LUI_A1,LUI_A2,LUI_B1,LUI_B2)
str(b)
summary(b)
# all management predictors have moire than 10% of NA
# yield: more than a half is missing... deleting that column

b <- b[,-5]
b <- b[complete.cases(b),]
# 94 samples left
summary(b)
# PH_2 missing now

# checking for managements that occur in at least 10% of dataset
apply(b>0,2,sum)    
0.1*length(b[,1])
# all ok
```

```{r S in regions}
###### species richness in regions: histograms and density plots #####

frame<-b
levels(frame$Region)

# create colors for regions: red/orange for PH, blue/bluegreen for VN, 
colour<-c("red", "orangered3", "royalblue1", "lightseagreen", "midnightblue", "lightblue")
# "orangered",
palette(colour)

par(mfrow=c(2,2))
hist(frame$S, breaks=20)
plot(density(frame$S, na.rm=TRUE), main="density only PH", col="grey", lwd=3, xlab="species richness", ylab="density", ylim=c(0,0.15))
for (i in 1:2) { 
  lines(density(frame[which(frame$Region==levels(factor(frame$Region))[i]),]$S, na.rm=TRUE), col=colour[i], lwd=3) }
legend("topright", lwd=3, col=c("grey", colour[1:2]), legend=c("all regions", levels(factor(frame$Region))[1:3]))

plot(density(frame$S, na.rm=TRUE), main="density only VN", col="grey", lwd=3, xlab="species richness", ylab="density", ylim=c(0,0.15))
for (i in 3:6) { lines(density(frame[which(frame$Region==levels(factor(frame$Region))[i]),]$S, na.rm=TRUE), col=colour[i], lwd=3)  }
legend("topright", lwd=3, col=c("grey", colour[3:6]), legend=c("all regions", levels(frame$Region)[4:7])) 

boxplot(frame$S~factor(frame$Region), whisklty=1, ylim=c(0,70), main= "boxplot with sample sizes")
for (i in 1:7) { 
  text( i ,68,labels=length(frame[frame$Region==levels(factor(frame$Region))[i],]$Region)) 
       }
# failing normal distribution
# distinct distribution with dominating values 20-40 species
# ok for PH_1, PH_3 
# VN_1, VN_2, VN_4 borderline... risk of artefact in modeling
# hard risk of artefact in VN_3 - is it really necessary???

# rest ok
# variance in regions not equal, distinct behavior of VN_4 
```

```{r preds in regions}
# management: histograms and distr. in areas #######
names(b)
frame <- b

par(mfrow=c(3,3)) 
for (i in c(4:10,12:15)) {
  hist(frame[,i], main=names(frame)[i]) }
# grazing 0/1
# LUI_A2 and B2 negative values

# preds in regions
par(mfrow=c(3,4)) 
for (i in c(4:10,12:15)) {
boxplot(frame[,i]~factor(frame$Region), whisklty=1, main=names(frame)[i])
  abline(h=0)
}
```

strange differences in areas with risk of bias
herbicide much higher in species-poor vn_4
because VN_4 is much lower in species richness, risk of dominating combination of specific VN_4 predictors
remove soil in only two regions, also mulching only seldom

```{r management fixed combinations and artefact risk}
# plots for fixed distinct combinations of management in regions

###### crops per year ######
testp <- frame$CropsPerYear
par(mfrow=c(2,4), mar=c(3,2,3,1))
for (i in 5:10) 
  { plot(jitter(testp, factor=0.6),jitter(frame[,i], factor=0.6), pch=19, col=factor(frame$Region), cex=1, main = names(frame)[i]) }
plot(testp,frame[,1], type="n", axes=F, xlab="", ylab="")
legend("left", legend=levels(factor(frame$Region)), cex=1.2, pch=19, col=colour)

# strong bindings to:
# Remove soil
# regions with specific managements:
# VN_3, VN_4, PH_3
# only VN_4 with 3 crops/season

###### Herbicide ######
testp <- frame$Herbicide
par(mfrow=c(2,4), mar=c(3,2,3,1))
for (i in c(4,6:10)) 
  { plot(jitter(testp, factor=0.6),jitter(frame[,i], factor=0.6), pch=19, col=factor(frame$Region), cex=1, main = names(frame)[i]) }
plot(testp,frame[,1], type="n", axes=F, xlab="", ylab="")
legend("left", legend=levels(factor(frame$Region)), cex=1.2, pch=19, col=colour)

# strong bindings to
# Remove soil, (weaker: crops per year)
# good scatter with others
# regions with specific management combination
# PH_1, PH_3

###### Cutting ######
testp <- frame$Cutting
par(mfrow=c(2,4), mar=c(3,2,3,1))
for (i in c(4:5,7:10)) 
  { plot(jitter(testp, factor=0.6),jitter(frame[,i], factor=0.6), pch=19, col=factor(frame$Region), cex=1, main = names(frame)[i]) }
plot(testp,frame[,1], type="n", axes=F, xlab="", ylab="")
legend("left", legend=levels(factor(frame$Region)), cex=1.2, pch=19, col=colour)

# no strong bindings, only weak: crops per year
# good scatter with others
# no regions with specific, unique management combination
# good predictor

###### Handweeding ######
testp <- frame$Handweeding
par(mfrow=c(2,4), mar=c(3,2,3,1))
for (i in c(4:6,8:10)) 
  { plot(jitter(testp, factor=0.6),jitter(frame[,i], factor=0.6), pch=19, col=factor(frame$Region), cex=1, main = names(frame)[i]) }
plot(testp,frame[,1], type="n", axes=F, xlab="", ylab="")
legend("left", legend=levels(factor(frame$Region)), cex=1.2, pch=19, col=colour)

# no strong bindings
# good scatter 
# regions with specific, unique management combination:
# PH_1, PH_3, VN_3

###### Remove soil ######
testp <- frame$RemoveSoil
par(mfrow=c(2,4), mar=c(3,2,3,1))
for (i in c(4:7,9:10)) 
  { plot(jitter(testp, factor=0.6),jitter(frame[,i], factor=0.6), pch=19, col=factor(frame$Region), cex=1, main = names(frame)[i]) }
plot(testp,frame[,1], type="n", axes=F, xlab="", ylab="")
legend("left", legend=levels(factor(frame$Region)), cex=1.2, pch=19, col=colour)

# strong bindings to:
# crops per year, weaker to herbicide and mulching
# regions with specific, unique management combination:
# all regions with more or less specific combination patterns
# no good predictor


###### Mulching ######
testp <- frame$Mulching
par(mfrow=c(2,4), mar=c(3,2,3,1))
for (i in c(4:8,10)) 
  { plot(jitter(testp, factor=0.6),jitter(frame[,i], factor=0.6), pch=19, col=factor(frame$Region), cex=1, main = names(frame)[i]) }
plot(testp,frame[,1], type="n", axes=F, xlab="", ylab="")
legend("left", legend=levels(factor(frame$Region)), cex=1.2, pch=19, col=colour)

# bindings to: crops per year, remove soil and grazing
# good scatter to cutting, handweeding
# regions with specific, unique management combination:
# PH_1, VN_3, 
# VN_4 completely without mulching - risk of artefact


###### Grazing ######
testp <- frame$Grazing
par(mfrow=c(2,4), mar=c(3,2,3,1))
for (i in 4:9) 
  { plot(jitter(testp, factor=0.6),jitter(frame[,i], factor=0.6), pch=19, col=factor(frame$Region), cex=1, main = names(frame)[i]) }
plot(testp,frame[,1], type="n", axes=F, xlab="", ylab="")
legend("left", legend=levels(factor(frame$Region)), cex=1.2, pch=19, col=colour)

# only 0/1, and therefore strong bindings to : 
# crops per year, herbicide, remove soil and mulching
# all regions with more or less specific combination patterns
# no good predictor
```

summary of qualtiy of preds and fixed combinations
not good as predictor because of high artefact risk
- remove soil, mulching, grazing


```{r collinearity of managements}
pairs(~S+CropsPerYear+Herbicide+Cutting+Handweeding, data=frame,upper.panel=panel.cor.pearson, diag.panel=panel.hist, lower.panel=panel.smooth, main="colinearity between predictors land use structures", na.action=na.omit)
```
 
cutting and handweeding correlated; not both 

no need for testing non-linearity - there is no evidence for nonlinear relation

```{r test of model technique}
# tested with crops per year and Handweeding
model<-m5
plotConventionalResiduals(model) 
simR<-simulateResiduals(model)
plotSimulatedResiduals(simR)    



frame <- b
m1<-lmer(S~CropsPerYear + (1|Region/Landscape), data=frame, REML=F)
m2<-glmer(S~CropsPerYear + (1|Region/Landscape), data=frame, family=poisson)
overdisp_fun(m2)
m3<-glmer.nb(S~CropsPerYear + (1|Region/Landscape), data=frame)
# lmer is best

m4<-lmer(S~Handweeding + (1|Region/Landscape), data=frame, REML=F)
m5<-glmer(S~Handweeding + (1|Region/Landscape), data=frame, family=poisson)
overdisp_fun(m5)
m6<-glmer.nb(S~Handweeding + (1|Region/Landscape), data=frame)
# lmer is best

# taking lmer for modeling
```

```{r model selection}
m0<-lmer(S~1 + (1|Region/Landscape), data=frame, REML=F)

m7<-lmer(S~Herbicide + (1|Region/Landscape), data=frame, REML=F)
m8<-lmer(S~Cutting + (1|Region/Landscape), data=frame, REML=F)

anova(m0,m1,m4,m7,m8)
# with crops per year significantly better

m9<-lmer(S~scale(CropsPerYear) + scale(Handweeding) + (1|Region/Landscape), data=frame, REML=F)
m10<-lmer(S~scale(CropsPerYear) + scale(Herbicide) + (1|Region/Landscape), data=frame, REML=F)
m11<-lmer(S~scale(CropsPerYear) + scale(Cutting) + (1|Region/Landscape), data=frame, REML=F)

anova(m1,m9,m10,m11)
# nothing is better

# interaction tests
m12<-lmer(S~scale(CropsPerYear) * scale(Handweeding) + (1|Region/Landscape), data=frame, REML=F)
m13<-lmer(S~scale(CropsPerYear) * scale(Herbicide) + (1|Region/Landscape), data=frame, REML=F)
m14<-lmer(S~scale(CropsPerYear) * scale(Cutting) + (1|Region/Landscape), data=frame, REML=F)

anova(m1,m12,m13,m14)
anova(m1,m13)
# not significantly better

# final model is m1 with only crops per year

m15<-lmer(S ~ Region + offset(scale(CropsPerYear)) + (1|Landscape), data=frame, REML=F)
m16<-lmer(S ~ Region + (1| Landscape), data=frame, REML=F)

anova(m1,m15,m16)
# regional impact is much more higher than explanation power of crops/year!
# interesting; impact is so low, that including crops/year as offset does not have any influence on model!!!

```

```{r model diagnostics}
model<-m1
plotConventionalResiduals(model) 
simR<-simulateResiduals(model)
plotSimulatedResiduals(simR)  
# ok

summary(model)

sjp.lmer(model, type="fe", showIntercept=FALSE, title="")
# huge negative effect
```

```{r generalized predictions}
frame <- b
# check min/max, variation... for setting up general gradients
sort(unique(frame$CropsPerYear)) 

# gradients
crops<-seq(1,3,1)

94*3    # sch?tzung der Gesamtdatens?tze f?r alle Kombinationen

# Create new data.frame from original (fixed) combinations of Region/Landscape 
subset<-data.frame(Region=frame$Region, Landscape=frame$Landscape)
subset<-subset[sample(1:nrow(subset), 300, replace=TRUE),]    # resample 

# Create data.frame with synthetic gradient
new<-data.frame(CropsPerYear=sample(crops, 300, replace=T))
general<-data.frame(subset, new)   # zusammenf?gen
# generalized model prediction
general$pred<-predict(model, newdata=general, type="response", allow.new.levels=T) 
hist(general$pred)
hist(frame$S)
# very good

```

```{r plot }
# region in original ranges 
par(mfrow=c(1,1))
plot(general$CropsPerYear, general$pred, col=factor(general$Region), type="n")
for (i in 1:7) {
  points(jitter(general[factor(general$Region)==levels(factor(general$Region))[i] &
         general$CropsPerYear <= max(frame[frame$Region==levels(factor(general$Region))[i],]$CropsPerYear) &
         general$CropsPerYear >= min(frame[frame$Region==levels(factor(general$Region))[i],]$CropsPerYear),]$CropsPerYear, factor=0.8), 
         jitter(general[factor(general$Region)==levels(factor(general$Region))[i] &
         general$CropsPerYear <= max(frame[frame$Region==levels(factor(general$Region))[i],]$CropsPerYear) &
                             general$CropsPerYear >= min(frame[frame$Region==levels(factor(general$Region))[i],]$CropsPerYear),]$pred, factor=0.8), col=colour[i], pch=19)
}
legend("topright", legend=levels(factor(general$Region)), pch=19, col=colour)




 # general effect from generalized predictions

boxplot(general$pred~general$CropsPerYear)
```

There is a clear influence of crops / year as land use factor on species richness.
BUT:
results +/- meaningless:
- regions completely distinct in crops per year already in original data
- explanation character for S is much better with Region as pred

If you want to use it anyway, add significance group differences tests (t-test, wilcox) to underline result
but mention that basic regional difference is maybe crucially inlfuencing this result; 
correlation is given, but causality not necessarily!!!


######### Land use intensity indices LUI ############

same dataset as above

```{r LUI plots}
# LUI histograms and distr. in areas #######
names(b)
frame <- b

par(mfrow=c(2,2)) 
for (i in c(12:15)) {
  hist(frame[,i], main=names(frame)[i]) }
# LUI_A2 and B2 negative values

# preds in regions
par(mfrow=c(2,2)) 
for (i in c(12:15)) {
boxplot(frame[,i]~factor(frame$Region), whisklty=1, main=names(frame)[i])
  abline(h=0)
}

pairs(~S+LUI_A1+LUI_A2+LUI_B1+LUI_B2, data=frame,upper.panel=panel.cor.pearson, diag.panel=panel.hist, lower.panel=panel.smooth, main="colinearity between predictors land use structures", na.action=na.omit)

# all correlated; only one is to be used
```


```{r model selection}
# m0 like above, same dataset, so it can be used

m17<-lmer(S~ LUI_A1 + (1|Region/Landscape), data=frame, REML=F)
m18<-lmer(S~ LUI_A2 + (1|Region/Landscape), data=frame, REML=F)
m19<-lmer(S~ LUI_B1 + (1|Region/Landscape), data=frame, REML=F)
m20<-lmer(S~ LUI_B2 + (1|Region/Landscape), data=frame, REML=F)

anova(m0,m17,m18,m19,m20)
# they are all not better than nullmodel!
# that means they cannot explain variation and S is completely independent from preds; no patterns found.

```

No patterns.
We should try life form pattern dependency on land use intensity. (As proposed last meeting; classify life form pattern groups and test them for differences in land unse intensity).
