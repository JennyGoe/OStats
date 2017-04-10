---
title: "# 4: Species richness is a function of environmental variables, e.g. MeanTemp, pH,..."
author: "Dipl.-Biol. Jenny Schellenberg"
date: "21 März 2017"
output: html_document
---

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


####### function to test for overdispersion (from http://glmm.wikidot.com/faq) #######
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

```

```{r data}
### env data "bhead"
bhead<-read.csv("bheadall.txt", header=T, sep="\t")
#Inspect the results
names(bhead)
str(bhead)
summary(bhead)

# select relevant columns
b<-bhead[,c(1:3, 8:15, 54)]    # only choosing predictors that have not more than ~10% of NA
b<-b[complete.cases(b),]
summary(b)
```





################ Overview plots ##################################
```{r overviews}
###### species richness in regions: histograms and density plots #####

levels(b$Region)
# create colors for regions: red/orange for PH, blue/bluegreen for VN, 
colour<-c("red", "orangered", "orangered3", "royalblue1", "lightseagreen", "midnightblue", "lightblue")

par(mfrow=c(2,2))
hist(bhead$S, breaks=20)
plot(density(b$S, na.rm=TRUE), main="density only PH", col="grey", lwd=3, xlab="species richness", ylab="density", ylim=c(0,0.15))
for (i in 1:3) { 
  lines(density(b[which(b$Region==levels(b$Region)[i]),]$S, na.rm=TRUE), col=colour[i], lwd=3) }
legend("topright", lwd=3, col=c("grey", colour[1:3]), legend=c("all regions", levels(b$Region)[1:3]))

plot(density(b$S, na.rm=TRUE), main="density only VN", col="grey", lwd=3, xlab="species richness", ylab="density", ylim=c(0,0.15))
for (i in 4:7) { lines(density(b[which(b$Region==levels(b$Region)[i]),]$S, na.rm=TRUE), col=colour[i], lwd=3)  }
legend("topright", lwd=3, col=c("grey", colour[4:7]), legend=c("all regions", levels(b$Region)[4:7])) 

boxplot(b$S~factor(b$Region), whisklty=1, ylim=c(0,70), main= "boxplot with sample sizes")
for (i in 1:7) { 
  text( i ,68,labels=length(b[b$Region==levels(b$Region)[i],]$Region)) 
       }

# failing normal distribution, but acceptable for all Regions, PH and VN_1, 
# VN_2 and VN_4 borderline... risk of artefact in modeling
# hard risk of artefact in VN_3 - is it really necessary???
```

######### influences of environmental variables on S ##################

```{r predictor check}
# histogram for all predictors
par(mfrow=c(2,4)) 
for (i in 4:11) {
  hist(b[,i], main=names(b)[i]) }

# preds in regions
par(mfrow=c(2,4)) 
for (i in 4:11) {
boxplot(b[,i]~factor(b$Region), whisklty=1, main=names(b)[i])
}

# PH_3 and VN_3 mountain region with much higher elevation and so completely different temperature.... not really comparable. I suppose to do analysis without them.

# Rock_Frag_Perc/MaxRockDiam only in PH_3, therefore not to be recognized
# other preds ok
```

only non-mountain regions.
 
 Predictors to test: 
 Elevation, pH, EC, CN, Humus and MeanTemp

```{r dataframe correction}
names(b)
blow<-b
# blow<-b[(b$Region!="PH_3"), c(1:4,7:12)]
hist(blow$S)

names(blow)
par(mfrow=c(2,4)) 
for (i in 4:10) {
boxplot(blow[,i]~factor(blow$Region), whisklty=1, main=names(blow)[i])
}
```

```{r collinearity of predictors}
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

##### colinearity between preds #####
pairs(~Elevation+pH+EC+C.N+Humus+MeanTemp+S, data=blow,upper.panel=panel.cor.pearson, diag.panel=panel.hist, lower.panel=panel.smooth, main="colinearity between predictors")
```

******Sehr interessant*********

in cor panel ist die beste erklärende Variable für S MeanTemp

nur lineare Beziehungen - kein test für nicht-lineare Beziehungen nötig 

da zwischen den >0.5 korr. Variablen kein kausaler Zusammenhang besteht, akzeptiere ich Collinearität bis 0.6

nur EC und Humus sind stärker korreliert und dürfen nicht zusammen in ein Modell

```{r model technique}
# setting up different models to test best fitting technique for response
# testing with pH and Elevation
model<-m2
plotConventionalResiduals(model) 
simR<-simulateResiduals(model)
plotSimulatedResiduals(simR)    

m1<-lmer(S~Elevation + (1|Region/Landscape), data=blow, REML=F)
m2<-glmer(S~Elevation + (1|Region/Landscape), data=blow, REML=F, family=poisson)
overdisp_fun(m2)   
m3<-glmer.nb(S~Elevation + (1|Region/Landscape), data=blow, REML=F)
overdisp_fun(m3)
# m1 not good
# m2 and m3 not better

m4<-lmer(S~pH + (1|Region/Landscape), data=blow, REML=F)
m5<-glmer(S~pH + (1|Region/Landscape), data=blow, REML=F, family=poisson)
overdisp_fun(m5)   
m6<-glmer.nb(S~pH + (1|Region/Landscape), data=blow, REML=F)
overdisp_fun(m6)
# m4 ok, m5 and m6 not better

# choosing lmer as model technique
```

```{r forward selection}
# starting with nullmodel
m0 <-lmer(S~1 + (1|Region/Landscape), data=blow, REML=F)
# m1 is model with Elevenation
# m4 is model with pH
m7<-lmer(S~EC + (1|Region/Landscape), data=blow, REML=F)
m8<-lmer(S~C.N + (1|Region/Landscape), data=blow, REML=F)
m9<-lmer(S~Humus + (1|Region/Landscape), data=blow, REML=F)
m10<-lmer(S~MeanTemp + (1|Region/Landscape), data=blow, REML=F)

anova(m0,m1,m4,m7,m8,m9,m10)
anova(m0,m4,m10)
# m4 is best explaining model

# testing additive effects on non-colinear predictors next
m11<-lmer(S~scale(pH) + scale(Elevation)+(1|Region/Landscape), data=blow, REML=F)
m12<-lmer(S~scale(pH) + scale(EC)+(1|Region/Landscape), data=blow, REML=F)
m13<-lmer(S~scale(pH) + scale(C.N) +(1|Region/Landscape), data=blow, REML=F)
m14<-lmer(S~scale(pH) + scale(Humus)+ (1|Region/Landscape), data=blow, REML=F)
m15<-lmer(S~scale(pH) + scale(MeanTemp) + (1|Region/Landscape), data=blow, REML=F)

anova(m4,m11,m12,m13,m14,m15)
anova(m4,m15)
# nothing is significantly better except MeanTemp: but this is only because of mountain region and not to be recognized

# testing interaction effects not done; high risk of artefact with interaction effect

# nothing is better

# final model is m4 pH
```


```{r model diagnostics}
model<-m4
plotConventionalResiduals(model) 
simR<-simulateResiduals(model)
plotSimulatedResiduals(simR)  
# quite good

sjp.lmer(model, type="fe", showIntercept=FALSE, title="")
# pH mit stark positiven Effekt auf S (rd 4 Arten Zunahme bei Steigen pH um 1)

```

#################   prediction plots    ######################
```{r generalized predictions}
sort(unique(blow$pH)) # check pH, min/max, variation... range from 4.05 to 7.49
length(sort(unique(blow$pH)))
# pH gradient for continuous line
pH <- seq(min(blow$pH), max(blow$pH),0.05)

123*69    # schätzung der Gesamtdatensätze für alle Kombinationen: 96 obs., 100 pH-levels  = ~9000

# Create new data.frame from original (fixed) combinations of Region/Landscape 
subset<-data.frame(Region=blow$Region, Landscape=blow$Landscape)
subset<-subset[sample(1:nrow(subset), 9000, replace=TRUE),]    # resample 

# Create data.frame with artificial pH gradient
new<-data.frame(pH=sample(pH, 9000, replace=TRUE))
general<-data.frame(subset, new)   # zusammenfügen
# generalized model prediction
general$pred<-predict(model, newdata=general, type="response", allow.new.levels=T) 
hist(general$pred)
# very good, without artefact in prediciton range



###### plot ########

# plot predicted species richness of the regions, depending on pH only, whole pH-range

colour<-c( "red", "orangered3", "royalblue1", "lightseagreen", "midnightblue", "green", "black") 
par(mfrow=c(1,1))
plot(general$pH, general$pred, col=factor(general$Region), type="n")
for (i in 1:7) {
  lines(names(tapply(general[factor(general$Region)==levels(factor(general$Region))[i],]$pred, general[factor(general$Region)==levels(factor(general$Region))[i],]$pH, mean)),
        tapply(general[factor(general$Region)==levels(factor(general$Region))[i],]$pred, general[factor(general$Region)==levels(factor(general$Region))[i],]$pH,mean), col=colour[i]) 
}
legend("topright", legend=levels(factor(general$Region)), lwd=1, col=colour)

# only range of real existing pH and MeanTemp range of regions
# setting dataframe 
gen<-general[general$Region==levels(factor(general$Region))[1] &  
             general$pH<=max(blow[blow$Region==levels(factor(blow$Region))[1],]$pH) &
             general$pH>=min(blow[blow$Region==levels(factor(blow$Region))[1],]$pH) |
             general$Region==levels(factor(general$Region))[2] &  
             general$pH<=max(blow[blow$Region==levels(factor(blow$Region))[2],]$pH) &
             general$pH>=min(blow[blow$Region==levels(factor(blow$Region))[2],]$pH) |
             general$Region==levels(factor(general$Region))[3] &  
             general$pH<=max(blow[blow$Region==levels(factor(blow$Region))[3],]$pH) &
             general$pH>=min(blow[blow$Region==levels(factor(blow$Region))[3],]$pH) |
             general$Region==levels(factor(general$Region))[4] &  
             general$pH<=max(blow[blow$Region==levels(factor(blow$Region))[4],]$pH) &
             general$pH>=min(blow[blow$Region==levels(factor(blow$Region))[4],]$pH) |   
             general$Region==levels(factor(general$Region))[5] &  
             general$pH<=max(blow[blow$Region==levels(factor(blow$Region))[5],]$pH) &
             general$pH>=min(blow[blow$Region==levels(factor(blow$Region))[5],]$pH) |
             general$Region==levels(factor(general$Region))[6] &  
             general$pH<=max(blow[blow$Region==levels(factor(blow$Region))[6],]$pH) &
             general$pH>=min(blow[blow$Region==levels(factor(blow$Region))[6],]$pH)|
             general$Region==levels(factor(general$Region))[7] &  
             general$pH<=max(blow[blow$Region==levels(factor(blow$Region))[7],]$pH) &
             general$pH>=min(blow[blow$Region==levels(factor(blow$Region))[7],]$pH),]

summary(gen)
# control
range(blow[blow$Region=="PH_1",]$pH)
range(gen[gen$Region=="PH_1",]$pH)

range(blow[blow$Region=="PH_2",]$pH)
range(gen[gen$Region=="PH_2",]$pH)

range(blow[blow$Region=="VN_1",]$pH)
range(gen[gen$Region=="VN_1",]$pH)
# ok

# plot predicted species richness of the regions, only real existing pH- and Temp-ranges
colour <- c(rep("black",3), rep("grey60",4))
linetype <- c(1,2,3,1,2,3,4)

png("spec_div_bunds_result2.png", width=800, height=800)
par(mfrow=c(1,1))
plot(gen$pH, gen$pred, col=factor(gen$Region), type="n", xlim=c(4,8), ylab="species", xlab="soil pH")
for (i in 1:7) {
  lines(loess.smooth(gen[factor(gen$Region)==levels(factor(gen$Region))[i],]$pH, gen[factor(gen$Region)==levels(factor(gen$Region))[i],]$pred), col=colour[i], lty=linetype[i], lwd=2)
}
legend("topright", legend=levels(factor(gen$Region)), lty=linetype, lwd=2, col=colour)
dev.off()

```


#####################   general pattern of species diversity over all Regions  ########################

```{r effect size of region vs pH*Temp}
# 5: The patterns found are independent from the study region

# for testeing this, model with Regions only should be not better than the one with pH


m30<-lmer(S~Region + (1|Landscape), data=blow, REML=F)
m31<-lmer(S~Region + offset(pH) + (1|Landscape), data=blow, REML=F)

anova(m4, m30, m31)
# others are significantly better

# is to be rejected.

```

# 5: The patterns found are independent from the study region
is to be rejected.


######### thinking of ignoring Region as random term #########
```{r}
par(mfrow=c(1,1), mar=c(10,2,2,2))
bymedian<-with(b, reorder(factor(Landscape),-S,median))
boxplot(S~bymedian, data=b, names=substr(levels(bymedian),1,4), las=2)
```

makes no sense, as regional impact is known and is not to be ignored. 

################# final boxplots of regional differences in species richness ###############

```{r final boxplots and tests}

names(blow)

pairwise.wilcox.test(blow$S,blow$Region, p.adjust="bonf")
# ph_2~all     vn1~vn2     vn4~all
pnotes <- c("cd","a","cd","c","d","cd","b")

pairwise.wilcox.test(gen$pred,gen$Region)
# all
ppnotes <- c("a","b","c","d","e","f","g")

par(mfrow=c(1,2))
boxplot(blow$S~blow$Region, col="grey", whisklty=1, pch=19, cex=0.6, main="original data", las=2, ylab="species richness", ylim=c(0,max(blow$S)*1.1))
for ( i in 1:7) {
  text(i,max(blow$S)*1.1, labels=pnotes[i], cex=0.8)
}

boxplot(gen$pred~gen$Region, col="grey", whisklty=1, pch=19, cex=0.6, main="model data", las=2,ylim=c(0,max(gen$pred)*1.1))
for ( i in 1:7) {
  text(i,max(gen$pred)*1.1, labels=ppnotes[i], cex=0.8)
}



```



