# load required packages

.libPaths("C:/R/library")
#setwd("//env/cllmm/Files/CLL Restoration/MA 6_Research and Monitoring/Projects/VP2015/Veg/code")

setwd("C:/VP2015")

library(extrafont)    # fonts for plotting (if using for first time need to load scripts - see help...)
library(reshape2)     # longfile to flatfile and vice versa
library(RODBC)        # import data from various ODBC datasources
library(FD)           # functional diversity indices
library(vegan)        # community analyses
library(lme4)         # multilevel modelling
library(arm)          # convenience functions for regression in R
library(sjPlot)       # plots effects of merMod objects, which were fitted using the lmer function of the lme4 package
library(labdsv)       # indval
library(lmerTest)
library(ggplot2)
library(plyr)




require(extrafont)    # fonts for plotting (if using for first time need to load scripts - see help...)
require(reshape2)     # longfile to flatfile and vice versa
require(RODBC)        # import data from various ODBC datasources
require(FD)           # functional diversity indices
require(vegan)        # community analyses
require(lme4)         # multilevel modelling
require(arm)          # convenience functions for regression in R
require(sjPlot)       # plots effects of merMod objects, which were fitted using the lmer function of the lme4 package
require(labdsv)       # indval
require(lmerTest)
require(ggplot2)
require(plyr)

# set working directory


setwd("//env/cllmm/Files/CLL Restoration/MA 6_Research and Monitoring/Projects/VP2015/Veg/rel")

dat <- read.csv("rel/q00200.csv", header = T)
dat2<- read.csv("rel/q00120.csv", header = T)
env<- read.csv("rel/q00300.csv", header = T)

names(dat) <- c("quad","spp","cov")
head(dat)
datmat <- dcast(dat, quad ~ spp, mean, value = "cov", fill = 0)
rownames(datmat) <- datmat[,1]
datmat <- datmat[,-1]
datmat[1:20,1:10]

dat2[1:10,1:10]
rownames(dat2) <- dat2[,1]
dat2 <- dat2[,-1]
dat2[1:10,1:10]

# #dat2 holds the species by trait matrix
# setwd("//env/cllmm/Files/CLL Restoration/MA 6_Research and Monitoring/Projects/VP2015/Veg/data")
# dat3<- read.csv("VP data Primer2.csv", header=T, row.names=1)

dat2<- as.matrix(dat2)
dat3<- as.matrix(datmat)

datpa<- decostand(dat3,method="pa")

is.matrix(dat2)
is.matrix(dat3)
is.matrix(datpa)

sitebytrait = dat2 %*% dat3
sitebytraitpa <-dat2 %*% datpa #not working need to check that all species are correctly named and present to do this

write.csv(dat2, file = "dat2.csv") 
write.csv(datpa, file = "datpa.csv")
dim(MacroT2)

storage.mode(trait2)
