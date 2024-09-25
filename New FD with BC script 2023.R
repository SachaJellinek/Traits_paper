#setwd("C:/VP2015")
setwd("~/uomShare/wergProj/W12 - Revegetation/CLLMM work/Traits paper/Code and data")

library(extrafont)    # fonts for plotting (if using for first time need to load scripts - see help...)
library(reshape2)     # longfile to flatfile and vice versa
#library(RODBC)        # import data from various ODBC datasources
library(FD)           # functional diversity indices
library(vegan)        # community analyses
library(lme4)         # multilevel modelling
library(arm)          # convenience functions for regression in R
library(sjPlot)       # plots effects of merMod objects, which were fitted using the lmer function of the lme4 package
library(labdsv)       # indval
library(lmerTest)
library(ggplot2)
library(plyr)
library(tidyverse)
#library(ecodist)
library(lme4)
library(dplyr)
library(tidyr)
library(corrr)
library(ggbiplot)
library("FactoMineR")
library(ggfortify)
library(cluster)


dat <- read.csv("q00200.csv", header = T)
dat2<- read.csv("q00120.csv", header = T)
env<- read.csv("q00300.csv", header = T)


names(dat) <- c("quad","spp","cov")
head(dat)
dat$ecosystem <- env$'iEcosystemID'[match(dat$'quad', env$'iQuadratID')]
dat$treat <- env$'iTreatID'[match(dat$'quad', env$'iQuadratID')]
dat$WptID <- env$'iWptID'[match(dat$'quad', env$'iQuadratID')]
dat$year <- env$'iPlantYear'[match(dat$'quad', env$'iQuadratID')]
dat2 <- filter(dat, ecosystem != '9')
plants$WptID <- as.factor(plants$WptID)
plants$treat <- as.factor(plants$treat)
plants$ecosystem <- as.factor(plants$ecosystem)
plants$quad <- as.factor(plants$quad)
plants$year <- as.factor(plants$year)
#plants$cov <- log(plants$cov)
#plants$cov <- as.numeric(plants$cov)
#plants_log <- mutate_if(plants$cov, is.numeric, log)
datmat <- dcast(dat2, quad + WptID+treat+ecosystem+year ~ spp, mean, value.var="cov", fill = 0)

names(datmat)[3] <- "ecosys"
names_ecosystem <- datmat[,3:3]

datmat2 <- datmat[,6:446]
#datmat2 <- subset(datmat, select = -c(quad))
str(datmat)

# NMDS
NMDS.count <- metaMDS(datmat2, distance = "bray", k = 2,trymax=100, autotransform = T)
stressplot(NMDS.count)
plot(NMDS.count)
plot(NMDS.count$points)
plot(NMDS.count, type = "n")
orditorp(NMDS.count, display = "sites", labels = F, pch = 15, col = c("green", "blue", "red","black") [as.factor(names_ecosystem)], cex = 1)

data.scores.gg.count = as.data.frame(scores(NMDS.count, "sites"))
data.scores.gg.count$treat = data.group2.count$sites

#ggplot Hulls
grp.a.c <- data.scores.gg.count[data.scores.gg.count$treat == "1/7", ][chull(data.scores.gg.count[data.scores.gg.count$treat == 
                                                                                                    "1/7", c("NMDS1", "NMDS2")]), ]  
grp.b.c <- data.scores.gg.count[data.scores.gg.count$treat == "7/7", ][chull(data.scores.gg.count[data.scores.gg.count$treat == 
                                                                                                    "7/7", c("NMDS1", "NMDS2")]), ]  
grp.c.c <- data.scores.gg.count[data.scores.gg.count$treat == "1/14", ][chull(data.scores.gg.count[data.scores.gg.count$treat == 
                                                                                                     "1/14", c("NMDS1", "NMDS2")]), ]  
grp.d.c <- data.scores.gg.count[data.scores.gg.count$treat == "7/14", ][chull(data.scores.gg.count[data.scores.gg.count$treat == 
                                                                                                     "7/14", c("NMDS1", "NMDS2")]), ]  
grp.e.c <- data.scores.gg.count[data.scores.gg.count$treat == "1/21", ][chull(data.scores.gg.count[data.scores.gg.count$treat == 
                                                                                                     "1/21", c("NMDS1", "NMDS2")]), ]  
grp.f.c <- data.scores.gg.count[data.scores.gg.count$treat == "7/21", ][chull(data.scores.gg.count[data.scores.gg.count$treat == 
                                                                                                     "7/21", c("NMDS1", "NMDS2")]), ]  
grp.g.c <- data.scores.gg.count[data.scores.gg.count$treat == "Dry", ][chull(data.scores.gg.count[data.scores.gg.count$treat == 
                                                                                                    "Dry", c("NMDS1", "NMDS2")]), ]  
grp.h.c <- data.scores.gg.count[data.scores.gg.count$treat == "Flood", ][chull(data.scores.gg.count[data.scores.gg.count$treat == 
                                                                                                      "Flood", c("NMDS1", "NMDS2")]), ]  


hull.data.count <- rbind(grp.a.c, grp.b.c, grp.c.c, grp.d.c, grp.e.c, grp.f.c, grp.g.c, grp.h.c)

hull.data.count

species.scores.count <- as.data.frame(scores(NMDS.count, "species"))  
species.scores.count$species <- rownames(species.scores.count)  

# Explore NMDS
sig.species <- envfit(NMDS.count, all.sp2.count, permutations = 999)
head(sig.species)
species.scores.count <- cbind(species.scores.count, pval = sig.species$vectors$pvals)
sig.spp.scrs <- subset(species.scores.count, pval<=0.05)
head(sig.spp.scrs)

pca_res <- prcomp(datmat, scale. = TRUE)
autoplot(pca_res)
autoplot(pca_res, data = datmat, colour = 'treat', frame = TRUE)

summary(pca_res)
round(pca_res$sdev^2,2)
screeplot(pca_res)

rownames(datmat) <- datmat[,1]
datmat <- datmat[,-1]
datmat[1:20,1:10]

dat2[1:10,1:10]
rownames(dat2) <- dat2[,1]
dat2 <- dat2[,-1]
dat2[1:10,1:10]

head(env)
env$iWptID <- factor(env$iWptID)
env$iQuadratID <- factor(env$iQuadratID)
env$iEcosystemID <- factor(env$iEcosystemID)
env$iTreatID <- factor(env$iTreatID)
env$Age <- mapply(function(x) if(x == "Remnant") "Remnant" else (2015-as.numeric(x))
                  , x = as.vector(env$iPlantYear)
)

env$iPlantYear <- factor(env$iPlantYear)

env$iPlantYear <- ordered(env$iPlantYear, levels=c("Remnant","2012","2013","2014","2015"))

samples <- rep(seq(1:9),100)
env <- env[order(env$iWptID,env$iQuadratID),]
env$sample <- as.factor(samples)

#
# Check species lists are the same between site*species matrix and trait matrix  
#

dat.species <- levels(as.factor(dat$spp)) 
dat2.species <- rownames(dat2)

length(levels(dat.species)) == length(levels(dat2.species))

setdiff(levels(dat.species), levels(dat2.species)) # null means no difference in species between dat.sp and dat2.sp
setdiff(levels(dat2.species), levels(dat.species))# 

################################################################Bray Curtis analysis##############################################################

##datmat<-dplyr::mutate(datmat(x^1/4)) - trying to apply overall transformation!

sim.method = 'bray'
dat.dist <- vegdist(datmat, method=sim.method)
head(dat.dist)

Bcdist <- as.matrix(dat.dist)
as.data.frame(Bcdist)
head(Bcdist)
A<-melt(Bcdist)
A$Rowyear<-env$iPlantYear[A$Var1]
A$Colyear<-env$iPlantYear[A$Var2]
#A$Rowyear<-env$Age[A$Var1]
#A$Colyear<-env$Age[A$Var2]

A$RowSite<-env$iWptID[A$Var1]
A$ColSite<-env$iWptID[A$Var2]
A$Rowecosystem<-env$iEcosystemID[A$Var1] #datmat$iEcosystemID !=9
A$Colecosystem<-env$iEcosystemID[A$Var2]
A<-A[A$Rowyear=="Remnant",]  #only compare against remnant
A<-A[A$Colyear!="Remnant",] #not remnant against remnant
A<-A[A$Rowecosystem==A$Colecosystem,] #only same ecosystems
A<- A[A$Rowecosystem !=9,]
A<- A[A$Colecosystem !=9,]
A<-na.omit(A)

head(A)
str(A)

levels(A$Colyear)
levels(A$Colecosystem)


#Plot the rowecosystem, the col year as it will always be against remnant. 

hist(A$value) #check for normality
A$trans_value <-asin(A$value) # tried many transformations but none have improved the distribution much. need to move to glm
hist(A$trans_value)

#A<-as.factor(A$ColSite)

A%>%
  filter(!is.na(value))%>%
  group_by(Rowecosystem)%>%
  ggplot(aes(x = Colyear, y = value)) + geom_boxplot() + facet_wrap( ~ Rowecosystem)

########################################GLM for Bray Curtis ###################################  

#glmBC <- lmer(A$value ~ as.numeric(A$Colyear) + A$Colecosystem, family = Gamma(link ="log"))
#summary(glmBC)
#plot


#glmBC2<-(lmer(A$value ~ as.numeric(Colyear) + Colecosystem + (1 | ColSite), data = A[A$Colyear != "9",], family =poisson))
#summary(glmBC)

BCGLM<- lmer(value ~ Colyear+ Colecosystem + (1 | ColSite), A)  
summary(BCGLM)
plot(BCGLM)

BCGLM2<- lmer(value ~ Colyear + (1 | ColSite), A)  
summary(BCGLM2)
plot(BCGLM2)

levels(A$Colyear)
levels(A$ColSite)
levels(A$Colecosystem)

#################################### FD analyses######################################################################  

dat.gowdis <- gowdis(datmat[-1,])
dat.functcomp <- functcomp(dat2, as.matrix(datmat))
dat.dbFD <- dbFD(dat2
                 , as.matrix(datmat)
                 , print.pco = TRUE
                 , corr = "lingoes"
)

dat.dbFD.2 <- dbFD(dat2, as.matrix(datmat)
                   , corr = "cailliez"
                   , calc.FGR = TRUE
                   , clust.type = "kmeans"
                   #, km.sup.gr = 10
)

# create dataframe of FDis results

FDis.result <- data.frame(dat.dbFD$FDis) # create column FDis
FDis.result <- data.frame(iQuadratID=rownames(FDis.result),FDis.result)
names(FDis.result) <- c("iQuadratID","datFDis")

nbsp.result <- data.frame(dat.dbFD$nbsp) #nbsp is vector listing the number of species in each community
nbsp.result <- data.frame(iQuadratID=rownames(nbsp.result),nbsp.result)
names(nbsp.result) <- c("iQuadratID","nbsp")

singsp.result <- data.frame(dat.dbFD$sing.sp) #sind.sp is vector listing the number of functionally singular species in each community. If all species are functionally different, sing.sp will be identical to nbsp.
singsp.result <- data.frame(iQuadratID=rownames(singsp.result),singsp.result)
names(singsp.result) <- c("iQuadratID","singsp")

FRic.result <- data.frame(dat.dbFD$FRic) # vector specifying functional eveness of each community
FRic.result <- data.frame(iQuadratID=rownames(FRic.result),FRic.result)
names(FRic.result) <- c("iQuadratID","FRic") 
head(FRic.result)

Func.result <- merge(env,FDis.result)
Func.result <- merge(Func.result,nbsp.result)
Func.result <- merge(Func.result,singsp.result)
Func.result <- merge(Func.result,FRic.result)
head(Func.result)

Func.result$iPlantYear <- relevel(Func.result$iPlantYear, ref = 5)

Func.result <- Func.result[order(Func.result$iWptID,Func.result$iQuadratID),]


head(Func.result)

#Func.result <- within(Func.result, sample <- factor(iWptID:iQuadratID))

#setwd("//env/cllmm/Files/CLL Restoration/MA 6_Research and Monitoring/Projects/VP2015/Veg")

write.csv(Func.result, file = "Funcdiv3.csv")

##################################  Remove ecosystem 9  ############################# 

Func.result <- Func.result[Func.result$iEcosystemID !=9,]

Func.result$iEcosystemID <- factor(Func.result$iEcosystemID)

hist(Func.result$datFDis)

# examine nature of FDis  

png(file=file.path("fig",paste("fig", "FDis_Histogram2.png",sep=""))
    , width=3000, height=3000, res=300
    , family = "Segoe UI"
    , bg = "transparent"
)

hist(Func.result$datFDis
     , main = "Histogram of FDis"
     , xlab = "Value of FDis"
)

dev.off()


hist(Func.result$datFDis) # may actually need a transformation but check with final data

#Analyse Fdis

dat.lmer.ref <- lmer(datFDis ~ iEcosystemID + (1 | iWptID), data = Func.result[Func.result$iPlantYear == "Remnant",])
display(dat.lmer.ref)
summary(dat.lmer.ref)
plot(dat.lmer.ref)



# including or excluding remnant (sub excludes remnant)

Func.result.sub <- Func.result[Func.result$iPlantYear != "Remnant",] #removal of remnant? check logic??
Func.result.sub$Age <- as.numeric(Func.result.sub$Age)


dat.lmer <- lmer(datFDis ~ Age*iEcosystemID + (1 | iWptID), data = Func.result.sub)
display(dat.lmer)
summary(dat.lmer)
plot(dat.lmer)

dat.lmer2 <- lmer(datFDis ~ Age + iEcosystemID + (1 | iWptID), data = Func.result.sub)
display(dat.lmer)
summary(dat.lmer)
plot(dat.lmer)

#plot FDis using dplyr for preparation
AveFDis3 <-ddply(Func.result, .(iEcosystemID, iPlantYear), summarise,  datFDis=mean(datFDis))
SDFDIs<-ddply(Func.result, .(iEcosystemID, iPlantYear), summarise, SDdatFDis=sd(datFDis))

AveFDis4<- merge(AveFDis3,SDFDIs)

ggplot(AveFDis4, aes(x=iPlantYear,iEcosystemID, y=datFDis)) + geom_errorbar(ymin = AveFDis4$datFDis+AveFDis4$SDdatFDis,  ymax = AveFDis4$datFDis-AveFDis4$SDdatFDis)

