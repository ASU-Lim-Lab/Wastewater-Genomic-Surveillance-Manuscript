setwd('~/wastewater/')
library(RColorBrewer)
library(gplots)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(png)
library(vegan)
library(tidyverse)
library(tidyr)
library(nlme)
library(phyloseq)
library(ape)
library(corrplot)
library(taxonomizr)
library("qvalue")
library(glmnet)
library(Rcpp)
library(decontam)
library(knitr)
library(broom)
library(lmerTest)

#AlphaDiversity
d.1<-read.csv("VSP_RPM_vegan_in_final.csv", row.names = 1, header = TRUE)
data_vegan<-diversity(d.1, index = 'shannon')
write.table(data_vegan,"VSP_RPM_AlphaDiversity_112923.txt", sep = '\t')

#BetaDiversity
m.1 <- read.csv("VSP_metadata_vegan_in_final.csv", header=TRUE)
m.1$Site<-as.character(m.1$Site)
d.1<-read.csv("VSP_RPM_vegan_in_final_relativeAbundance.csv", row.names = 1, header = TRUE)
dataTransposed<-t(d.1)
#m.2<-m.1[!m.1$Richness=="#N/A",]
#d.2<-d.1[,colnames(d.1) %in% m.2$Sample]
beta_dist <- vegdist(t(dataTransposed),index = "bray")
distMatrix <- as.matrix(beta_dist)
write.csv(distMatrix,"weightedBeta_VSP_RPM_112923.csv")

library(reshape2)
df <- melt(as.matrix(beta_dist))
p    <- t(apply(df[,c(1,2)],1,FUN=sort))
rmv1 <- which(p[,1] == p[,2])
p    <- paste(p[,1],p[,2],sep="|")
rmv2 <- which(duplicated(p))
df   <- df[-c(rmv1,rmv2),] 
write.csv(df,"paired_WeightedBeta_120523.csv")

#k-means
library(factoextra)
library(cluster)
#calculate optimal number of clusters (sum of squares)
WUDM<-read.csv("weightedBeta_VSP_RPM_112923.csv",row.names=1,header = T)
a<-fviz_nbclust(WUDM, pam, method = "wss")
a
b<-fviz_nbclust(WUDM, pam, method = "silhouette")
b
#number of clusters vs. gap statistic
gap_stat <- clusGap(WUDM,
                    FUN = pam,
                    K.max = 20, #max clusters to consider
                    B = 50) #total bootstrapped iterations
fviz_gap_stat(gap_stat)

#kmeans clusters
set.seed(123)
km.res3<-kmeans(WUDM, 3, iter.max = 10, nstart = 25)
km.res3
k3<-as.table(km.res3$cluster)
write.csv(k3,"k3.csv")
set.seed(123)
km.res4<-kmeans(WUDM, 4, iter.max = 10, nstart = 25)
km.res4
k4<-as.table(km.res4$cluster)
write.csv(k4,"k4.csv")
set.seed(123)
km.res5<-kmeans(WUDM, 5, iter.max = 10, nstart = 25)
km.res5
k5<-as.table(km.res5$cluster)
write.csv(k5,"k5.csv")
set.seed(123)


#multinomial kmeans4
#2021 winter
library(mclogit)
m.1 <- read.csv("VSP_metadata_vegan_in_final.csv", header=TRUE)
m.1$Site<-as.factor(m.1$Site)
m.1$Meteorological_Season<-as.factor(m.1$Meteorological_Season)
m.1$kmeans3<-as.factor(m.1$kmeans3)
m.1$kmeans4<-as.factor(m.1$kmeans4)
m.1$kmeans5<-as.factor(m.1$kmeans5)

Subset1<-subset(m.1,kmeans4 =='2'|kmeans4 =='3'|kmeans4 =='4')
Subset2<-subset(m.1,kmeans4 =='3'|kmeans4 =='4')

set.seed(123)
(comm.mblogit <- mblogit(kmeans4~Meteorological_Season+Week,  data = m.1, random=~1|Site))
summary(comm.mblogit)
set.seed(123)
(comm.mblogit.subset1 <- mblogit(kmeans4~Meteorological_Season+Week,  data = Subset1, random=~1|Site))
summary(comm.mblogit.subset1)
set.seed(123)
(comm.mblogit.subset2 <- mblogit(kmeans4~Meteorological_Season+Week,  data = Subset2, random=~1|Site))
summary(comm.mblogit.subset2)

library(mclogit)
#2022_fall control
m.1 <- read.csv("VSP_metadata_vegan_in_final.csv", header=TRUE)
m.1$Site<-as.factor(m.1$Site)
m.1$Meteorological_Season2.<-m.1$Meteorological_Season
m.1$kmeans3<-as.factor(m.1$kmeans3)
m.1$kmeans4<-as.factor(m.1$kmeans4)
m.1$kmeans5<-as.factor(m.1$kmeans5)
m.1$Meteorological_Season2.[m.1$Meteorological_Season2.=='2021_Winter']<-'2023_Winter'
m.1$Meteorological_Season<-as.factor(m.1$Meteorological_Season)
m.1$Meteorological_Season2.<-as.factor(m.1$Meteorological_Season2.)

Subset1<-subset(m.1,kmeans4 =='2'|kmeans4 =='3'|kmeans4 =='4')
Subset2<-subset(m.1,kmeans4 =='3'|kmeans4 =='4')

set.seed(123)
(comm.mblogit <- mblogit(kmeans4~Meteorological_Season2.+Week,  data = m.1, random=~1|Site))
summary(comm.mblogit)
set.seed(123)
(comm.mblogit.subset1 <- mblogit(kmeans4~Meteorological_Season2.+Week,  data = Subset1, random=~1|Site))
summary(comm.mblogit.subset1)
set.seed(123)
(comm.mblogit.subset2 <- mblogit(kmeans4~Meteorological_Season2.+Week,  data = Subset2, random=~1|Site))
summary(comm.mblogit.subset2)

#2022_spring control
m.1 <- read.csv("VSP_metadata_vegan_in_final.csv", header=TRUE)
m.1$Site<-as.factor(m.1$Site)
m.1$Meteorological_Season2.<-m.1$Meteorological_Season
m.1$kmeans3<-as.factor(m.1$kmeans3)
m.1$kmeans4<-as.factor(m.1$kmeans4)
m.1$kmeans5<-as.factor(m.1$kmeans5)
m.1$Meteorological_Season2.[m.1$Meteorological_Season2.=='2021_Winter']<-'2023_winter'
m.1$Meteorological_Season2.[m.1$Meteorological_Season2.=='2022_Fall']<-'2023_fall'
m.1$Meteorological_Season<-as.factor(m.1$Meteorological_Season)
m.1$Meteorological_Season2.<-as.factor(m.1$Meteorological_Season2.)

Subset1<-subset(m.1,kmeans4 =='2'|kmeans4 =='3'|kmeans4 =='4')
Subset2<-subset(m.1,kmeans4 =='3'|kmeans4 =='4')

set.seed(123)
(comm.mblogit <- mblogit(kmeans4~Meteorological_Season2.+Week,  data = m.1, random=~1|Site))
summary(comm.mblogit)
set.seed(123)
(comm.mblogit.subset1 <- mblogit(kmeans4~Meteorological_Season2.+Week,  data = Subset1, random=~1|Site))
summary(comm.mblogit.subset1)
set.seed(123)
(comm.mblogit.subset2 <- mblogit(kmeans4~Meteorological_Season2.+Week,  data = Subset2, random=~1|Site))
summary(comm.mblogit.subset2)

#2022_summer control
m.1 <- read.csv("VSP_metadata_vegan_in_final.csv", header=TRUE)
m.1$Site<-as.factor(m.1$Site)
m.1$Meteorological_Season2.<-m.1$Meteorological_Season
m.1$kmeans3<-as.factor(m.1$kmeans3)
m.1$kmeans4<-as.factor(m.1$kmeans4)
m.1$kmeans5<-as.factor(m.1$kmeans5)
m.1$Meteorological_Season2.[m.1$Meteorological_Season2.=='2021_Winter']<-'2023_winter'
m.1$Meteorological_Season2.[m.1$Meteorological_Season2.=='2022_Fall']<-'2023_fall'
m.1$Meteorological_Season2.[m.1$Meteorological_Season2.=='2022_Spring']<-'2023_spring'
m.1$Meteorological_Season<-as.factor(m.1$Meteorological_Season)
m.1$Meteorological_Season2.<-as.factor(m.1$Meteorological_Season2.)

Subset1<-subset(m.1,kmeans4 =='2'|kmeans4 =='3'|kmeans4 =='4')
Subset2<-subset(m.1,kmeans4 =='3'|kmeans4 =='4')

set.seed(123)
(comm.mblogit <- mblogit(kmeans4~Meteorological_Season2.+Week,  data = m.1, random=~1|Site))
summary(comm.mblogit)
set.seed(123)
(comm.mblogit.subset1 <- mblogit(kmeans4~Meteorological_Season2.+Week,  data = Subset1, random=~1|Site))
summary(comm.mblogit.subset1)
set.seed(123)
(comm.mblogit.subset2 <- mblogit(kmeans4~Meteorological_Season2.+Week,  data = Subset2, random=~1|Site))
summary(comm.mblogit.subset2)

#2022_winter control
m.1 <- read.csv("VSP_metadata_vegan_in_final.csv", header=TRUE)
m.1$Site<-as.factor(m.1$Site)
m.1$Meteorological_Season2.<-m.1$Meteorological_Season
m.1$kmeans3<-as.factor(m.1$kmeans3)
m.1$kmeans4<-as.factor(m.1$kmeans4)
m.1$kmeans5<-as.factor(m.1$kmeans5)
m.1$Meteorological_Season2.[m.1$Meteorological_Season2.=='2021_Winter']<-'2023_winter'
m.1$Meteorological_Season2.[m.1$Meteorological_Season2.=='2022_Fall']<-'2023_fall'
m.1$Meteorological_Season2.[m.1$Meteorological_Season2.=='2022_Spring']<-'2023_spring'
m.1$Meteorological_Season2.[m.1$Meteorological_Season2.=='2022_Summer']<-'2023_summer'
m.1$Meteorological_Season<-as.factor(m.1$Meteorological_Season)
m.1$Meteorological_Season2.<-as.factor(m.1$Meteorological_Season2.)

Subset1<-subset(m.1,kmeans4 =='2'|kmeans4 =='3'|kmeans4 =='4')
Subset2<-subset(m.1,kmeans4 =='3'|kmeans4 =='4')

set.seed(123)
(comm.mblogit <- mblogit(kmeans4~Meteorological_Season2.+Week,  data = m.1, random=~1|Site))
summary(comm.mblogit)
set.seed(123)
(comm.mblogit.subset1 <- mblogit(kmeans4~Meteorological_Season2.+Week,  data = Subset1, random=~1|Site))
summary(comm.mblogit.subset1)
set.seed(123)
(comm.mblogit.subset2 <- mblogit(kmeans4~Meteorological_Season2.+Week,  data = Subset2, random=~1|Site))
summary(comm.mblogit.subset2)


################Alpha diversity################
setwd("/Users/rmaqsood/Documents/wastewater/")
m.1 <- read.csv("VSP_metadata_vegan_in_final.csv", header=TRUE)
m.1$Site<-as.character(m.1$Site)
m.1$Meteorological_Season<-as.character(m.1$Meteorological_Season)
#season significance
set.seed(100)
allrichness<-lme(Richness~  Season*Week , random=~1|Site, data=m.1)
summary(allrichness)
set.seed(100)
allalphaDiversity<-lme(AlphaDiversity~ Season*Week , random=~1|Site, data=m.1)
summary(allalphaDiversity)

#season significance
set.seed(100)
allrichness<-lme(Richness~  Meteorological_Season , random=~1|Site, data=m.1)
summary(allrichness)
set.seed(100)
allalphaDiversity<-lme(AlphaDiversity~ Meteorological_Season , random=~1|Site, data=m.1)
summary(allalphaDiversity)

#plots seasons
loessPlotRichness <- ggplot(m.1, aes(x = Week, y = Richness, color = Season, fill = Season)) +
  geom_point() +scale_colour_manual(values = c("lightblue1", "springgreen2","yellow","darkorange3","darkblue"))+
  stat_smooth(method="loess", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(m.1$Week), max(m.1$Week), by =1 ),1)) + labs (x= "Week",y="Virome Richness", title = "Richness by seasons" )
loessPlotRichness


loessPlotShannon <- ggplot(m.1, aes(x = Week, y = AlphaDiversity, color = Season, fill = Season)) +
  geom_point() +scale_colour_manual(values = c("lightblue1", "springgreen2","yellow","darkorange3","darkblue"))+
  stat_smooth(method="loess", se=TRUE, span=0.5, level=0.95)+
  theme_bw()+ scale_x_continuous(breaks = round(seq(min(m.1$Week), max(m.1$Week), by =1 ),1)) + labs (x= "Week",y="Virome Alpha Diversity", title = "Alpha diversity by seasons" )
loessPlotShannon


library(lme4)
library(lmerTest)
set.seed(123)
model1 <- lmer(AlphaDiversity ~ Meteorological_Season + Week+ (1|Site), data = m.1, REML = FALSE)
model1
set.seed(123)
model1 %>% anova()

library(lme4)
library(lmerTest)
set.seed(123)
model2 <- lmer(Richness ~ Meteorological_Season + Week+ (1|Site), data = m.1, REML = FALSE)
model2
set.seed(123)
model2 %>% anova()

library(readxl)
VSP_Analysis_112923 <- read_excel("VSP_Analysis_112923.xlsx", sheet = "weigtbetaPaired")
#all
scatter<-ggplot(VSP_Analysis_112923, aes(x=`Absolute Weeks Apart`, y=`Bray-Curtis distance`)) + 
  geom_point(color="lightpink", alpha=0)+geom_smooth(colour="black")+theme_classic()+labs(title="All sites")
scatter
#samesite
samesite<-VSP_Analysis_112923[which(VSP_Analysis_112923$samesite=='TRUE'),]
c<-samesite[which(samesite$site1=='C'),]
scatter.1<-ggplot(c, aes(x=`Absolute Weeks Apart`, y=`Bray-Curtis distance`)) + 
  geom_point(color="lightpink")+ylim(0, 1)+xlim(0, 51)+geom_smooth(colour="black")+theme_classic()+labs(title="Site C")
scatter.1
DT<-samesite[which(samesite$site1=='DT'),]
scatter.2<-ggplot(DT, aes(x=`Absolute Weeks Apart`, y=`Bray-Curtis distance`)) + 
  geom_point(color="lightpink")+ylim(0, 1)+xlim(0, 51)+geom_smooth(colour="black")+theme_classic()+labs(title="Site DT")
scatter.2
GUAD<-samesite[which(samesite$site1=='GUAD'),]
scatter.3<-ggplot(GUAD, aes(x=`Absolute Weeks Apart`, y=`Bray-Curtis distance`)) + 
  geom_point(color="lightpink")+ylim(0, 1)+xlim(0, 51)+geom_smooth(colour="black")+theme_classic()+labs(title="Site GUAD")
scatter.3
MEO1<-samesite[which(samesite$site1=='MEO1'),]
scatter.4<-ggplot(MEO1, aes(x=`Absolute Weeks Apart`, y=`Bray-Curtis distance`)) + 
  geom_point(color="lightpink")+ylim(0, 1)+xlim(0, 51)+geom_smooth(colour="black")+theme_classic()+labs(title="Site MEO1")
scatter.4
MEO2<-samesite[which(samesite$site1=='MEO2'),]
scatter.5<-ggplot(MEO2, aes(x=`Absolute Weeks Apart`, y=`Bray-Curtis distance`)) + 
  geom_point(color="lightpink")+ylim(0, 1)+xlim(0, 51)+geom_smooth(colour="black")+theme_classic()+labs(title="Site MEO2")
scatter.5
R<-samesite[which(samesite$site1=='R'),]
scatter.6<-ggplot(R, aes(x=`Absolute Weeks Apart`, y=`Bray-Curtis distance`)) + 
  geom_point(color="lightpink")+ylim(0, 1)+xlim(0, 51)+geom_smooth(colour="black")+theme_classic()+labs(title="Site R")
scatter.6
STLUKES<-samesite[which(samesite$site1=='STLUKES'),]
scatter.7<-ggplot(STLUKES, aes(x=`Absolute Weeks Apart`, y=`Bray-Curtis distance`)) + 
  geom_point(color="lightpink")+ylim(0, 1)+xlim(0, 51)+geom_smooth(colour="black")+theme_classic()+labs(title="Site STLUKES")
scatter.7
TP01<-samesite[which(samesite$site1=='TP01'),]
scatter.8<-ggplot(TP01, aes(x=`Absolute Weeks Apart`, y=`Bray-Curtis distance`)) + 
  geom_point(color="lightpink")+ylim(0, 1)+xlim(0, 51)+geom_smooth(colour="black")+theme_classic()+labs(title="Site TP01")
scatter.8
TP02<-samesite[which(samesite$site1=='TP02'),]
scatter.9<-ggplot(TP02, aes(x=`Absolute Weeks Apart`, y=`Bray-Curtis distance`)) + 
  geom_point(color="lightpink")+ylim(0, 1)+xlim(0, 51)+geom_smooth(colour="black")+theme_classic()+labs(title="Site TP02")
scatter.9
TP03<-samesite[which(samesite$site1=='TP03'),]
scatter.10<-ggplot(TP03, aes(x=`Absolute Weeks Apart`, y=`Bray-Curtis distance`)) + 
  geom_point(color="lightpink")+ylim(0, 1)+xlim(0, 51)+geom_smooth(colour="black")+theme_classic()+labs(title="Site TP03")
scatter.10
TP04<-samesite[which(samesite$site1=='TP04'),]
scatter.11<-ggplot(TP04, aes(x=`Absolute Weeks Apart`, y=`Bray-Curtis distance`)) + 
  geom_point(color="lightpink")+ylim(0, 1)+xlim(0, 51)+geom_smooth(colour="black")+theme_classic()+labs(title="Site TP04")
scatter.11
TP05<-samesite[which(samesite$site1=='TP05'),]
scatter.12<-ggplot(TP05, aes(x=`Absolute Weeks Apart`, y=`Bray-Curtis distance`)) + 
  geom_point(color="lightpink")+ylim(0, 1)+xlim(0, 51)+geom_smooth(colour="black")+theme_classic()+labs(title="Site TP05")
scatter.12

###PCoA and PERMANOVA
set.seed(100)

d.1 <-  read_csv("VSP_RPM_vegan_in_final.csv") %>% 
  column_to_rownames("Sequence_ID")

d.1 <- d.1[colSums(abs(d.1), na.rm = TRUE) > 0]

m.1 <- read.csv("VSP_metadata_vegan_in_final.csv", header=TRUE) %>% 
  column_to_rownames(var="Sequence_ID")

feature_table <- rotate_df(d.1, cn = F)
feature_table <- as.data.frame(feature_table)
OTU = otu_table(feature_table,taxa_are_rows = TRUE)

taxtable <- read.csv("taxtable.csv") 
taxtable <- taxtable %>%
  distinct() %>% 
  remove_rownames %>%
  column_to_rownames(var="X")
taxmat <- as.matrix(taxtable)
TAX = tax_table(taxmat)

sampledata = sample_data(m.1) 

physeq1_RPM = merge_phyloseq(OTU, TAX, sampledata) 
OTU_relativeAbundance = transform_sample_counts(OTU, function(x) x / sum(x))
physeq2_relativeAbundance = merge_phyloseq(OTU_relativeAbundance, TAX, sampledata) 
OTU_filter1 =  filter_taxa(physeq2_relativeAbundance, function(x) sum(x) > 0.005, prune = TRUE)
physeq3_nondetectsFiltered = merge_phyloseq(OTU_filter1, TAX, sampledata)
physeq4_speciesAgglom = tax_glom(physeq3_nondetectsFiltered, "species")
sample_data(physeq4_speciesAgglom)$Meteorological_Season <- factor(sample_data(physeq4_speciesAgglom)$Meteorological_Season,
                                                                   levels = c("2021_Winter", "2022_Spring", "2022_Summer","2022_Fall","2022_Winter"))
#Fig2B PCoA
ord2 <- ordinate(physeq4_speciesAgglom,
                 method = "PCoA",
                 distance = "bray")
PCoA_Season_species  <- plot_ordination(physeq = physeq4_speciesAgglom,
                                        ordination = ord2,
                                        color = "Meteorological_Season", # metadata variable
                                        #shape = "Site", # metadata variable
                                        axes = c(1,2),
                                        title='Weighted Bray-Curtis; Species-level taxonomic agglomeration') +
  theme_bw() +
  theme(text=element_text(size=20)) +
  scale_shape_manual(values = (1:12)) +
  stat_ellipse(type = "t") +
  geom_point(size = 2) +
  coord_equal()

PCoA_Season_species

m.1 <- m.1 %>%
  rownames_to_column(var="Sequence_ID")

ord2_PCoA_eigenvectors.df <- as.data.frame(ord2$vectors[,1:2]) %>% rownames_to_column(., "Sequence_ID") %>% left_join(.,m.1, by = "Sequence_ID")
dim(ord2_PCoA_eigenvectors.df)
write.csv(x = ord2_PCoA_eigenvectors.df, file = "PCoA_eigenvectors.csv")

#Fig2B PERMANOVA
relAbundance <- as.matrix(physeq4_speciesAgglom@otu_table)
relAbundance <- as.data.frame(relAbundance)
relAbundance <- rotate_df(relAbundance, cn = F)
brayDist <- vegdist(relAbundance,
                    method = 'bray')
brayDist.mat <- as.matrix(brayDist)
brayDist.df <- as.data.frame(brayDist.mat)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = m.1$Site))
adonis2(brayDist ~  Week + Meteorological_Season, data = m.1, permutations = perm, by = "margin")
