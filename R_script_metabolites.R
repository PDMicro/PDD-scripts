### METABOLITES
library(survival)
library(ggpubr)
library(PairedData)

library(reshape)
library(ggplot2)
library(ape)
library(vegan)
library(dplyr)
library(factoextra)
library(RColorBrewer)
library(cowplot)
library(tibble)

library(philentropy)
library(ggrepel)
library(pheatmap)
library(psych)
library(tidyr)
library(rstatix)
library(ggforce)
library(dendextend)
library(compositions)

# - SCFA
SCFA<- read.delim("~/R/MGSv2/metabolites/SCFA1.txt",row.names = 1, sep = "\t", header = TRUE, dec = ".")
summary(SCFA)
SCFA=(SCFA[, !colnames(SCFA) %in% c("Individual", "Pair", "Time", "Group", "Visite")])
SCFA <-as.matrix(t(SCFA))
#exclude T2 samples and dropout samples
SCFA=(SCFA[, !colnames(SCFA) %in% c("PD10T0", "C10T0","PD10T1", "C10T1", "PD4T2", "PD6T2")])

# for univariate tests see R_script_univariate_tests

#compute effect size
library(rstatix)

effPD=c()
effCO=c()

magPD=c()
magCO=c()
for(i in 1:dim(SCFA)[1]) {
  
  if(require("coin")){
    sel=Status=="PD"
    temp_dat = data.frame(metab=SCFA[i,sel], `Status`=Status[sel], `Time`=Time[sel]);
    EFPD=wilcox_effsize(temp_dat, metab ~ Time, paired = TRUE)
    effPD[i]=EFPD$effsize
    magPD[i]=EFPD$magnitude
    sel=Status=="CO"
    temp_dat = data.frame(metab=SCFA[i,sel], `Status`=Status[sel], `Time`=Time[sel]);
    EFCO=wilcox_effsize(temp_dat, metab ~ Time, paired = TRUE)
    effCO[i]=EFCO$effsize
    magCO[i]=EFCO$magnitude
  }
}

names(effCO)=rownames(SCFA)
names(effPD)=rownames(SCFA)

names(magCO)=rownames(SCFA)
names(magPD)=rownames(SCFA)

#boxplot paired samples 
library(ggplot2)
library(ggpubr)

SCFA1=as.data.frame(t(SCFA))
SCFA1 <-SCFA1[order(row.names(SCFA1)), ]

Pair=Meta1$Pair
SCFA2=cbind(SCFA1, Status, Time, Pair) 


SCFA2$status_time = paste(SCFA2$Status, SCFA2$Time, sep = "_")


p<-ggpaired(SCFA2, x = "Time", y = "Acetate", color="black", fill = "status_time", facet.by = "Status", palette = col1, title= "Acetate",             
            id="Pair",line.color = "gray", line.size = 0.4) +
  theme_classic() + ylab("Relative abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"))


p1<-ggpaired(SCFA2, x = "Time", y = "Propionate",  color="black", fill = "status_time",  facet.by = "Status", palette = col1, title= "Propionate",                
             id="Pair",line.color = "gray", line.size = 0.4) +
  theme_classic() + ylab("Relative abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


p2<-ggpaired(SCFA2, x = "Time", y = "Butyrate",  color="black", fill = "status_time",  facet.by = "Status", palette = col1, title= "Butyrate",       
             id="Pair",line.color = "gray", line.size = 0.4) +
  theme_classic() + ylab("Relative abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


p3<-ggpaired(SCFA2, x = "Time", y = "iso.Butyrate",  color="black", fill = "status_time",  facet.by = "Status", palette = col1, title= "Iso-Butyrate",       
             id="Pair",line.color = "gray", line.size = 0.4) +
  theme_classic() + ylab("Relative abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


p4<-ggpaired(SCFA2, x = "Time", y = "Valerate",  color="black", fill = "status_time", facet.by = "Status", palette = col1, title= "Valerate",             
             id="Pair",line.color = "gray", line.size = 0.4)+
  theme_classic() + ylab("Relative abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


p5<-ggpaired(SCFA2, x = "Time", y = "iso.Valerate",  color="black", fill = "status_time", facet.by = "Status", palette = col1, title= "Iso-Valerate",       
             id="Pair",line.color = "gray", line.size = 0.4)+
  stat_compare_means(label = "p.format",paired = TRUE) +
  theme_classic() + ylab("Relative abundance") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))



p.all <-ggarrange(p,p1,p2,p3,p4,p5,  ncol = 3, nrow = 2,common.legend=TRUE) 


#metabolites feces NMR
feces<- read.delim2("~/R/MGSv2/metabolites/feces.txt",row.names=1,sep = "\t", header = TRUE, dec = ".")
names(feces)[13] = "C9T0"

#exclude T2 samples and dropout samples
feces=(feces[, !colnames(feces) %in% c("PD10T0", "C10T0","PD10T1", "C10T1", "PD4T2", "PD6T2")])
feces=as.matrix(feces)

#remove unknowns and "Medikament"
feces_clean = feces[-(grep("NA_|UK_|Medikament", row.names(feces), perl = T, value = F)),]

#log transformation

feces_clean= log(feces_clean +1)

# for univariate tests see R_script_univariate_tests

#plot trends
Trend_pre=which(prepvals<0.05)
col=rep(rgb(0,0,0,1), length(postpvals))
col[Trend_pre]=rgb(0,0,0.6,0.5)
names(col)=names(postpvals)
pch=rep(20, length(postpvals))
pch[Trend_pre]=19

Trend_post=which(postpvals<0.05)
col[Trend_post]=rgb(0.6,1,0.8,0.5)
pch[Trend_post]=17

sel=unique(names(col[c(Trend_pre,Trend_post)]))


windows()
plot(postpvals, prepvals, col=col, pch= pch, cex=2.0, log="xy", xlim=c(0.005, 0.9), ylim=c(0.005, 0.9), 
     main = "Faecal SCFA equalize between PD and CO after prebiotics", bty="L"
     , xlab = "after prebiotics (p-value)", ylab = "before prebiotics (p-value)")
abline(h=0.05, v=0.05, lt=2, col="darkgrey")
text(postpvals[sel],prepvals[sel], names(col[sel]), cex=0.8, pos=1, col="black")

legend('bottom',c("Trend before (p<0.05, q>0.1) = SCFA","Trend after (p<0.05, q>0.1) = Amino acids"),   
       col=c(rgb(0,0,0.6,0.5), rgb(0.6,1,0.8,0.5)), cex=0.8,
       pch= c(19,17), pt.cex=1.5, bg='white',box.lty=0)

#metabolites - urine NMR
urine <- read.delim2("~/R/MGSv2/metabolites/urine.txt",row.names=1,sep = "\t", header = TRUE, dec = ".")

#remove T2 samples and negative values
urine=(urine[, !colnames(urine) %in% c("PD10T0", "C10T0","PD10T1", "C10T1", "PD4T2", "PD6T2")])

urine_filter <- urine                    
urine_filter[urine_filter < 0] <- NA      
urine_filter2 <- na.omit(urine_filter)         

negatives=which(urine_filter2[,]<0)

urine_clean = urine_filter2[-grep("NA_|UK", row.names(urine_filter2)),]
urine_clean=as.matrix(urine_clean)

name_match = match(row.names(metadata2), colnames(urine_clean))
urine_clean = urine_clean[,name_match]

# log
urine_clean = log(urine_clean +1)

# for univariate tests see R_script_univariate_tests

#dbRDA with conditioning out family/household
dbRDA7=dbrda(t(urine_clean)~  Condition(family), metadata2, distance = "euclidean", sqrt.dist =FALSE ,add = FALSE, dfun = vegdist, metaMDSdist = FALSE,na.action = na.fail)

sppscores(dbRDA7) <- t(urine_clean)
dbRDA7 


#visualize
#extract scores
dbRDA7_scores = scores(dbRDA7, scaling = 3, correlation = TRUE)

dbRDA7_site = data.frame(dbRDA7_scores$sites)
name_match = match(row.names(dbRDA7_site), row.names(metadata2))
dbRDA7_site[,c(3:4)] =  metadata2[name_match,c(2,3)]
dbRDA7_site$status_time = as.factor(paste(dbRDA7_site$status, dbRDA7_site$time, sep = "_"))

#extract species
dbRDA7_species = data.frame(dbRDA7_scores$species)
dbRDA7_species$length = round(sqrt(dbRDA7_species[,1]^2 + dbRDA7_species[,2]^2), digits = 3)
name_match = match(row.names(dbRDA7_species), row.names(t(urine_clean)))
dbRDA7_species = dbRDA7_species %>% arrange(desc(length))


#plot
windows()
ggplot(dbRDA7_site, aes(x = MDS1, y = MDS2, color=status_time, shape=status_time)) +
  geom_point(size = 3.5, alpha = 0.5)  +
  scale_color_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  scale_shape_manual(values = c(19,19,17,17),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  stat_ellipse() +
  facet_wrap(~time) +
  theme_minimal_grid() +
  ggtitle("Dissimilarity in urine NMR metabolites \n between PD and CO before and after prebiotics")+
  theme(plot.title = element_text(hjust = 0.5))



#get data
dbRDA7_site[,6]=metadata2[,9]
colnames(dbRDA7_site)[colnames(dbRDA7_site) == "V6"] <- "family"
dbRDA7_site$family = factor(dbRDA7_site$family, levels = unique(dbRDA7_site$family))


#test for differences
library(permute)
adonis2(vegdist(t(urine_clean), method = "euclidean") ~ time + status + family , data = dbRDA7_site, permutations = 1000, parallel = 8)


#block family
set.seed(1)
perm <- how(nperm = 1000)
setBlocks(perm) <- with(dbRDA7_site, family)
adonis2(vegdist(t(urine_clean), method = "euclidean") ~ status, data = dbRDA7_site, permutations = perm, parallel = 8)


set.seed(1)
perm <- how(nperm = 1000)
setBlocks(perm) <- with(dbRDA7_site, family)
adonis2(vegdist(t(urine_clean), method = "euclidean") ~ time, data = dbRDA7_site, permutations = perm, parallel = 8)

#test for different timepoint per group
#pre
sel=dbRDA7_site$time=="0"
set.seed(1)
perm <- how(nperm = 1000)
setBlocks(perm) <- with(dbRDA7_site[sel,], family)
adonis2(vegdist(t(urine_clean[,sel]), method = "euclidean") ~ status, data = dbRDA7_site[sel,], permutations = perm, parallel = 8)

#post
sel=dbRDA7_site$time=="28"
set.seed(1)
perm <- how(nperm = 1000)
setBlocks(perm) <- with(dbRDA7_site[sel,], family)
adonis2(vegdist(t(urine_clean[,sel]), method = "euclidean") ~ status, data = dbRDA7_site[sel,], permutations = perm, parallel = 8)


#plot metabolites of interest
urine=as.data.frame(t(urine_clean))

MI=urine$"metabolite_of_interest"
`

windows()
boxplot(MI~factor(paste(Status, Time)),col=c("palegreen3","seagreen","lightblue","steelblue"), xlab = "Timepoint", ylab = "metabolite_of_interest", main="Urine NMR pre/post")
mtext(side = 3, line = 0.25, at = 1, adj = -1, subtitle)
stripchart(MI~factor(paste(Status, Time)), method = "jitter", vertical=TRUE, pch = 19, cex=1.5, add = TRUE, col = c("mediumseagreen", "darkgreen","lightsteelblue3","deepskyblue4"))

