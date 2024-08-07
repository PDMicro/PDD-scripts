#univariate tests

#groups differences (PD vs CO)
# load an prepare input data 

#Genus level
M <- read.delim("~/R/MGSv2/specI.Genus",sep = "\t", header = TRUE, row.names = 1)
#remove T2 samples and PD/CO10 (dropouts)
M=M[, !colnames(M) %in% c("PD4T2", "PD6T2","PD10T1", "C10T1", "PD10T0", "C10T0")]
#remove"?" unknowns
M=M[!rownames(M) %in% c("?"),]
M=as.matrix(M)
#normalize
M1 = sweep(M,2,colSums(M),"/")

# Filter taxa present in less than 4 samples and taxa with relative abundance less than 0.001%
M1= M1[rowMeans(M1)>0.001 & rowSums(M1>0)>4, ]

#metadata subsets/grouping variables
Time=as.factor(metadata2$time)
Status=as.factor(metadata2$status)

### Group differences with wilcox test pre and post diet

prepvals=prestats=c()  
postpvals=poststats=c() 

for(i in 1:dim(M1)[1]) {
  temp_dat = data.frame(taxon=M1[i,], `Status`=Status, `Time`=Time);
  sel=Time=="0"
  wt=wilcox.test(M1[i, sel]~ Status[sel], paired=TRUE)
  prepvals[i]= wt$p.value
  prestats[i]=  wt$statistic
  sel=Time=="28"
  wt=wilcox.test(M1[i, sel]~ Status[sel], paired=TRUE)
  postpvals[i]= wt$p.value
  poststats[i]=  wt$statistic
} 

names(prepvals)=rownames(M1)
names(postpvals)=rownames(M1)

preqvals <- p.adjust(prepvals,"BH")
postqvals <- p.adjust(postpvals,"BH")

#make subsets of interesting pvals/qvals

sum(preqvals<0.1, na.rm=TRUE)
sum(postqvals<0.1, na.rm=TRUE)

idx1=which(preqvals<0.1 & prepvals<0.05) 

idx2=which(postqvals<0.1 & postpvals<0.05) 


#visualite results incl. trends

idx2=which(postqvals<0.1 & postpvals<0.05)
col=rep(rgb(0,0,0,1), length(postpvals))
col[idx2]=rgb(0,0,0.6,0.5)
names(col)=names(postpvals)
pch=rep(20, length(postpvals))
pch[idx2]=19

idx1a=which(prepvals<0.05 & preqvals>0.1)
col[idx1a]=rgb(0.6,1,0.8,0.5)

idx1b=which(postpvals<0.05 & postqvals>0.1)
col[idx1b]=rgb(0.4,0.6,0.8,0.5)

sel=unique(names(col[c(idx2,idx1a,idx1b)]))


windows()
plot(postpvals, prepvals, col=col, pch= pch, cex=2.0, log="xy", xlim=c(0.0001, 1), ylim=c(0.0001, 1), main = "PD-dysbiosis persists after prebiotics (Genus level)", bty="L"
     , xlab = "PD vs. CO after prebiotics (p-value)", ylab = "PD vs. CO before prebiotics (p-value)")
abline(h=0.05, v=0.05, lt=2, col="darkgrey")
text(postpvals[sel],prepvals[sel], gsub(".*;","", names(col[sel])), cex=0.7, pos=3, col="black")

legend('bottomleft',c("p<0.05, q<0.1 after","p<0.05, q<0.1 after with trend before (q>0.1)", "Trend before (q>0.1)","Trend after (q>0.1)" ),   
       col=c(rgb(0,0,0.4,0.5),rgb(0.6,1,0.8,0.5),rgb(0.6,1,0.8,0.5),rgb(0.4,0.6,0.8,0.5)), cex=0.7,
       pch=list(19,19,20,20), pt.cex=1.5, bg='white',box.lty=0)

#save Set of different Genera for visualization

GENUS_group=M1[idx2, ]
GENUS_group_pre=M1[idx1a, ]
GENUS_group_post=M1[idx1b, ]



# within group differences (PD before vs. after/CO before vs. after)

#load MGS table
M <- read.delim("~/R/MGSv2/specIMGSv2.mat", sep = "\t", header = TRUE, row.names = 1)
#remove T2 samples and PD/CO10 (dropouts)
M=M[, !colnames(M) %in% c("PD4T2", "PD6T2","PD10T1", "C10T1", "PD10T0", "C10T0")]
#remove"?" unknowns
M=M[!rownames(M) %in% c("?"),]
M=as.matrix(M)
#normalize
M1 = sweep(M,2,colSums(M),"/") 

# Filter taxa present in less than 4 samples and taxa with relative abundance less than 0.001%
M1= M1[rowMeans(M1)>0.001 & rowSums(M1>0)>4, ]

##perform wt paired with all samples, then select all q<0.1 and make post hoc test (wt paired for PD or CO)

#wt paired
pvals2=c()
for(i in 1:dim(M1)[1]) {
  temp_dat = data.frame(taxon=M1[i,], `Status`=Status, `Time`=Time);
  
  wt = wilcox.test(taxon ~ Time, data=temp_dat, paired=TRUE);
  
  pval2 <- wt$p.value
  pvals2[i]=pval2
  
}


names(pvals2)=rownames(M1)
qvals2 <- p.adjust(pvals2,"BH")

sum(qvals2<0.1, na.rm=TRUE)
idx3=which(qvals2<0.1 & pvals2<0.05) 


#select subset with q<0.1
subset1=M1[idx3, ]


#post hoc test
COpvals=COstats=c()
PDpvals=PDstats=c()
for(i in 1:dim(subset1)[1]) {
  sel=Status=="CO"
  wtCO=wilcox.test(subset1[i, sel]~ Time[sel], paired=TRUE)
  COpvals[i]= wtCO$p.value
  COstats[i]=  wtCO$statistic
  sel=Status=="PD"
  wtPD=wilcox.test(subset1[i, sel]~ Time[sel], paired=TRUE)
  PDpvals[i]= wtPD$p.value
  PDstats[i]=  wtPD$statistic
}
#load taxonomy
Tax <- read.delim("~/R/MGSv2/specI.tax",sep = "\t", header = TRUE, row.names = 1)
TaxNames1=Tax[rownames(subset1),]

names(PDpvals)=TaxNames1$Species


#select sign. taxa
idx4=which(PDpvals <0.05 & COpvals<0.05)
idx5=which(PDpvals <0.05 & COpvals>0.05)
idx6=which(COpvals<0.05 & PDpvals>0.05)


# for Figure other Version to include all taxa

COpvals=COstats=c()
PDpvals=PDstats=c()
for(i in 1:dim(M1)[1]) {
  sel=Status=="CO"
  wtCO=wilcox.test(M1[i, sel]~ Time[sel], paired=TRUE)
  COpvals[i]= wtCO$p.value
  COstats[i]=  wtCO$statistic
  sel=Status=="PD"
  wtPD=wilcox.test(M1[i, sel]~ Time[sel], paired=TRUE)
  PDpvals[i]= wtPD$p.value
  PDstats[i]=  wtPD$statistic
}

#load taxonomy
Tax <- read.delim("~/R/MGSv2/specI.tax",sep = "\t", header = TRUE, row.names = 1)
TaxNames1=Tax[rownames(M1),]

names(PDpvals)=TaxNames1$Species
names(COpvals)=TaxNames1$Species

sel1=c("Bifidobacterium adolescentis","Bifidobacterium longum","Bifidobacterium pseudocatenulatum","UMGS1975 sp900546685","Bifidobacterium breve")
sel2=c("Streptococcus thermophilus","Bifidobacterium angulatum","Bifidobacterium ruminantium")
sel3=c("Bifidobacterium animalis", "Bifidobacterium bifidum", "Bifidobacterium catenulatum", "Bifidobacterium dentium", "Bifidobacterium gallinarum")

col3=rep(rgb(0,0,0,1), length(PDpvals))
names(col3)=names(PDpvals)
col3[sel1]=rgb(0.2,0.4,0.6, 0.5)
col3[sel2]=rgb(0.4,0.8,1,0.5)
col3[sel3]=rgb(0.2,0.6,1,0.5)
pch=rep(20, length(PDpvals))
names(pch)=names(PDpvals)
pch[c(sel1, sel2)]=19



#visualize
windows()

plot(PDpvals, COpvals, col=  col3, pch= pch, cex=2.0,  log="xy", xlim=c(0.002, 1), ylim=c(0.002, 1), main="Several Bifidobacteria spp. increase after prebiotics \n (Species level)", 
     bty="L", xlab="PD before vs. after (p-value)", ylab = "CO before vs. after (p-value)")
abline(h=0.05, v=0.05, lt=2, col="darkgrey")
text(PDpvals[c(sel1,sel2,sel3)],COpvals[c(sel1,sel2,sel3)], c(sel1,sel2,sel3), cex=0.7, pos=4, col="black")

legend('bottomleft',c("sign. diff. in PD and CO", "sign. diff. uniquely in PD", "non sign. Bifidobacteria"),  col= c(rgb(0.2,0.4,0.6, 0.5),rgb(0.4,0.8,1,0.5),rgb(0.2,0.6,1,0.5)), cex=0.7,
       pch=list(19,19,20) , pt.cex = 1.5, bg='white',box.lty=0)

#Save set of different MGS for visualization
MGS_pre_post=as.data.frame(t(subset1))


###DATA VISULAZATION

library(ggplot2)
library(ggpubr)

#pre/post differences MGS
#include metadata
Meta <- read.delim("~/R/MGSv2/metaPDD",sep = "\t", header = TRUE, dec = ",", row.names = 1)
Meta <- Meta[ order(row.names(Meta)), ]

#exclude T2 samples and dropout samples at T1
Meta1 = Meta[!row.names(Meta) %in% c( "PD10T1", "C10T1", "PD4T2","PD6T2","PD10T0", "C10T0"), ]
Meta1 <- Meta1[ order(row.names(Meta1)), ]
Pair=as.factor(Meta1$Pair)


MGS_pre_post=cbind(MGS_pre_post, Status, Time, Pair)


MGS_pre_post$status_time = paste(MGS_pre_post$Status, MGS_pre_post$Time, sep = "_")

p<-ggpaired(MGS_pre_post, x = "Time", y = "MB2bin135",fill="status_time", facet.by = "Status", title= "B. adolescentis",id="Pair",line.color = "gray", line.size = 0.4) +
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  theme_classic() + ylab("Relative abundance") + xlab("Timepoint") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


p1<-ggpaired(MGS_pre_post, x = "Time", y = "MB2bin141", fill="status_time", facet.by = "Status", title= "B. longum", id="Pair",line.color = "gray", line.size = 0.4) +
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  theme_classic() + ylab("Relative abundance") + xlab("Timepoint") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


p2<-ggpaired(MGS_pre_post, x = "Time", y = "MB2bin145",fill="status_time", facet.by = "Status",title= "Strept. thermophilus",id="Pair",line.color = "gray", line.size = 0.4) +
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  theme_classic() + ylab("Relative abundance") + xlab("Timepoint") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


p3<-ggpaired(MGS_pre_post, x = "Time", y = "MB2bin265",fill="status_time", facet.by = "Status", title= "B. pseudocatenulatum",id="Pair",line.color = "gray", line.size = 0.4) +
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  theme_classic() + ylab("Relative abundance") + xlab("Timepoint") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


p4<-ggpaired(MGS_pre_post, x = "Time", y = "MB2bin272", fill="status_time", facet.by = "Status", title= "UMGS1975 (Christensenellales)",id="Pair",line.color = "gray", line.size = 0.4) +
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  theme_classic() + ylab("Relative abundance") + xlab("Timepoint") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


p5<-ggpaired(MGS_pre_post, x = "Time", y = "MB2bin288",fill="status_time", facet.by = "Status", title= "B. angulatum", id="Pair",line.color = "gray", line.size = 0.4) +
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  theme_classic() + ylab("Relative abundance") + xlab("Timepoint") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))



p6<-ggpaired(MGS_pre_post, x = "Time", y = "specI_v3_Cluster1098",fill="status_time", facet.by = "Status", title= "B. breve", id="Pair",line.color = "gray", line.size = 0.4) +
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  theme_classic() + ylab("Relative abundance") + xlab("Timepoint") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


p7<-ggpaired(MGS_pre_post, x = "Time", y = "specI_v3_Cluster2702", fill="status_time", facet.by = "Status", title= "B. ruminantium", id="Pair",line.color = "gray", line.size = 0.4) +
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  theme_classic() + ylab("Relative abundance") + xlab("Timepoint") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


p.all <- ggarrange(p,p1, p3, p5,p6,p7, p2, p4, ncol=4, nrow = 2, common.legend= TRUE)
p.all

#Group differences
#core differences pre and post
GENUS_group=as.data.frame(t(GENUS_group))
GENUS_group=cbind(GENUS_group, Status, Time)
GENUS_group$status_time = paste(GENUS_group$Status, GENUS_group$Time, sep = "_")

S1 <- ggplot(GENUS_group, aes(Status,`Bacteria;Firmicutes_A;Clostridia;Lachnospirales;Lachnospiraceae;Blautia_A`, fill = status_time))+
  geom_boxplot(outlier.colour = NA) + 
  geom_point(position=position_jitter(0.1), size = 2, alpha = 0.5) +
  labs(title = "Blautia A") +
  facet_wrap(Time) +   
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + ylab("Relative abundance") + xlab("Group")


S2 <- ggplot(GENUS_group, aes(Status,`Bacteria;Firmicutes_A;Clostridia;Lachnospirales;Lachnospirales;Dorea`, fill = status_time))+
  geom_boxplot(outlier.colour = NA) + 
  geom_point(position=position_jitter(0.1), size = 2, alpha = 0.5) +
  labs(title = "Dorea")+
  facet_wrap(Time) +
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + ylab("Relative abundance") + xlab("Group")



S3 <- ggplot(GENUS_group, aes(Status,`Bacteria;Firmicutes_A;Clostridia;Peptostreptococcales;Anaerovoracaceae;UBA1191`, fill = status_time))+ 
  geom_boxplot(outlier.colour = NA) + 
  geom_point(position=position_jitter(0.1), size = 2,alpha = 0.5) +
  labs(title = "UBA1191 (Anaerovoracaceae)")+
  facet_wrap(Time) +
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))  + ylab("Relative abundance") + xlab("Group")



S4 <- ggplot(GENUS_group, aes(Status,`Bacteria;Firmicutes_A;Clostridia_A;Christensenellales;QAND01;UMGS1975`, fill =status_time))+
  geom_boxplot(outlier.colour = NA) + 
  geom_point(position=position_jitter(0.1), size = 2,alpha = 0.5) +
  labs(title = "UMGS1975 (Christensenellales)")+
  facet_wrap(Time) +
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))  + ylab("Relative abundance") + xlab("Group")


S5 <- ggplot(GENUS_group, aes(Status,`Bacteria;Firmicutes_I;Bacilli_A;Erysipelotrichales;Erysipelatoclostridiaceae;Erysipelatoclostridium`, fill =status_time))+ 
  geom_boxplot(outlier.colour = NA) + 
  geom_point(position=position_jitter(0.1), size = 2,alpha = 0.5) +
  labs(title = "Erysipelatoclostridium")+
  facet_wrap(Time) +
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"))  + ylab("Relative abundance") + xlab("Group")

S6 <- ggplot(GENUS_group, aes(Status,`Bacteria;Firmicutes_A;Clostridia;Lachnospirales;Lachnospiraceae;Eubacterium_I`, fill =status_time))+ 
  geom_boxplot(outlier.colour = NA) + 
  geom_point(position=position_jitter(0.1), size = 2,alpha = 0.5) +
  labs(title = "Eubacterium I")+
  facet_wrap(Time) +
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"))  + ylab("Relative abundance") + xlab("Group")

#combine  in one plot  
windows()
ggarrange(S1, S2, S3, S4, S5, S6,  ncol = 3, nrow = 2,common.legend=TRUE)
