library(tidyverse)
library(ggpubr)
library(glue)
library(ggrepel)
library(vegan)

es_composition2 <- read_delim("apply_es_results/h.tsv")
es_composition2=es_composition2[-c(1,22,35,38), ]
rownames(es_composition2) <- es_composition2$sample

es_composition2$sample==sample_md2$sample


Status=sample_md2$condition
Time=sample_md2$time_point
es_composition2=t(es_composition2)
es_composition2=es_composition2[-1,]

pvals2=c()
for(i in 1:dim(es_composition2[,])[1]) {
  temp_dat = data.frame(es=as.numeric(es_composition2[i,]), `Status`=Status, `Time`=Time);
  
  wt = wilcox.test(es ~ Time, data=temp_dat, paired=TRUE);
  
  pval2 <- wt$p.value
  pvals2[i]=pval2
  
}


names(pvals2)=rownames(es_composition2)
qvals2 <- p.adjust(pvals2,"BH")

sum(qvals2<0.1, na.rm=TRUE)
idx3=which(qvals2<0.1 & pvals2<0.05) 


#select subset with q<0.1
subset1=es_composition2[idx3, ]


#post hoc test
COpvals=COstats=c()
PDpvals=PDstats=c()
for(i in 1:dim(subset1)[1]) {
  sel=Status=="Control"
  wtCO=wilcox.test(as.numeric(subset1[i, sel])~ Time[sel], paired=TRUE)
  COpvals[i]= wtCO$p.value
  COstats[i]=  wtCO$statistic
  sel=Status=="Parkinsons Disease"
  wtPD=wilcox.test(as.numeric(subset1[i, sel])~ Time[sel], paired=TRUE)
  PDpvals[i]= wtPD$p.value
  PDstats[i]=  wtPD$statistic
}

names(COpvals)=rownames(subset1)

names(PDpvals)=rownames(subset1)

COpvals
PDpvals

# Use this to subset the model fit
#remove unpaired PD/C10 sample 
es_composition=es_composition[-c(1:10), ]


es_composition %>%
  arrange(participant_pair, condition, time_point) %>%
  ggplot(aes(x = time_point, y = abundance, fill = status_time)) +
  scale_fill_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  geom_boxplot() +
  geom_line(aes(group = participant_id), linewidth = 0.2) +
  geom_point() +
  facet_grid(rows = vars(condition), cols = vars(es)) +
  #stat_compare_means(
  #  method = "wilcox.test", 
  # paired = TRUE,
  # Comparisons plots some nice brackets rather than a label
  #comparisons = list(c("0", "1"))
  # ) +
  # Adjust the scale to make a bit of room to plot braces
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  xlab("Timepoint") +
  ylab("ES Abundance") +
  theme_classic() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  ggtitle("Enterosignature shifts after prebiotics") +
  theme(plot.title = element_text(hjust = 0.5))


