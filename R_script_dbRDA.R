#dbrDA and household effect on taxonomy

# packages 
library(vegan)
library(dplyr)
library(tidyr)
library(reshape)
library(ape)
library(tibble)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(cowplot)

#load input data
M <- read.delim("~/R/MGSv2/specIMGSv2.mat", sep = "\t", header = TRUE, row.names = 1)

#remove T2 samples and PD/CO10 (dropouts)
M=M[, !colnames(M) %in% c("PD4T2", "PD6T2","PD10T1", "C10T1", "PD10T0", "C10T0")]

#remove"?"/unknowns
M=M[!rownames(M) %in% c("?"),]
M=as.matrix(M)

#normalize (TSS)
M1 = sweep(M,2,colSums(M),"/") 

# load tax-table
taxa1<- read.delim("~/R/MGSv2/specI.tax", sep = "\t", header = TRUE)

#remove"?"/unknowns
taxa1=taxa1[-1,]

rownames(taxa1)<- taxa1$MGS
#set taxa as factor
taxa1<- taxa1 %>% 
  mutate_if(is.character, as.factor)
str(taxa1)

taxa1<- taxa1 %>% 
  select(-MGS)


# Metadata 
#exclude T2 samples and dropout samples at T1 and format environmental variables
metadata2 = read.delim("~/R/MGSv2/metaPDD_mean.txt", sep = "\t", header = TRUE, row.names = 1)
metadata2 = metadata2[!row.names(metadata2) %in% c("PD4T2", "PD6T2","PD10T1", "C10T1", "PD10T0", "C10T0"),]
metadata2 = data.frame(row.names = row.names(metadata2),
                       Individual=as.factor(metadata2$Individual),
                       time = as.factor(metadata2$Time),
                       status = as.factor(metadata2$Group),
                       bloating = as.numeric(metadata2$Bloating),
                       UPDRS = as.numeric(metadata2$UPDRSIII),
                       Age = as.numeric(metadata2$Age),
                       GSRS = as.numeric(metadata2$GSRS),
                       StoolFreq = as.numeric(metadata2$StoolFrequencyScore),
                       name = row.names(metadata2))

#create new variable "family" based on spouses/pairs of PD/CO
metadata2 = metadata2 %>% separate(col = name, into = c("time2", "family", "time_point"), sep = "PD|C|T")
metadata2 = metadata2[,c(1:8,10)]

#perform partialing out of family/household using an unconstrained dbRDA
dbRDA2=dbrda(t(M1)~  Condition(family), metadata2, distance = "bray", sqrt.dist =FALSE ,add = FALSE, dfun = vegdist, metaMDSdist = FALSE,na.action = na.fail)
sppscores(dbRDA2) <- t(M1)
dbRDA2 # family explains 46.4% of the data variability

windows()
screeplot(dbRDA2)

summary(dbRDA2)

#test dissimilarity for the different timepoints separately
#get data
dbRDA2_site[,6]=metadata2[,9]
colnames(dbRDA2_site)[colnames(dbRDA2_site) == "V6"] <- "family"
dbRDA2_site$family = factor(dbRDA2_site$family, levels = unique(dbRDA2_site$family))


#test all 
set.seed(1)
adonis2(vegdist(t(M1), method = "bray") ~ time + status + family, data = dbRDA2_site, permutations = 1000, parallel = 8, by = "margin")

#block family
set.seed(1)
perm <- how(nperm = 1000)
setBlocks(perm) <- with(dbRDA2_site, family)
adonis2(vegdist(t(M1), method = "bray") ~ time + status, data = dbRDA2_site, permutations = perm, parallel = 8)

#test for different timepoints per group
#before
sel=dbRDA2_site$time=="0"
set.seed(1)
perm <- how(nperm = 1000)
setBlocks(perm) <- with(dbRDA2_site[sel,], family)
adonis2(vegdist(t(M1[,sel]), method = "bray") ~ status, data = dbRDA2_site[sel,], permutations = perm, parallel = 8)

#after
sel=dbRDA2_site$time=="28"
set.seed(1)
perm <- how(nperm = 1000)
setBlocks(perm) <- with(dbRDA2_site[sel,], family)
adonis2(vegdist(t(M1[,sel]), method = "bray") ~status, data = dbRDA2_site[sel,], permutations = perm, parallel = 8)



#visualize

#extract scores
dbRDA2_scores = scores(dbRDA2, scaling = 3, correlation = TRUE)
dbRDA2_site = data.frame(dbRDA2_scores$sites)
name_match = match(row.names(dbRDA2_site), row.names(metadata2))
dbRDA2_site[,c(3:4)] =  metadata2[name_match,c(2,3)]
dbRDA2_site$status_time = as.factor(paste(dbRDA2_site$status, dbRDA2_site$time, sep = "_"))


#extract species
dbRDA2_species = data.frame(dbRDA2_scores$species)
dbRDA2_species$length = round(sqrt(dbRDA2_species[,1]^2 + dbRDA2_species[,2]^2), digits = 3)
name_match = match(row.names(dbRDA2_species), row.names(taxa1))
dbRDA2_species$genus = as.vector(taxa1[name_match,6])
dbRDA2_species$species = as.vector(taxa1[name_match,7])
dbRDA2_species = dbRDA2_species %>% arrange(desc(length))
dbRDA2_species$genus_species = dbRDA2_species$species
dbRDA2_species$genus_species[which(dbRDA2_species$genus_species=="?")] = dbRDA2_species$genus[which(dbRDA2_species$genus_species=="?")]
dbRDA2_species[,c(3:6)]

#plot
windows()
ggplot(dbRDA2_site, aes(x = MDS1, y = MDS2)) +
  geom_point(data=dbRDA2_site,aes(color=status_time, shape=status_time),size = 3.5, alpha = 0.5) +
  stat_ellipse(data=dbRDA2_site,aes(color=status_time)) +
  scale_color_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  scale_shape_manual(values = c(19,19,17,17),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  facet_wrap(~time) +
  ylab("MDS2") + xlab("MDS1") +
  theme_minimal_grid() +
  background_grid() +
  ggtitle("Dissimilarity bewteen PD and CO before and after prebiotics")+
  theme(plot.title = element_text(hjust = 0.5))
 
#plot with species (top 15)

#extract top 15 
dbRDA2_species_filter = dbRDA2_species %>% 
  arrange(desc(length)) %>%
  filter(between(row_number(), 1,15))

#enfit with environmental variables (Status and Time)
env.factors=dbRDA2_site[,3:4]

set.seed(1)
en = envfit(dbRDA2, env.factors, permutations = 999, na.rm = TRUE)

en

#extract centroids
en_centroid = rownames_to_column(data.frame(en$factors$centroids))
en_centroid$length = round(sqrt(en_centroid$MDS1^2 + en_centroid$MDS2^2), digits = 3)


data.scores = as.data.frame(scores(dbRDA2)$sites)

#plot top 15 species with env variables
windows()
ggplot(data=dbRDA2_site, aes(x = MDS1, y = MDS2)) +
  geom_point(data=dbRDA2_site,aes(color=status_time, shape=status_time),size = 3.5, alpha = 0.5)  +
  scale_color_manual(values = c("palegreen3","seagreen","lightblue","steelblue"),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  scale_shape_manual(values = c(19,19,17,17),name="Group", labels=c("CO before", "CO after", "PD before", "PD after")) +
  geom_segment(data = en_centroid, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), arrow = arrow(length = unit(0.2, "cm")), color = "grey10", lwd = 0.3) + 
  geom_label_repel(data = en_centroid, aes(label = rowname), box.padding = 0.3, size = 3, color = "grey10", max.overlaps = Inf,force = 1) +
  geom_point(data=dbRDA2_species_filter, aes(x = MDS1, y = MDS2), size = 2.5, alpha=0.5) +
  geom_text_repel(data=dbRDA2_species_filter, aes(label = genus_species), size =3) +
  theme_minimal_grid() +
  ggtitle("Bacterial composition conditioned by household \n constrained by disease and treatment") +
  theme(plot.title = element_text(hjust = 0.5))



#partial out Status and Time
dbRDA3=dbrda(t(M1)~ family + Condition(status+time), metadata2, distance = "bray", sqrt.dist =FALSE ,add = FALSE, dfun = vegdist, metaMDSdist = FALSE,na.action = na.fail)
sppscores(dbRDA3) <- t(M1)
dbRDA3 # status+time explain 9.8% of the data variability

windows()
screeplot(dbRDA3)

#test
set.seed(1)
perm <- how(nperm = 1000)
setBlocks(perm) <- with(dbRDA3_site, time, status)
adonis2(vegdist(t(M1), method = "bray") ~ family , data = dbRDA3_site, permutations = perm, parallel = 8)


#visualize
#extract scores
dbRDA3_scores = scores(dbRDA3, scaling = 3, correction = TRUE)
dbRDA3_site = data.frame(dbRDA3_scores$sites)
name_match = match(row.names(dbRDA3_site), row.names(metadata2))
dbRDA3_site[,c(3:4)] = metadata2[name_match, c(2,3)]
dbRDA3_site$status_time = paste(dbRDA3_site$status, dbRDA3_site$time, sep = "_")
dbRDA3_site$family = metadata2[name_match, 9]

#extract centroids
dbRDA3_centroid = rownames_to_column(data.frame(dbRDA3_scores$centroids))
dbRDA3_centroid$length = round(sqrt(dbRDA3_centroid$dbRDA1^2 + dbRDA3_centroid$dbRDA2^2), digit = 3)
dbRDA3_centroid %>% arrange(desc(length))


#extract species
dbRDA3_species = data.frame(dbRDA3_scores$species)
dbRDA3_species$length = round(sqrt(dbRDA3_species[,1]^2 + dbRDA3_species[,2]^2), digits = 3)
name_match = match(row.names(dbRDA3_species), row.names(taxa1))
dbRDA3_species$genus = as.vector(taxa1[name_match,6])
dbRDA3_species$species = as.vector(taxa1[name_match,7])
dbRDA3_species = dbRDA3_species %>% arrange(desc(length))
dbRDA3_species$genus_species = dbRDA3_species$species
dbRDA3_species$genus_species[which(dbRDA3_species$genus_species=="?")] = dbRDA3_species$genus[which(dbRDA3_species$genus_species=="?")]
dbRDA3_species[,c(3:6)]

#plot with spiders
library(ggordiplots)

plot<- gg_ordiplot(dbRDA3,groups = dbRDA3_site$family, label=TRUE, spiders = TRUE, ellipse = FALSE, pt.size = 4, plot= FALSE) 

windows()
plot$plot+
  theme_cowplot() +
  labs(color="Family", x="dbRDA1 (25.3%)", y="dbRDA2 (17.1%)", title="Bacterial species composition \n in different housedolds") +
  theme(plot.title = element_text(hjust = 0.5))

