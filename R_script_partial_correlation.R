library(ppcor)
library(corrplot)


#load data

#MGS
M <- read.delim("~/R/MGSv2/specIMGSv2.mat",sep = "\t", header = TRUE, row.names = 1)
#remove T2 samples and PD/CO10 (dropouts)
#remove"?" unknowns
M=M[!rownames(M) %in% c("?"),]
M=as.matrix(M)
#normalize
M1 = sweep(M,2,colSums(M),"/") 
# Filter taxa present in less than 4 samples and taxa with relative abundance less than 0.001%
M1= M1[rowMeans(M1)>0.001 & rowSums(M1>0)>4, ]

#SCFA
SCFA<- read.delim("~/R/MGSv2/metabolites/SCFA1.txt",row.names = 1, sep = "\t", header = TRUE, dec = ".")
SCFA1=(SCFA[, !colnames(SCFA) %in% c("Individual", "Pair", "Time", "Group", "Visite")])
SCFA1 <-as.data.frame(t(SCFA1))

#KEGG
K <- read.delim("~/R/MGSv2/funcs/modules/modules/KEGG.mat",sep = "\t", header = TRUE, row.names = 1)

#brain relevant KEGG
K <- read.delim("~/R/MGSv2/funcs/KEGG_brain.txt", sep = "\t", header = TRUE, row.names = 1)
H<- read.delim("~/R/MGSv2/funcs/KEGG_brain_hiera.txt", sep = "\t", header = TRUE, row.names = 1)

#remove T2 samples and PD/CO10 (dropouts)
K=K[, !colnames(K) %in% c("PD4T2", "PD6T2","PD10T1", "C10T1", "PD10T0", "C10T0")]
K=as.matrix(K)
#normalize
K1 = sweep(K,2,colSums(K),"/") 


#correlate SCFA to clinical scores in PD subset

#select clinical scores
Clinmeta=metadata2[,c(4,5,7,8)]
Clinmeta$bloating=as.numeric(as.factor(Clinmeta$bloating))

Indv= as.numeric(as.factor(Meta2$Individual))

Clinmeta=t(Clinmeta)

#sel only PD

sel=Status=="PD"

output.adj <- data.frame(row.names = rownames(Clinmeta))
output2.adj <- data.frame(row.names = rownames(Clinmeta))


for (j in 1:dim(Clinmeta)[1]) {
  for (i in 1:dim(SCFA)[1]) {
    a <- pcor.test(unlist(Clinmeta[j,sel]), unlist(SCFA[i,sel]), z=Indv, method = "spearman") 
    output.adj[j,i] <- a$p.value
    output2.adj[j,i] <- a$estimate
    
    colnames(output.adj)[i] <- rownames(SCFA)[i]
    colnames(output2.adj)[i] <- rownames(SCFA)[i]
    
  }
}



output.p.adj= output.adj
for (i in 1:dim(SCFA)[1]) {
  output.p.adj[, i]<- p.adjust(output.adj[,i], "BH")
}


colnames(output2.adj)=colnames(output2.adj)
colnames(output.p.adj)=colnames(output.p.adj)


#visualize all with corrmatrix
windows()
corrplot(as.matrix(output2.adj), method= "color",addCoef.col =TRUE, p.mat = as.matrix(output.p.adj), tl.cex = 0.9, number.cex=0.7,
         addgrid.col = "white", tl.col ="black",tl.srt=70,  insig = "label_sig", sig.level = c(0.1,0.01, 0.001),
         pch.cex = 1.5, pch.col = "darkgray", is.corr = FALSE,
         title = "Correlation of SCFA to clinical scores in PD subset", cex.main=0.8, mar=c(0,0,2,0))


# corr coeff and sign. labels overlap, change

trace(corrplot, edit=TRUE)
#change position of sig.locs with +0.25 (place_points = function(sig.locs, point) { text(pos.pNew[, 1][sig.locs], pos.pNew[, 2][sig.locs], 
     #labels = point, col = pch.col, cex = pch.cex, 
     #lwd = 2)
#}) 

#redo plot, improve fonts
windows()
corrplot(as.matrix(output2.adj), method= "color",addCoef.col ="lightskyblue4", p.mat = as.matrix(output.p.adj), tl.cex = 1, number.cex=1, cl.pos = "b",
         addgrid.col = "white", tl.col ="black",tl.srt=70,  insig = "label_sig", sig.level = 0.1,
         pch.cex = 1.5, pch.col = "lightskyblue4", is.corr = FALSE,
         title = "Correlation of SCFA to clinical scores in PD subset", cex.main=1.2, mar=c(0,0,2,0))


# correlate taxa with functional metabolic potential

#define K1 and M1 matrices based on univariate testing


#correlation corrected for the individual
Indv= as.numeric(as.factor(Meta2$Individual))


output.adj <- data.frame(row.names = rownames(M1))
output2.adj <- data.frame(row.names = rownames(M1))

for (j in 1:dim(M1)[1]) {
  for (i in 1:dim(K1)[1]) {
    a <- pcor.test(unlist(M1[j, ]), unlist(K1[i, ]), z=Indv, method = "spearman") 
    output.adj[j,i] <- a$p.value
    output2.adj[j,i] <- a$estimate
    
    colnames(output.adj)[i] <- rownames(K1)[i]
    colnames(output2.adj)[i] <- rownames(K1)[i]
    
  }
}

#p.adjust
output.p.adj= output.adj
for (i in 1:dim(K1)[1]) {
  output.p.adj[, i]<- p.adjust(output.adj[,i], "BH")
}


colnames(output2.adj)=colnames(output2.adj)
colnames(output.p.adj)=colnames(output.p.adj)


# make long format
# all sign. raw p-values
output.adj= as.matrix(output.adj)
subset_KEGG = data.frame(Row=rownames(output.adj)[row(output.adj)], Col=colnames(output.adj)[col(output.adj)], 
                         Value=as.vector(output.adj))

a<-subset_KEGG[ subset_KEGG$Value < 0.05,]

# all sign. q-values   

output.p.adj= as.matrix(output.p.adj)
subset_KEGG1 = data.frame(Row=rownames(output.p.adj)[row(output.p.adj)],Col=colnames(output.p.adj)[col(output.p.adj)], 
                          Value=as.vector(output.p.adj))

b<-subset_KEGG1[ subset_KEGG1$Value < 0.1,]

# all corr-coeff (r)
output2.adj= as.matrix(output2.adj)
subset_KEGG2 = data.frame(Row=rownames(output2.adj)[row(output2.adj)],Col=colnames(output2.adj)[col(output2.adj)], 
                          Value=as.vector(output2.adj))


# generate final dataframe with q-vals and r-values

Corr_KEGG=cbind(subset_KEGG, subset_KEGG1$Value,subset_KEGG2$Value) 

colnames(Corr_KEGG)[3]<-"p-Value"
colnames(Corr_KEGG)[4]<-"q-Value"
colnames(Corr_KEGG)[5]<-"corr-coeff"

# proof qvals again/ select only sign. correlations
c<-Corr_KEGG[Corr_KEGG$`p-Value`<0.05 & Corr_KEGG$`q-Value` < 0.1,]
