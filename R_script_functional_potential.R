# load an prepare input data 
M <- read.delim("~/R/MGSv2/funcs/modules/modules/KEGG.mat",sep = "\t", header = TRUE, row.names = 1)
H<- read.delim("~/R/MGSv2/funcs/modules/modules/KEGG.descr", sep = "\t", header = TRUE, row.names = 1)

#or load

M <- read.delim("~/R/MGSv2/funcs/KEGG_brain.txt", sep = "\t", header = TRUE, row.names = 1)
H<- read.delim("~/R/MGSv2/funcs/KEGG_brain_hiera.txt", sep = "\t", header = TRUE, row.names = 1)

#remove T2 samples and PD/CO10 (dropouts)
M=M[, !colnames(M) %in% c("PD4T2", "PD6T2","PD10T1", "C10T1", "PD10T0", "C10T0")]
M=as.matrix(M) 
#normalize
M1 = sweep(M,2,colSums(M),"/") 

# for univariate tests see R_script_univariate_tests

#plot trends for brain relevant KEGG

# for Figure other Version

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
Tax=H
TaxNames1=Tax[rownames(M1),]

names(PDpvals)=TaxNames1$X.2
names(COpvals)=TaxNames1$X.2

sel1=c("GABA synthesis III", "ClpB (ATP-dependent chaperone protein)")
sel2=c("Acetate synthesis II", "Acetate synthesis III", "p-Cresol synthesis","Inositol synthesis", "Inositol degradation", "Glutamate synthesis I")
sel3=c("Acetate synthesis I")

col3=rep(rgb(0,0,0,1), length(PDpvals))
names(col3)=names(PDpvals)
col3[sel1]=rgb(0.2,0.4,0.6, 0.5)
col3[sel2]=rgb(0.4,0.8,1,0.5)
col3[sel3]=rgb(0.2,0.6,1,0.5)
pch=rep(20, length(PDpvals))
names(pch)=names(PDpvals)
pch[c(sel1, sel2)]=19



#plot 
windows()

plot(PDpvals, COpvals, col=  col3, pch= pch, cex=2.0,  log="xy", xlim=c(0.002, 1), ylim=c(0.002, 1), main="Brain???relevant bacterial metabolic functional potential \n improves after prebiotics", 
     bty="L", xlab="PD before vs. after (p-value)", ylab = "CO before vs. after (p-value)")
abline(h=0.05, v=0.05, lt=2, col="darkgrey")
text(PDpvals[c(sel1,sel2,sel3)],COpvals[c(sel1,sel2,sel3)], c(sel1,sel2,sel3), cex=0.7, pos=4, col="black")

legend('bottomleft',c("sign. diff. in PD and CO", "sign. diff. uniquely in PD", "sign. diff. uniquely in CO"),  col= c(rgb(0.2,0.4,0.6, 0.5),rgb(0.4,0.8,1,0.5),rgb(0.2,0.6,1,0.5)), cex=0.7,
       pch=list(19,19,20) , pt.cex = 1.5, bg='white',box.lty=0)

#plot modules of interest
M1=as.data.frame(t(M1))
MI=M1$"module_of_interest"
`

windows()
boxplot(MI~factor(paste(Status, Time)),col=c("palegreen3","seagreen","lightblue","steelblue"), xlab = "Timepoint", ylab = "module_of_interest", main="Urine NMR pre/post")
mtext(side = 3, line = 0.25, at = 1, adj = -1, subtitle)
stripchart(MI~factor(paste(Status, Time)), method = "jitter", vertical=TRUE, pch = 19, cex=1.5, add = TRUE, col = c("mediumseagreen", "darkgreen","lightsteelblue3","deepskyblue4"))


#fort plotting of paired data see _script_univariate_tests
# for dbRDA with envfitted top 15 modules  see R_script_dbRDA

