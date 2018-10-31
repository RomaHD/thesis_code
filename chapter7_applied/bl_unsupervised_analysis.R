# written in 2015-2018, modified on 30.10.18 
# Chapter 7 of the thesis. 
# Burkitt lymphoma drug sensitivity screen analysis:
# Unsupervised analysis:
# single screen: calculation of IC50 and AUC from raw viability data, plotting cor-cor plot (AUC) and heatmap with raw data
# comb. screen: calculating of IC50, AUC, Combination Index (CI) values, plotting heatmap with CI values

library(drc)
library(flux)
library(gplots)
library(corrplot)
library(RColorBrewer)

path <- "/data/kurilov/genestack/phd/work_2018/thesis_code/"
setwd(file.path(path, "chapter7_patient_applied"))

### 1. calculating AUC and IC50 avlues from single screen raw data ###
data <- read.csv(file.path(path, "data/Single_drug.csv"), check.names = F)

rownames(data) <- data$Compound
data <- data[,-1]  

table_auc <- matrix(NA, 34, 32)
table_ic50 <- matrix(NA, 34, 32)
rownames(table_auc) <- rownames(data)[-1]
rownames(table_ic50) <- rownames(data)[-1]
seq <- seq(1, 320, 10)
colnames(table_auc) <- colnames(data)[seq]
colnames(table_ic50) <- colnames(data)[seq]

# calculating AUC values from raw data
for (i in 1:32)
{
  conc <- unlist(data[1,((i-1)*10+1):(i*10)])
  
  table_auc[,i] <- apply(data[-1,((i-1)*10+1):(i*10)], 1, function(x) {
    a=auc(conc,x)/100 
    print(a)
    a
  })
  
}  

# calculating IC50 values from raw data
# (some cell lines-drugs combinations produces an error)
for (i in 1:32)
{
  conc <- unlist(data[1,((i-1)*10+1):(i*10)])
  table_ic50[,i] <- apply(data[-c(1),((i-1)*10+1):(i*10)], 1, function(x) {
    t="NA"
    try({
      fit <- drm(x ~ conc, fct = LL.4())
      t <- ED(fit,50, lref=0,uref=100, type="absolute", display=FALSE)[1] 
    })
    if (is.na(t)) {
      if (mean(x)>50) t <- 30
      if ((x[1]<50) & (mean(x)<50)) t <- 0.0015
    }
    as.numeric(t)
  })
}

write.table(table_auc, file="table_auc.txt")
write.table(table_ic50, file="table_ic50.txt")

heatmap(table_auc, scale="none",main="AUC, single drug screen")
heatmap(table_ic50, main="IC50, single drug screen")

### 2.drug-drug correlation plot and heatmap with raw values ###

# 1) drug-drug cor plot
#let's now set maximal ic50 or auc value to 30uM
auc_table <- as.matrix(table_auc)
n <- which(auc_table>30, arr.ind =T)
auc_table[n] <- 30

#let's change some drug names
drug_names <- read.csv(file.path(path, "data/drug-names.csv"))
drug_names <- drug_names[,1]
colnames(auc_table) <- drug_names[c(1:8,33,9,34,10,35,11,12,36,13,14,37,38,15,39:41,16:20,42,21,22)]

pdf(file="fig_30.pdf")
M <- cor(auc_table)
corrplot(M, method="color", tl.cex=0.8, order = "hclust")
dev.off()

# 2) heatmap with raw values

#median subsrtuction
data_ms <- data
data_ms[-1,] <- apply(data[-1,],2, function(x) {
  m <- median(x, na.rm=TRUE)
  x-m
}) 

col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(n = 206) 

col_breaks2 = c( seq(-111,-100.1,length=1),
                 seq(-100,-1.1,length=100), 
                 seq(-1,1,length=5), # for white
                 seq(1.1,100,length=100), 
                 seq(100.1,171,length=1)) 


par(cex.main=1)
par(mar=c(5,5,5,5))

pdf(file="fig_31.pdf")
heatmap.2(as.matrix(data_ms[-c(1:6,29, 31),]), col=col2,breaks=col_breaks2, key=T, keysize=1.5,margin=c(4, 8),
          density.info="none", trace="none",cexCol=0.07,cexRow=0.8, main="raw values (mean subtracted), single drug screen")
dev.off()

###. 3. Combination screen: AUC, IC50, CI ###
data <- read.csv(file.path(path, "data/Comb_drugs2.csv"), check.names = F)
rownames(data) <- paste0(data$Compound, "_", data$`Combination compound`)

# AUC calculation
table2_auc <- matrix(NA, 78, 32)
rownames(table2_auc) <- rownames(data)[-1]
table2_auc <- cbind(data[-1,1:2],table2_auc)
seq <- seq(1, 320, 10)
seq<- seq+3
colnames(table2_auc) <- colnames(data)[c(1,2,seq)]


for (i in 1:32)
{
  conc <- unlist(data[1,(3+(i-1)*10+1):(3+i*10)])
  print(conc)
  table2_auc[,2+i] <- apply(data[-1,(3+(i-1)*10+1):(3+i*10)], 1, function(x) {
    auc(conc,x)/100
  })
  
}  

write.table(table2_auc, file="comb_auc_new.txt")


# IC50 calculation
table2_ic50 <- matrix(NA, 78, 32)
rownames(table2_ic50) <- rownames(data)[-1]
table2_ic50 <- cbind(data[-1,1:2],table2_ic50)
seq <- seq(1, 320, 10)
seq<- seq+3
colnames(table2_ic50) <- colnames(data)[c(1,2,seq)]

for (i in 1:32)
{
  conc <- unlist(data[1,(3+(i-1)*10+1):(3+i*10)])
  table2_ic50[,2+i] <- apply(data[-1,(3+(i-1)*10+1):(3+i*10)], 1, function(x) {
    t="NA"
    try({fit <- drm(x ~ conc, fct = LL.4())
    t <- ED(fit,50, lref=0,uref=100, type="absolute", display=FALSE)[1] 
    })
    if (is.na(t)) {
      if (mean(x)>50) t <- 30
      if ((x[1]<50) & (mean(x)<50)) t <- 0.0015
    }
    as.numeric(t)
  })
}

table2_ic50[70,16] <- 5.2
nn <- which(is.na(table2_ic50), arr.ind =T)
table2_ic50[nn] <- 30


#let's now set maximal ic50 or auc value to 30uM
auc_table <- table2_auc[,3:34]
n <- which(auc_table>30, arr.ind =T)
auc_table[n] <- 30

ic50_table <- table2_ic50[,3:34]
n <- which(ic50_table>30, arr.ind =T)
ic50_table[n] <- 30

# saving IC50, AUC tables
write.table(auc_table, file="auc_table_comb.txt")
write.table(ic50_table, file="ic50_table_comb.txt")

# calculating Combination Index (CI)
CI <- matrix(NA, 54, 32)
dmso_num <- seq(1,69,4)
colnames(CI) <- colnames(ic50_table)
rownames(CI) <- rownames(ic50_table[1:72,])[-dmso_num]
head(CI)

for (i in 1:18)
{
  for (j in 1:32)
  {
    if (ic50_table[(4*i-3),1] >0.2) {
      CI[(3*i-2),j]=(0.2/ic50_table[(4*i-3),1])+(ic50_table[(4*i-2),j]/ic50_table[(4*i-3),j])
    } else CI[(3*i-2),j] <- NA
    
    if (ic50_table[(4*i-3),2] >0.2) {
      CI[(3*i-1),j]=(0.2/ic50_table[(4*i-3),2])+(ic50_table[(4*i-1),j]/ic50_table[(4*i-3),j])
    } else CI[(3*i-1),j] <- NA
    
    if (ic50_table[(4*i-3),12] >0.2) {
      CI[(3*i-0),j]=(0.2/ic50_table[(4*i-3),12])+(ic50_table[(4*i-0),j]/ic50_table[(4*i-3),j])
    } else CI[(3*i-0),j] <- NA
  }
}

write.table(CI, file="comb_index.txt")

###. 4. Combination screen: plotting CI values heatmap ###
ci <- CI

ci_mod <- cbind(ci[seq(1,54, by=3),], ci[seq(2,54, by=3),], ci[seq(3,54, by=3),])
colnames(ci_mod) <- c(paste0("PCI-32765 + ", colnames(ci_mod)[1:32]), paste0("CAL-101 + ", colnames(ci_mod)[1:32]), paste0("OTX15 + ", colnames(ci_mod)[1:32]))
rownames(ci_mod) <- unlist(strsplit(rownames(ci_mod), split="_"))[seq(1,36, by=2)]
colnames(ci_mod)[c(29,61,93)] <- c("PCI-32765 + PRT062607", "CAL-101 + PRT062607","OTX15 + PRT062607" )
n <- which(ci_mod>2, arr.ind=TRUE)
ci_mod[n] <- 2

# making ALL values for BL-7,BL-41, Gumbus equal to NA (asked by Kasia)
ci_mod[c(1,10,17),] <- NA

col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(n = 206) 

col_breaks2 = c( seq(0,0.01,length=1),
                 seq(0.02,0.84,length=100), 
                 seq(0.85,1.15,length=5), # for white
                 seq(1.16,1.98,length=100), 
                 seq(1.99,2,length=1)) 

par(cex.main=1)
par(oma=c(7,0,0,1))
pdf(file="fig_33.pdf")
heatmap.2(as.matrix(ci_mod[-c(1,10,17),]), col=col2,breaks=col_breaks2, key=T, keysize=1.5, na.color = "grey",
          density.info="none", trace="none",cexCol=0.7,cexRow=1, main="Combination Index (CI)")
dev.off()