# written on 1.08.206, modified on 29.10.18 
# Chapter 6 of the thesis. 
# Patient treatment outcome prediction using classification models trained on cell lines:
# Data preparation

source("http://bioconductor.org/biocLite.R")
biocLite("sva")
library(devtools)
install_github("vqv/ggbiplot")
library(PharmacoGx)
library(Biobase)
library(plyr)
library(ggplot2)
library(ggbiplot)
library(sva)

path <- "/data/kurilov/genestack/phd/work_2018/thesis_code/"
setwd(file.path(path, "chapter6_patient_response"))

source("prep_functions/binarization.R")
source("prep_functions/binarization_nibr.R")

setwd(file.path(path, "data"))

load("exp_nibr.RData")
load("bestresp_nibr.RData")
load("bestavgresp_nibr.RData")
load("tissue_nibr.RData")
bestresp_nibr <- as.data.frame(bestresp_nibr)
bestavgresp_nibr <- as.data.frame(bestavgresp_nibr)

# getting data from PharmacoGx
CCLE <- downloadPSet("CCLE")
CTRP <- downloadPSet("CTRPv2")
GDSC <- downloadPSet("GDSC1000")
gCSI <- downloadPSet("gCSI")

#exp
ccle_exp <- summarizeMolecularProfiles(pSet=CCLE, mDataType="rna", cell.lines=cellNames(CCLE), summary.stat = 'mean', fill.missing = TRUE, verbose=TRUE)
ccle_exp <- t(exprs(ccle_exp))

gdsc_exp <- summarizeMolecularProfiles(pSet=GDSC, mDataType="rna", cell.lines=cellNames(GDSC), summary.stat = 'mean', fill.missing = TRUE, verbose=TRUE)
gdsc_exp <- t(exprs(gdsc_exp))

gcsi_exp <- summarizeMolecularProfiles(gCSI, mDataType = "rnaseq", cell.lines=cellNames(gCSI), summary.stat = 'mean', fill.missing = TRUE, verbose=TRUE)
gcsi_exp <- t(exprs(gcsi_exp))

# ensemble -> gene symbol conversion for ccle
gene_order <- match(colnames(ccle_exp), featureInfo(CCLE, "rna")[,"Probe"])
colnames(ccle_exp) <- featureInfo(CCLE, "rna")[gene_order,"Symbol"]
na_names <- which(duplicated(colnames(ccle_exp)))
ccle_exp <- ccle_exp[,-na_names]

# ensemble -> gene symbol conversion for gdsc
gene_order <- match(colnames(gdsc_exp), featureInfo(GDSC, "rna")[,"Probe"])
colnames(gdsc_exp) <- featureInfo(GDSC, "rna")[gene_order,"Symbol"]
na_names <- which(duplicated(colnames(gdsc_exp)))
gdsc_exp <- gdsc_exp[,-na_names]

# gene Symbols for gCSI
colnames(gcsi_exp) <- featureInfo(gCSI, "rnaseq")$Symbol

#prefixes
colnames(ccle_exp) <- paste0("exp_", colnames(ccle_exp))
colnames(gdsc_exp) <- paste0("exp_", colnames(gdsc_exp))
colnames(gcsi_exp) <- paste0("exp_", colnames(gcsi_exp))

# renaming
exp_ctrp <- ccle_exp
exp_gdsc <- gdsc_exp
exp_gcsi <- gcsi_exp

# removing rows that have only NAs
na_ctrp <- which(is.na(apply(exp_ctrp,1,function(x){sd(x,na.rm=T)})))
na_gdsc <- which(is.na(apply(exp_gdsc,1,function(x){sd(x,na.rm=T)})))
na_gcsi <- which(is.na(apply(exp_gcsi,1,function(x){sd(x,na.rm=T)})))

exp_ctrp <- exp_ctrp[-na_ctrp,]
exp_gdsc <- exp_gdsc[-na_gdsc,]
exp_gcsi <- exp_gcsi[-na_gcsi,]


# drug_response
drugNames(GDSC) <- GDSC@drug$DRUG.NAME
drugNames(CTRP) <- CTRP@drug$cpd_name
#drugNames(gCSI) <- gCSI@drug$originalDrugId

# AUC
auc_ctrp <- t(summarizeSensitivityProfiles(CTRP, sensitivity.measure="auc_recomputed", drugs=c("erlotinib", "bortezomib","docetaxel")))
auc_gdsc <- t(summarizeSensitivityProfiles(GDSC, sensitivity.measure="auc_recomputed", drugs=c("Erlotinib", "Bortezomib","Docetaxel")))
auc_gcsi <- t(summarizeSensitivityProfiles(gCSI ,sensitivity.measure = "auc_recomputed", drugs=c("Erlotinib", "Bortezomib","Docetaxel")))

# 1-AUC correction for AUC data
auc_ctrp <- 1-auc_ctrp
auc_gdsc <- 1-auc_gdsc
auc_gcsi <- 1-auc_gcsi

# IC50
ic50_ctrp <- t(summarizeSensitivityProfiles(CTRP, sensitivity.measure="ic50_recomputed", drugs=c("erlotinib", "bortezomib","docetaxel")))
ic50_gdsc <- t(summarizeSensitivityProfiles(GDSC, sensitivity.measure="ic50_recomputed", drugs=c("Erlotinib", "Bortezomib","Docetaxel")))
ic50_gcsi <- t(summarizeSensitivityProfiles(gCSI ,sensitivity.measure = "ic50_recomputed", drugs=c("Erlotinib", "Bortezomib","Docetaxel")))

outlier_capping <- function(table){
  for (i in 1:ncol(table)){
    fin <- which(is.finite(table[,i]))
    qv <- quantile(table[fin,i], probs=0.85, na.rm=T)
    out <- which(table[,i]>qv)
    table[out,i] <- qv
  }
  return(table)
}

ic50_ctrp <- outlier_capping(ic50_ctrp)
ic50_gdsc <- outlier_capping(ic50_gdsc)
ic50_gcsi <- outlier_capping(ic50_gcsi)

# viability
load("viab_gdsc1000.RData")
viab_gdsc <- viab_gdsc1000
load("viab_gCSI.RData")
viab_gcsi <- viab_gCSI
load("viab_ctrp.RData")


# renaming rownames for viability data
temp <- read.csv("viab_names_diff.txt")
temp <- temp[,2:3]
ord_ctrp <- match(temp[,1], rownames(viab_ctrp))
rownames(viab_ctrp)[ord_ctrp] <- as.character(temp[,2])
viab_gdsc <- viab_gdsc[intersect(rownames(viab_gdsc), cellInfo(GDSC)$Sample.Name),]
rownames(viab_gdsc) <- cellInfo(GDSC)$cellid[match(rownames(viab_gdsc), cellInfo(GDSC)$Sample.Name)]
rownames(viab_gcsi) <- rownames(auc_gcsi)[-grep("SR", rownames(auc_gcsi))[2]]

#converting drug. resp data into data frames
auc_ctrp <- data.frame(auc_ctrp)
auc_gdsc <- data.frame(auc_gdsc)
auc_gcsi <- data.frame(auc_gcsi)

ic50_ctrp <- data.frame(ic50_ctrp)
ic50_gdsc <- data.frame(ic50_gdsc)
ic50_gcsi <- data.frame(ic50_gcsi)

viab_ctrp <- data.frame(viab_ctrp)
viab_gdsc <- data.frame(viab_gdsc)
viab_gcsi <- data.frame(viab_gcsi)


# preparinh tissue files: gdsc, ctrp, gcsi

tissue_ctrp <- cellInfo(CTRP)[,c("cellid","tissueid")]
tissue_gdsc <- cellInfo(GDSC)[,c("cellid","tissueid")]
tissue_gcsi <- cellInfo(gCSI)[,c("unique.id","tissueid")]

colnames(tissue_ctrp) <-c("sample","tissue")
colnames(tissue_gdsc) <-c("sample","tissue")
colnames(tissue_gcsi) <-c("sample","tissue")

# asessing overlap
length(intersect(rownames(exp_gdsc), rownames(viab_gdsc)))
length(intersect(rownames(exp_gcsi), rownames(viab_gcsi)))
length(intersect(rownames(exp_ctrp), rownames(viab_ctrp)))

# binarized (with and without grey_zone) auc, ic50, viab labels and tissue inf. for cell line data
data_prep <- function(dataset){
  exp <- get(paste0("exp_",dataset))
  tissue <- get(paste0("tissue_",dataset))
  ord_t <- match(rownames(exp), tissue$sample)
  tissue <- as.character(tissue$tissue[ord_t])
  batch <- rep(dataset, dim(exp)[1])
  
  if (dataset != "nibr") {
    auc <- get(paste0("auc_",dataset))
    ic50 <- get(paste0("ic50_",dataset))
    viab <- get(paste0("viab_",dataset))
    
    colnames(auc) <- tolower(colnames(auc))
    colnames(ic50) <- tolower(colnames(ic50))
    colnames(viab) <- tolower(colnames(viab))
    
    ord <- match(rownames(exp), rownames(auc))
    
    erlotinib <- binarization(auc$erlotinib[ord])
    bortezomib <- binarization(auc$bortezomib[ord])
    docetaxel <- binarization(auc$docetaxel[ord])
    erlotinib_g <- binarization(auc$erlotinib[ord], grey_zone = T)
    bortezomib_g <- binarization(auc$bortezomib[ord], grey_zone = T)
    docetaxel_g <- binarization(auc$docetaxel[ord], grey_zone = T)
    auc <- list(erlotinib=erlotinib, bortezomib=bortezomib, docetaxel=docetaxel,
                erlotinib_g=erlotinib_g, bortezomib_g=bortezomib_g, docetaxel_g=docetaxel_g)
    
    erlotinib <- binarization(ic50$erlotinib[ord])
    bortezomib <- binarization(ic50$bortezomib[ord])
    docetaxel <- binarization(ic50$docetaxel[ord])
    erlotinib_g <- binarization(ic50$erlotinib[ord], grey_zone = T)
    bortezomib_g <- binarization(ic50$bortezomib[ord], grey_zone = T)
    docetaxel_g <- binarization(ic50$docetaxel[ord], grey_zone = T)
    ic50 <- list(erlotinib=erlotinib, bortezomib=bortezomib, docetaxel=docetaxel,
                 erlotinib_g=erlotinib_g, bortezomib_g=bortezomib_g, docetaxel_g=docetaxel_g)
    
    erlotinib <- binarization(viab$erlotinib[ord])
    bortezomib <- binarization(viab$bortezomib[ord])
    docetaxel <- binarization(viab$docetaxel[ord])
    erlotinib_g <- binarization(viab$erlotinib[ord], grey_zone = T)
    bortezomib_g <- binarization(viab$bortezomib[ord], grey_zone = T)
    docetaxel_g <- binarization(viab$docetaxel[ord], grey_zone = T)
    viab <- list(erlotinib=erlotinib, bortezomib=bortezomib, docetaxel=docetaxel,
                 erlotinib_g=erlotinib_g, bortezomib_g=bortezomib_g, docetaxel_g=docetaxel_g)
    
    result <- list(exp=exp, auc=auc, ic50=ic50, viab=viab, tissue=tissue, batch=batch)
  } else {
    bestresp <- get(paste0("bestresp_",dataset))
    bestavgresp <- get(paste0("bestavgresp_",dataset))
    
    ord <- match(rownames(exp), rownames(bestresp))
    erlotinib <- binarization_nibr(bestresp$erlotinib[ord])
    erlotinib_g <- binarization_nibr(bestresp$erlotinib[ord], grey_zone = T)
    bestresp <- list(erlotinib=erlotinib, erlotinib_g=erlotinib_g)
    
    erlotinib <- binarization_nibr(bestavgresp$erlotinib[ord])
    erlotinib_g <- binarization_nibr(bestavgresp$erlotinib[ord], grey_zone = T)
    bestavgresp <- list(erlotinib=erlotinib, erlotinib_g=erlotinib_g)
    
    result <- list(exp=exp, bestresp=bestresp, bestavgresp=bestavgresp, tissue=tissue, batch=batch)
  }
  return(result)
}

ctrp_data <- data_prep("ctrp")
gdsc_data <- data_prep("gdsc")
gcsi_data <- data_prep("gcsi")
nibr_data <- data_prep("nibr")


# let's z-transform cell lines expression
ctrp_data$exp <- scale(ctrp_data$exp)
gdsc_data$exp <- scale(gdsc_data$exp)
gcsi_data$exp <- scale(gcsi_data$exp)
nibr_data$exp <- scale(nibr_data$exp)

# fix for columns with only NA values after this transform
nn <- which(is.na(apply(nibr_data$exp, 2, function(x) {mean(x, na.rm=T)})))
nibr_data$exp <- nibr_data$exp[,-nn]

#preparing patient data for 3 drugs (3 cohorts) from BHK paper data

# Bortezamib
load("bortezomib.patient.RData")
pt_bortezomib_exp <- t(bortezomib.patient_ComBat)
pt_bortezomib_resp <- revalue(as.character(binaryResponse), c("0"="resist", "1"="sens"))
# let's also make z-transformation (to have data in the same format as erlotinib and docetaxel)
pt_bortezomib_exp <- scale(pt_bortezomib_exp)

# Erlotinib
load("erlotinib_data.RData")
pt_erlotinib_exp <- erlotinib$patient
pt_erlotinib_resp <- revalue(as.character(erlotinib.labels$patient), c("0"="resist", "1"="sens"))

# Docetaxel
load("pp.RData")
pt_docetaxel_exp <- docetaxel.patient
pt_docetaxel_resp <- revalue(as.character(pp.ground_truth), c("0"="resist", "1"="sens"))

colnames(pt_bortezomib_exp) <- paste0("exp_", colnames(pt_bortezomib_exp))
colnames(pt_erlotinib_exp) <- paste0("exp_", colnames(pt_erlotinib_exp))
colnames(pt_docetaxel_exp) <- paste0("exp_", colnames(pt_docetaxel_exp))

# building new tables and batch vectors

combine <- function(cell_set_name, drug_name){
  cell_data <- get(paste0(cell_set_name,"_data"))
  cell_exp_set <- cell_data$exp
  cell_tissue <- cell_data$tissue
  pt_exp_set <- get(paste0("pt_",drug_name,"_exp"))
  
  genes <- intersect(colnames(cell_exp_set), colnames(pt_exp_set))
  exp <- rbind(cell_exp_set[,genes], pt_exp_set[,genes])
  tissue <- c(cell_tissue, rep("pt", dim(pt_exp_set)[1]))
  batch <- c(cell_data$batch, rep(paste0("pt_",drug_name), dim(pt_exp_set)[1]))
  
  pt_labels <- get(paste0("pt_", drug_name,"_resp"))
  
  #labels
  if (cell_set_name != "nibr") {
    auc <- get(drug_name, cell_data$auc)
    auc_g <- get(paste0(drug_name,"_g"), cell_data$auc)
    ic50 <- get(drug_name, cell_data$ic50)
    ic50_g <- get(paste0(drug_name,"_g"), cell_data$ic50)
    viab <- get(drug_name, cell_data$viab)
    viab_g <- get(paste0(drug_name,"_g"), cell_data$viab)
    
    auc <- c(auc, pt_labels)
    auc_g <- c(auc_g, pt_labels)
    ic50 <- c(ic50, pt_labels)
    ic50_g <- c(ic50_g, pt_labels)
    viab <- c(viab, pt_labels)
    viab_g <- c(viab_g, pt_labels)
  } else {
    
    bestresp <- get(drug_name, cell_data$bestresp)
    bestresp_g <- get(paste0(drug_name,"_g"), cell_data$bestresp)
    bestavgresp <- get(drug_name, cell_data$bestavgresp)
    bestavgresp_g <- get(paste0(drug_name,"_g"), cell_data$bestavgresp)
    
    bestresp <- c(bestresp, pt_labels)
    bestresp_g <- c(bestresp_g, pt_labels)
    bestavgresp <- c(bestavgresp, pt_labels)
    bestavgresp_g <- c(bestavgresp_g, pt_labels)
    
  }
  
  if (cell_set_name != "nibr") {
    
    n <- which(is.na(auc))
    
    auc <- auc[-n]
    auc_g <- auc_g[-n]
    ic50 <- ic50[-n]
    ic50_g <- ic50_g[-n]
    viab <- viab[-n]
    viab_g <- viab_g[-n]
    
    labels <- list(auc=auc,auc_g=auc_g,ic50=ic50,ic50_g=ic50_g,viab=viab,viab_g=viab_g)
  } else {
    n <- which(is.na(bestresp))
    
    bestresp <- bestresp[-n]
    bestresp_g <- bestresp_g[-n]
    bestavgresp <- bestavgresp[-n]
    bestavgresp_g <- bestavgresp_g[-n]
    
    labels=list(bestresp=bestresp, bestresp_g=bestresp_g,
                bestavgresp=bestavgresp, bestavgresp_g=bestavgresp_g)
  }
  
  exp <- exp[-n,]
  batch <- batch[-n]
  tissue <- tissue[-n]
  
  comb <- list(exp=exp,labels=labels,batch=batch,tissue=tissue)
  return(comb)
}

comb_ctrp_erlotinib <-combine("ctrp","erlotinib")
comb_ctrp_bortezomib <-combine("ctrp","bortezomib")
comb_ctrp_docetaxel <-combine("ctrp","docetaxel")

comb_gdsc_erlotinib <-combine("gdsc","erlotinib")
comb_gdsc_bortezomib <-combine("gdsc","bortezomib")
comb_gdsc_docetaxel <-combine("gdsc","docetaxel")

comb_gcsi_erlotinib <-combine("gcsi","erlotinib")
comb_gcsi_bortezomib <-combine("gcsi","bortezomib")
comb_gcsi_docetaxel <-combine("gcsi","docetaxel")

comb_nibr_erlotinib <-combine("nibr","erlotinib")

# now let's make tables with more than one cell line set
# function to create cell line data list for multiple sets (combining single cell line data)
data_prep_mult <- function(data1, data2) {
  genes <- intersect(colnames(data1$exp), colnames(data2$exp))
  exp <- rbind(data1$exp[,genes], data2$exp[,genes])
  rownames(exp) <- make.names(rownames(exp), unique = T)
  
  tissue <- c(data1$tissue, data2$tissue)
  batch <- c(data1$batch, data2$batch)
  
  if (deparse(substitute(data2))!="nibr_data")
  {
    auc <- as.list(rbind(as.data.frame(data1$auc), as.data.frame(data2$auc)))
    ic50 <- as.list(rbind(as.data.frame(data1$ic50), as.data.frame(data2$ic50)))
    viab <- as.list(rbind(as.data.frame(data1$viab), as.data.frame(data2$viab)))
    
    #let's convert factor vectors to the character vectors to make it consistent with already created objects of this type
    auc <-  lapply(auc, as.character)
    ic50 <-  lapply(ic50, as.character)
    viab <-  lapply(viab, as.character)
    
  } else {
    auc <- list(erlotinib=c(data1$auc$erlotinib, data2$bestresp$erlotinib),
                erlotinib_g=c(data1$auc$erlotinib_g, data2$bestresp$erlotinib_g))
    ic50 <- list(erlotinib=c(data1$ic50$erlotinib, data2$bestresp$erlotinib),
                 erlotinib_g=c(data1$ic50$erlotinib_g, data2$bestresp$erlotinib_g))
    viab <- list(erlotinib=c(data1$viab$erlotinib, data2$bestresp$erlotinib),
                 erlotinib_g=c(data1$viab$erlotinib_g, data2$bestresp$erlotinib_g))
  }
  
  result <- list(exp=exp, auc=auc, ic50=ic50, viab=viab, tissue=tissue, batch=batch)
}

#ctrp_gdsc
ctrp_gdsc_data <- data_prep_mult(ctrp_data, gdsc_data)
#ctrp_gdsc_gcsi
ctrp_gdsc_gcsi_data <- data_prep_mult(ctrp_gdsc_data, gcsi_data)
#ctrp_gdsc_gcsi_nibr (only for erlotinib)
ctrp_gdsc_gcsi_nibr_data <- data_prep_mult(ctrp_gdsc_gcsi_data, nibr_data)


#adding patient data 
comb_ctrp_gdsc_erlotinib <-combine("ctrp_gdsc","erlotinib")
comb_ctrp_gdsc_bortezomib <-combine("ctrp_gdsc","bortezomib")
comb_ctrp_gdsc_docetaxel <-combine("ctrp_gdsc","docetaxel")

comb_ctrp_gdsc_gcsi_erlotinib <-combine("ctrp_gdsc_gcsi","erlotinib")
comb_ctrp_gdsc_gcsi_bortezomib <-combine("ctrp_gdsc_gcsi","bortezomib")
comb_ctrp_gdsc_gcsi_docetaxel <-combine("ctrp_gdsc_gcsi","docetaxel")

comb_ctrp_gdsc_gcsi_nibr_erlotinib <-combine("ctrp_gdsc_gcsi_nibr","erlotinib")





comb_tables_list <- c("comb_ctrp_bortezomib", "comb_ctrp_docetaxel", "comb_ctrp_erlotinib",
                      "comb_gdsc_bortezomib", "comb_gdsc_docetaxel", "comb_gdsc_erlotinib",
                      "comb_gcsi_bortezomib", "comb_gcsi_docetaxel", "comb_gcsi_erlotinib",
                      "comb_nibr_erlotinib",
                      "comb_ctrp_gdsc_bortezomib", "comb_ctrp_gdsc_docetaxel", "comb_ctrp_gdsc_erlotinib",
                      "comb_ctrp_gdsc_gcsi_bortezomib", "comb_ctrp_gdsc_gcsi_docetaxel", "comb_ctrp_gdsc_gcsi_erlotinib",
                      "comb_ctrp_gdsc_gcsi_nibr_erlotinib")

#removing some objects from memory
rm(list=setdiff(ls(), c(comb_tables_list,"comb_tables_list")))

# Assessing principal components (on expression data) before ComBat
show_pca <- function(data, ber=F) {
  if (ber==T){
    pca <- prcomp(data$exp_ber)}
  else {pca <- prcomp(data$exp)}
  
  g <- ggbiplot(pca, obs.scale = 1, choices = c(1,2), var.scale = 1, 
                groups = as.factor(data$batch), ellipse = TRUE, 
                circle = TRUE, varname.size=0, var.axes = FALSE, 
                labels.size = 10)
  return(g)
}

pca_before <- lapply(comb_tables_list, function(x){
  show_pca(get(x))
})
#save(pca_before, file="pca_before.RData")

# applying ComBat
# we need to transpose exp matrix
combating <- function(comb_data){
  exp <- t(comb_data$exp)
  batch <- as.factor(comb_data$batch)
  if("auc" %in% names(comb_data$labels)) {
    response <- as.factor(comb_data$labels$auc) 
  } else { response <- as.factor(comb_data$labels$bestresp) }
  
  modcombat = model.matrix(~1, data=response)
  combat_edata = ComBat(dat=exp, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=F, BPPARAM = MulticoreParam(1))
  return(t(combat_edata))
}

#running combat on each of 17 comb. sets
for (i in 1:17)
{
  x <- comb_tables_list[i]
  command <- paste0(x,"$exp_ber <- combating(get(x))")
  eval(parse(text=command))
}

#comb_ctrp_docetaxel$exp_ber <- combating(comb_ctrp_docetaxel)
pca_after <- lapply(comb_tables_list, function(x){
  show_pca(get(x),ber=T)
})

#plot before and after
library(gridExtra)
erl1 <- arrangeGrob(grobs=pca_before[c(3,6,9,10,13,16,17)], ncol=1, nrow=7)
erl2 <- arrangeGrob(grobs=pca_after[c(3,6,9,10,13,16,17)], ncol=1, nrow=7) 
bor1 <- arrangeGrob(grobs=pca_before[c(1,4,7,11,14)], ncol=1, nrow=5)
bor2 <- arrangeGrob(grobs=pca_after[c(1,4,7,11,14)], ncol=1, nrow=5) 
doc1 <- arrangeGrob(grobs=pca_before[c(2,5,8,12,15)], ncol=1, nrow=5)
doc2 <- arrangeGrob(grobs=pca_after[c(2,5,8,12,15)], ncol=1, nrow=5) 
grid.arrange(erl1,erl2,nrow=1,ncol=2,top="Erlotinib")
grid.arrange(bor1,bor2,nrow=1,ncol=2,top="Bortezomib")
grid.arrange(doc1,doc2,nrow=1,ncol=2,top="Docetaxel")

#saving all comb. data
dir.create("comb_data")
for (i in 1:17)
{
  x <- comb_tables_list[i]
  command <- paste0('save(',x,', file="comb_data/',x,'.RData")')
  eval(parse(text=command))
}
save(comb_tables_list, file="comb_data/comb_tables_list.RData")
