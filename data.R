# 25.10.18 getting data for analyses in the chapters 2,3,6,7 of the thesis. (c) Roman Kurilov

source("http://bioconductor.org/biocLite.R")
biocLite("PharmacoGx")
library(PharmacoGx)
library(Biobase)
library(glmnet)

path <- "/data/kurilov/genestack/phd/work_2018/thesis_code/"
setwd(file.path(path, "data"))
source(file.path(path, "completeMatrixFunction.R"))

### I. data for chapters 2,3,7 ###

CCLE <- downloadPSet("CCLE")
CTRP <- downloadPSet("CTRPv2")
GDSC <- downloadPSet("GDSC1000")
gCSI <- downloadPSet("gCSI")

#CCLE/CTRP
ccle_exp <- summarizeMolecularProfiles(pSet=CCLE, mDataType="rna", cell.lines=cellNames(CCLE), summary.stat = 'mean', fill.missing = TRUE, verbose=TRUE)
ccle_cnv <- summarizeMolecularProfiles(pSet=CCLE, mDataType="cnv", cell.lines=cellNames(CCLE), summary.stat = 'mean', fill.missing = TRUE, verbose=TRUE)
ccle_mut <- summarizeMolecularProfiles(pSet=CCLE, mDataType="mutation", cell.lines=cellNames(CCLE), summary.stat = 'or', fill.missing = TRUE, verbose=TRUE)

ccle_exp <- t(exprs(ccle_exp))
ccle_cnv <- t(exprs(ccle_cnv))
ccle_mut <- t(exprs(ccle_mut))

#GDSC
gdsc_exp <- summarizeMolecularProfiles(pSet=GDSC, mDataType="rna", cell.lines=cellNames(GDSC), summary.stat = 'mean', fill.missing = TRUE, verbose=TRUE)
# cnv, mut, meth data
load("cn_mut_meth_gdsc1000.RData")

gdsc_exp <- t(exprs(gdsc_exp))
gdsc_cnv <- t(exprs(gdsc_cnv))
gdsc_mut <- t(exprs(gdsc_mut))

# gCSI
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
colnames(ccle_cnv) <- paste0("cnv_", colnames(ccle_cnv))
colnames(gdsc_cnv) <- paste0("cnv_", colnames(gdsc_cnv))
colnames(ccle_mut) <- paste0("mut_", colnames(ccle_mut))
colnames(gdsc_mut) <- paste0("mut_", colnames(gdsc_mut))

# drug_response
drugNames(GDSC) <- GDSC@drug$DRUG.NAME
drugNames(CTRP) <- CTRP@drug$cpd_name

ctrp_auc <- t(summarizeSensitivityProfiles(CTRP, sensitivity.measure="auc_recomputed", drugs=drugNames(CTRP)))
gdsc_auc <- t(summarizeSensitivityProfiles(GDSC, sensitivity.measure="auc_recomputed", drugs=drugNames(GDSC)))
gcsi_auc <- t(summarizeSensitivityProfiles(gCSI ,sensitivity.measure = "auc_recomputed", drugs=drugNames(gCSI)))


# 1-AUC correction for AUC data
ctrp_auc <- 1-ctrp_auc
gdsc_auc <- 1-gdsc_auc
gcsi_auc <- 1-gcsi_auc

#combining data
lines_ctrp <- intersect(cellInfo(CCLE)$cellid, cellInfo(CTRP)$cellid)
ctrp_table <- cbind(ccle_exp[lines_ctrp,], ccle_cnv[lines_ctrp,], ccle_mut[lines_ctrp,], ctrp_auc[lines_ctrp,])
gcsi_table <- cbind(gcsi_exp, gcsi_auc)
 # gdsc case
rownames(gdsc_exp) <- GDSC@cell$Sample.Name
rownames(gdsc_auc) <- GDSC@cell$Sample.Name
com_lines_gdsc <- intersect(rownames(cn_mut_meth_gdsc1000), rownames(gdsc_auc))
gdsc_table <- cbind(gdsc_exp[com_lines_gdsc,], cn_mut_meth_gdsc1000[com_lines_gdsc,], gdsc_auc[com_lines_gdsc,])


#removing redundant cols 
x <- apply(ctrp_table,2,function(x){sd(x, na.rm=T)})
sdz <- which(x==0)
ctrp_table <- ctrp_table[,-sdz]

# x <- apply(gdsc_table,2,function(x){sd(x, na.rm=T)})
# sdz <- which(x==0)
# sdn <- which(is.na(x))
# gdsc_table <- gdsc_table[,-c(sdz,sdn[1:16])]

##
y1 <- apply(ctrp_table[,1:45987],1,function(x){length(which(is.na(x)))==length(x)})
y2 <- apply(ctrp_table[,45988:46008],1,function(x){length(which(is.na(x)))==length(x)})
sdn <- which(y1)
sdn2 <- which(y2)
ctrp_table <- ctrp_table[-c(sdn, sdn2),]

# y1 <- apply(gdsc_table[,1:36721],1,function(x){length(which(is.na(x)))==length(x)})
# y2 <- apply(gdsc_table[,36722:36742],1,function(x){length(which(is.na(x)))==length(x)})
# sdn <- which(y1)
# sdn2 <- which(y2)
# gdsc_table <- gdsc_table[-c(sdn, sdn2),]


ctrp_table <- as.matrix(ctrp_table)
gdsc_table <- as.matrix(gdsc_table)
gcsi_table <- as.matrix(gcsi_table)
mode(ctrp_table) <- "numeric"
mode(gdsc_table) <- "numeric"
mode(gcsi_table) <- "numeric"

# imputing missing AUC values using the completeMatrix function from https://github.com/paulgeeleher/glds/blob/master/glds

gen_last_gdsc <- max(grep("_", colnames(gdsc_table)))
gen_last_ctrp <- max(grep("_", colnames(ctrp_table)))
gen_last_gcsi <- max(grep("_", colnames(gcsi_table)))

auc_gdsc_imp <- completeMatrix(gdsc_table[,(gen_last_gdsc+1):ncol(gdsc_table)], nPerms=15)
auc_ctrp_imp <- completeMatrix(ctrp_table[,(gen_last_ctrp+1):ncol(ctrp_table)], nPerms=15)
auc_gcsi_imp <- completeMatrix(gcsi_table[,(gen_last_gcsi+1):ncol(gcsi_table)], nPerms=15)

gdsc_table_imp <- cbind(gdsc_table[,1:gen_last_gdsc], auc_gdsc_imp)
ctrp_table_imp <- cbind(ctrp_table[,1:gen_last_ctrp], auc_ctrp_imp)
gcsi_table_imp <- cbind(gcsi_table[,1:gen_last_gcsi], auc_gcsi_imp)

# saving tables
save(ctrp_table, file="ctrp_table.RData")
save(gdsc_table, file="gdsc_table.RData")
save(gcsi_table, file="gcsi_table.RData")

save(ctrp_table_imp, file="ctrp_table_imp.RData")
save(gdsc_table_imp, file="gdsc_table_imp.RData")
save(gcsi_table_imp, file="gcsi_table_imp.RData")


