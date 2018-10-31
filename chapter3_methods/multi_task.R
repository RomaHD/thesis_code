# written in February 2017, modified on 27.10.18 
# Chapter 3 of the thesis. 
# Methods for improving accuracy of drug response prediction: 
# 1. Multi-task modelling: all tissues, tissue-specific, plotting the results
#
# datasets:CTRP, GDSC

source("https://bioconductor.org/biocLite.R")
biocLite("PharmacoGx")
library(PharmacoGx)
library(glmnet)
library(survcomp)
library(ggplot2)
library(tidyr)
library(gridExtra)


path <- "/data/kurilov/genestack/phd/work_2018/thesis_code/"
setwd(file.path(path, "chapter3_methods"))

# load CTRP, GDSC (imputed drug resp data because of glmnet can't handle missing data)
load(file.path(path, "data/ctrp_table_imp.RData"))
load(file.path(path, "data/gdsc_table_imp.RData"))

# PharmacoGx annotation data
setwd(file.path(path, "data"))
CTRP <- downloadPSet("CTRPv2")
GDSC <- downloadPSet("GDSC1000")
setwd(file.path(path, "chapter3_methods"))

names <- cbind(c("MS-275","Tubastatin A", "Belinostat", "Vorinostat", "Lapatinib", "Erlotinib", "Afatinib", "Gefitinib",
                 "Trametinib", "selumetinib", "PLX4720", "Dabrafenib", "SNX-2112", "17-AAG", "OSI-027", 
                 "Rapamycin", "BEZ235", "Temsirolimus", "AZD8055"),
               c("entinostat","tubastatin A", "belinostat","vorinostat", "lapatinib", "erlotinib", "afatinib", "gefitinib",
                 "trametinib","selumetinib", "PLX-4720","dabrafenib", "SNX-2112","tanespimycin", "OSI-027", 
                 "sirolimus", "NVP-BEZ235", "temsirolimus", "AZD8055"),
               c(rep("hdac", 4), rep("egfr", 4), rep("mek", 2), rep("braf", 2), rep("hsp90", 2), rep("mtor", 5)))
colnames(names) <- c("gdsc", "ctrp", "class")

gdsc_classes <- list(hdac=drugInfo(GDSC)$DRUG.NAME[grep("HDAC",drugInfo(GDSC)$TARGET)],
                     egfr=drugInfo(GDSC)$DRUG.NAME[(grep("EGFR",drugInfo(GDSC)$TARGET))[c(3:5,8:10,12,14)]],
                     mek=drugInfo(GDSC)$DRUG.NAME[grep("MEK1, MEK2",drugInfo(GDSC)$TARGET)],
                     braf=drugInfo(GDSC)$DRUG.NAME[grep("BRAF",drugInfo(GDSC)$TARGET)],
                     hsp90=drugInfo(GDSC)$DRUG.NAME[grep("HSP90",drugInfo(GDSC)$TARGET)],
                     mtor=drugInfo(GDSC)$DRUG.NAME[grep("mTOR",drugInfo(GDSC)$TARGET)])

ctrp_classes <- list(hdac=drugInfo(CTRP)$cpd_name[grep("HDAC",drugInfo(CTRP)$gene_symbol_of_protein_target)][c(1:17,23:26)],
                     egfr=drugInfo(CTRP)$cpd_name[grep("EGFR",drugInfo(CTRP)$gene_symbol_of_protein_target)][c(1:6,8:11)],
                     mek=drugInfo(CTRP)$cpd_name[grep("MAP2K1;MAP2K2",drugInfo(CTRP)$gene_symbol_of_protein_target)][c(11,12,14)] ,
                     braf=drugInfo(CTRP)$cpd_name[grep("BRAF",drugInfo(CTRP)$gene_symbol_of_protein_target)][c(1:5,10,13:14)],
                     hsp90=drugInfo(CTRP)$cpd_name[grep("HSP90",drugInfo(CTRP)$gene_symbol_of_protein_target)][c(1,2,7)],
                     mtor=drugInfo(CTRP)$cpd_name[grep("MTOR",drugInfo(CTRP)$gene_symbol_of_protein_target)][-5])

drug_classes <- list(gdsc=gdsc_classes, ctrp=ctrp_classes)

gen_last_gdsc <- max(grep("_", colnames(gdsc_table_imp)))
gen_last_ctrp <- max(grep("_", colnames(ctrp_table_imp)))
f_set <- intersect(colnames(gdsc_table_imp)[1:gen_last_gdsc], colnames(ctrp_table_imp)[1:gen_last_ctrp])
# 16126 exp + 146 mut 

#let's remove rows with NA values
na_gdsc <- which(is.na(gdsc_table_imp[,1]))
na_ctrp <- c(which(is.na(ctrp_table_imp[,1])), which(is.na(ctrp_table_imp[,"mut_MAP3K1"])))

table_gdsc <- cbind(gdsc_table_imp[-na_gdsc,f_set], gdsc_table_imp[-na_gdsc,(1+gen_last_gdsc):ncol(gdsc_table_imp)])
table_ctrp <- cbind(ctrp_table_imp[-na_ctrp,f_set], ctrp_table_imp[-na_ctrp,(1+gen_last_ctrp):ncol(ctrp_table_imp)])


r2 <- function(v1, v2) {
  value <- cor(v1, v2, use="na.or.complete")
  return(value^2)
}

# concordance.index 
ci <- function(v1, v2) {
  value <- concordance.index(x=v1, surv.time=v2, surv.event=rep(1,length(v1)))
  return(1-value$c.index)
}


require(doMC)
registerDoMC(cores=4)

### 1. analysis for all tissues ###


res_table <- matrix(NA, nrow=12, ncol=38)
res_table_ci <- matrix(NA, nrow=12, ncol=38)

alpha_vect=c(0, 0.5, 1)

for (j in 1:3) {
  alpha  <- alpha_vect[j]
  
  gdsc_hdac <- cv.glmnet(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$hdac], family = "mgaussian", alpha=alpha)
  gdsc_egfr <- cv.glmnet(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$egfr], family = "mgaussian", alpha=alpha)
  gdsc_mek <- cv.glmnet(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$mek], family = "mgaussian", alpha=alpha)
  gdsc_braf <- cv.glmnet(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$braf], family = "mgaussian", alpha=alpha)
  gdsc_hsp90 <- cv.glmnet(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$hsp90], family = "mgaussian", alpha=alpha)
  gdsc_mtor <- cv.glmnet(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$mtor], family = "mgaussian", alpha=alpha)
  
  ctrp_hdac <- cv.glmnet(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$hdac], family = "mgaussian", alpha=alpha)
  ctrp_egfr <- cv.glmnet(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$egfr], family = "mgaussian", alpha=alpha)
  ctrp_mek <- cv.glmnet(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$mek], family = "mgaussian", alpha=alpha)
  ctrp_braf <- cv.glmnet(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$braf], family = "mgaussian", alpha=alpha)
  ctrp_hsp90 <- cv.glmnet(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$hsp90], family = "mgaussian", alpha=alpha)
  ctrp_mtor <- cv.glmnet(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$mtor], family = "mgaussian", alpha=alpha)
  
  for (i in 1:19) {
    print(i)
    drug_gdsc <- names[i,1]
    drug_ctrp <- names[i,2]
    class <- names[i,3]
    
    # GDSC -> CTRP
    gdsc_m = get(paste0("gdsc_", class))
    gdsc_s = cv.glmnet(x=table_gdsc[,f_set], y=table_gdsc[, drug_gdsc], family = "gaussian", alpha=alpha)
    
    pred_gdsc_m_gdsc <- predict.cv.glmnet(gdsc_m, newx=table_gdsc[,f_set], s="lambda.min")
    pred_gdsc_s_gdsc <- predict.cv.glmnet(gdsc_s, newx=table_gdsc[,f_set], s="lambda.min")
    pred_gdsc_m_ctrp <- predict.cv.glmnet(gdsc_m, newx=table_ctrp[,f_set], s="lambda.min")
    pred_gdsc_s_ctrp <- predict.cv.glmnet(gdsc_s, newx=table_ctrp[,f_set], s="lambda.min")
    
    # CTRP -> GDSC
    ctrp_m = get(paste0("ctrp_", class))
    ctrp_s = cv.glmnet(x=table_ctrp[,f_set], y=table_ctrp[, drug_ctrp], family = "gaussian", alpha=alpha)
    
    pred_ctrp_m_ctrp <- predict.cv.glmnet(ctrp_m, newx=table_ctrp[,f_set], s="lambda.min")
    pred_ctrp_s_ctrp <- predict.cv.glmnet(ctrp_s, newx=table_ctrp[,f_set], s="lambda.min")
    pred_ctrp_m_gdsc <- predict.cv.glmnet(ctrp_m, newx=table_gdsc[,f_set], s="lambda.min")
    pred_ctrp_s_gdsc <- predict.cv.glmnet(ctrp_s, newx=table_gdsc[,f_set], s="lambda.min")
    
    # R2
    r2_1 <- r2(pred_gdsc_m_gdsc[,drug_gdsc,1], table_gdsc[,drug_gdsc])
    r2_2 <- r2(pred_gdsc_m_ctrp[,drug_gdsc,1], table_ctrp[,drug_ctrp])
    r2_3 <- r2(pred_ctrp_m_ctrp[,drug_ctrp,1], table_ctrp[,drug_ctrp])
    r2_4 <- r2(pred_ctrp_m_gdsc[,drug_ctrp,1], table_gdsc[,drug_gdsc])
    r2_5 <- r2(pred_gdsc_s_gdsc, table_gdsc[,drug_gdsc])
    r2_6 <- r2(pred_gdsc_s_ctrp, table_ctrp[,drug_ctrp])
    r2_7 <- r2(pred_ctrp_s_ctrp, table_ctrp[,drug_ctrp])
    r2_8 <- r2(pred_ctrp_s_gdsc, table_gdsc[,drug_gdsc])
    
    #saving res. 
    col <- i*2-1
    row <- j*4-3
    res_table[row:(row+3),col:(col+1)] <- signif(cbind(c(r2_1,r2_2,r2_3,r2_4), c(r2_5,r2_6,r2_7,r2_8)), digits=3)
    res_table_ci[row:(row+3),col:(col+1)] <- signif(cbind(c(r2_1,r2_2,r2_3,r2_4), c(r2_5,r2_6,r2_7,r2_8)), digits=3)
  } 
}

write.table(res_table, file="multitask_glmnet_res.txt")


### 2. tissue-specific analysis ###

# tables for tissue subsetting
table_gdsc0 <- table_gdsc
table_ctrp0 <- table_ctrp

# tissue information
t_int <- intersect(cellInfo(GDSC)$tissueid, cellInfo(CTRP)$tissueid)
tissue_numbers <- cbind(summary(as.factor(cellInfo(GDSC)$tissueid))[t_int], summary(as.factor(cellInfo(CTRP)$tissueid))[t_int])
tissue_numbers <- tissue_numbers[order(tissue_numbers[,1]+tissue_numbers[,2], decreasing = T),]


#
require(doMC)
registerDoMC(cores=2)

r2_table <- matrix(NA, nrow=36, ncol=38)
ci_table <- matrix(NA, nrow=36, ncol=38)

alpha=1
for (j in 1:18) {
  tissue  <- rownames(tissue_numbers)[j]
  print(tissue)
  gdsc_lines <- intersect(rownames(table_gdsc0), cellInfo(GDSC)$Sample.Name[which(cellInfo(GDSC)$tissueid==tissue)])
  ctrp_lines <- intersect(rownames(table_ctrp0), cellInfo(CTRP)$cellid[which(cellInfo(CTRP)$tissueid==tissue)])
  
  table_gdsc <- table_gdsc0[gdsc_lines,]
  table_ctrp <- table_ctrp0[ctrp_lines,]
  
  gdsc_hdac <- cv.glmnet(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$hdac], family = "mgaussian", alpha=alpha)
  gdsc_egfr <- cv.glmnet(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$egfr], family = "mgaussian", alpha=alpha)
  gdsc_mek <- cv.glmnet(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$mek], family = "mgaussian", alpha=alpha)
  gdsc_braf <- cv.glmnet(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$braf], family = "mgaussian", alpha=alpha)
  gdsc_hsp90 <- cv.glmnet(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$hsp90], family = "mgaussian", alpha=alpha)
  gdsc_mtor <- cv.glmnet(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$mtor], family = "mgaussian", alpha=alpha)
  
  ctrp_hdac <- cv.glmnet(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$hdac], family = "mgaussian", alpha=alpha)
  ctrp_egfr <- cv.glmnet(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$egfr], family = "mgaussian", alpha=alpha)
  ctrp_mek <- cv.glmnet(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$mek], family = "mgaussian", alpha=alpha)
  ctrp_braf <- cv.glmnet(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$braf], family = "mgaussian", alpha=alpha)
  ctrp_hsp90 <- cv.glmnet(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$hsp90], family = "mgaussian", alpha=alpha)
  ctrp_mtor <- cv.glmnet(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$mtor], family = "mgaussian", alpha=alpha)
  
  for (i in 1:19) {
    print(i)
    drug_gdsc <- names[i,1]
    drug_ctrp <- names[i,2]
    class <- names[i,3]
    
    # GDSC -> CTRP
    gdsc_m = get(paste0("gdsc_", class))
    gdsc_s = cv.glmnet(x=table_gdsc[,f_set], y=table_gdsc[, drug_gdsc], family = "gaussian", alpha=alpha)
    
    # pred_gdsc_m_gdsc <- predict.cv.glmnet(gdsc_m, newx=table_gdsc[,f_set], s="lambda.min")
    # pred_gdsc_s_gdsc <- predict.cv.glmnet(gdsc_s, newx=table_gdsc[,f_set], s="lambda.min")
    pred_gdsc_m_ctrp <- predict.cv.glmnet(gdsc_m, newx=table_ctrp[,f_set], s="lambda.min")
    pred_gdsc_s_ctrp <- predict.cv.glmnet(gdsc_s, newx=table_ctrp[,f_set], s="lambda.min")
    
    # CTRP -> GDSC
    ctrp_m = get(paste0("ctrp_", class))
    ctrp_s = cv.glmnet(x=table_ctrp[,f_set], y=table_ctrp[, drug_ctrp], family = "gaussian", alpha=alpha)
    
    # pred_ctrp_m_ctrp <- predict.cv.glmnet(ctrp_m, newx=table_ctrp[,f_set], s="lambda.min")
    # pred_ctrp_s_ctrp <- predict.cv.glmnet(ctrp_s, newx=table_ctrp[,f_set], s="lambda.min")
    pred_ctrp_m_gdsc <- predict.cv.glmnet(ctrp_m, newx=table_gdsc[,f_set], s="lambda.min")
    pred_ctrp_s_gdsc <- predict.cv.glmnet(ctrp_s, newx=table_gdsc[,f_set], s="lambda.min")
    
    # R2
    r2_1 <- r2(pred_gdsc_m_ctrp[,drug_gdsc,1], table_ctrp[,drug_ctrp])
    r2_2 <- r2(pred_ctrp_m_gdsc[,drug_ctrp,1], table_gdsc[,drug_gdsc])
    r2_3 <- r2(pred_gdsc_s_ctrp, table_ctrp[,drug_ctrp])
    r2_4 <- r2(pred_ctrp_s_gdsc, table_gdsc[,drug_gdsc])
    
    # Concordance index
    ci_1 <- ci(pred_gdsc_m_ctrp[,drug_gdsc,1], table_ctrp[,drug_ctrp])
    ci_2 <- ci(pred_ctrp_m_gdsc[,drug_ctrp,1], table_gdsc[,drug_gdsc])
    ci_3 <- ci(pred_gdsc_s_ctrp, table_ctrp[,drug_ctrp])
    ci_4 <- ci(pred_ctrp_s_gdsc, table_gdsc[,drug_gdsc])
    
    #saving res. 
    col <- i*2-1
    row <- j*2-1
    #res_table[row:(row+3),col:(col+1)] <- signif(cbind(c(r2_1,r2_2,r2_3,r2_4), c(r2_5,r2_6,r2_7,r2_8)), digits=3)
    r2_table[row:(row+1),col:(col+1)] <- signif(cbind(c(r2_1,r2_2), c(r2_3,r2_4)), digits=3)
    ci_table[row:(row+1),col:(col+1)] <- signif(cbind(c(ci_1,ci_2), c(ci_3,ci_4)), digits=3)
  } 
}

write.table(r2_table, file="r2_res_glmnet_tissue_spec.txt")
write.table(ci_table, file="ci_res_glmnet_tissue_spec.txt")


### 3. plotting the results ###
library(tidyr)
library(gridExtra)
library(ggplot2)

# plotting R2 data for multi-task glmnet tests

res_table <- read.table("multitask_glmnet_res.txt")

title_list <- c("Ridge", "Elastic net", "Lasso")
#col_list <- list(c(2,4), c(6,8), c(10,12))
col_list <- list(c(4), c(8), c(12))
for (i in 1:3) {

  table1 <- res_table[col_list[[i]],]
  colnames(table1) <- NULL
  rownames(table1) <- NULL
  table2 <- matrix(c(table1[,seq(1,38, by=2)], table1[,seq(2,38, by=2)]), nrow=2, byrow = T)
  #rownames(table2) <- c("m_gdsc->ctrp", "m_ctrp->gdsc", "s_gdsc->ctrp", "s_ctrp->gdsc")
  rownames(table2) <- c("multi-task", "single output")
  colnames(table2) <- c("entinostat","tubastatin A", "belinostat","vorinostat", "lapatinib", "erlotinib", "afatinib", "gefitinib",
                        "trametinib","selumetinib", "PLX-4720","dabrafenib", "SNX-2112","tanespimycin", "OSI-027", 
                        "sirolimus", "NVP-BEZ235", "temsirolimus", "AZD8055")
  table3 <- cbind(rownames(table2), table2)
  colnames(table3)[1] <- "model"
  
  table3 <- as.data.frame(table3)
  table4 <- table3 %>% 
    gather(drug, r2, -model) 
  table4$r2<- as.numeric(as.character(table4$r2))
  table4$drug <- factor(table4$drug, levels = colnames(table2) )
  table4$model <- as.factor(unlist(table4$model))
  
  p <- ggplot(table4, aes(y=r2, x=drug)) +
    geom_bar(aes(fill=model), position="dodge", stat="identity") +
    scale_fill_manual(values=c("#0072B2", "#D55E00", "#9999CC", "pink1")) +
    labs(title = title_list[i])
  
  assign(paste0("s",i), p)
}

grid.arrange(s1,s2,s3, nrow=3, ncol=1)
ggsave("fig_13.pdf")  

#grid.arrange(p1,s1,p2,s2,p3,s3, nrow=3, ncol=2, top="left: training on GDSC, right: training on CTRP")


# plotting R2/ci from glmnet tisuue-specific tests
# tissue information
t_int <- intersect(cellInfo(GDSC)$tissueid, cellInfo(CTRP)$tissueid)
tissue_numbers <- cbind(summary(as.factor(cellInfo(GDSC)$tissueid))[t_int], summary(as.factor(cellInfo(CTRP)$tissueid))[t_int])
tissue_numbers <- tissue_numbers[order(tissue_numbers[,1]+tissue_numbers[,2], decreasing = T),]

r2_table <- read.table("r2_res_glmnet_tissue_spec.txt")

table_comb <- matrix(NA,18,38)
for (i in 1:18) {
  table_comb[i,] <- apply(r2_table[(2*i-1):(2*i),], 2, function(x) {mean(x, na.rm=T)})
}

table_comb_m <- table_comb[,seq(1,38, by=2)]
table_comb_s <- table_comb[,seq(2,38, by=2)]
table_comb_m <- cbind(rownames(tissue_numbers)[1:18], table_comb_m)
colnames(table_comb_m) <- c("tissue","entinostat","tubastatin A", "belinostat","vorinostat", "lapatinib", "erlotinib", "afatinib", "gefitinib",
                            "trametinib","selumetinib", "PLX-4720","dabrafenib", "SNX-2112","tanespimycin", "OSI-027", 
                            "sirolimus", "NVP-BEZ235", "temsirolimus", "AZD8055")
table_comb_s <- cbind(table_comb_m[,1],table_comb_s)
colnames(table_comb_s) <- colnames(table_comb_m)

table_comb_m <- as.data.frame(table_comb_m)
table_comb_m2 <- table_comb_m %>% 
  gather(drug, r2, -tissue) 

table_comb_s <- as.data.frame(table_comb_s)
table_comb_s2 <- table_comb_s %>% 
  gather(drug, r2, -tissue) 

table_comb_tidy <- cbind(rbind(table_comb_m2, table_comb_s2), c(rep("multi-task", 342), rep("single-task", 342)))
colnames(table_comb_tidy)[4] <- "model"
table_comb_tidy$drug <- as.factor(table_comb_tidy$drug)
table_comb_tidy$r2 <- as.numeric(table_comb_tidy$r2)
table_comb_tidy$model <- as.factor(table_comb_tidy$model)

theme_set(theme_bw(base_size = 12))
for (i in 1:18)
{
  tissue=rownames(tissue_numbers)[i]
  table <- table_comb_tidy[which(table_comb_tidy$tissue==tissue),]
  
  p <- ggplot(table, aes(y=r2, x=drug, fill=model)) +
    #geom_bar(position="dodge", stat="identity") + scale_y_continuous(limits=c(0,0.85),oob = rescale_none) +
    geom_bar(position="dodge", stat="identity") + scale_y_continuous(limits=c(0,0.85)) +
    #scale_fill_manual(values=c("#0072B2", "#D55E00")) +
    labs(title = tissue) + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())
  assign(paste0("p",i), p)
}

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,
             p10,p11,p12,p13,p14,p15,p16,p17,p18, nrow=9, ncol=2, top="R2 for Lasso models")

grid.arrange(p1,p2,p3,p4,p5,p6, nrow=6, ncol=1, top="R2 for Lasso models")
ggsave("fig_14.pdf") 
