# written in January 2017, modified on 26.10.18 
# Chapter 2 of the thesis. 
# Methods for improving drug response consistency: I. Cell line filtering, II. GLDS correction 

source("https://bioconductor.org/biocLite.R")
biocLite("PharmacoGx")
library(PharmacoGx)

path <- "/data/kurilov/genestack/phd/work_2018/thesis_code/"
setwd(file.path(path, "chapter2_inconsistency"))

# load CTRP, GDSC data
load(file.path(path, "data/ctrp_table.RData"))
load(file.path(path, "data/gdsc_table.RData"))


### I. cell line filtering ###

#tables for common drugs on my data. Targets: HDAC, EGFR, MEK1/MEK2, BRAF, HSP90, mTOR

com_lines <- intersect(rownames(gdsc_table), rownames(ctrp_table))
table_hdac <- cbind(gdsc_table[com_lines,c("MS-275", "Belinostat","Tubastatin A","Vorinostat")], ctrp_table[com_lines, c("entinostat", "belinostat","tubastatin A", "vorinostat")])
table_egfr <- cbind(gdsc_table[com_lines,c("Lapatinib", "Erlotinib", "Afatinib", "Gefitinib")], ctrp_table[com_lines, c("lapatinib", "erlotinib", "afatinib", "gefitinib")])
table_mek <- cbind(gdsc_table[com_lines,c("Trametinib", "selumetinib")], ctrp_table[com_lines, c("trametinib","selumetinib")])
table_braf <- cbind(gdsc_table[com_lines,c("PLX4720", "Dabrafenib")], ctrp_table[com_lines, c("PLX-4720","dabrafenib")])
table_hsp90 <- cbind(gdsc_table[com_lines,c("SNX-2112", "17-AAG")], ctrp_table[com_lines, c("SNX-2112","tanespimycin")])
table_mtor <- cbind(gdsc_table[com_lines,c("OSI-027", "Rapamycin", "BEZ235", "Temsirolimus", "AZD8055")], ctrp_table[com_lines, c("OSI-027", "sirolimus", "NVP-BEZ235", "temsirolimus", "AZD8055")])

# cell line filtering method, steps:
# 1) corr before and #of lines, 2) scaling, 3) subtract from each column dataset's mean column 
# 4) calc. corr for rows < threshold and # of lines 5) calculating average cor. for number of random subsets (10 subs.)
filtering <- function(table, percent) {
  n <- ncol(table)/2
  percentile <- percent/2
  
  table_scaled <- scale(table, center=F)
  mean1 <- apply(table_scaled[,1:n], 1, function(x) {mean(x, na.rm=T)})
  mean2 <- apply(table_scaled[,(n+1):(2*n)], 1, function(x) {mean(x, na.rm=T)})
  table_mean_subt <- cbind((table_scaled[,1:n]-mean1), (table_scaled[,(n+1):(2*n)]-mean2))
  #print(summary(table_mean_subt))
  
  for (i in 1:n) {
    print(" ")
    print(colnames(table)[i])
    cor <- cor(table[,i], table[,i+n], method="spearman", use="na.or.complete")
    print(paste0("corr. before: ", signif(cor)))
    lines_num <- nrow(na.omit(cbind(table[,i], table[,i+n])))
    print(paste0("number of common lines: ",lines_num))
    
    table_new <- cbind(table[,i], table[,i+n])
    threshold11 <- quantile(table_mean_subt[,i], percentile, na.rm=T)
    threshold12 <- quantile(table_mean_subt[,i], 1-percentile, na.rm=T)
    threshold21 <- quantile(table_mean_subt[,i+n], percentile, na.rm=T)
    threshold22 <- quantile(table_mean_subt[,i+n], 1-percentile, na.rm=T)
    
    ex1 <- c(which(table_mean_subt[,i]<threshold11), which(table_mean_subt[,i]>threshold12))
    ex2 <- c(which(table_mean_subt[,i+n]<threshold21), which(table_mean_subt[,i+n]>threshold22))
    ex <- union(ex1, ex2)
    table_new <- table_new[-ex,]
    
    cor_new <- cor(table_new[,1], table_new[,2], method="spearman", use="na.or.complete")
    print(paste0("corr. after: ", signif(cor_new)))
    lines_num_new <- nrow(na.omit(table_new))
    print(paste0("number of common lines after filtering: ",lines_num_new))
    
    table_new2 <- cbind(table[,i], table[,i+n])
    cor_random_subset <- vector()
    for (j in 1:10){
      a <- which(!(is.na(table_new2[,1])))
      b <- which(!(is.na(table_new2[,2])))
      a_ex <- a[sample(length(a)*percent)]
      b_ex <- b[sample(length(b)*percent)]
      ab_ex <- union(a_ex, b_ex)
      cor_subs <- cor(table_new2[-ab_ex,1], table_new2[-ab_ex,2], method="spearman", use="na.or.complete")
      cor_random_subset[j] <- cor_subs
    }
    print(paste0("average corr. after 10 random subsetting: ",signif(mean(cor_random_subset))))
    
    
  }
}

#producing filtering results
filtering(table_hdac, 0.1)
filtering(table_egfr, 0.1)
filtering(table_mek, 0.1)
filtering(table_braf, 0.1)
filtering(table_hsp90, 0.1)
filtering(table_mtor, 0.1)



### II. GLDS correction ###

# getting PharmacoGx annotation information
setwd(file.path(path, "data"))
CTRP <- downloadPSet("CTRPv2")
GDSC <- downloadPSet("GDSC1000")
setwd(file.path(path, "chapter2_inconsistency"))

# getting lists of unrelated drugs for each component
summary(as.factor(drugInfo(GDSC)$TARGET.PATHWAY))
drugInfo(GDSC)$TARGET[which(drugInfo(GDSC)$TARGET.PATHWAY=="TOR signaling")]

# HDAC
gdsc_hdac_drugs <- drugInfo(GDSC)$DRUG.NAME[which(drugInfo(GDSC)$TARGET.PATHWAY=="chromain  histone acetylation")]
ctrp_hdac_drugs <- drugInfo(CTRP)$cpd_name[grep("HDAC", drugInfo(CTRP)$gene_symbol_of_protein_target, value=F)]
# EGFR
gdsc_egfr_drugs <- drugInfo(GDSC)$DRUG.NAME[which(drugInfo(GDSC)$TARGET.PATHWAY=="EGFR signaling")]
ctrp_egfr_drugs <- c(drugInfo(CTRP)$cpd_name[grep("EGFR", drugInfo(CTRP)$gene_symbol_of_protein_target, value=F)],
                     drugInfo(CTRP)$cpd_name[grep("ERBB2", drugInfo(CTRP)$gene_symbol_of_protein_target, value=F)])


# MEK1, MEK2 and BRAF
gdsc_mek_braf_drugs <- drugInfo(GDSC)$DRUG.NAME[which(drugInfo(GDSC)$TARGET.PATHWAY=="ERK MAPK signaling")]
ctrp_mek_braf_drugs <- c(drugInfo(CTRP)$cpd_name[grep("ERK", drugInfo(CTRP)$gene_symbol_of_protein_target, value=F)],
                         drugInfo(CTRP)$cpd_name[grep("CRAF", drugInfo(CTRP)$gene_symbol_of_protein_target, value=F)],
                         drugInfo(CTRP)$cpd_name[grep("RSK", drugInfo(CTRP)$gene_symbol_of_protein_target, value=F)],
                         drugInfo(CTRP)$cpd_name[grep("BRAF", drugInfo(CTRP)$gene_symbol_of_protein_target, value=F)],
                         drugInfo(CTRP)$cpd_name[grep("MAP2K", drugInfo(CTRP)$gene_symbol_of_protein_target, value=F)])

gdsc_mek_drugs <- gdsc_mek_braf_drugs
gdsc_braf_drugs <- gdsc_mek_braf_drugs
ctrp_mek_drugs <- ctrp_mek_braf_drugs
ctrp_braf_drugs <- ctrp_mek_braf_drugs

# HSP90
gdsc_hsp90_drugs <- drugInfo(GDSC)$DRUG.NAME[grep("HSP", drugInfo(GDSC)$TARGET, value=F)]
ctrp_hsp90_drugs <- drugInfo(CTRP)$cpd_name[grep("HSP", drugInfo(CTRP)$gene_symbol_of_protein_target, value=F)]
#mTOR
gdsc_mtor_drugs <- drugInfo(GDSC)$DRUG.NAME[which(drugInfo(GDSC)$TARGET.PATHWAY=="TOR signaling")]
ctrp_mtor_drugs <- c(drugInfo(CTRP)$cpd_name[grep("MTOR", drugInfo(CTRP)$gene_symbol_of_protein_target, value=F)],
                     drugInfo(CTRP)$cpd_name[grep("p70", drugInfo(CTRP)$gene_symbol_of_protein_target, value=F)])


auc_gdsc_c <- gdsc_table[com_lines,18483:18733]
auc_ctrp_c <- ctrp_table[com_lines,45988:46532]

unrel_drugs <- function(drug, data) {
  cor <- apply(data, 2, function(x) {
    cor(data[,drug],x, use="na.or.complete", method="spearman")
  })
  unrel_drugs <- colnames(data)[which(cor<0.15)]
  return(unrel_drugs)
}

hdac_unrel <- list(unrel_drugs("MS-275", auc_gdsc_c), unrel_drugs("Belinostat", auc_gdsc_c), unrel_drugs("Tubastatin A", auc_gdsc_c), unrel_drugs("Vorinostat", auc_gdsc_c),
                   unrel_drugs("entinostat", auc_ctrp_c), unrel_drugs("belinostat", auc_ctrp_c), unrel_drugs("tubastatin A", auc_ctrp_c), unrel_drugs("vorinostat", auc_ctrp_c))

egfr_unrel <- list(unrel_drugs("Lapatinib", auc_gdsc_c), unrel_drugs("Erlotinib", auc_gdsc_c), unrel_drugs("Afatinib", auc_gdsc_c), unrel_drugs("Gefitinib", auc_gdsc_c),
                   unrel_drugs("lapatinib", auc_ctrp_c), unrel_drugs("erlotinib", auc_ctrp_c), unrel_drugs("afatinib", auc_ctrp_c), unrel_drugs("gefitinib", auc_ctrp_c))

mek_unrel <- list(unrel_drugs("Trametinib", auc_gdsc_c), unrel_drugs("selumetinib", auc_gdsc_c), 
                  unrel_drugs("trametinib", auc_ctrp_c), unrel_drugs("selumetinib", auc_ctrp_c))

braf_unrel <- list(unrel_drugs("PLX4720", auc_gdsc_c), unrel_drugs("Dabrafenib", auc_gdsc_c), 
                   unrel_drugs("PLX-4720", auc_ctrp_c), unrel_drugs("dabrafenib", auc_ctrp_c))

hsp90_unrel <- list(unrel_drugs("SNX-2112", auc_gdsc_c), unrel_drugs("17-AAG", auc_gdsc_c), 
                    unrel_drugs("SNX-2112", auc_ctrp_c), unrel_drugs("tanespimycin", auc_ctrp_c))

mtor_unrel <- list(unrel_drugs("OSI-027", auc_gdsc_c), unrel_drugs("Rapamycin", auc_gdsc_c), unrel_drugs("BEZ235", auc_gdsc_c), unrel_drugs("Temsirolimus", auc_gdsc_c), unrel_drugs("AZD8055", auc_gdsc_c),
                   unrel_drugs("OSI-027", auc_ctrp_c), unrel_drugs("sirolimus", auc_ctrp_c), unrel_drugs("NVP-BEZ235", auc_ctrp_c), unrel_drugs("temsirolimus", auc_ctrp_c), unrel_drugs("AZD8055", auc_ctrp_c)) 


# GLDS filtering steps: removing drugs from the same class, subtracting mean AUC from unrelated drugs and assessing the consistency 
glds_correction <- function(table, unrel_drug_data, rel_drug_gdsc, rel_drug_ctrp) {
  n <- ncol(table)/2
  
  for (i in 1:n) {
    print(" ")
    print(colnames(table)[i])
    cor <- cor(table[,i], table[,i+n], method="spearman", use="na.or.complete")
    print(paste0("corr. before: ", signif(cor)))
    lines_num <- nrow(na.omit(cbind(table[,i], table[,i+n])))
    print(paste0("number of common lines: ",lines_num))
    
    unrel_drugs_gdsc <- setdiff(unrel_drug_data[[i]], rel_drug_gdsc)
    unrel_drugs_ctrp <- setdiff(unrel_drug_data[[i+n]], rel_drug_ctrp)
    
    drug1_new <- table[,i]-rowMeans(auc_gdsc_c[,unrel_drugs_gdsc], na.rm=T)
    drug2_new <- table[,i+n]-rowMeans(auc_ctrp_c[,unrel_drugs_ctrp], na.rm=T)
    cor_new <- cor(drug1_new, drug2_new, method="spearman", use="na.or.complete")
    print(paste0("corr. after: ", signif(cor_new)))
  }
}

# producing GLDS filtering results
glds_correction(table_hdac, hdac_unrel, gdsc_hdac_drugs, ctrp_hdac_drugs)
glds_correction(table_egfr, egfr_unrel, gdsc_egfr_drugs, ctrp_egfr_drugs)
glds_correction(table_mek, mek_unrel, gdsc_mek_braf_drugs, ctrp_mek_braf_drugs)
glds_correction(table_braf, braf_unrel, gdsc_mek_braf_drugs, ctrp_mek_braf_drugs)
glds_correction(table_hsp90, hsp90_unrel, gdsc_hsp90_drugs, ctrp_hsp90_drugs)
glds_correction(table_mtor, mtor_unrel, gdsc_mtor_drugs, ctrp_mtor_drugs)


