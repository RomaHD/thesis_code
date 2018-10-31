# written in February 2017, modified on 26.10.18 
# Chapter 2 of the thesis. 
# Methods for improving biomarker's consistency: GLDS correction 
# datasets:CTRP, GDSC1000, gCSI and NIBR PDXE

source("https://bioconductor.org/biocLite.R")
biocLite("PharmacoGx")
library(PharmacoGx)

path <- "/data/kurilov/genestack/phd/work_2018/thesis_code/"
setwd(file.path(path, "chapter2_inconsistency"))

# load CTRP, GDSC, gCSI and NIBR PDXE data
load(file.path(path, "data/ctrp_table_imp.RData"))
load(file.path(path, "data/gdsc_table_imp.RData"))
load(file.path(path, "data/gcsi_table_imp.RData"))

# NIBR PDXE data was extracted from NIBR PDXE paper (Gao H. et al., Nature Medicine 2015) supplementary data 
# https://media.nature.com/original/nature-assets/nm/journal/v21/n11/extref/nm.3954-S2.xlsx and combined into a single table
load(file.path(path, "data/table_bestresp_nibr.RData"))

# column numbers of the last column with genomic data in each dataset
gen_last_gdsc <- max(grep("_", colnames(gdsc_table_imp)))
gen_last_ctrp <- max(grep("_", colnames(ctrp_table_imp)))
gen_last_nibr <- max(grep("_", colnames(table_bestresp_nibr))) 
gen_last_gcsi <- max(grep("_", colnames(gcsi_table_imp)))

#imputing drug. resp values in nibr dataset
bestresp_nibr_imp <- completeMatrix(table_bestresp_nibr[,(gen_last_nibr+1):ncol(table_bestresp_nibr)], nPerms=15)
nibr_table_imp <- cbind(table_bestresp_nibr[,1:gen_last_nibr], bestresp_nibr_imp)

# PharmacoGx annotation data
setwd(file.path(path, "data"))
CTRP <- downloadPSet("CTRPv2")
GDSC <- downloadPSet("GDSC1000")
gCSI <- downloadPSet("gCSI")
setwd(file.path(path, "chapter2_inconsistency"))

## annotation information and main functions

# drug annotations
names <- cbind(c("5-Fluorouracil", "Bortezomib", "Crizotinib", "Docetaxel","Doxorubicin","Erlotinib", "GDC0941", "Gemcitabine","Lapatinib","MS-275","Paclitaxel","Rapamycin","Vorinostat","Tamoxiofen", "Trametinib"),
                c("fluorouracil", "bortezomib", "crizotinib", "docetaxel","doxorubicin","erlotinib", "GDC-0941", "gemcitabine","lapatinib","entinostat","paclitaxel","sirolimus","vorinostat","tamoxiofen", "trametinib"),
                c("NA",           "Bortezomib", "Crizotinib", "Docetaxel","Doxorubicin","Erlotinib", "GDC-0941", "Gemcitabine","lapatinib","MS-275","paclitaxel","Rapamycin","Vorinostat","NA",        "NA"),
               c("5FU",           "NA",         "NA",         "NA",       "NA",         "erlotinib", "NA",       "gemcitabine-50mpk", "NA",      "NA",    "paclitaxel","NA",       "NA",       "tamoxiofen", "trametinib"),
                c("dna_rep","bortezomib", "rtk", "cytoskeleton","dna_rep","egfr","pi3k","dna_rep","egfr","hdac","cytoskeleton","mtor","hdac","tamoxifen","trametinib"))
colnames(names) <- c("gdsc", "ctrp", "gcsi","nibr", "class")
names <- as.data.frame(names, stringsAsFactors=F)

rel_drugs_gdsc <- list("hdac"=drugInfo(GDSC)$DRUG.NAME[which(drugInfo(GDSC)$TARGET.PATHWAY == "chromain  histone acetylation")],
                       "dna_rep"=drugInfo(GDSC)$DRUG.NAME[which(drugInfo(GDSC)$TARGET.PATHWAY == "DNA replication")],
                       "rtk"=drugInfo(GDSC)$DRUG.NAME[which(drugInfo(GDSC)$TARGET.PATHWAY == "RTK signaling")],
                       "egfr"=drugInfo(GDSC)$DRUG.NAME[which(drugInfo(GDSC)$TARGET.PATHWAY == "EGFR signaling")],
                       "bortezomib"="MG-132",
                       "mtor"=drugInfo(GDSC)$DRUG.NAME[which(drugInfo(GDSC)$TARGET.PATHWAY == "TOR signaling")],
                       "cytoskeleton"=drugInfo(GDSC)$DRUG.NAME[which(drugInfo(GDSC)$TARGET.PATHWAY == "cytoskeleton")],
                       "pi3k"=drugInfo(GDSC)$DRUG.NAME[which(drugInfo(GDSC)$TARGET.PATHWAY == "PI3K signaling")],
                       "tamoxifen"=NULL,
                       "trametinib"=c("VX-11e",  "TL-2-105", "FMK", "FR-180204", "HG-6-64-1", "CMK", "AZ628", "SL 0101-1", "Trametinib", "RDEA119", "PLX4720", "selumetinib", "CI-1040",  
                                      "dabrafenib", "PD-0325901","SB590885"))

rel_drugs_ctrp <- list("hdac"=drugInfo(CTRP)$cpd_name[grep("HDAC", drugInfo(CTRP)$gene_symbol_of_protein_target, value=F)],
                       "dna_rep"=c("fluorouracil","tanespimycin:gemcitabine (1:1 mol/mol)", "navitoclax:gemcitabine (1:1 mol/mol)", "etoposide", "mitomycin", "doxorubicin", "doxorubicin:navitoclax (2:1 mol/mol)",
                                   "temozolomide", "methotrexate", "SN-38", "carboplatin", "procarbazine", "cytarabine hydrochloride", "austocystin D", "bleomycin A2", 
                                   "vorinostat:carboplatin (1:1 mol/mol)", "tretinoin:carboplatin (2:1 mol/mol)", "carboplatin:etoposide (40:17 mol/mol)", "JQ-1:carboplatin (1:1 mol/mol)",
                                   "BRD-A02303741:carboplatin (1:1 mol/mol)", "clofarabine",  "carboplatin:UNC0638 (2:1 mol/mol)", "decitabine:carboplatin (1:1 mol/mol)", 
                                   "teniposide", "topotecan", "nelarabine", "gemcitabine"),
                       "rtk"=c("NVP-TAE684", "SU11274", "SGX-523", "quizartinib", "NVP-BSK805", "TG-101348", "tivantinib", "crizotinib:PLX-4032 (2:1 mol/mol)", "crizotinib", "foretinib", "MGCD-265", "SGX-523", "cabozantinib","tivantinib"),
                       "egfr"=c("neratinib", "afatinib", "canertinib", "lapatinib", "WZ8040", "erlotinib:PLX-4032 (2:1 mol/mol)", "gefitinib", "PD 153035", "WZ4002"),
                       "bortezomib"=c("MLN2238","sirolimus:bortezomib (250:1 mol/mol)", "ISOX:bortezomib (250:1 mol/mol)", "tanespimycin:bortezomib (250:1 mol/mol)","SNX-2112:bortezomib (250:1 mol/mol)"),
                       "mtor"=drugInfo(CTRP)$cpd_name[grep("MTOR", drugInfo(CTRP)$gene_symbol_of_protein_target, value=F)],
                       "cytoskeleton"=c("paclitaxel","parbendazole", "docetaxel", "nakiterpiosin", "darinaparsin", "tivantinib", "docetaxel:tanespimycin (2:1 mol/mol)", "PF-3758309", "tubastatin A", "CHM-1", "NSC23766"),
                       "pi3k"=drugInfo(CTRP)$cpd_name[grep("PI3K", drugInfo(CTRP)$target_or_activity_of_compound, value=F)],
                       "tamoxifen"=c("fulvestrant"),
                       "trametinib"=c("selumetinib:PLX-4032 (8:1 mol/mol)", "selumetinib:GDC-0941 (4:1 mol/mol)", "selumetinib:tretinoin (2:1 mol/mol)",
                                      "selumetinib:vorinostat (8:1 mol/mol)", "selumetinib:BRD-A02303741 (4:1 mol/mol)", "selumetinib:MK-2206 (8:1 mol/mol)",
                                      "selumetinib:piperlongumine (8:1 mol/mol)", "selumetinib:navitoclax (8:1 mol/mol)", "selumetinib:decitabine (4:1 mol/mol)",
                                      "selumetinib:UNC0638 (4:1 mol/mol)", "PD318088", "selumetinib:JQ-1 (4:1 mol/mol)", "PLX-4720", "selumetinib", "dabrafenib" ))

rel_drugs_gcsi <- list("hdac"=c("MS-275", "Vorinostat"),
                       "dna_rep"=c("Doxorubicin", "Gemcitabine"),
                       "rtk"=NULL,
                       "egfr"=c("Erlotinib", "Lapatinib"),
                       "bortezomib"=NULL,
                       "mtor"=NULL,
                       "cytoskeleton"=c("Docetaxel", "paclitaxel"),
                       "pi3k"=NULL)

rel_drugs_nibr <- list("5fu"=c("abraxane + gemcitabine", "gemcitabine-50mpk"),
                       "egfr"=c("cetuximab", "cetuximab + encorafenib", "BYL719 + cetuximab" , "BYL719 + cetuximab + encorafenib"),
                       "dna_rep"=c("5FU", "abraxane + gemcitabine"),
                       "cytoskeleton"=c("LCL161 + paclitaxel"),
                       "tamoxifen"=NULL,
                       "trametinib"=NULL)

rel_drugs <- list(gdsc=rel_drugs_gdsc, ctrp=rel_drugs_ctrp, gcsi=rel_drugs_gcsi, nibr=rel_drugs_nibr)

#function for generating the list of unrelated drugs (based on drug resp. correlation and exclusion of known related drugs)
unrel_drugs <- function(drug, dataset_name) {
  if (dataset_name=="gdsc") {data <- gdsc_table_imp[,(gen_last_gdsc+1):ncol(gdsc_table_imp)]}
  if (dataset_name=="ctrp") {data <- ctrp_table_imp[,(gen_last_ctrp+1):ncol(ctrp_table_imp)]}
  if (dataset_name=="gcsi") {data <- gcsi_table_imp[,(gen_last_gcsi+1):ncol(gcsi_table_imp)]}
  if (dataset_name=="nibr") {data <- nibr_table_imp[,(gen_last_nibr+1):ncol(nibr_table_imp)]}
  print(drug)
  cor <- apply(data, 2, function(x) {
    cor(data[,drug],x, use="na.or.complete", method="spearman")
  })
  unrel_drugs <- colnames(data)[which(cor<0.15)]
  rel_drugs1 <- (rel_drugs[[dataset_name]])[[names$class[which(names[,dataset_name]==drug)]]]
  unrel_drugs <- setdiff(unrel_drugs, rel_drugs1)
  
  return(unrel_drugs)
}

# testing
unrel_drugs(names[12,3], "gcsi")

# main analysis function !! check number of prin. comp.
check_all_biomarkers <- function(drug, dataset_name, feature_overlap) {
  n_pc=10
  if (dataset_name=="gdsc") {data <- gdsc_table_imp}
  if (dataset_name=="ctrp") {data <- ctrp_table_imp}
  if (dataset_name=="nibr") {data <- nibr_table_imp}
  if (dataset_name=="gcsi") {data <- gcsi_table_imp
  n_pc=1
  }
  a <- scale(data[,drug], center=F)
  unrel_drugs <- unrel_drugs(drug, dataset_name)
  controlPCsAll <- prcomp(data[, unrel_drugs])$x
  
  res <- sapply(feature_overlap, function(x){
    
    b <- scale(data[,x], center=F)
    model_1 <- try(lm(a ~ b))
    if(class(model_1)!="try-error") {
      model_2 <- lm(a ~ b + scale(controlPCsAll[,1:n_pc], center=F))
      res <- signif(c(model_1$coefficients[2], model_2$coefficients[2], anova(model_1)$"Pr(>F)"[1], anova(model_2)$"Pr(>F)"[1]))
    } else {res <- rep(NA, 4)}
    return(res)
  })
  
  res <- t(res)
  colnames(res) <- c("reg_coef", "reg_coef_new", "p","p_new")
  return(res)
}

# checking biomarker's overlap
check_overlap <- function(res_table, n_top) {
  mode(res_table) <- "numeric"
  a1 <- order(abs(res_table[,2]), decreasing = T)[1:n_top]
  b1 <- order(abs(res_table[,6]), decreasing = T)[1:n_top]
  c1 <- order(abs(res_table[,10]), decreasing = T)[1:n_top]
  
  a2 <- order(abs(res_table[,3]), decreasing = T)[1:n_top]
  b2 <- order(abs(res_table[,7]), decreasing = T)[1:n_top]
  c2 <- order(abs(res_table[,11]), decreasing = T)[1:n_top]
  
  int1 <- length(intersect(intersect(a1,b1), c1))
  int2 <- length(intersect(intersect(a2,b2), c2))
  
  top_f1 <- unique(c(a1,b1,c1))
  top_f2 <- unique(c(a2,b2,c2))
  
  res_table1 <- res_table[top_f1,]
  res_table2 <- res_table[top_f2,]
  cor_ab1 <- cor(res_table1[,2], res_table1[,6], use="na.or.complete")
  cor_ac1 <- cor(res_table1[,2], res_table1[,10], use="na.or.complete")
  cor_bc1 <- cor(res_table1[,6], res_table1[,10], use="na.or.complete")
  
  cor_ab2 <- cor(res_table2[,3], res_table2[,7], use="na.or.complete")
  cor_ac2 <- cor(res_table2[,3], res_table2[,11], use="na.or.complete")
  cor_bc2 <- cor(res_table2[,7], res_table2[,11], use="na.or.complete")
  
  inf <- rbind(signif(c(int1, cor_ab1, cor_ac1, cor_bc1)), signif(c(int2, cor_ab2, cor_ac2, cor_bc2)))
  rownames(inf) <- c("before", "after")
  colnames(inf) <- c(paste0("overlap within top ",n_top," predictors"), "cor GDSC-CTRP", "cor GDSC-NIBR", "cor CTRP-NIBR")
  return(inf)
}

## 1) GDSC-CTRP-NIBR analysis

# features_overlap
fov1 <- intersect(colnames(gdsc_table_imp)[1:gen_last_gdsc], colnames(ctrp_table_imp)[1:gen_last_ctrp])
fov2 <- intersect(fov1, colnames(table_bestresp_nibr)[1:gen_last_nibr])
# overlap: 15749 features: 15607 expr, 142 mut

res_fluorouracil <- cbind(fov2, check_all_biomarkers(names[1,1], "gdsc", fov2), check_all_biomarkers(names[1,2], "ctrp", fov2), check_all_biomarkers(names[1,4], "nibr", fov2))
res_erlotinib <- cbind(fov2, check_all_biomarkers(names[6,1], "gdsc", fov2), check_all_biomarkers(names[6,2], "ctrp", fov2), check_all_biomarkers(names[6,4], "nibr", fov2))
res_gemcitabine <- cbind(fov2, check_all_biomarkers(names[8,1], "gdsc", fov2), check_all_biomarkers(names[8,2], "ctrp", fov2), check_all_biomarkers(names[8,4], "nibr", fov2))
res_paclitaxel <- cbind(fov2, check_all_biomarkers(names[11,1], "gdsc", fov2), check_all_biomarkers(names[11,2], "ctrp", fov2), check_all_biomarkers(names[11,4], "nibr", fov2))
res_tamoxifen <- cbind(fov2, check_all_biomarkers(names[14,1], "gdsc", fov2), check_all_biomarkers(names[14,2], "ctrp", fov2), check_all_biomarkers(names[14,4], "nibr", fov2))
res_trametinib <- cbind(fov2, check_all_biomarkers(names[15,1], "gdsc", fov2), check_all_biomarkers(names[15,2], "ctrp", fov2), check_all_biomarkers(names[15,4], "nibr", fov2))



res_check1 <- rbind(check_overlap(res_fluorouracil, 100),
                    check_overlap(res_erlotinib, 100),
                    check_overlap(res_gemcitabine, 100),
                    check_overlap(res_paclitaxel, 100),
                    check_overlap(res_tamoxifen, 100),
                    check_overlap(res_trametinib, 100))

res_check2 <- rbind(check_overlap(res_fluorouracil, 500),
                    check_overlap(res_erlotinib, 500),
                    check_overlap(res_gemcitabine, 500),
                    check_overlap(res_paclitaxel, 500),
                    check_overlap(res_tamoxifen, 500),
                    check_overlap(res_trametinib, 500))
res_check <-cbind(res_check1, res_check2)

# saving
res_6drugs_glds <- list(fluorouracil=res_fluorouracil, erlotinib=res_erlotinib, gemcitabine=res_gemcitabine, 
                        paclitaxel=res_paclitaxel, tamoxifen=res_tamoxifen, trametinib=res_trametinib)
save(res_6drugs_glds, file="res_6drugs_glds.RData")
write.table(res_check, file="res_check_6drugs_glds.txt")


## 2) GDSC-CTRP-gCSI analysis

# features_overlap
fov1 <- intersect(colnames(gdsc_table_imp)[1:gen_last_gdsc], colnames(ctrp_table_imp)[1:gen_last_ctrp])
fov_gcsi <- intersect(fov1, colnames(gcsi_table_imp)[1:gen_last_gcsi])
# overlap: 15595 features: all exp.

for (i in 2:13) {
  res <- cbind(fov_gcsi, check_all_biomarkers(names[i,1], "gdsc",fov_gcsi), check_all_biomarkers(names[i,2], "ctrp",fov_gcsi), check_all_biomarkers(names[i,3], "gcsi",fov_gcsi))
  assign(paste0("res_", names[i,2]), res)
}

res_check_gcsi_100 <- rbind(check_overlap(res_bortezomib, 100),
                    check_overlap(res_crizotinib, 100),
                    check_overlap(res_docetaxel, 100),
                    check_overlap(res_doxorubicin, 100),
                    check_overlap(res_erlotinib, 100),
                    check_overlap(`res_GDC-0941`, 100),
                    check_overlap(res_gemcitabine, 100),
                    check_overlap(res_lapatinib, 100),
                    check_overlap(res_entinostat, 100),
                    check_overlap(res_paclitaxel, 100),
                    check_overlap(res_sirolimus, 100),
                    check_overlap(res_vorinostat, 100))

res_check_gcsi_100 <- rbind(check_overlap(res_bortezomib, 500),
                            check_overlap(res_crizotinib, 500),
                            check_overlap(res_docetaxel, 500),
                            check_overlap(res_doxorubicin, 500),
                            check_overlap(res_erlotinib, 500),
                            check_overlap(`res_GDC-0941`, 500),
                            check_overlap(res_gemcitabine, 500),
                            check_overlap(res_lapatinib, 500),
                            check_overlap(res_entinostat, 500),
                            check_overlap(res_paclitaxel, 500),
                            check_overlap(res_sirolimus, 500),
                            check_overlap(res_vorinostat, 500))

res_check_gcsi <-cbind(res_check_gcsi_100, res_check_gcsi_500)

# saving
res_12drugs_glds <- list(bortezomib=res_bortezomib, crizotinib=res_crizotinib, docetaxel=res_docetaxel, 
                        doxorubicin=res_doxorubicin, erlotinib=res_erlotinib, "GDC-0941"=`res_GDC-0941`,
                        gemcitabine=res_gemcitabine, lapatinib=res_lapatinib, entinostat=res_entinostat,
                        paclitaxel=res_paclitaxel, sirolimus=res_sirolimus, vorinostat=res_vorinostat)
save(res_12drugs_glds, file="res_12drugs_glds.RData")
write.table(res_check_gcsi, file="res_check_12drugs_glds_1pc.txt")

