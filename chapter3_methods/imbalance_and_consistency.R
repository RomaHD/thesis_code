# written in March 2017, modified on 28.10.18 
# Chapter 3 of the thesis. 
# Methods for improving accuracy of drug response prediction: 
# 5. Testing influence of class imbalance and cross-set inconsistency on accuracy
#
# datasets:CTRP, GDSC

path <- "/data/kurilov/genestack/phd/work_2018/thesis_code/"
setwd(file.path(path, "chapter3_methods"))

# load CTRP, GDSC (imputed drug resp data because of glmnet can't handle missing data)
load(file.path(path, "data/ctrp_table_imp.RData"))
load(file.path(path, "data/gdsc_table_imp.RData"))

# getting single-task R2 results from single- vs.multitask glmnet tests 
# (produced by the script multi_task.R)
res <- read.table("multitask_glmnet_res.txt")
res_lasso <- res[c(10,12),seq(2,38,by=2)]
colnames(res_lasso) <-  c("entinostat","tubastatin A", "belinostat","vorinostat", "lapatinib", "erlotinib", "afatinib", "gefitinib",
                          "trametinib","selumetinib", "PLX-4720","dabrafenib", "SNX-2112","tanespimycin", "OSI-027", 
                          "sirolimus", "NVP-BEZ235", "temsirolimus", "AZD8055")
rownames(res_lasso) <- c("GDSC","CTRP")


names <- cbind(c("MS-275","Tubastatin A", "Belinostat", "Vorinostat", "Lapatinib", "Erlotinib", "Afatinib", "Gefitinib",
                 "Trametinib", "selumetinib", "PLX4720", "Dabrafenib", "SNX-2112", "17-AAG", "OSI-027", 
                 "Rapamycin", "BEZ235", "Temsirolimus", "AZD8055"),
               c("entinostat","tubastatin A", "belinostat","vorinostat", "lapatinib", "erlotinib", "afatinib", "gefitinib",
                 "trametinib","selumetinib", "PLX-4720","dabrafenib", "SNX-2112","tanespimycin", "OSI-027", 
                 "sirolimus", "NVP-BEZ235", "temsirolimus", "AZD8055"))
colnames(names) <- c("gdsc", "ctrp")

# calculating imbalance
imb <- function(x) {
  sens <- length(which(x<0.5))
  resist  <- length(which(x>=0.5))
  value <- abs(0.5-(sens/(sens+resist)))
  return(value)
}

imb_gdsc <- sapply(names[,1], function(x) {
  imb(gdsc_table_imp[,x])
})
imb_ctrp <- sapply(names[,2], function(x) {
  imb(ctrp_table_imp[,x])
})

# calculating drug response correlation
com_lines <- intersect(rownames(gdsc_table_imp),rownames(ctrp_table_imp))
cor_auc <- apply(names,1, function(x) {
  cor(gdsc_table_imp[com_lines,x[1]], ctrp_table_imp[com_lines,x[2]])
})

# testing correlation between [accuracy and imbalance] and [accuracy and drug. resp. consistency]
cor(unlist(res_lasso[1,]), as.vector(imb_gdsc), method="spearman")
cor(unlist(res_lasso[2,]), as.vector(imb_ctrp), method="spearman")
cor(unlist(res_lasso[1,]), cor_auc, method="spearman")
cor(unlist(res_lasso[2,]), cor_auc, method="spearman")