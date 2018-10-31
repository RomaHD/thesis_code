# written in March 2017, modified on 28.10.18 
# Chapter 3 of the thesis. 
# Methods for improving accuracy of drug response prediction: 
# 4. Modelling with sample weights (for imbalanced data)
#
# datasets:CTRP, GDSC

library(glmnet)
library(survcomp)

path <- "/data/kurilov/genestack/phd/work_2018/thesis_code/"
setwd(file.path(path, "chapter3_methods"))

# load CTRP, GDSC (imputed drug resp data because of glmnet can't handle missing data)
load(file.path(path, "data/ctrp_table_imp.RData"))
load(file.path(path, "data/gdsc_table_imp.RData"))

### 1. Annotation and functions ###

names <- cbind(c("MS-275","Tubastatin A", "Belinostat", "Vorinostat", "Lapatinib", "Erlotinib", "Afatinib", "Gefitinib",
                 "Trametinib", "selumetinib", "PLX4720", "Dabrafenib", "SNX-2112", "17-AAG", "OSI-027", 
                 "Rapamycin", "BEZ235", "Temsirolimus", "AZD8055"),
               c("entinostat","tubastatin A", "belinostat","vorinostat", "lapatinib", "erlotinib", "afatinib", "gefitinib",
                 "trametinib","selumetinib", "PLX-4720","dabrafenib", "SNX-2112","tanespimycin", "OSI-027", 
                 "sirolimus", "NVP-BEZ235", "temsirolimus", "AZD8055"))
colnames(names) <- c("gdsc", "ctrp")

gen_last_gdsc <- max(grep("_", colnames(gdsc_table_imp)))
gen_last_ctrp <- max(grep("_", colnames(ctrp_table_imp)))
f_set <- intersect(colnames(gdsc_table_imp)[1:gen_last_gdsc], colnames(ctrp_table_imp)[1:gen_last_ctrp])
# 16126 exp + 24960 cn+ 146 mut 

table_gdsc <- cbind(gdsc_table_imp[,f_set], gdsc_table_imp[,(1+gen_last_gdsc):ncol(gdsc_table_imp)])
table_ctrp <- cbind(ctrp_table_imp[,f_set], ctrp_table_imp[,(1+gen_last_ctrp):ncol(ctrp_table_imp)])

gen_last_gdsc <- max(grep("_", colnames(table_gdsc)))
gen_last_ctrp <- max(grep("_", colnames(table_ctrp)))

imp_values_gdsc <- apply(table_gdsc[,1:gen_last_gdsc],2, function(x){
  na_n <- which(is.na(x))
  x[na_n] <- median(x, na.rm=T)
  x
})

table_gdsc <- cbind(imp_values_gdsc, table_gdsc[,(gen_last_gdsc+1):ncol(table_gdsc)])

imp_values_ctrp <- apply(table_ctrp[,1:gen_last_ctrp],2, function(x){
  na_n <- which(is.na(x))
  x[na_n] <- median(x, na.rm=T)
  x
})

table_ctrp <- cbind(imp_values_ctrp, table_ctrp[,(gen_last_ctrp+1):ncol(table_ctrp)])


r2 <- function(v1, v2) {
  value <- cor(v1, v2, use="na.or.complete")
  return(value^2)
}

# concordance.index 
ci <- function(v1, v2) {
  value <- concordance.index(x=v1, surv.time=v2, surv.event=rep(1,length(v1)))
  return(value$c.index)
}


require(doMC)
registerDoMC(cores=2)

### 2. Tests ###

r2_gdsc <- matrix(NA, nrow=8, ncol=19)
r2_ctrp <- matrix(NA, nrow=8, ncol=19)
ci_gdsc <- matrix(NA, nrow=8, ncol=19)
ci_ctrp <- matrix(NA, nrow=8, ncol=19)

# weighting schemes:
## for all test samples
# 1: no wieghting
# 2: w=1/AUC
# 3: w=1, AUC>0.5 and w=2, AUC<0.5
# 4: w=1, AUC>0.5 and w=10, AUC<0.5
## for test samples with AUC<0.5
# 5: no wieghting
# 6: w=1/AUC
# 7: w=1, AUC>0.5 and w=2, AUC<0.5
# 8: w=1, AUC>0.5 and w=10, AUC<0.5


# j==1 GDSC -> CTRP
# j==2 CTRP -> GDSC
for (j in 1:2) {
  
  if (j==1) { table <- table_gdsc
  table_test <- table_ctrp }
  if (j==2) { table <- table_ctrp
  table_test <- table_gdsc }  
  
  for (i in 1:19) {
    print(i)
    
    if (j==1) {drug <- names[i,1]
    drug_test <- names[i,2] }
    if (j==2) {drug <- names[i,2]
    drug_test <- names[i,1] } 
    
    #weights
    weights2 <- 1/table[, drug]
    weights34 <- t(sapply(table[, drug], function(x) {
      w3=1
      w4=1
      if(x<0.5) {
        w3=2
        w4=10}
      c(w3,w4)
    }))
    weights3 <- weights34[,1]
    weights4 <- weights34[,2]
    
    
    m_1 = cv.glmnet(x=table[,f_set], y=table[, drug], family = "gaussian")
    m_2 = cv.glmnet(x=table[,f_set], y=table[, drug], family = "gaussian", weights=weights2)
    m_3 = cv.glmnet(x=table[,f_set], y=table[, drug], family = "gaussian", weights=weights3)
    m_4 = cv.glmnet(x=table[,f_set], y=table[, drug], family = "gaussian", weights=weights4)
    
    pred_1 <- predict.cv.glmnet(m_1, newx=table_test[,f_set], s="lambda.min")
    pred_2 <- predict.cv.glmnet(m_2, newx=table_test[,f_set], s="lambda.min")
    pred_3 <- predict.cv.glmnet(m_3, newx=table_test[,f_set], s="lambda.min")
    pred_4 <- predict.cv.glmnet(m_4, newx=table_test[,f_set], s="lambda.min")
    
    sens_samples <- which(table_test[,drug_test]<median(table_test[,drug_test]))
    pred_5 <- predict.cv.glmnet(m_1, newx=table_test[sens_samples,f_set], s="lambda.min")
    pred_6 <- predict.cv.glmnet(m_2, newx=table_test[sens_samples,f_set], s="lambda.min")
    pred_7 <- predict.cv.glmnet(m_3, newx=table_test[sens_samples,f_set], s="lambda.min")
    pred_8 <- predict.cv.glmnet(m_4, newx=table_test[sens_samples,f_set], s="lambda.min")
    
    
    # R2
    r2_1 <- r2(pred_1, table_test[,drug_test])
    r2_2 <- r2(pred_2, table_test[,drug_test])
    r2_3 <- r2(pred_3, table_test[,drug_test])
    r2_4 <- r2(pred_4, table_test[,drug_test])
    r2_5 <- r2(pred_5, table_test[sens_samples,drug_test])
    r2_6 <- r2(pred_6, table_test[sens_samples,drug_test])
    r2_7 <- r2(pred_7, table_test[sens_samples,drug_test])
    r2_8 <- r2(pred_8, table_test[sens_samples,drug_test])
    
    # ci
    ci_1 <- ci(pred_1, table_test[,drug_test])
    ci_2 <- ci(pred_2, table_test[,drug_test])
    ci_3 <- ci(pred_3, table_test[,drug_test])
    ci_4 <- ci(pred_4, table_test[,drug_test])
    ci_5 <- ci(pred_5, table_test[sens_samples,drug_test])
    ci_6 <- ci(pred_6, table_test[sens_samples,drug_test])
    ci_7 <- ci(pred_7, table_test[sens_samples,drug_test])
    ci_8 <- ci(pred_8, table_test[sens_samples,drug_test])
    
    
    #saving res. 
    if (j==1) {
      r2_gdsc[,i] <- c(r2_1, r2_2, r2_3, r2_4, r2_5, r2_6, r2_7, r2_8)
      ci_gdsc[,i] <- c(ci_1, ci_2, ci_3, ci_4, ci_5, ci_6, ci_7, ci_8) }
    
    if (j==2) {
      r2_ctrp[,i] <- c(r2_1, r2_2, r2_3, r2_4, r2_5, r2_6, r2_7, r2_8)
      ci_ctrp[,i] <- c(ci_1, ci_2, ci_3, ci_4, ci_5, ci_6, ci_7, ci_8) }
  } 
}

write.table(r2_gdsc, file="glmnet_weights_r2_gdsc.txt")
write.table(r2_ctrp, file="glmnet_weights_r2_ctrp.txt")
write.table(ci_gdsc, file="glmnet_weights_ci_gdsc.txt")
write.table(ci_ctrp, file="glmnet_weights_ci_ctrp.txt")