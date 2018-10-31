# written in March 2017, modified on 28.10.18 
# Chapter 3 of the thesis. 
# Methods for improving accuracy of drug response prediction: 
# 3. Modelling using feature interactions: glmnet tests, RF tests, plotting glmnet results
#
# datasets:CTRP, GDSC

library(glmnet)
library(survcomp)
library(caret)
library(ggplot2)

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
                 "sirolimus", "NVP-BEZ235", "temsirolimus", "AZD8055"),
               c(rep("hdac", 4), rep("egfr", 4), rep("mek", 2), rep("braf", 2), rep("hsp90", 2), rep("mtor", 5)))
colnames(names) <- c("gdsc", "ctrp", "class")

gen_last_gdsc <- max(grep("_", colnames(gdsc_table_imp)))
gen_last_ctrp <- max(grep("_", colnames(ctrp_table_imp)))
f_set <- intersect(colnames(gdsc_table_imp)[1:gen_last_gdsc], colnames(ctrp_table_imp)[1:gen_last_ctrp])
# 16126 exp + 24960 cn + 146 mut 

table_gdsc <- cbind(gdsc_table_imp[,f_set], gdsc_table_imp[,gen_last_gdsc:ncol(gdsc_table_imp)])
table_ctrp <- cbind(ctrp_table_imp[,f_set], ctrp_table_imp[,gen_last_ctrp:ncol(ctrp_table_imp)])

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

## feature engineering functions

# function for calc. gene multiplication features
mult_features <- function(table) {
  
  table_names <- matrix(NA, nrow=ncol(table), ncol=ncol(table)) 
  for (i in 1:ncol(table)){
    for (j in 1:ncol(table)){
      table_names[i,j] <- paste0(colnames(table)[i],"_",colnames(table)[j])
    }
  }
  
  result <- t(apply(table, 1,function(x){
    res <- sapply(x, function(z){
      x*z
    })
    res[lower.tri(res, diag = FALSE)]
  }))
  colnames(result) <- table_names[lower.tri(table_names, diag = FALSE)]
  
  return(result)
}

# function for calc. binary gene pairs features

bgp_features <- function(table) {
  
  table_names <- matrix(NA, nrow=ncol(table), ncol=ncol(table)) 
  for (i in 1:ncol(table)){
    for (j in 1:ncol(table)){
      table_names[i,j] <- paste0(colnames(table)[i],"_",colnames(table)[j])
    }
  }
  
  result <- t(apply(table, 1,function(x){
    res <- matrix(0, nrow=length(x), ncol=length(x))
    for (i in 1:length(x)){
      for (j in 1:length(x)){
        if (x[i]>x[j]) {res[i,j]<-1}
      }
    }
    res[lower.tri(res, diag = FALSE)]
  }))
  colnames(result) <- table_names[lower.tri(table_names, diag = FALSE)]
  return(result) 
}

# caret regression
simple_regression <- function(y, x, method="rf") {
  ctrl <- trainControl(method = "cv", number = 5, search = "random")
  set.seed(100)
  
  print(dim(x))
  Tune <- train(x, y,
                method = method,
                tuneLength = 15,
                trControl = ctrl,
                preProc = c("center","scale","medianImpute"))
  
  return(Tune)
}

r2 <- function(v1, v2) {
  value <- cor(v1, v2, use="na.or.complete")
  return(value^2)
}


# concordance.index 
ci <- function(v1, v2) {
  value <- concordance.index(x=v1, surv.time=v2, surv.event=rep(1,length(v1)))
  return(value$c.index)
}


### 2. Tests with glmnet ###

require(doMC)
registerDoMC(cores=2)

r2_gdsc <- matrix(NA, nrow=11, ncol=19)
r2_ctrp <- matrix(NA, nrow=11, ncol=19)
ci_gdsc <- matrix(NA, nrow=11, ncol=19)
ci_ctrp <- matrix(NA, nrow=11, ncol=19)

# models:
# 1: all features
# 2: 200 top correlation
# 3: 200 top variance
# 4: BGP (200 top cor) all
# 5: BGP (200 top cor) 200 top corr
# 6: Multipl (200 top cor) all
# 7: Multipl (200 top cor) 200 top corr
# 8: BGP (200 top var) all
# 9: BGP (200 top var) 200 top corr
# 10: Multipl (200 top var) all
# 11: Multipl (200 top var) 200 top corr

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
    
    m_1 = cv.glmnet(x=table[,f_set], y=table[, drug], family = "gaussian")
    
    cor_data <- apply(table[,f_set], 2, function(z) cor.test(z, table[, drug], method="spearman", use="na.or.complete")$p.value)
    features_2 <- names(sort(cor_data))[1:200]
    m_2 = cv.glmnet(x=table[,features_2], y=table[, drug], family = "gaussian")
    
    sd_data <- apply(table[,f_set], 2, function(z) sd(z, na.rm=T))
    features_3 <- names(sort(sd_data, decreasing = T))[1:200]
    m_3 = cv.glmnet(x=table[,features_3], y=table[, drug], family = "gaussian")
    
    bgp200cor <- bgp_features(table[,features_2])
    m_4 = cv.glmnet(x=bgp200cor, y=table[, drug], family = "gaussian")
    
    cor_data <- apply(bgp200cor, 2, function(z) cor.test(z, table[, drug], method="spearman", use="na.or.complete")$p.value)
    features_5 <- names(sort(cor_data))[1:200]
    m_5 = cv.glmnet(x=bgp200cor[,features_5], y=table[, drug], family = "gaussian")
    
    mult200cor <- mult_features(table[,features_2])
    m_6 = cv.glmnet(x=mult200cor, y=table[, drug], family = "gaussian")
    
    cor_data <- apply(mult200cor, 2, function(z) cor.test(z, table[, drug], method="spearman", use="na.or.complete")$p.value)
    features_7 <- names(sort(cor_data))[1:200]
    m_7 = cv.glmnet(x=mult200cor[,features_7], y=table[, drug], family = "gaussian")
    
    bgp200var <- bgp_features(table[,features_3])
    m_8 = cv.glmnet(x=bgp200var, y=table[, drug], family = "gaussian")
    
    cor_data <- apply(bgp200var, 2, function(z) cor.test(z, table[, drug], method="spearman", use="na.or.complete")$p.value)
    features_9 <- names(sort(cor_data))[1:200]
    m_9 = cv.glmnet(x=bgp200var[,features_9], y=table[, drug], family = "gaussian")
    
    mult200var <- mult_features(table[,features_3])
    m_10 = cv.glmnet(x=mult200var, y=table[, drug], family = "gaussian")
    
    cor_data <- apply(mult200var, 2, function(z) cor.test(z, table[, drug], method="spearman", use="na.or.complete")$p.value)
    features_11 <- names(sort(cor_data))[1:200]
    m_11 = cv.glmnet(x=mult200var[,features_11], y=table[, drug], family = "gaussian")
    
    
    pred_1 <- predict.cv.glmnet(m_1, newx=table_test[,f_set], s="lambda.min")
    pred_2 <- predict.cv.glmnet(m_2, newx=table_test[,features_2], s="lambda.min")
    pred_3 <- predict.cv.glmnet(m_3, newx=table_test[,features_3], s="lambda.min")
    
    bgp200cor_test <- bgp_features(table_test[,features_2])
    pred_4 <- predict.cv.glmnet(m_4, newx=bgp200cor_test, s="lambda.min")
    pred_5 <- predict.cv.glmnet(m_5, newx=bgp200cor_test[,features_5], s="lambda.min")
    
    mult200cor_test <- mult_features(table_test[,features_2])
    pred_6 <- predict.cv.glmnet(m_6, newx=mult200cor_test, s="lambda.min")
    pred_7 <- predict.cv.glmnet(m_7, newx=mult200cor_test[,features_7], s="lambda.min")
    
    bgp200var_test <- bgp_features(table_test[,features_3])
    pred_8 <- predict.cv.glmnet(m_8, newx=bgp200var_test, s="lambda.min")
    pred_9 <- predict.cv.glmnet(m_9, newx=bgp200var_test[,features_9], s="lambda.min")
    
    mult200var_test <- mult_features(table_test[,features_3])
    pred_10 <- predict.cv.glmnet(m_10, newx=mult200var_test, s="lambda.min")
    pred_11 <- predict.cv.glmnet(m_11, newx=mult200var_test[,features_11], s="lambda.min")
    
    # R2
    r2_1 <- r2(pred_1, table_test[,drug_test])
    r2_2 <- r2(pred_2, table_test[,drug_test])
    r2_3 <- r2(pred_3, table_test[,drug_test])
    r2_4 <- r2(pred_4, table_test[,drug_test])
    r2_5 <- r2(pred_5, table_test[,drug_test])
    r2_6 <- r2(pred_6, table_test[,drug_test])
    r2_7 <- r2(pred_7, table_test[,drug_test])
    r2_8 <- r2(pred_8, table_test[,drug_test])
    r2_9 <- r2(pred_9, table_test[,drug_test])
    r2_10 <- r2(pred_10, table_test[,drug_test])
    r2_11 <- r2(pred_11, table_test[,drug_test])
    
    # ci
    ci_1 <- ci(pred_1, table_test[,drug_test])
    ci_2 <- ci(pred_2, table_test[,drug_test])
    ci_3 <- ci(pred_3, table_test[,drug_test])
    ci_4 <- ci(pred_4, table_test[,drug_test])
    ci_5 <- ci(pred_5, table_test[,drug_test])
    ci_6 <- ci(pred_6, table_test[,drug_test])
    ci_7 <- ci(pred_7, table_test[,drug_test])
    ci_8 <- ci(pred_8, table_test[,drug_test])
    ci_9 <- ci(pred_9, table_test[,drug_test])
    ci_10 <- ci(pred_10, table_test[,drug_test])
    ci_11 <- ci(pred_11, table_test[,drug_test])
    
    
    #saving res. 
    if (j==1) {
      r2_gdsc[,i] <- c(r2_1, r2_2, r2_3, r2_4, r2_5, r2_6, r2_7, r2_8, r2_9, r2_10, r2_11)
      ci_gdsc[,i] <- c(ci_1, ci_2, ci_3, ci_4, ci_5, ci_6, ci_7, ci_8, ci_9, ci_10, ci_11) }
    
    if (j==2) {
      r2_ctrp[,i] <- c(r2_1, r2_2, r2_3, r2_4, r2_5, r2_6, r2_7, r2_8, r2_9, r2_10, r2_11)
      ci_ctrp[,i] <- c(ci_1, ci_2, ci_3, ci_4, ci_5, ci_6, ci_7, ci_8, ci_9, ci_10, ci_11) }
  } 
}

write.table(r2_gdsc, file="feature_interaction_r2_gdsc.txt")
write.table(r2_ctrp, file="feature_interaction_r2_ctrp.txt")
write.table(ci_gdsc, file="feature_interaction_ci_gdsc.txt")
write.table(ci_ctrp, file="feature_interaction_ci_ctrp.txt")

### 2. Tests with Random Forest ###

require(doMC)
registerDoMC(cores=2)

r2_gdsc <- matrix(NA, nrow=4, ncol=19)
r2_ctrp <- matrix(NA, nrow=4, ncol=19)
ci_gdsc <- matrix(NA, nrow=4, ncol=19)
ci_ctrp <- matrix(NA, nrow=4, ncol=19)

# models:
# 2: 200 top correlation
# 3: 200 top variance
# 5: BGP (200 top cor) 200 top corr
# 7: Multipl (200 top cor) 200 top corr


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
    
    cor_data <- apply(table[,f_set], 2, function(z) cor.test(z, table[, drug], method="spearman", use="na.or.complete")$p.value)
    features_2 <- names(sort(cor_data))[1:200]
    m_2 = simple_regression(x=table[,features_2], y=table[, drug])
    
    sd_data <- apply(table[,f_set], 2, function(z) sd(z, na.rm=T))
    features_3 <- names(sort(sd_data, decreasing = T))[1:200]
    m_3 = simple_regression(x=table[,features_3], y=table[, drug])
    
    bgp200cor <- bgp_features(table[,features_2])
    cor_data <- apply(bgp200cor, 2, function(z) cor.test(z, table[, drug], method="spearman", use="na.or.complete")$p.value)
    features_5 <- names(sort(cor_data))[1:200]
    m_5 = simple_regression(x=bgp200cor[,features_5], y=table[, drug])
    
    mult200cor <- mult_features(table[,features_2])
    cor_data <- apply(mult200cor, 2, function(z) cor.test(z, table[, drug], method="spearman", use="na.or.complete")$p.value)
    features_7 <- names(sort(cor_data))[1:200]
    m_7 = simple_regression(x=mult200cor[,features_7], y=table[, drug])
    
    # bgp200var <- bgp_features(table[,features_3])
    # cor_data <- apply(bgp200var, 2, function(z) cor.test(z, table[, drug], method="spearman", use="na.or.complete")$p.value)
    # features_9 <- names(sort(cor_data))[1:200]
    # m_9 = simple_regression(x=bgp200var[,features_9], y=table[, drug])
    # 
    # mult200var <- mult_features(table[,features_3])
    # cor_data <- apply(mult200var, 2, function(z) cor.test(z, table[, drug], method="spearman", use="na.or.complete")$p.value)
    # features_11 <- names(sort(cor_data))[1:200]
    # m_11 = simple_regression(x=mult200var[,features_11], y=table[, drug])
    
    pred_2 <- predict.train(m_2, newdata=table_test[,features_2], s="lambda.min")
    pred_3 <- predict.train(m_3, newdata=table_test[,features_3], s="lambda.min")
    
    bgp200cor_test <- bgp_features(table_test[,features_2])
    pred_5 <- predict.train(m_5, newdata=bgp200cor_test[,features_5], s="lambda.min")
    
    mult200cor_test <- mult_features(table_test[,features_2])
    pred_7 <- predict.train(m_7, newdata=mult200cor_test[,features_7], s="lambda.min")
    
    # bgp200var_test <- bgp_features(table_test[,features_3])
    # pred_9 <- predict.train(m_9, newdata=bgp200var_test[,features_9], s="lambda.min")
    # 
    # mult200var_test <- mult_features(table_test[,features_3])
    # pred_11 <- predict.train(m_11, newdata=mult200var_test[,features_11], s="lambda.min")
    
    # R2
    r2_2 <- r2(pred_2, table_test[,drug_test])
    r2_3 <- r2(pred_3, table_test[,drug_test])
    r2_5 <- r2(pred_5, table_test[,drug_test])
    r2_7 <- r2(pred_7, table_test[,drug_test])
    # r2_9 <- r2(pred_9, table_test[,drug_test])
    # r2_11 <- r2(pred_11, table_test[,drug_test])
    
    # ci
    ci_2 <- ci(pred_2, table_test[,drug_test])
    ci_3 <- ci(pred_3, table_test[,drug_test])
    ci_5 <- ci(pred_5, table_test[,drug_test])
    ci_7 <- ci(pred_7, table_test[,drug_test])
    # ci_9 <- ci(pred_9, table_test[,drug_test])
    # ci_11 <- ci(pred_11, table_test[,drug_test])
    
    
    #saving res. 
    r2_res <- c(r2_2, r2_3, r2_5, r2_7)
    ci_res <- c(ci_2, ci_3, ci_5, ci_7)
    if (j==1) {
      r2_gdsc[,i] <- r2_res
      ci_gdsc[,i] <- ci_res }
    
    if (j==2) {
      r2_ctrp[,i] <- r2_res
      ci_ctrp[,i] <- ci_res }
  } 
}

write.table(r2_gdsc, file="feature_interaction_rf_r2_gdsc.txt")
write.table(r2_ctrp, file="feature_interaction_rf_r2_ctrp.txt")
write.table(ci_gdsc, file="feature_interaction_rf_ci_gdsc.txt")
write.table(ci_ctrp, file="feature_interaction_rf_ci_ctrp.txt")


### 3. Plotting glmnet results ###

r2_gdsc <- read.table("feature_interaction_r2_gdsc.txt")

f_list <- c("1) all features", "2) 200 top correlation","3) 200 top variance",
            "4) BGP (200 top cor) all", "5) BGP (200 top cor) 200 top corr",
            "6) Multipl (200 top cor) all", "7) Multipl (200 top cor) 200 top corr",
            "8) BGP (200 top var) all", "9) BGP (200 top var) 200 top corr", 
            "10) Multipl (200 top var) all", "11) Multipl (200 top var) 200 top corr")

table <- cbind(f_list, apply(r2_gdsc, 1, mean),apply(r2_gdsc, 1, sd), c(rep("original",3), rep(c("BGP","BGP" ,"multiplication", "multiplication"),2)))
colnames(table) <- c("features", "R2","sd", "feature_type")
table <- as.data.frame(table)
table$R2 <- as.numeric(as.character(table$R2))
table$sd <- as.numeric(as.character(table$sd))
table$features <- factor(table$features, levels = table$features[order(table$R2, decreasing=F)])
table$feature_type <- factor(table$feature_type, levels = c("original","BGP", "multiplication"))

p <- ggplot(table, aes(y=R2, x=features, fill=feature_type)) +
  geom_bar(position="dodge", stat="identity") +  geom_errorbar(aes(ymin=R2-sd, ymax=R2+sd)) +
  labs(title = "average R2 across 19 drugs") + 
  theme(axis.text=element_text(size=12), axis.text.y = element_text(hjust = 0)) +
  coord_flip() 
p
ggsave("fig_17.pdf")
