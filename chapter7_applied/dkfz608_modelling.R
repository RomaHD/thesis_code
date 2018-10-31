# written on 2.08.2018, modified on 30.10.18 
# Chapter 7 of the thesis. 
# Drug response prediction model for DKFZ-608 compound
# building DFKZ-608 model using local IC50 data and genomics data from GDSC 

library(PharmacoGx)
library(caret)

path <- "/data/kurilov/genestack/phd/work_2018/thesis_code/"

# getting data from PharmacoGx
setwd(file.path(path, "data"))
GDSC <- downloadPSet("GDSC1000")

load("gdsc_table.RData")
ic50 <- read.csv("ic50_updated.csv", header=F)

setwd(file.path(path, "chapter7_applied"))

### 1. identifying features which are correlated with IC50 values ###

#adding tissue variables to the genomic table
tissue_gdsc <- cellInfo(GDSC)[,c("cellid","tissueid")]
colnames(tissue_gdsc) <-c("sample","tissue")
order <- match(ic50$V1, tissue_gdsc$sample)
ic50  <- ic50[-which(is.na(order)),]
order <- na.omit(order)

dummy_test <- model.matrix(~tissue_gdsc[order,"tissue"])
dummy_test <- dummy_test[,-which(apply(dummy_test,2,sd)==0)]
colnames(dummy_test) <- substr(colnames(dummy_test),33,50)

# gdsc_table <- as.matrix(gdsc_table)
# mode(gdsc_table) <- "numeric"

gen_last_gdsc <- max(grep("_", colnames(gdsc_table)))
table_gdsc <- cbind(gdsc_table[as.character(ic50$V1),1:gen_last_gdsc], dummy_test)

# function for assessing the correlation
cor_fs <- function(x, y) {
  
  p_vect <- apply(x,2, function(z) {cor.test(z, y, metho="spearman", use="na.or.complete")$p.value})
  return(p_vect)
}

f_set <- cor_fs(x=table_gdsc, y=ic50$V2)
f_set_adj <- p.adjust(f_set, method="BH")
f_set_adj2 <- sort(f_set_adj[which(f_set_adj<0.001)])
cor_coef <- apply(table_gdsc[,names(f_set_adj2)],2, function(z) {cor.test(z, y=ic50$V2, method="spearman", use="na.or.complete")$estimate})
write.table(cbind(names(f_set_adj2), f_set_adj2, cor_coef), file="associations_aug18.txt")

### 2. identifying drugs with similar drug resp. (IC50) profiles ###

# getting IC50 data from PharmacoGx
ic50_gdsc <- t(summarizeSensitivityProfiles(GDSC, sensitivity.measure="ic50_recomputed", drugs=drugNames(GDSC)))

# function for truncating IC50 values
outlier_capping <- function(table){
  for (i in 1:ncol(table)){
    fin <- which(is.finite(table[,i]))
    qv <- quantile(table[fin,i], probs=0.85, na.rm=T)
    out <- which(table[,i]>qv)
    table[out,i] <- qv
  }
  return(table)
}

ic50_gdsc <- outlier_capping(ic50_gdsc)
ic50_gdsc_subset <- ic50_gdsc[as.character(ic50$V1),]
n <- which(is.na(apply(ic50_gdsc_subset,2,function(x){sd(x,na.rm=T)})))
ic50_gdsc_subset <- ic50_gdsc_subset[,-n]

drug_sim <- apply(ic50_gdsc_subset,2, function(z) {cor.test(z, ic50$V2, method="spearman", use="na.or.complete")$estimate})
drug_sim2 <- head(sort(drug_sim, decreasing =T), n=20)
write.table(cbind(names(drug_sim2), drug_sim2, drugInfo(GDSC)[match(names(drug_sim2), drugInfo(GDSC)$DRUG.NAME),c(4,5)]), file="drug_similarity_aug18.txt")

### 3. creating a model for drug response 
### and generating predictions for those GDSC cell lines that were not tested with the compound ###

simple_regression <- function(y, x, method) {
  ctrl <- trainControl(method = "cv", number = 10, search = "random")
  na <- which(is.na(y) | is.infinite(y))
  if(length(na)>0) {
    x <- as.matrix(x[-na,])
    y <- y[-na]
  } else {  x <- as.matrix(x)}
  
  print(dim(x))
  Tune <- train(x, y,
                method = method,
                tuneLength = 30,
                trControl = ctrl,
                preProc = c("medianImpute", "center","scale"))
  
  return(Tune)
}

model_rf <- simple_regression(x=table_gdsc[,names(f_set_adj2)], y=ic50$V2, method="rf")
model_svm <- simple_regression(x=table_gdsc[,names(f_set_adj2)], y=ic50$V2, method="svmRadial")
models <- list(model_rf=model_rf,model_svm=model_svm)
save(models, file="models_aug18.RData")
# best R2 (rf)=0.43, best R2 (svm)=0.62 

test_lines <- setdiff(rownames(gdsc_table), as.character(ic50$V1))
leukemia <- lung_NSCLC <- rep(0, nrow(gdsc_table))
leukemia[which(tissue_gdsc$tissue[match(rownames(gdsc_table), tissue_gdsc$sample)]=="leukemia")] <- 1
lung_NSCLC[which(tissue_gdsc$tissue[match(rownames(gdsc_table), tissue_gdsc$sample)]=="lung_NSCLC")] <- 1
gdsc_table <- cbind(gdsc_table, leukemia, lung_NSCLC)

pred_rf <- predict.train(model_rf, gdsc_table[test_lines,names(f_set_adj2)])
pred_svm <- predict.train(model_svm, gdsc_table[test_lines,names(f_set_adj2)])

# ZR-75 is the last cell line with expression data, predictions for cell lines w/o exp. doesn't make sense
pred <- cbind(pred_rf, pred_svm)[1:grep("ZR-75", test_lines),]
ord <- match(rownames(pred), cellInfo(GDSC)$Sample.Name)
pred <- cbind(pred, cellInfo(GDSC)[ord,c(8,9)], (pred[,1]+pred[,2])/2)
write.table(pred, file="predicted_ic50_gdsc_aug18.txt")