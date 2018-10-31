# written in February 2017, modified on 27.10.18 
# Chapter 3 of the thesis. 
# Methods for improving accuracy of drug response prediction: 
# 2. Modelling on aggregated data: analysis, plotting the results
#
# datasets:CTRP, GDSC

source("https://bioconductor.org/biocLite.R")
biocLite("PharmacoGx")
library(PharmacoGx)
library(caret)
library(survcomp)
library(ggplot2)


path <- "/data/kurilov/genestack/phd/work_2018/thesis_code/"
setwd(file.path(path, "chapter3_methods"))

# load CTRP, GDSC (w/o imputation)
load(file.path(path, "data/ctrp_table.RData"))
load(file.path(path, "data/gdsc_table.RData"))

# PharmacoGx annotation data
setwd(file.path(path, "data"))
CTRP <- downloadPSet("CTRPv2")
GDSC <- downloadPSet("GDSC1000")
setwd(file.path(path, "chapter3_methods"))


### 1. Analysis ###

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

gen_last_gdsc <- max(grep("_", colnames(gdsc_table)))
gen_last_ctrp <- max(grep("_", colnames(ctrp_table)))
f_set <- intersect(colnames(gdsc_table)[1:gen_last_gdsc], colnames(ctrp_table)[1:gen_last_ctrp])
# 16126 exp + 146 mut 

# regr. function
simple_regression <- function(y, x, method, aggregate=F) {
  ctrl <- trainControl(method = "cv", number = 10, search = "random")
  set.seed(100)
  
  # table for regression on aggregated data (from several drugs)
  if (aggregate==T){
    drug_n <- ncol(y)
    text <- paste0(rep("x", drug_n), collapse=",")
    text <- paste0("rbind(", text ,")")
    x <- eval(parse(text=text))
    y <- as.vector(y)
  }
  
  na <- which(is.na(y))
  if(length(na)>0) {
    x <- as.matrix(x[-na,])
    y <- y[-na]
  } else {  x <- as.matrix(x)}
  
  cor_data <- apply(x, 2, function(z) cor.test(z, y, method="spearman", use="na.or.complete")$p.value)
  features <- names(sort(cor_data))[1:200]
  print(dim(x))
  Tune <- train(x[,features], y,
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

#concordance.index 
ci <- function(v1, v2) {
  value <- concordance.index(x=v1, surv.time=v2, surv.event=rep(1,length(v1)), na.rm=TRUE)
  value2 <- 1-value$c.index
  return(value2)
}

require(doMC)
registerDoMC(cores=2)

r2_table <- matrix(NA, nrow=4, ncol=38)
ci_table <- matrix(NA, nrow=4, ncol=38)
method_list <- c("rf", "svmRadial")

#renanaming tables (to use table names from the original script)
table_gdsc <- gdsc_table
table_ctrp <- ctrp_table

for (j in 1:2) {
  method <- method_list[j]
  
  # aggregated models  
  gdsc_hdac <- simple_regression(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$hdac],method=method ,aggregate=T)
  gdsc_egfr <- simple_regression(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$egfr], method=method,aggregate=T)
  gdsc_mek <- simple_regression(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$mek], method=method,aggregate=T)
  gdsc_braf <- simple_regression(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$braf],method=method ,aggregate=T)
  gdsc_hsp90 <- simple_regression(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$hsp90], method=method,aggregate=T)
  gdsc_mtor <- simple_regression(x=table_gdsc[,f_set], y=table_gdsc[, drug_classes$gdsc$mtor], method=method,aggregate=T)
  
  ctrp_hdac <- simple_regression(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$hdac],method=method, aggregate=T)
  ctrp_egfr <- simple_regression(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$egfr],method=method ,aggregate=T)
  ctrp_mek <- simple_regression(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$mek],method=method ,aggregate=T)
  ctrp_braf <- simple_regression(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$braf],method=method ,aggregate=T)
  ctrp_hsp90 <- simple_regression(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$hsp90], method=method,aggregate=T)
  ctrp_mtor <- simple_regression(x=table_ctrp[,f_set], y=table_ctrp[, drug_classes$ctrp$mtor], method=method,aggregate=T)
  
  for (i in 1:19) {
    print(i)
    drug_gdsc <- names[i,1]
    drug_ctrp <- names[i,2]
    class <- names[i,3]
    
    # GDSC -> CTRP
    gdsc_m = get(paste0("gdsc_", class))
    gdsc_s = simple_regression(x=table_gdsc[,f_set], y=table_gdsc[, drug_gdsc], method=method)
    
    # pred_gdsc_m_gdsc <- predict.train(gdsc_m, newdata=table_gdsc[,colnames(gdsc_m$trainingData)[1:200]])
    # pred_gdsc_s_gdsc <- predict.train(gdsc_s, newdata=table_gdsc[,colnames(gdsc_s$trainingData)[1:200]])
    pred_gdsc_m_ctrp <- predict.train(gdsc_m, newdata=table_ctrp[,colnames(gdsc_m$trainingData)[1:200]])
    pred_gdsc_s_ctrp <- predict.train(gdsc_s, newdata=table_ctrp[,colnames(gdsc_s$trainingData)[1:200]])
    
    # CTRP -> GDSC
    ctrp_m = get(paste0("ctrp_", class))
    ctrp_s = simple_regression(x=table_ctrp[,f_set], y=table_ctrp[, drug_ctrp], method=method)
    
    # pred_ctrp_m_ctrp <- predict.cv.glmnet(ctrp_m, newdata=table_ctrp[,colnames(m_test$trainingData)[1:200]])
    # pred_ctrp_s_ctrp <- predict.cv.glmnet(ctrp_s, newdata=table_ctrp[,colnames(m_test$trainingData)[1:200]])
    pred_ctrp_m_gdsc <- predict.train(ctrp_m, newdata=table_gdsc[,colnames(ctrp_m$trainingData)[1:200]])
    pred_ctrp_s_gdsc <- predict.train(ctrp_s, newdata=table_gdsc[,colnames(ctrp_s$trainingData)[1:200]])
    
    # R2
    r2_1 <- r2(pred_gdsc_m_ctrp, table_ctrp[,drug_ctrp])
    r2_2 <- r2(pred_ctrp_m_gdsc, table_gdsc[,drug_gdsc])
    r2_3 <- r2(pred_gdsc_s_ctrp, table_ctrp[,drug_ctrp])
    r2_4 <- r2(pred_ctrp_s_gdsc, table_gdsc[,drug_gdsc])
    
    # concordance index
    ci_1 <- ci(pred_gdsc_m_ctrp, table_ctrp[,drug_ctrp])
    ci_2 <- ci(pred_ctrp_m_gdsc, table_gdsc[,drug_gdsc])
    ci_3 <- ci(pred_gdsc_s_ctrp, table_ctrp[,drug_ctrp])
    ci_4 <- ci(pred_ctrp_s_gdsc, table_gdsc[,drug_gdsc])
    
    #saving res. 
    col <- i*2-1
    row <- j*2-1
    r2_table[row:(row+1),col:(col+1)] <- signif(cbind(c(r2_1,r2_2), c(r2_3,r2_4)), digits=3)
    ci_table[row:(row+1),col:(col+1)] <- signif(cbind(c(ci_1,ci_2), c(ci_3,ci_4)), digits=3)
    
  } }

#save results

write.table(r2_table, file="svm_rf_aggr_r2.txt")
write.table(ci_table, file="svm_rf_aggr_ci.txt")


### 2. Plotting the results ###

r2_table <- read.table("svm_rf_aggr_r2.txt")

t1 <- r2_table[,seq(1,38,by=2)]
t2 <- r2_table[,seq(2,38,by=2)]
colnames(t1) <- colnames(t2) <- NULL
table1 <- rbind(as.matrix(t1), as.matrix(t2))
table1 <- cbind(rep(c("rf", "rf", "svm", "svm"), 2),rep(c("GDSC", "CTRP"),4), c(rep("aggregated",4), rep("simple",4)), table1)
colnames(table1) <- c("method","training","model","entinostat","tubastatin A", "belinostat","vorinostat", "lapatinib", "erlotinib", "afatinib", "gefitinib",
                      "trametinib","selumetinib", "PLX-4720","dabrafenib", "SNX-2112","tanespimycin", "OSI-027", 
                      "sirolimus", "NVP-BEZ235", "temsirolimus", "AZD8055")

table2 <- as.data.frame(table1)
table2 <- table2 %>% 
  gather(drug, r2, -c(method, training, model)) 

table2$r2 <- as.numeric(as.character(table2$r2))
table2$training <- as.factor(table2$training)
table2$model <- as.factor(table2$model)
table2$drug <- factor(table2$drug, levels = c("entinostat","tubastatin A", "belinostat","vorinostat", "lapatinib", "erlotinib", "afatinib", "gefitinib",
                                              "trametinib","selumetinib", "PLX-4720","dabrafenib", "SNX-2112","tanespimycin", "OSI-027", 
                                              "sirolimus", "NVP-BEZ235", "temsirolimus", "AZD8055"))

method_list <- c("rf","rf", "svm", "svm")
training_list <- c("GDSC", "CTRP", "GDSC", "CTRP")
for (i in 1:4) {
  table3 <- table2[which(table2$method==method_list[i] & table2$training==training_list[i]),]
  
  p <- ggplot(table3, aes(y=r2, x=drug, fill=model)) +
    geom_bar(position="dodge", stat="identity") + #scale_y_continuous(limits=c(0.5,0.85),oob = rescale_none) +
    #scale_fill_manual(values=c("#0072B2", "#D55E00")) +
    labs(title = paste0(method_list[i],", training on ", training_list[i])) +# theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())
    theme(axis.text=element_text(size=8))
  assign(paste0("p",i), p)
}

grid.arrange(p1,p2,p3,p4, nrow=2,ncol=2)
ggsave("fig_16.pdf") 