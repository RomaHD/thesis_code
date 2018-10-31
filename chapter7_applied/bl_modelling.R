# written on 20.10.2016, modified on 30.10.18 
# Chapter 7 of the thesis. 
# Burkitt lymphoma drug sensitivity screen analysis:
# modelling using GDSC data

# training models on pan-cancer cohort (1000 samples)/ lymphoma lines (50 samples) 
# and predicting on lines from GDSC matching those from Burkitt lymphoma screen (15 samples). 
# Methods: SVM (FS: corr), elastic net  

library(PharmacoGx)
library(caret)
library(doMC)
registerDoMC(cores = 1)
library(e1071)

path <- "/data/kurilov/genestack/phd/work_2018/thesis_code/"

# getting data from PharmacoGx
setwd(file.path(path, "data"))
GDSC <- downloadPSet("GDSC1000")

load("gdsc_table.RData")

setwd(file.path(path, "chapter7_applied"))

### 1.  Modelling example: SVM, all genomic data, 200 variables in each model (selected using Spearman cor.)###

# lists of drugs and lines 

# lines from cell line screen
Kasia_lines <- c("MEC-1", "NAMALWA", "Raji", "Ramos-2G6-4C10", "CA46", "DG-75", "BL-41", "K-562",
                 "SU-DHL-5", "HT", "GRANTA-519", "JEKO-1", "JVM-2", "OPM-2", "RPMI-8226", "U-266")
Kasia_lines <- intersect(Kasia_lines, rownames(gdsc_table))
drugs <- c("CAL-101", "Olaparib", "Vorinostat", "Nutlin-3a (-)", "selumetinib",
           "Afatinib", "AZD7762", "MK-2206", "NU-7441", "Thapsigargin", "YM155", 
           "Saracatinib", "Dasatinib", "Gefitinib", "Lapatinib", "Ruxolitinib",
           "Obatoclax Mesylate", "Doxorubicin", "JQ1", "Trametinib")

# feature selection: Spearman corrlation
cor_table <- function(data_table, drug_number){
  
  data_table <- as.matrix(data_table)
  mode(data_table) <- "numeric"
  features_num <- dim(data_table)[2]-drug_number
  cor_table2=matrix(NA,drug_number,200)
  rownames(cor_table2) <- colnames(data_table)[(features_num+1):dim(data_table)[2]]
  
  for (i in 1:drug_number) {
    na_num <- which(is.na(data_table[,(features_num+i)]))
    if (length(na_num)<1) {data <- apply(data_table[,1:features_num],2, function(x) cor.test(x, data_table[,(features_num+i)], method="spearman", use="na.or.complete")$p.value)
    } else { data <- apply(data_table[-na_num,1:features_num],2, function(x) cor.test(x, data_table[-na_num,(features_num+i)], method="spearman", use="na.or.complete")$p.value)  }
    cor_table2[i,]<- names(sort(data))[1:200]
  }
  print("Correlation table:")
  print(cor_table2[1:5,1:5])
  return(cor_table2)
}

pan_samples <- setdiff(rownames(gdsc_table), Kasia_lines)
lymph_samples <- setdiff(cellInfo(GDSC)$Sample.Name[which(cellInfo(GDSC)$GDSC.Tissue.descriptor.1=="lymphoma")], Kasia_lines)
lymph_samples <- intersect(lymph_samples, rownames(gdsc_table))
mode(pan_samples) <- "character"
mode(lymph_samples) <- "character"

gen_last_gdsc <- max(grep("_", colnames(gdsc_table)))
cor_table_pan <- cor_table(cbind(gdsc_table[pan_samples,1:gen_last_gdsc], gdsc_table[pan_samples,drugs]),20)
cor_table_lymph<- cor_table(cbind(gdsc_table[lymph_samples,1:gen_last_gdsc], gdsc_table[lymph_samples,drugs]),20)

# save/load FS results
save(cor_table_pan, file="cor_table_pan.RData")
save(cor_table_lymph, file="cor_table_lymph.RData")

load("cor_table_pan.RData")
load("cor_table_lymph.RData")

# regression function
simple_regression <- function(y, x) {
  ctrl <- trainControl(method = "cv", number = 10, search = "random")
  set.seed(100)
  na <- which(is.na(y))
  if(length(na)>0) {
    x <- as.matrix(x[-na,])
    y <- y[-na]
  } else {  x <- as.matrix(x)}
  
  print(dim(x))
  Tune <- train(x, y,
                method = "svmRadial",
                tuneLength = 30,
                trControl = ctrl,
                preProc = c("center","scale","medianImpute"))
  
  return(Tune)
}

#for results
results_table <- matrix(NA, 20, 4)
rownames(results_table) <- drugs
colnames(results_table) <- c("R2_cor_pan", "RMSE_cor_pan", "R2_cor_lymph", "RMSE_cor_lymph")

# model's trainig
for (i in 1:20) {
  
  print(drugs[i])
  
  y_pan <- gdsc_table[pan_samples,drugs[i]]
  y_lymph <- gdsc_table[lymph_samples,drugs[i]]
  
  x_pan_cor <- gdsc_table[pan_samples,cor_table_pan[i,]]
  x_lymph_cor <- gdsc_table[lymph_samples,cor_table_lymph[i,]]
  
  model1 <- try(simple_regression(y_pan, x_pan_cor))
  model2 <- try(simple_regression(y_lymph, x_lymph_cor))
  
  new_data <- gdsc_table[Kasia_lines,]
  y_new_data <- new_data[,drugs[i]]
  y_model1 <- predict.train(model1, newdata=new_data[,cor_table_pan[i,]])
  y_model2 <- predict.train(model2, newdata=new_data[,cor_table_lymph[i,]])
  
  r2_1 <- R2(pred=y_model1, obs = y_new_data, na.rm=T)
  rmse_1 <- RMSE(pred=y_model1, obs = y_new_data, na.rm=T)
  r2_2 <- R2(pred=y_model2, obs = y_new_data, na.rm=T)
  rmse_2 <- RMSE(pred=y_model2, obs = y_new_data, na.rm=T)
  
  results_table[i,1:4] <- c(r2_1, rmse_1,r2_2, rmse_2) 
}

save(results_table, file="results_table_kasia_lines_modelling.RData")

### 2. Plotting average R2 for all results (SVM/enet, 200/1000 var, all genomic data) ###
load(file.path(path, "data/kasia_mod_results_all.RData"))

titles <- c("FS: 200 var. by cor; modelling: SVM", "FS: 200 var. by mRMR; modelling: SVM", 
            "FS: 1000 var. by cor; modelling: elastic net","FS: 200 var. by cor; modelling: elastic net", 
            "FS: 200 var. by mRMR; modelling: elastic net", "FS: 1000 var. by cor; modelling: SVM")


# avereging
pan <- cbind(kasia_mod_results[[1]][,1], kasia_mod_results[[2]][,1], kasia_mod_results[[3]][,1], kasia_mod_results[[4]][,1], kasia_mod_results[[5]][,1], kasia_mod_results[[6]][,1])
lymph <- cbind(kasia_mod_results[[1]][,3], kasia_mod_results[[2]][,3], kasia_mod_results[[3]][,3], kasia_mod_results[[4]][,3], kasia_mod_results[[5]][,3], kasia_mod_results[[6]][,3])
average_pan <- apply(pan,1, mean)
average_lymph <- apply(lymph,1, mean)
sd_pan <- apply(pan,1, sd)
sd_lymph <- apply(lymph,1, sd)
avg <- (average_pan+average_lymph)/2
table_avg <- cbind(c(average_pan, average_lymph),c(sd_pan, sd_lymph), 
                   rep(rownames(pan), 2), 
                   c(rep("pan-cancer", 20), rep("lymphoma", 20)))
colnames(table_avg) <- c("R2","sd", "drugs", "training")
table_avg <- as.data.frame(table_avg)
table_avg$R2 <- as.numeric(as.character(table_avg$R2)) 
table_avg$sd <- as.numeric(as.character(table_avg$sd)) 
table_avg$training <- factor(table_avg$training, levels = c("pan-cancer", "lymphoma"))
table_avg$drugs <- factor(table_avg$drugs, levels=table_avg$drugs[order(avg, decreasing = T)])

# modified on 20.09.2018
g <- ggplot(table_avg, aes(x=drugs, y=R2, fill=training))
g2 <- g + geom_bar(position = position_dodge(), stat="identity") + labs(title = "Average R2", ylab="R2") +geom_errorbar(aes(ymax=R2+sd, ymin=pmax(R2-sd,0)), position = position_dodge(0.9)) 
g2 <- g2 +theme(axis.text=element_text(size=9), legend.position = "bottom")
g2
ggsave("Fig_34.pdf")
