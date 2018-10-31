# written on 1.08.2016, modified on 30.10.18 
# Chapter 6 of the thesis. 
# Patient treatment outcome prediction using classification models trained on cell lines:
# Main analysis


library(caret)
library(e1071)
library(pROC)
library(doMC)
registerDoMC(cores = 1)

path <- "/data/kurilov/genestack/phd/work_2018/thesis_code/"

### 1. loading data ###
setwd(file.path(path, "data/comb_data"))
load("comb_tables_list.RData")

for (i in 1:17)
{
  load(paste0(comb_tables_list[i],".RData"))
}
setwd(file.path(path, "chapter6_patient_response"))


###. 2. functions ###
# twoClassSummary function modified for returning the balanced accuracy

twoClassSummary_mod <- function (data, lev = NULL, model = NULL)
{
  if(length(levels(data$obs)) > 2)
    stop(paste("Your outcome has", length(levels(data$obs)),
               "levels. The twoClassSummary() function isn't appropriate."))
  caret:::requireNamespaceQuietStop('pROC')
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"])))
    stop("levels of observed and predicted data do not match")
  rocObject <- try(pROC::roc(data$obs, data[, lev[1]], direction = ">"), silent = TRUE)
  rocAUC <- if(class(rocObject)[1] == "try-error") NA else rocObject$auc
  out <- c(rocAUC,
           sensitivity(data[, "pred"], data[, "obs"], lev[1]),
           specificity(data[, "pred"], data[, "obs"], lev[2]))
  out <- c(out, (out[2]+out[3])/2)
  names(out) <- c("ROC", "Sens", "Spec", "Balanced acc")
  out
}

# SVM modelling func.
simple_classification <- function(y, x) {
  
  ctrl <- trainControl(method = "cv", 
                       number = 10, 
                       search = "random", 
                       classProbs = T,
                       summaryFunction = twoClassSummary_mod)
  set.seed(100)
  na <- which(is.na(y) | y=="NA")
  if(length(na)>0) {
    x <- as.matrix(x[-na,])
    y <- y[-na]
  } else {  x <- as.matrix(x)}
  
  print(dim(x))
  mode(x) <- "numeric"
  
  print(summary(as.factor(y)))
  y <- factor(y, levels=c("sens","resist"))
  
  Tune <- train(x, y,
                method = "svmRadial",
                tuneLength = 30,
                trControl = ctrl,
                metric="Balanced acc",
                maximize=T,
                preProc = c("center","scale","medianImpute"))
  
  return(Tune)
}

# main analysis function 
analysis <- function(data_obj_name) {
  
  var_num <- 200
  dataset <- get(data_obj_name)
  data_obj_name_digest <- unlist(strsplit(data_obj_name ,split = "_"))
  dnd_len <- length(data_obj_name_digest)
  dataset_name <- paste(data_obj_name_digest[-c(1,dnd_len)], collapse="_")
  drug_name <- data_obj_name_digest[dnd_len]
  print(paste0("dataset name: ",dataset_name))
  print(paste0("drug name: ",drug_name))
  
  parameters <- matrix(c("ic50","no_g","pan","ic50","no_g","spec","ic50","g","pan","ic50","g","spec",
                         "auc","no_g","pan","auc","no_g","spec","auc","g","pan","auc","g","spec",
                         "viab","no_g","pan","viab","no_g","spec","viab","g","pan","viab","g","spec"),byrow = T,nrow=12)
  colnames(parameters) <- c("metric","grey_zone","tissue_spec")
  parameters <- as.data.frame(parameters)
  
  if (dataset_name == "nibr"){
    parameters <- parameters[1:8,]
    parameters$metric <- c(rep("bestresp",4),rep("bestavgresp",4))
    parameters$metric <- as.factor(parameters$metric)
  }
  
  print(parameters)
  
  tissue_type <- list(blood=c("haematopoietic_and_lymphoid_tissue","leukemia","lymphoma","myeloma","Blood","Lymphoid"),
                      lung=c("lung","lung_SCLC","lung_NSCLC","Lung","NSCLC"),
                      breast=c("breast","Breast","BRCA"))
  if (drug_name=="erlotinib") {tissue_type <- tissue_type$lung}
  if (drug_name=="bortezomib") {tissue_type <- tissue_type$blood}
  if (drug_name=="docetaxel") {tissue_type <- tissue_type$breast}
  
  # data split
  table <- dataset$exp_ber
  n <- dim(table)[1]
  # when we split just cell line data, we first sort data by drug response, here we basically don't
  test_n <- which(dataset$tissue=="pt")
  
  final_result <- matrix(NA, dim(parameters)[1], 11)
  colnames(final_result) <- c("drug","dataset","metric","grey_zone", "tissue_spec","train_size", "tr_bal_acc","tr_bal_accSD","accuracy","balanced acc" ,"auroc")
  
  for (i in 1:dim(parameters)[1]) {
    
    train_n <- setdiff(1:n, test_n)
    if(parameters$tissue_spec[i]=="spec") {
      tissue_n <- which(dataset$tissue %in% tissue_type)
      #         if(length(tissue_n)==0) {
      #           result <- c(drug_name,dataset_name,
      #                       as.character(parameters$metric[i]),as.character(parameters$grey_zone[i]), as.character(parameters$tissue_spec[i]),
      #                       0, NA, NA)
      #           final_result[i,] <-result
      #           break
      #        }
      train_n <- intersect(train_n, tissue_n)}
    
    if (parameters$grey_zone[i]=="no_g"){
      metric_full <- parameters$metric[i]
    } else {
      metric_full <- paste0(parameters$metric[i],"_g")
    }
    print(metric_full)
    drug_resp <- get(as.character(metric_full), dataset$labels)
    train <- cbind(table[train_n,], drug_resp[train_n])
    test <-  cbind(table[test_n,], drug_resp[test_n])
    
    print(dim(train))
    print(dim(test))
    print(metric_full)
    sd1 <- apply(train[,-dim(train)[2]],2, function(x) {sd(x, na.rm=TRUE)})
    sd2 <- apply(test[,-dim(test)[2]],2, function(x) {sd(x, na.rm=TRUE)})
    sd<- sd1+sd2
    sd_na <- which(is.na(sd) | sd==0)
    if (length(sd_na)>0) {
      train <- train[,-sd_na]
      test <- test[,-sd_na]
    }
    
    #### FS + model fit
    y <- train[,dim(train)[2]]
    
    if (length(unique(y))==1) {
      
      rmse <- NA
      r2 <- NA
      train_size <- 0
      
    } else {
      p_vect <- apply(train[,-dim(train)[2]],2, function(x) {anova(lm(x ~ y), test = "F")[1, "Pr(>F)"]})
      print(summary(p_vect))
      
      
      f_num <- which(p_vect<=(sort(p_vect, decreasing = FALSE)[var_num]))
      f_num <- f_num[1:var_num]
      x <- train[,f_num]
      
      model1 <- try(simple_classification(y, x))
      
      if(class(model1)!="try-error" & length(train_n)!=0) {
        
        newdata=test[,f_num]
        mode(newdata) <- "numeric"
        
        obs <- factor(test[,dim(test)[2]], levels = c("sens","resist"))
        pred <- predict.train(object=model1, newdata=newdata, type="raw")
        sens <- predict.train(object=model1, newdata=newdata, type="prob")[,"sens"]
        resist <- predict.train(object=model1, newdata=newdata, type="prob")[,"resist"]
        df <- data.frame(obs = obs, pred=pred, sens=sens, resist=resist)
        results <- twoClassSummary_mod(df, lev=c("sens","resist"))
        auroc <- results[1]
        bal_acc <- results[4]
        
        #auroc <- auc(roc(response=real, predictor=preds, levels=c("sens","resist")))
        accuracy <- postResample(pred=factor(pred, levels=c("sens","resist")),
                                 obs=factor(obs, levels=c("sens","resist")))[1]
        
        train_size <- length(y) - length(which(is.na(y)))
        tr_ba <- model1$resample$`Balanced acc`
        tr_bal_acc <- mean(tr_ba, na.rm=T)
        tr_bal_accSD <- sd(tr_ba, na.rm=T)
        
        
      } else {
        
        train_size <- 0
        tr_bal_acc <- NA
        tr_bal_accSD <- NA
        accuracy <- NA
        bal_acc <- NA
        auroc <- NA
      }
    }
    
    result <- c(drug_name,dataset_name,
                as.character(parameters$metric[i]),as.character(parameters$grey_zone[i]), as.character(parameters$tissue_spec[i]),
                train_size, tr_bal_acc, tr_bal_accSD, 
                accuracy,bal_acc ,auroc)
    final_result[i,] <-result
  }
  
  
  print(final_result)
  #final_result <- apply(result, 2, function(x) {mean(x, na.rm=TRUE)})
  return(final_result)
}



### 3. results collection and saving ###
validation_res <- matrix(NA,nrow=0,ncol=11)
colnames(validation_res) <- c("drug","dataset","metric","grey_zone", "tissue_spec","train_size","tr_bal_acc","tr_bal_accSD", "accuracy","balanced acc" ,"auroc")
for (i in 1:17)
{
  print(i)
  x <- comb_tables_list[i]
  print(x)
  
  results_part <- analysis(x)
  validation_res <- rbind(validation_res, results_part)
}
save(validation_res, file="validation_res.RData")

# test run
analysis("comb_nibr_erlotinib")