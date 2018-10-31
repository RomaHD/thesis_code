# written on 14.07.2017, modified on 31.10.18 
# Chapter 7 of the thesis. 
# Shiny application for complex drug response visualisation:
# preparing files with raw viabilities from CTRP and GDSC datasets

# raw data from GDSC and CTRP projects should be downloaded and put in the script directory:
# GDSC:
# ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/v17a_public_raw_data.csv
# https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS1F.xlsx -- converted to "Screened_Compounds_GDSC1000.csv"
#
# CTRP:
# files from ftp://anonymous:guest@caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/
  
path <- "/data/kurilov/genestack/phd/work_2018/thesis_code/"
setwd(file.path(path, "chapter7_applied"))

### 1. GDSC ###

gdsc_data <- read.csv("v17a_public_raw_data.csv")
gdsc_drugs <- read.csv("Screened_Compounds_GDSC1000.csv")

pos_control <- apply(gdsc_data[,23:70], 1, function(x) {
  mean(as.numeric(as.character(unlist(x))), na.rm=TRUE)
})

neg_control <- apply(gdsc_data[,71:102], 1, function(x) {
  mean(as.numeric(as.character(unlist(x))), na.rm=TRUE)
})

pos_control_norm <- pos_control-neg_control
response <- matrix(as.numeric(as.character(unlist(gdsc_data[,14:22]))), 225384,9, byrow = FALSE)
response_norm <- response-neg_control
response_norm <- (response_norm/pos_control_norm)*100

drug_vect <- gdsc_drugs$DRUG.NAME[match(gdsc_data$DRUG_ID, gdsc_drugs$DRUG.ID)]
table_gdsc_ind_conc <- cbind(drug_vect, gdsc_data[,c(10,12,13)], response_norm)
colnames(table_gdsc_ind_conc)[1] <- "DRUG_NAME"
save(table_gdsc_ind_conc, file="table_gdsc_ind_conc.RData")


### 1. CTRP ###

ctrp_data <- read.table("CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.per_cpd_pre_qc.txt", header=TRUE)
experiment_data <- read.table("CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt", header=TRUE)
cell_data <- read.table("CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt", header=TRUE, sep="\t")
drug_data <- read.csv("CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt", header=TRUE, sep="\t")

sizes <- c(5:17, 25, 29)
s_numbers <- sapply(sizes, function(x) {
  which(ctrp_data$num_conc_pts == x)
})

for (i in 1:15)
{
  data <- ctrp_data[s_numbers[[i]],]
  response <- matrix(data$cpd_avg_pv, ncol=sizes[i], byrow=TRUE) 
  conc <- matrix(data$cpd_conc_umol, ncol=sizes[i], byrow=TRUE) 
  
  conc <- t(apply(conc, 1, function(x){
    rev(as.numeric(as.character(x)))
  }))
  response <- t(apply(response,1, function(x){
    rev(as.numeric(as.character(x)))
  }))
  response <- response*100
  
  select <- seq(1, length(data$experiment_id), sizes[i])
  experiment_id <- data$experiment_id[select]
  master_cpd_id <- data$master_cpd_id[select]
  
  cell_line <- sapply(experiment_id, function(x) {
    cell_data$ccl_name[which(cell_data$master_ccl_id==experiment_data$master_ccl_id[which(experiment_data$experiment_id==x)])]
  })
  drug <- sapply(master_cpd_id, function(x) {
    drug_data$cpd_name[which(drug_data$master_cpd_id==x)]
  })
  
  
  if (ncol(response)<29) {
    response <- cbind(response, matrix(NA, nrow(response), (29-ncol(response))))
    conc <- cbind(conc, matrix(NA, nrow(conc), (29-ncol(conc))))
  }
  table_chunk <- cbind(as.character(cell_line), as.character(drug), response, conc)
  if (i!=1)
  {
    main_table <- rbind(main_table, table_chunk)
  } else {main_table <- table_chunk}
  
  
}

table_ctrp_ind_conc <- as.data.frame(main_table)
colnames(table_ctrp_ind_conc) <- c("cell_line", "drug", c(1:29), paste0("c",1:29))
save(table_ctrp_ind_conc, file="table_ctrp_ind_conc.RData")

