# written on 26.07.2017, modified on 31.10.18 
# Chapter 7 of the thesis. 
# Shiny application for complex drug response visualisation:
# app.R file

# "table_gdsc_ind_conc.RData" and "table_ctrp_ind_conc.RData" are produced by the script shiny_data.R
# "ann_data.RData" can be found in /data
# "table_auc_gdsc1000_pgx_2_cn.RData" and "table_auc_ctrp_pgx_2.RData" can be substituted with "gdsc_table.RData" and "ctrp_table.RData" from /data folder


library(shiny)
library(ggplot2)
library(dplyr)

# load gdsc and ctrp data

load("table_gdsc_ind_conc.RData")
load("table_ctrp_ind_conc.RData")
load("ann_data.RData")

load("table_auc_gdsc1000_pgx_2_cn.RData")
load("table_auc_ctrp_pgx_2.RData")

ui <- fluidPage(
  titlePanel("Drug response, viability for certain concentration", windowTitle = "Drug response, AUC"),
  sidebarLayout(
    sidebarPanel(selectInput("datasetInput", "dataset",
                             c("GDSC", "CTRP"),
                             selected = "GDSC"),
                 uiOutput("tissue_Output"),
                 uiOutput("mutation_Output"),
                 hr(),
                 "1st drug",
                 uiOutput("target1_Output"),
                 uiOutput("drug1_Output"),
                 uiOutput("conc1_Output"),
                 hr(),
                 "2nd drug", 
                 uiOutput("target2_Output"),
                 uiOutput("drug2_Output"),
                 uiOutput("conc2_Output")),
    mainPanel("the results will go here", 
              checkboxGroupInput("moltypes_Input", "types of mol. data", c("exp", "cn", "mut"), selected=c("exp", "cn", "mut")),
              plotOutput("coolplot", brush="plot_brush"),
              actionButton("cluster1", "define cluster1"),
              actionButton("cluster2", "define cluster2"),
              hr(),
              numericInput("pval_Input", "p-value cutoff",value=0.01 ,min=0, max=1),
              actionButton("showres_Input", "show mol. differences between clusters"),
              tableOutput("results")))
    #mainPanel("the results will go here", plotOutput("coolplot"), br(), br(), tableOutput("results")))
)

server <- function(input, output) {
  
  mut_gdsc <- grep("mut_", colnames(table_auc_gdsc1000_pgx_2), value=T) 
  mut_ctrp <- grep("mut_", colnames(table_auc_ctrp_pgx_2), value=T) 
  
  mutations <- reactive({
    if (input$datasetInput=="GDSC") {
      mut <- sort(mut_gdsc)
    }
    if (input$datasetInput=="CTRP") {
      mut <- sort(mut_ctrp)}
    return(mut)
  })
  
  tissues <- reactive({
    if (input$datasetInput=="GDSC") {
      t <- sort(unique(ann_data$cellInfo_gdsc$tissueid))
    }
    if (input$datasetInput=="CTRP") {
      t <- sort(unique(ann_data$cellInfo_ctrp$tissueid))}
    return(t)
  })
  
  targets <- reactive({
    if (input$datasetInput=="GDSC") {
      t <- sort(unique(unlist(strsplit(ann_data$drugInfo_gdsc$TARGET, split=", "))))
    }
    if (input$datasetInput=="CTRP") {
     t <- sort(unique(unlist(strsplit(ann_data$drugInfo_ctrp$gene_symbol_of_protein_target, split=";"))))}
    return(t)
  })
  
  drugs <- reactive({
    if (input$datasetInput=="GDSC") {
      d <- sort(ann_data$drugInfo_gdsc$DRUG.NAME)
    }
    if (input$datasetInput=="CTRP") {
      d <- sort(ann_data$drugInfo_ctrp$cpd_name)
    }
    return(d)
  })
  
  drugs1 <- reactive({
    d1 <- drugs()

    if (input$target1_Input!="all_targets") {
      if (input$datasetInput=="GDSC") {
      d1 <- ann_data$drugInfo_gdsc$DRUG.NAME[grep(input$target1_Input, ann_data$drugInfo_gdsc$TARGET)]}
      if (input$datasetInput=="CTRP") {
      d1 <- ann_data$drugInfo_ctrp$cpd_name[grep(input$target1_Input, ann_data$drugInfo_ctrp$gene_symbol_of_protein_target)]}
    } 
    return(sort(d1))
  })
  
  drugs2 <- reactive({
    d2 <- drugs()
    
    if (input$target2_Input!="all_targets") {
      if (input$datasetInput=="GDSC") {
        d2 <- ann_data$drugInfo_gdsc$DRUG.NAME[grep(input$target2_Input, ann_data$drugInfo_gdsc$TARGET)]}
      if (input$datasetInput=="CTRP") {
        d2 <- ann_data$drugInfo_ctrp$cpd_name[grep(input$target2_Input, ann_data$drugInfo_ctrp$gene_symbol_of_protein_target)]}
    } 
    return(sort(d2))
  })
  
  concentrations1 <- reactive({

      if (input$datasetInput=="GDSC") {
        record <- as.vector(table_gdsc_ind_conc[which(table_gdsc_ind_conc[,1]==input$drug1_Input)[1],])
        num <- 9 - length(which(is.na(record)))
        conc1 <- record[3]
 
        for (i in 2:num)
        {
          c <- conc1[i-1]/record[4]
          conc1 <- c(conc1, c)
        }
     
      }
      if (input$datasetInput=="CTRP") {
  
          record <- table_ctrp_ind_conc[which(table_ctrp_ind_conc[,2]==input$drug1_Input)[1],]
          #record <- as.numeric(as.character(record))
          num <- 29 - length(which(is.na(record[3:31])))
          conc1 <- record[32:(32+num-1)]
      }
   
    conc1 <- unname(conc1)
    return(unlist(conc1))
  })
  
 
  
  concentrations2 <- reactive({
    
    if (input$datasetInput=="GDSC") {
      record <- as.vector(table_gdsc_ind_conc[which(table_gdsc_ind_conc[,1]==input$drug2_Input)[1],])
      num <- 9 - length(which(is.na(record)))
      conc2 <- record[3]
      
      for (i in 2:num)
      {
        c <- conc2[i-1]/record[4]
        conc2 <- c(conc2, c)
      }
   
    }
    if (input$datasetInput=="CTRP") {
      
      record <- table_ctrp_ind_conc[which(table_ctrp_ind_conc[,2]==input$drug2_Input)[1],]
      #record <- as.numeric(as.character(record))
      num <- 29 - length(which(is.na(record[3:31])))
      conc2 <- record[32:(32+num-1)]
    }

    conc2 <- unname(conc2)
    return(unlist(conc2))
  })
  
  output$mutation_Output <- renderUI({
    selectizeInput("mutation_Input", "show mutation", c(mutations(), "not_selceted"), selected = "not_selected", options = list(maxOptions=2000))
  })
  
  output$tissue_Output <- renderUI({
    selectInput("tissue_Input", "tissue", c(tissues(), "all_tissues"), selected = "all_tissues")
  })
  
  output$target1_Output <- renderUI({
      selectInput("target1_Input", "target 1", c(targets(), "all_targets"), selected = "all_targets")
  })
  
  output$drug1_Output <- renderUI({
    selectInput("drug1_Input", "drug 1", drugs1())
  })
  
  output$conc1_Output <- renderUI({
    selectInput("conc1_Input", "concentr. drug1", concentrations1())
  })
  
  output$target2_Output <- renderUI({
    selectInput("target2_Input", "target 2", c(targets(), "all_targets"), selected = "all_targets")
  })
  
  output$drug2_Output <- renderUI({
    selectInput("drug2_Input", "drug 2", drugs2())
  })
  
  output$conc2_Output <- renderUI({
    selectInput("conc2_Input", "concentr. drug2", concentrations2())
  })

  # stuff for plot and analysis
  gdsc_lines <- unique(table_gdsc_ind_conc[,2])
  ctrp_lines <- unique(table_ctrp_ind_conc[,1])
  
  lines <- reactive({
    if (input$datasetInput=="GDSC") {
      lines <- ann_data$cellInfo_gdsc$Sample.Name[which(ann_data$cellInfo_gdsc$tissueid==input$tissue_Input)]
      if (input$tissue_Input!="all_tissues"){
        lines2 <- intersect(gdsc_lines, lines)
      } else lines2 <- gdsc_lines
      
    }
    if (input$datasetInput=="CTRP") {
      # for matching with raw viabilities data we use "ccl_name" identifiers for cell lines, when we work with AUC data we use "cellid"
      lines <- ann_data$cellInfo_ctrp$ccl_name[which(ann_data$cellInfo_ctrp$tissueid==input$tissue_Input)]
      if (input$tissue_Input!="all_tissues"){
        lines2 <- intersect(ctrp_lines, lines)
      } else lines2 <- ctrp_lines
    }
    return(lines2)
  })
 
  table_xy <- reactive({
    if (input$datasetInput=="GDSC") {
      table0 <- table_gdsc_ind_conc[which(table_gdsc_ind_conc[,2] %in% lines()),]
      table_x <- table0[which(table0[,1]==input$drug1_Input),c(2,(which(concentrations1()==input$conc1_Input)+4))]
      table_y <- table0[which(table0[,1]==input$drug2_Input),c(2,(which(concentrations2()==input$conc2_Input)+4))]
      table <- merge(x = table_x, y = table_y, by = "CELL_LINE_NAME", all = TRUE)
      
      ord_type <- match(table[,1], ann_data$cellInfo_gdsc$Sample.Name)
      ord_mut <- match(table[,1], rownames(table_auc_gdsc1000_pgx_2))
      table <- cbind(table, ann_data$cellInfo_gdsc$GDSC.Tissue.descriptor.2[ord_type], table_auc_gdsc1000_pgx_2[ord_mut,input$mutation_Input])
    }
    
    if (input$datasetInput=="CTRP") {
      table0 <- table_ctrp_ind_conc[which(table_ctrp_ind_conc[,1] %in% lines()),]
      table_x <- table0[which(table0[,2]==input$drug1_Input),c(1,(which(concentrations1()==input$conc1_Input)+2))]
      table_y <- table0[which(table0[,2]==input$drug2_Input),c(1,(which(concentrations2()==input$conc2_Input)+2))]
      table <- merge(x = table_x, y = table_y, by = "cell_line", all = TRUE)
      
      ord_type <- match(table[,1], ann_data$cellInfo_ctrp$ccl_name)
      ord_mut <- match(ann_data$cellInfo_ctrp$cellid[ord_type], rownames(table_auc_ctrp_pgx_2))
      table <- cbind(table, ann_data$cellInfo_ctrp$ccle_hist_subtype_1[ord_type], table_auc_ctrp_pgx_2[ord_mut,input$mutation_Input])
      
    }
    #table <- table[,-1]
    table <- as.data.frame(table)
    table <- na.omit(table)
    colnames(table) <- c("cell_line", "drug1", "drug2", "type", input$mutation_Input)
    table[,2] <- as.numeric(as.character(table[,2]))
    table[,3] <- as.numeric(as.character(table[,3]))
    table[,5] <- as.factor(table[,5])
    #print(head(table))
    return(table)
  })
  
  # for manual clustering
  clust_data <- reactiveValues(
    clust1 = vector(),
    clust2 = vector()
  )
  
  observeEvent(input$cluster1, {
    res <- brushedPoints(table_xy(),input$plot_brush, allRows = FALSE)
    clust_data$clust1 <- res[,1]
  })
  
  observeEvent(input$cluster2, {
    res <- brushedPoints(table_xy(),input$plot_brush, allRows = FALSE)
    clust_data$clust2 <- res[,1]
  })
  
  # km_clusters <- reactive({
  #   clusters <- kmeans(na.omit(table_xy()), 2, nstart = 30)
  # })
  # 
  # diag_clusters <- reactive({
  #   table <- table_xy()
  #   dist <- abs(-table[,1]+table[,2])/sqrt(2)
  #   a <- rownames(table)[which(dist>0.03 & (table[,1]>table[,2]))]
  #   b <- rownames(table)[which(dist>0.03 & (table[,2]>table[,1]))]
  #   cluster <- rep(1, nrow(table))
  #   cluster[which(table[,1]>table[,2])] <- 2
  #   return(list(a=a, b=b, cluster=cluster))
  # })
  
  t_test <- eventReactive(input$showres_Input, {

      if (input$datasetInput=="GDSC") {
      table <- table_auc_gdsc1000_pgx_2
      a <- as.character(clust_data$clust1)
      b <- as.character(clust_data$clust2)
      #lines <- lines()
    }
      
      if (input$datasetInput=="CTRP") {
        table <- table_auc_ctrp_pgx_2
        a <- ann_data$cellInfo_ctrp$cellid[match(as.character(clust_data$clust1), ann_data$cellInfo_ctrp$ccl_name)]
        b <- ann_data$cellInfo_ctrp$cellid[match(as.character(clust_data$clust2), ann_data$cellInfo_ctrp$ccl_name)]
        #lines <- ann_data$cellInfo_ctrp$cellid[match(lines(), ann_data$cellInfo_ctrp$ccl_name)]
      }
      print(a)
      print(b)
 
      # mol types
      prefixes <- paste0(input$moltypes_Input, "_")
      features <- vector()
      for (i in 1:length(prefixes)){
        features <- c(features, grep(prefixes[i], colnames(table), value=F))
      }
      

      #gen_last <- max(grep("_", colnames(table)))
    
      
      tt_res <- apply(table[c(a,b),features],2, function(x){
        tt <- t.test(x[a], x[b])
        c(tt$p.value, tt$estimate)})
      
      tt_res <- t(tt_res)
      
      res_table <- cbind(rownames(tt_res), tt_res)
      res_table <- res_table[order(res_table[,2]),]
      
      num <- which(res_table[,2]<input$pval_Input)
      if(length(num)>100) {num <- num[1:100]}
      res_table <- res_table[num,]
      #return(res_table) 

      return(res_table)

  })
  
  output$coolplot <- renderPlot({
    if (is.null(table_xy())) {
      return()}
    data0 <- table_xy()
    limits <- c(min(data0$drug1), max(data0$drug1), min(data0$drug2), max(data0$drug2))
    cl1 <- which(data0[,1] %in% clust_data$clust1)
    cl2 <- which(data0[,1] %in% clust_data$clust2)
    cl12 <- c(cl1,cl2)
    print(length(cl12))
    if (length(cl12)>0){
    data <- data0[-cl12,] } else data <- data0
    
    e <- ggplot(data, aes(drug1, drug2))
    e <- e + geom_point(aes(colour=factor(type), size=factor(data[,5]))) + labs(x = input$drug1_Input, y=input$drug2_Input) + xlim(limits[1], limits[2]) + ylim(limits[3], limits[4])
    if (length(cl1)>0){ e <- e + geom_point(data=data0[cl1,], color="black") }
    if (length(cl2)>0){ e <- e + geom_point(data=data0[cl2,], color="navyblue") }
    e
    
    # if (input$clustInput=="diagonal") {
    # plot(table_xy(), col=diag_clusters()$cluster, xlab=input$drug1_Input, ylab=input$drug2_Input, cex=1.5, cex.lab=1.5, xlim=c(0,100), ylim=c(0,100))
    # }
    # 
    # if (input$clustInput=="k-means") {
    #   plot(na.omit(table_xy()), col = km_clusters()$cluster,xlab=input$drug1_Input, ylab=input$drug2_Input, cex=1.5, cex.lab=1.5, xlim=c(0,100), ylim=c(0,100))
    # }
    #   abline(coef=c(0,1))
  })
  
  output$results <- renderTable({
    input$showres_Input
    isolate({
    res <- t_test()
    colnames(res) <- c("feature", "p-value", "mean cl. 1", "mean cl. 2")
    return(res) })
  })

}
shinyApp(ui = ui, server = server)