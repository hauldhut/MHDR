#v2: Change e.g., DiSimNet_OMIM to DiSimNetO, with O is subscript

####################
library(ROCR)
library(ggplot2)
library(Metrics)
library(hash)
library('cowplot')
setwd("~/Manuscripts/99MHDR/Code")

h.All = hash()

# Methods = c("M1", "M1_10", "M1_15", "M1_03") #For reviewer's comment
# Methods = c("M1_10", "M1_102", "M1_103", "M1_1023") #For reviewer's comment
# Methods = c("M1_15", "M1_152", "M1_153", "M1_1523") #For reviewer's comment
# Methods = c("M1_03", "M1_032", "M1_033", "M1_0323") #For reviewer's comment
# Methods = c("M123", "M1_1023", "M1_1523", "M1_0323") #For reviewer's comment

# Methods = c("M1", "M2", "M3", "M12", "M13", "M23", "M123") #Figure 3
# DrugSimNets = c("PREDICT")


# Methods = c("H1_10","MH1_102","MH1_103","MH1_1023") #For reviewer's comment 
# Methods = c("H1_15","MH1_152","MH1_153","MH1_1523") #For reviewer's comment
# Methods = c("H1_03","MH1_032","MH1_033","MH1_0323") #For reviewer's comment
# Methods = c("H1","H1_10","H1_15","H1_03") #For reviewer's comment --> S3
Methods = c("MH123","MH1_1023","MH1_1523","MH1_0323") #For reviewer's comment --> S4

# Methods = c("H1","H2","H3","MH12","MH13","MH23","MH123") #Figure 4
DrugSimNets = c("PREDICT")#PREDICT/CHEM

gamma = 0.5
delta = 0.5

# Methods = c("H1","H2","H3")
h.cap2AUPRC = hash()
h.cap2dis = hash()
for(DrugSimNet in DrugSimNets){
  cat("Analyzing for ", DrugSimNet, "\n")
  
  df.All = NULL
  
  for(Method in Methods){
    cat("\t==> Analyzing for ", Method, "\n")
    Result_by_Score_n_LabelFile = paste0("../Results/",Method,"_byTrial_",DrugSimNet,"_gamma_",gamma,"_delta_",delta,"_Score_n_Label.csv")
    print(Result_by_Score_n_LabelFile)
    df.ScoreLabel = read.csv(Result_by_Score_n_LabelFile)
    Scores = df.ScoreLabel$Scores
    Labels = df.ScoreLabel$Labels
    resultspred = prediction(Scores, Labels)
    # Calculate performance for Precision-Recall curve
    pr.perf <- performance(resultspred, "prec", "rec")
    
    # Extract Recall (x-axis) and Precision (y-axis) values for plotting
    Recall = pr.perf@x.values[[1]]    # Recall values
    Precision = pr.perf@y.values[[1]] # Precision values
    Precision[is.nan(Precision)] <- 0
    Precision[is.na(Precision)] <- 0
    
    
    # Calculate AUPRC (Area Under Precision-Recall Curve)
    aupr.perf = performance(resultspred, measure = "aucpr") # 'aucpr' is the measure for AUPR
    aupravgbyAll = round(aupr.perf@y.values[[1]], 3)       # Round to 3 decimal places
    
    if(DrugSimNet=="CHEM"){
      DrSNet = "C"
    }else{
      DrSNet = "P"
    }
    
      
    if(Method == "H1"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetO")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[O] * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
      # Network = rep(paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_OMIM] (AUC = ",aupravgbyAll,")"),length(Recall))  
    }else if(Method == "H2"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetH")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[H] * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
      # Network = rep(paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_HPO] (AUC = ",aupravgbyAll,")"),length(Recall))  
    }else if(Method == "H3"){ 
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetG")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[G] * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
      # Network = rep(paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_GeneNet] (AUC = ",aupravgbyAll,")"),length(Recall)) 
    }else if(Method == "MH12"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetOH")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[OH] * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
      # Network = rep(paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_OMIM_HPO] (AUC = ",aupravgbyAll,")"),length(Recall))  
    }else if(Method == "MH13"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetOG")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[OG] * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
      # Network = rep(paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_OMIM_GeneNet] (AUC = ",aupravgbyAll,")"),length(Recall))  
    }else if(Method == "MH23"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetHG")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[HG] * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
      # Network = rep(paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_HPO_GeneNet] (AUC = ",aupravgbyAll,")"),length(Recall))  
    }else if(Method == "MH123"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetOHG")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[OHG] * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
      # Network = rep(paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_OMIM_HPO_GeneNet] (AUC = ",aupravgbyAll,")"),length(Recall))
    }else if(Method == "M1"){
      cap = "DiSimNetO"
      dis = bquote("DiSimNet"[O] * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "M2"){
      cap = "DiSimNetH"
      dis = bquote("DiSimNet"[H] *  " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))
    }else if(Method == "M3"){
      cap = "DiSimNetG"
      dis = bquote("DiSimNet"[G] * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "M12"){
      cap = "DiSimNetOH"
      dis = bquote("DiSimNet"[OH] * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "M13"){
      cap = "DiSimNetOG"
      dis = bquote("DiSimNet"[OG] * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "M23"){
      cap = "DiSimNetHG"
      dis = bquote("DiSimNet"[HG] * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "M123"){
      cap = "DiSimNetOHG"
      dis = bquote("DiSimNet"[OHG] * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
      # Network = rep(paste0("DiSimNet_OMIM_HPO_GeneNet (AUC = ",aupravgbyAll,")"),length(Recall))  
    }else if(Method == "M1_10"){
      cap = "DiSimNetO (kLN=10)"
      dis = bquote("DiSimNet"[O] * " (kLN=10) " * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "M1_15"){
      cap = "DiSimNetO (kLN=15)"
      dis = bquote("DiSimNet"[O] * " (kLN=15) " * " (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "M1_03"){
      cap = "DiSimNetO (sim>=0.3)"
      dis = bquote("DiSimNet"[O] * " (sim>=0.3) " * "(AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUPRC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }
    
    #kLN=10
    if(Method == "M1_102"){
      cap = "DiSimNetOH (kLN=10)"
      dis = bquote("DiSimNet"[OH] * " (kLN=10) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "M1_103"){
      cap = "DiSimNetOG (kLN=10)"
      dis = bquote("DiSimNet"[OG] * " (kLN=10) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "M1_1023"){
      cap = "DiSimNetOHG (kLN=10)"
      dis = bquote("DiSimNet"[OHG] * " (kLN=10) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }
    
    #kLN=15
    if(Method == "M1_152"){
      cap = "DiSimNetOH (kLN=15)"
      dis = bquote("DiSimNet"[OH] * " (kLN=15) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "M1_153"){
      cap = "DiSimNetOG (kLN=15)"
      dis = bquote("DiSimNet"[OG] * " (kLN=15) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "M1_1523"){
      cap = "DiSimNetOHG (kLN=15)"
      dis = bquote("DiSimNet"[OHG] * " (kLN=15) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }
    
    #sim>=0.3
    if(Method == "M1_032"){
      cap = "DiSimNetOH (sim>=0.3)"
      dis = bquote("DiSimNet"[OH] * " (sim>=0.3) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "M1_033"){
      cap = "DiSimNetOG (sim>=0.3)"
      dis = bquote("DiSimNet"[OG] * " (sim>=0.3) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "M1_0323"){
      cap = "DiSimNetOHG (sim>=0.3)"
      dis = bquote("DiSimNet"[OHG] * " (sim>=0.3) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }
    
    #kLN=10
    if(Method == "H1_10"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetO (kLN=10)")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[O] * " (kLN=10) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "MH1_102"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetOH (kLN=10)")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[OH] * " (kLN=10) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "MH1_103"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetOG (kLN=10)")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[OG] * " (kLN=10) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "MH1_1023"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetOHG (kLN=10)")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[OHG] * " (kLN=10) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }
    
    #kLN=15
    if(Method == "H1_15"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetO (kLN=15)")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[O] * " (kLN=15) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "MH1_152"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetOH (kLN=15)")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[OH] * " (kLN=15) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "MH1_153"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetOG (kLN=15)")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[OG] * " (kLN=15) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "MH1_1523"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetOHG (kLN=15)")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[OHG] * " (kLN=15) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }
    
    #sim>=0.3
    #Methods = c("H1_03","MH1_032","MH1_033","MH1_0323") #For reviewer's comment 
    if(Method == "H1_03"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetO (sim>=0.3)")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[O] * " (sim>=0.3) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "MH1_032"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetOH (sim>=0.3)")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[OH] * " (sim>=0.3) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "MH1_033"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetOG (sim>=0.3)")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[OG] * " (sim>=0.3) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }else if(Method == "MH1_0323"){
      cap = paste0("DrSimNet",DrSNet,"-DiSimNetOHG (sim>=0.3)")
      dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[OHG] * " (sim>=0.3) (AUPRC = " * .(aupravgbyAll) * ")")
      h.cap2AUROC[[cap]] = aupravgbyAll
      h.cap2dis[[cap]] = dis
      Network = rep(cap,length(Recall))  
    }
    
    df.All = rbind(df.All, data.frame(Recall = Recall, Precision = Precision, Network = Network))
  }
  h.All[[DrugSimNet]] = df.All
}
length(h.All)

#Store
l.All = c(h.All, h.cap2AUPRC, h.cap2dis)
SummaryFile = paste0("../Results/Summary_AUPRC_",paste0(Methods,collapse = "-"),"_",paste0(DrugSimNets,collapse = "-"),"_gamma_",gamma,"_delta_",delta,".rdata")
saveRDS(l.All,SummaryFile)#readRDS


# #Load:
# # Methods = c("H1","H2","H3","MH12","MH13","MH23","MH123")
# # Methods = c("M1", "M2", "M3", "M12", "M13", "M23", "M123")
# # Methods = c("M1", "M1_10", "M1_15", "M1_03") #For reviewer's comment
# Methods = c("H1_03","MH1_032","MH1_033","MH1_0323") #For reviewer's comment
# 
# DrugSimNets = c("PREDICT")#PREDICT/CHEM
# #
# setwd("~/Manuscripts/99MHDR/Code")
# # DrugSimNetsSTR = paste0(DrugSimNets,collapse = "-")
# SummaryFile = paste0("../Results/Summary_AUPRC_",paste0(Methods,collapse = "-"),"_",paste0(DrugSimNets,collapse = "-"),"_gamma_",gamma,"_delta_",delta,".rdata")
# 
# l.All = readRDS(SummaryFile)
# h.All = l.All[[1]]
# h.cap2AUPRC = l.All[[2]]
# h.cap2dis = l.All[[3]]


######
FigureFile = paste0("../Figures/Figure_AUPRC_",paste0(Methods,collapse = "-"))

size = 15
l.plot = list()
pi=0
for(DrugSimNet in keys(h.All)){
  pi = pi+1
  print(paste0("Processing plot for ", DrugSimNet))
  df.O = h.All[[DrugSimNet]]
  df.O = df.O[(df.O$Recall!=0) & (df.O$Precision != 0),]
  
  if(nrow(df.O)>=10^6){
    df.D = df.O[sample(1:nrow(df.O),nrow(df.O)/100),]  
  }else{
    df.D = df.O[sample(1:nrow(df.O),nrow(df.O)/10),]
  }
  
  df.D.new = df.D 
  # Get Network levels in alphabetical order
  sorted_network_levels = sort(unique(h.All[[DrugSimNet]]$Network))
  
  # Set Network as factor with sorted levels
  df.D.new$Network <- factor(df.D.new$Network, levels = sorted_network_levels)
  
  # Create named list of captions in sorted order
  caption_list <- setNames(
    lapply(sorted_network_levels, function(cap) h.cap2dis[[cap]]),
    sorted_network_levels
  )
  
  # Extract captions in sorted order
  labels <- caption_list[sorted_network_levels]
  
  # Debugging: Print to verify mapping
  print("Network levels and their captions (alphabetically ordered):")
  print(data.frame(Network = sorted_network_levels, Caption = I(labels)))
  
  p <- ggplot(df.D.new, aes(x = Recall, y = Precision, shape = Network, color = Network)) +
    geom_line(size = 1) +
    xlab("\nRecall") +
    ylab("Precision\n") +
    theme_light() +
    theme(
      text = element_text(size = size),
      plot.title = element_text(hjust = 0.5),
      legend.text = element_text(size = size),
      legend.position = c(0.55, 0.8),
      legend.title = element_blank(),
      axis.text = element_text(size = size),
      axis.title = element_text(size = size)
    ) +
    scale_color_discrete(limits = sorted_network_levels, labels = labels) +
    scale_shape_discrete(limits = sorted_network_levels, labels = labels)  # Ensure shape aligns too
  
  l.plot[[pi]] = p
  saveRDS(p, paste0(FigureFile, "_", DrugSimNet, "_gamma_", gamma, "_delta_", delta, ".rdata"))
}

library(cowplot)
l.plot[[1]]
# plot_grid(l.plot[[1]], l.plot[[2]], labels=c("A", "B"), ncol = 2, nrow = 1)

ggsave(paste0(FigureFile,"_",paste0(DrugSimNets,collapse = "-"),"_gamma_",gamma,"_delta_",delta,".pdf"), width = 7, height = 7)





