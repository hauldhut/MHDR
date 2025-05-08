setwd("~/Manuscripts/99MHDR/Code")

####################
library(ROCR)
library(ggplot2)
library(Metrics)
library(hash)
library('cowplot')
setwd("~/Manuscripts/99MHDR/Code")


Methods = c("M123", "MH123")
DrugSimNets = c("PREDICT")
DrugSimNet = "PREDICT"

Params = c("gamma","delta")#,"delta"
h.config2dfAUC = hash()
h.cap2dis = hash()
for(Param in Params){
  config = paste0(DrugSimNet,"_",Param)
  df.AUC = NULL
  for(Method in Methods){
    for(pv in c(0.1, 0.3, 0.5, 0.7, 0.9)){
      cat("Analyzing for", Method, Param,"=",pv, "\n")
      if(Param == "gamma"){
        Result_bySLFile = paste0("../Results/",Method,"_byTrial_",DrugSimNet,"_gamma_",pv,"_delta_0.5_Score_n_Label.csv")  
      }else{
        Result_bySLFile = paste0("../Results/",Method,"_byTrial_",DrugSimNet,"_gamma_0.5_delta_",pv,"_Score_n_Label.csv")
      }
      
      print(Result_bySLFile)
      df.ScoreLabel = read.csv(Result_bySLFile)
      Scores = df.ScoreLabel$Scores
      Labels = df.ScoreLabel$Labels
      resultspred = prediction(Scores, Labels)
      roc.perf <- performance(resultspred,"tpr","fpr")
      FPR = roc.perf@x.values[[1]]
      TPR = roc.perf@y.values[[1]]
      
      auc.perf = performance(resultspred, measure = "auc")
      aucavgbyAll = round(auc.perf@y.values[[1]],3)
      
      if(DrugSimNet=="CHEM"){
        DrSNet = "C"
      }else{
        DrSNet = "P"
      }
      
      if(Method == "MH123"){
        cap = paste0("DrSimNet",DrSNet,"-DiSimNetOHG")
        dis = bquote("DrSimNet"[.(DrSNet)] * "-DiSimNet"[OHG])
        h.cap2dis[[cap]] = dis
      }else{
        cap = "DiSimNetOHG"
        dis = bquote("DiSimNet"[OHG])
        h.cap2dis[[cap]] = dis
      }
      df.AUC = rbind(df.AUC, data.frame(param = pv, AU = aucavgbyAll, Network = cap))
    }
  }
  #Store
  SummaryFile = paste0("../Results/Summary_AUROC_",paste0(Methods,collapse = "-"),"_",config,".rdata")
  saveRDS(df.AUC,SummaryFile)#readRDS
  
  h.config2dfAUC[[config]] = df.AUC
}
length(h.config2dfAUC)

h.config2dfAUC[["PREDICT_gamma"]]
h.config2dfAUC[["PREDICT_delta"]]

# # #Load:
# Methods = c("M123", "MH123")
# DrugSimNets = c("PREDICT")
# DrugSimNet = "PREDICT"
# 
# Params = c("gamma","delta")
# Param = "gamma"
# setwd("~/Manuscripts/99MHDR/Code")
# DrugSimNetsSTR = paste0(DrugSimNets,collapse = "-")
# 
# h.config2dfAUC = hash()
# for(Param in Params){
#   config = paste0(DrugSimNet,"_",Param)
#   SummaryFile = paste0("../Results/Summary_AUROC_",paste0(Methods,collapse = "-"),"_",config,".rdata")
#   print(SummaryFile)
# 
#   df.AUC = readRDS(SummaryFile)
# 
#   h.config2dfAUC[[config]] = df.AUC
# 
# }
# 
# h.cap2dis = hash()
# cap = paste0("DrSimNetP-DiSimNetOHG")
# dis = bquote("DrSimNet"[P] * "-DiSimNet"[OHG])
# h.cap2dis[[cap]] = dis
# 
# cap = "DiSimNetOHG"
# dis = bquote("DiSimNet"[OHG])
# h.cap2dis[[cap]] = dis
# ######


#Draw
h.config2AUplot = hash()
for (Param in Params) {
  for (DrugSimNet in DrugSimNets) {
    config = paste0(DrugSimNet, "_", Param)
    df.AU = h.config2dfAUC[[config]]
    
    # Ensure param is a factor with sorted levels
    df.AU$param = factor(df.AU$param, levels = sort(unique(df.AU$param)))
    
    # Set ParamSymbol using expression
    if (Param == "gamma") {
      ParamSymbol = expression("restart-probability (" * gamma * ")")
    } else {
      ParamSymbol = expression("jumping-probability (" * delta * ")")
    }
    
    # Get Network levels in alphabetical order
    sorted_network_levels = sort(unique(df.AU$Network))
    
    # Set Network as factor with sorted levels
    df.AU$Network <- factor(df.AU$Network, levels = sorted_network_levels)
    
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
    
    fontsize = 15
    p <- ggplot(data = df.AU, aes(x = param, y = AU, group = Network)) +
      geom_line(stat = "identity", position = position_dodge(width = 0.2), aes(color = Network)) +
      geom_point(aes(color = Network), position = position_dodge(width = 0.2)) +
      ylim(0, 1) +
      xlab(ParamSymbol) +
      ylab("AUROC\n") +
      theme_light() +
      theme(
        text = element_text(size = fontsize),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = fontsize),
        legend.title = element_blank(),
        legend.position = c(0.5, 0.11),
        axis.text = element_text(size = fontsize),
        axis.title = element_text(size = fontsize)
      ) +
      scale_color_discrete(limits = sorted_network_levels, labels = labels) +
      scale_shape_discrete(limits = sorted_network_levels, labels = labels)
    
    h.config2AUplot[[config]] = p
  }
}

library('cowplot')
plot_grid(h.config2AUplot[["PREDICT_gamma"]], h.config2AUplot[["PREDICT_delta"]], labels = c("(a)", "(b)"), ncol = 2, nrow = 1)

FigureFile = paste0("../Figures/Figure2.pdf")
ggsave(FigureFile, width = 14, height = 5)