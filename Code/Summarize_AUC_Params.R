library(ROCR)
library(ggplot2)
library(Metrics)
library(hash)

# Methods = c("M1", "M2", "M3", "M12", "M13", "M23", "M123")
# Methods = c("H1","H2","H3","MH12","MH13","MH23","MH123")
Method = "MH123"
DrugSimNet = "PREDICT"

setwd("~/Manuscripts/99HDR2/Code")

Result_bySLFile = paste0("../Results/",Method,"_byTrial_",DrugSimNet,"_Score_n_Label_KFold10.csv")
df.ScoreLabel = read.csv(Result_bySLFile)
Scores = df.ScoreLabel$Scores
Labels = df.ScoreLabel$Labels

resultspred = prediction(Scores, Labels)
auc.perf = performance(resultspred, measure = "auc")
aucavgbyAll = round(auc.perf@y.values[[1]],3)

# roc.perf <- performance(resultspred,"tpr","fpr")
# plot(roc.perf@x.values[[1]], roc.perf@y.values[[1]], xlab="FPR", ylab="TPR", main = as.character(aucavgbyAll), sub = "Hello")


Result_byTrialFile = paste0("../Results/",Method,"_byTrial_",DrugSimNet,"_ROC_KFold10.csv")
df.AUCbyTrial = read.csv(Result_byTrialFile)

aucavgbyTrial = round(mean(df.AUCbyTrial$auc),3)
aucavgbyTrial.sd = round(sd(df.AUCbyTrial$auc),3)

# Result_byDiseaseFile = paste0("../Results/",Method,"_byDisease_",DrugSimNet,"_ROC.csv")
# df.AUCbyDisease = read.csv(Result_byDiseaseFile)

df.AUCbyDisease <- aggregate(auc~disease, data=df.AUCbyTrial, FUN=function(x) c(mean=mean(x), count=length(x)))

aucavgbyDisease = round(mean(as.numeric(df.AUCbyDisease$auc[,1])),3)
aucavgbyDisease.sd = round(sd(as.numeric(df.AUCbyDisease$auc[,1])),3)


cat("Method=",Method,'DrugSimNet=',DrugSimNet,'aucavgbyAll=',aucavgbyAll,'aucavgbyTrial=',aucavgbyTrial,"(+-",aucavgbyTrial.sd,")","aucavgbyDisease=",aucavgbyDisease,"(+-",aucavgbyDisease.sd,")\n")

####################
library(ROCR)
library(ggplot2)
library(Metrics)
library(hash)
library('cowplot')
setwd("~/Manuscripts/99HDR2/Code")


Methods = c("M123", "MH123")
DrugSimNets = c("PREDICT")
DrugSimNet = "PREDICT"

Params = c("r","delta")#,"delta"
h.config2dfAUC = hash()
for(Param in Params){
  config = paste0(DrugSimNet,"_",Param)
  df.AUC = NULL
  for(Method in Methods){
    for(pv in c(0.1, 0.3, 0.5, 0.7, 0.9)){
      cat("Analyzing for", Param, Method, "\n")
      
      Result_bySLFile = paste0("../Results/",Method,"_byTrial_",DrugSimNet,"_",Param,pv,"_Score_n_Label.csv")
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
      
      if(Method == "MH123"){
        NetworkType = paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_OMIM_HPO_GeneNet]")
      }else{
        NetworkType = "DiSimNet_OMIM_HPO_GeneNet"  
      }
      df.AUC = rbind(df.AUC, data.frame(param = pv, AU = aucavgbyAll, NetworkType = NetworkType))
    }
  }
  #Store
  SummaryFile = paste0("../Results/Summary_",paste0(Methods,collapse = "-"),"_",config,".rdata")
  saveRDS(df.AUC,SummaryFile)#readRDS
  
  h.config2dfAUC[[config]] = df.AUC
}
length(h.config2dfAUC)

h.config2dfAUC[["PREDICT_r"]]
h.config2dfAUC[["PREDICT_delta"]]

# # #Load:
# Methods = c("M123", "MH123")
# 
# DrugSimNet = "PREDICT"
# 
# Params = c("r","delta")
# Param = "r"
# config = paste0(DrugSimNet,"_",Param)
# 
# setwd("~/Manuscripts/99HDR2/Code")
# DrugSimNetsSTR = paste0(DrugSimNets,collapse = "-")
# SummaryFile = paste0("../Results/Summary_",paste0(Methods,collapse = "-"),"_",config,".rdata")
# 
# df.AUC = readRDS(SummaryFile)
# h.config2dfAUC = hash()
# h.config2dfAUC[[config]] = df.AUC

######


#Draw
h.config2AUplot = hash()
for(Param in Params){
  for(DrugSimNet in DrugSimNets){
    config = paste0(DrugSimNet,"_",Param)
    df.AU = h.config2dfAUC[[config]]
    
    df.AU$param = as.factor(df.AU$param)
    
    if (Param=="r"){
      ParamSymbol = "restart-probability"
    }else{
      ParamSymbol = "jumping-probability"
    }
    
    fontsize = 14
    p = ggplot(data=df.AU, aes(x=param, y=AU, group = NetworkType)) +
      geom_line(stat="identity",position=position_dodge(),aes(color=NetworkType)) +
      geom_point(aes(color=NetworkType)) +
      ylim(0,1)
    
    p = p + labs(x=paste0("\n",ParamSymbol), y = "AUC\n") +
      scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
      # scale_fill_brewer(palette="Blues") +
      theme_light() +
      theme(
        text = element_text(size=fontsize),# All font sizes
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=fontsize-1),
        legend.title=element_blank(),#Remove legend title (Network)
        legend.position = c(0.5, 0.11),
        axis.text = element_text(size = fontsize),
        axis.title = element_text(size = fontsize)
        # axis.title.x = element_text(size = fontsize, family="Arial")
      )
    
    h.config2AUplot[[config]] = p
  }
}

library('cowplot')
plot_grid(h.config2AUplot[["PREDICT_r"]],h.config2AUplot[["PREDICT_delta"]], labels=c("(a)", "(b)"), ncol = 2, nrow = 1)

FigureFile = paste0("../Figures/Figure_Parameter.pdf")
ggsave(FigureFile, width = 14, height = 7)
# embed_fonts(FigureFile)
