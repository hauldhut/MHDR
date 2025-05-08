library(ROCR)
library(ggplot2)
library(Metrics)
library(hash)

# Methods = c("M1", "M2", "M3", "M12", "M13", "M23", "M123")
# Methods = c("H1","H2","H3","MH12","MH13","MH23","MH123")
Method = "M1"
DrugSimNet = "PREDICT"

setwd("~/Manuscripts/99MHDR/Code")

Result_by_Score_n_LabelFile = paste0("../Results/MH123_byTrial_PREDICT_Score_n_Label_KFold10.csv")
df.ScoreLabel = read.csv(Result_by_Score_n_LabelFile)
Scores = df.ScoreLabel$Scores
Labels = df.ScoreLabel$Labels

resultspred = prediction(Scores, Labels)
auc.perf = performance(resultspred, measure = "auc")
aucavgbyAll = round(auc.perf@y.values[[1]],3)

# roc.perf <- performance(resultspred,"tpr","fpr")
# plot(roc.perf@x.values[[1]], roc.perf@y.values[[1]], xlab="FPR", ylab="TPR", main = as.character(aucavgbyAll), sub = "Hello")

# Extract Recall (x-axis) and Precision (y-axis) values for plotting
pr.perf <- performance(resultspred, "prec", "rec")
# Recall = pr.perf@x.values[[1]]    # Recall values
# Precision = pr.perf@y.values[[1]] # Precision values
# Precision[is.nan(Precision)] <- 0
# Precision[is.na(Precision)] <- 0


# Calculate AUPRC (Area Under Precision-Recall Curve)
aupr.perf = performance(resultspred, measure = "aucpr") # 'aucpr' is the measure for AUPR
aupravgbyAll = round(aupr.perf@y.values[[1]], 3)       # Round to 3 decimal places

cat("aucavgbyAll",aucavgbyAll,"aupravgbyAll",aupravgbyAll,"\n")

# Result_byTrialFile = paste0("../Results/",Method,"_byTrial_",DrugSimNet,"_ROC.csv")
# df.AUCbyTrial = read.csv(Result_byTrialFile)
# 
# aucavgbyTrial = round(mean(df.AUCbyTrial$auc),3)
# aucavgbyTrial.sd = round(sd(df.AUCbyTrial$auc),3)
# 
# # Result_byDiseaseFile = paste0("../Results/",Method,"_byDisease_",DrugSimNet,"_ROC.csv")
# # df.AUCbyDisease = read.csv(Result_byDiseaseFile)
# 
# df.AUCbyDisease <- aggregate(auc~disease, data=df.AUCbyTrial, FUN=function(x) c(mean=mean(x), count=length(x)))
# 
# aucavgbyDisease = round(mean(as.numeric(df.AUCbyDisease$auc[,1])),3)
# aucavgbyDisease.sd = round(sd(as.numeric(df.AUCbyDisease$auc[,1])),3)
# 
# 
# cat("Method=",Method,'DrugSimNet=',DrugSimNet,'aucavgbyAll=',aucavgbyAll,'aucavgbyTrial=',aucavgbyTrial,"(+-",aucavgbyTrial.sd,")","aucavgbyDisease=",aucavgbyDisease,"(+-",aucavgbyDisease.sd,")\n")

# ####################
# library(ROCR)
# library(ggplot2)
# library(Metrics)
# library(hash)
# library('cowplot')
# setwd("~/Manuscripts/99MHDR/Code")
# 
# h.All = hash()
# 
# Methods = c("M1", "M2", "M3", "M12", "M13", "M23", "M123")
# DrugSimNets = c("PREDICT")
# 
# # Methods = c("H1","H2","H3","MH12","MH13","MH23","MH123")
# # DrugSimNets = c("PREDICT")#PREDICT/CHEM
# 
# gamma = 0.5
# delta = 0.5
# 
# # Methods = c("H1","H2","H3")
# for(DrugSimNet in DrugSimNets){
#   cat("Analyzing for ", DrugSimNet, "\n")
#   
#   df.All = NULL
#   
#   for(Method in Methods){
#     cat("\t==> Analyzing for ", Method, "\n")
#     Result_by_Score_n_LabelFile = paste0("../Results/",Method,"_byTrial_",DrugSimNet,"_gamma_",gamma,"_delta_",delta,"_Score_n_Label.csv")
#     print(Result_by_Score_n_LabelFile)
#     
#     df.ScoreLabel = read.csv(Result_by_Score_n_LabelFile)
#     Scores = df.ScoreLabel$Scores
#     Labels = df.ScoreLabel$Labels
#     resultspred = prediction(Scores, Labels)
#     roc.perf <- performance(resultspred,"tpr","fpr")
#     FPR = roc.perf@x.values[[1]]
#     TPR = roc.perf@y.values[[1]]
#     
#     auc.perf = performance(resultspred, measure = "auc")
#     aucavgbyAll = round(auc.perf@y.values[[1]],3)
#     
#     if(Method == "MH123"){
#       Network = rep(paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_OMIM_HPO_GeneNet] (AUC = ",aucavgbyAll,")"),length(FPR))  
#     }else if(Method == "MH12"){
#       Network = rep(paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_OMIM_HPO] (AUC = ",aucavgbyAll,")"),length(FPR))  
#     }else if(Method == "MH13"){
#       Network = rep(paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_OMIM_GeneNet] (AUC = ",aucavgbyAll,")"),length(FPR))  
#     }else if(Method == "MH23"){
#       Network = rep(paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_HPO_GeneNet] (AUC = ",aucavgbyAll,")"),length(FPR))  
#     }else if(Method == "H1"){
#       Network = rep(paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_OMIM] (AUC = ",aucavgbyAll,")"),length(FPR))  
#     }else if(Method == "H2"){
#       Network = rep(paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_HPO] (AUC = ",aucavgbyAll,")"),length(FPR))  
#     }else if(Method == "H3"){ 
#       Network = rep(paste0("[DrSimNet_",DrugSimNet,"]-[DiSimNet_GeneNet] (AUC = ",aucavgbyAll,")"),length(FPR))  
#     }else if(Method == "M123"){
#       Network = rep(paste0("DiSimNet_OMIM_HPO_GeneNet (AUC = ",aucavgbyAll,")"),length(FPR))  
#     }else if(Method == "M12"){
#       Network = rep(paste0("DiSimNet_OMIM_HPO (AUC = ",aucavgbyAll,")"),length(FPR))  
#     }else if(Method == "M13"){
#       Network = rep(paste0("DiSimNet_OMIM_GeneNet (AUC = ",aucavgbyAll,")"),length(FPR))  
#     }else if(Method == "M23"){
#       Network = rep(paste0("DiSimNet_HPO_GeneNet (AUC = ",aucavgbyAll,")"),length(FPR))  
#     }else if(Method == "M1"){
#       Network = rep(paste0("DiSimNet_OMIM (AUC = ",aucavgbyAll,")"),length(FPR))  
#     }else if(Method == "M2"){
#       Network = rep(paste0("DiSimNet_HPO (AUC = ",aucavgbyAll,")"),length(FPR))  
#     }else{ 
#       Network = rep(paste0("DiSimNet_GeneNet (AUC = ",aucavgbyAll,")"),length(FPR))  
#     }
#     
#     df.All = rbind(df.All, data.frame(FPR = FPR, TPR = TPR, Network = Network))
#   }
#   h.All[[DrugSimNet]] = df.All
# }
# length(h.All)
# 
# #Store
# SummaryFile = paste0("../Results/Summary_",paste0(Methods,collapse = "-"),"_",paste0(DrugSimNets,collapse = "-"),"_gamma_",gamma,"_delta_",delta,".rdata")
# saveRDS(h.All,SummaryFile)#readRDS
# 
# # # #Load: 
# # # Methods = c("H1","H2","H3","MH12","MH13","MH23","MH123")
# # Methods = c("M1", "M2", "M3", "M12", "M13", "M23", "M123")
# # DrugSimNets = c("PREDICT")#PREDICT/CHEM
# # 
# # setwd("~/Manuscripts/99HDR2/Code")
# # DrugSimNetsSTR = paste0(DrugSimNets,collapse = "-")
# # SummaryFile = paste0("../Results/Summary_",paste0(Methods,collapse = "-"),"_",DrugSimNetsSTR,".rdata")
# # 
# # h.All = readRDS(SummaryFile)
# 
# # h.All[[DrugSimNetsSTR]]$Network = gsub("Disease", "Di", h.All[[DrugSimNetsSTR]]$Network)
# # h.All[[DrugSimNetsSTR]]$Network = gsub("Drug", "Dr", h.All[[DrugSimNetsSTR]]$Network)
# 
# 
# # Legendvec = unique(h.All[["CHEM"]]$Network)
# 
# ######
# FigureFile = paste0("../Figures/Figure_",paste0(Methods,collapse = "-"))
# 
# size = 12
# l.plot = list()
# pi=0
# for(DrugSimNet in keys(h.All)){
#   pi = pi+1
#   print(paste0("Processing plot for ", DrugSimNet))
#   df.O = h.All[[DrugSimNet]]
#   if(nrow(df.O)>=10^6){
#     df.D = df.O[sample(1:nrow(df.O),nrow(df.O)/100),]  
#   }else{
#     df.D = df.O[sample(1:nrow(df.O),nrow(df.O)/10),]
#   }
#   
#   df.D.new = df.D 
#   df.D.new$Network <- factor(df.D.new$Network, unique(h.All[[DrugSimNet]]$Network))
#   
#   p <- ggplot(df.D.new, aes(x=FPR, y=TPR, shape=Network, color = Network)) +
#     # geom_point() +
#     geom_line(size=1) +
#     # ggtitle(DrugSimNet) + #Set plot title
#     xlab("\nFalse Positive Rate") + #Set x label
#     ylab("True Positive Rate\n") + #Set y label
#     
#     theme_light() +
#     theme(
#       # panel.background = element_rect(fill = "lightgray",colour = "lightgray",linetype = "solid"),
#       text = element_text(size=size),# All font sizes
#       plot.title = element_text(hjust = 0.5),
#       legend.text = element_text(size=size-1),
#       legend.position = c(0.5, 0.20),
#       legend.title=element_blank(),#Remove legend title (Network)
#       axis.text = element_text(size = size),
#       axis.title = element_text(size = size)
#     )
#   p = p + scale_fill_discrete(limits = unique(h.All[[DrugSimNet]]$Network))
#   l.plot[[pi]] = p
#   saveRDS(p,paste0(FigureFile,"_",DrugSimNet,"_gamma_",gamma,"_delta_",delta,".rdata"))#readRDS
# }
# library('cowplot')
# l.plot[[1]]
# # plot_grid(l.plot[[1]], l.plot[[2]], labels=c("A", "B"), ncol = 2, nrow = 1)
# 
# ggsave(paste0(FigureFile,"_",paste0(DrugSimNets,collapse = "-"),"_gamma_",gamma,"_delta_",delta,".pdf"), width = 7, height = 7)
# 
# 
# 
# 
# #### read stored plot and plot then save to pdf file (For MH and H Methods)
# library(ggplot2)
# library('cowplot')
# p1 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_H1-H2-H3-MH12-MH13-MH23-MH123_PREDICT_gamma_0.5_delta_0.5.rdata")
# 
# p2 = readRDS("~/Manuscripts/99MHDR/Figures/Figure_H1-H2-H3-MH12-MH13-MH23-MH123_CHEM_gamma_0.5_delta_0.5.rdata")
# 
# plot_grid(p1, p2, labels=c("(a)", "(b)"), ncol = 2, nrow = 1)
# ggsave("~/Manuscripts/99MHDR/Figures/Figure_H1-H2-H3-MH12-MH13-MH23-MH123_PREDICT_CHEM_gamma_0.5_delta_0.5.pdf", width = 14, height = 7)
# 
# # ####################################
# # #Compare the Chart between Origin (O) and Downsampled (D) data
# # df.O = h.All[["miR2Disease"]]#miR2Disease/HMDD
# # p.O = ggplot(df.O, aes(x=FPR, y=TPR, shape=Network, color = Network)) +
# #   # geom_point() +
# #   geom_line(size=1) +
# #   theme(
# #     legend.position = c(0.5, 0.15)
# #   )
# # 
# # # #Sol 1: Using downsample in groupdata2, or downSample in caret. Both down sampling to the number of sample in minority class (Network)
# # # library(groupdata2)
# # # df.D = downsample(df.O, cat_col = "Network")
# # # # library(caret)
# # # # df.D = downsample(df.O[,1:2], as.factor(df.O$Network))
# # # 
# # # p.D = ggplot(df.D, aes(x=FPR, y=TPR, shape=Network, color = Network)) +
# # #   geom_point() +
# # #   theme(
# # #     legend.position = c(0.5, 0.15)
# # #   )
# # 
# # #Sol 2: Use sample() function. More plexible
# # df.D = df.O[sample(1:nrow(df.O),nrow(df.O)/10),]
# # 
# # p.D = ggplot(df.D, aes(x=FPR, y=TPR, shape=Network, color = Network)) +
# #   # geom_point() +
# #   geom_line(size=1) +
# #   theme(
# #     legend.position = c(0.5, 0.15)
# #   )
# # 
# # plot_grid(p.O, p.D, labels=c("A", "B"), ncol = 2, nrow = 1)
# # 
