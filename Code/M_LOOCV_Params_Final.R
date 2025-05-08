# Methods = c("M1", "M2", "M3", "M12", "M13", "M23", "M123")
# Methods = c("M3", "M13", "M23")
Methods = c("M123")

Method = "M123"
delta = 0.5
gamma = 0.5#gamma
# for(delta in c(0.1,0.3, 0.7, 0.9)){#0.1,
for(gamma in c(0.1,0.3, 0.7, 0.9)){#0.1,
# for(Method in Methods){   
  cat("PROCESSING FOR Method=", Method, " and delta=", delta,"gamma=",gamma,"\n")
  # Method = "M2"#M123/M12/M13/M23/M1/M2/M3
  
  start_time <- Sys.time()
  
  library(RandomWalkRestartMH)
  library(igraph)
  library(ROCR)
  
  setwd("~/Manuscripts/99MHDR/Code")
  
  #DiseaseSimNet_OMIM.txt
  #DiseaseSimNet_HPO.sif
  #DiseaseSimNet_GeneNet.txt
  
  DiSimNet1 <- read.delim("../Data/DiseaseSimNet_OMIM.txt",header = FALSE)
  DiSimNet1.frame <- data.frame(DiSimNet1[[1]], DiSimNet1[[3]])
  DiSimNet1.g <- graph.data.frame(d = DiSimNet1.frame, directed = FALSE)
  DiSimNet1.weight = DiSimNet1[[2]]
  E(DiSimNet1.g)$weight <- DiSimNet1.weight
  
  
  DiSimNet2 <- read.delim("../Data/DiseaseSimNet_HPO.sif",header = FALSE)
  DiSimNet2.frame <- data.frame(DiSimNet2[[1]], DiSimNet2[[3]])
  DiSimNet2.g <- graph.data.frame(d = DiSimNet2.frame, directed = FALSE)
  DiSimNet2.weight = DiSimNet2[[2]]
  E(DiSimNet2.g)$weight <- DiSimNet2.weight
  
  DiSimNet3 <- read.delim("../Data/DiseaseSimNet_GeneNet.txt",header = FALSE)
  DiSimNet3.frame <- data.frame(DiSimNet3[[1]], DiSimNet3[[3]])
  DiSimNet3.g <- graph.data.frame(d = DiSimNet3.frame, directed = FALSE)
  DiSimNet3.weight = DiSimNet3[[2]]
  E(DiSimNet3.g)$weight <- DiSimNet3.weight
  
  if(Method == "M123"){
    disease_MultiplexObject <- create.multiplex(list(DiSimNet1.g,DiSimNet2.g,DiSimNet3.g),Layers_Name = c("DiSimNet1","DiSimNet2","DiSimNet3"))  
    tau1 = 1
    tau2 = 1
    tau3 = 1
    tau <- c(tau1, tau2, tau3)
  }else if(Method == "M12"){
    disease_MultiplexObject <- create.multiplex(list(DiSimNet1.g,DiSimNet2.g),Layers_Name = c("DiSimNet1","DiSimNet2"))
    tau1 = 1
    tau2 = 1
    tau <- c(tau1, tau2)
  }else if(Method == "M13"){
    disease_MultiplexObject <- create.multiplex(list(DiSimNet1.g,DiSimNet3.g),Layers_Name = c("DiSimNet1","DiSimNet3"))
    tau1 = 1
    tau2 = 1
    tau <- c(tau1, tau2)
  }else if(Method == "M23"){
    disease_MultiplexObject <- create.multiplex(list(DiSimNet2.g,DiSimNet3.g),Layers_Name = c("DiSimNet2","DiSimNet3"))
    tau1 = 1
    tau2 = 1
    tau <- c(tau1, tau2)
  }else if(Method == "M1"){
    disease_MultiplexObject <- create.multiplex(list(DiSimNet1.g),Layers_Name = c("DiSimNet1"))
    tau <- c(1)
  }else if(Method == "M2"){
    disease_MultiplexObject <- create.multiplex(list(DiSimNet2.g),Layers_Name = c("DiSimNet2"))
    tau <- c(1)
  }else{ 
    disease_MultiplexObject <- create.multiplex(list(DiSimNet3.g),Layers_Name = c("DiSimNet3"))
    tau <- c(1)
  }
  

  
  
  AdjMatrix_disease <- compute.adjacency.matrix(disease_MultiplexObject, delta = delta)
  AdjMatrixNorm_disease <- normalize.multiplex.adjacency(AdjMatrix_disease)
  
  #Add DiDrRelation
  BipartiteNet = "PREDICT"
  DiDr.frame <- read.csv("../Data/Drug2Disease_PREDICT_BinaryInteraction.csv", header = TRUE)
  DiDr.frame <- DiDr.frame[which(DiDr.frame$diseaseid %in% disease_MultiplexObject$Pool_of_Nodes),]
  
  #loop through
  res = NULL
  for (i in 1:length(DiDr.frame$diseaseid)) {
  
    prd_disease = DiDr.frame$diseaseid[i]
    seeddrug = DiDr.frame$drugid[i]
  
    cat("==>",i,"/",length(DiDr.frame$diseaseid),":",seeddrug,"\n")
    
    disease_relation = DiDr.frame[which(DiDr.frame$drugid==seeddrug),]
    Seeddisease = disease_relation$diseaseid[-c(which(disease_relation$diseaseid==prd_disease))]
    
  
    if (length(disease_relation$diseaseid)>=2) {
  
      Ranking_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_disease,
                                                            disease_MultiplexObject,
                                                            Seeddisease,tau = tau, r = gamma)
  
      tf = Ranking_Results$RWRM_Results
      #create labels for ranking results
      tf$labels <- ifelse(tf$NodeNames==prd_disease, 1, 0)
      
      # calculating AUC
      resultspred = prediction(tf$Score, tf$labels)
      
      pauc.perf = performance(resultspred, measure = "auc")
      
      DiDr.frame$auc[i] <- pauc.perf@y.values[[1]]
      
      res = rbind(res, data.frame(Scores=tf$Score, Labels=tf$labels))
    } else {
      DiDr.frame$auc[i] <- NA
    }
  }
  
  dim(res)
  
  #Store auc by Trial --> From this, auc by Disease can be obtain by aggregate function when summarizing AUCs
  DiDr.frame = DiDr.frame[!is.na(DiDr.frame$auc),]
  
  Result_byTrialFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_gamma_",gamma,"_delta_",delta,"_ROC.csv")
  cat(Result_byTrialFile,"\n")
  write.csv(DiDr.frame,Result_byTrialFile, row.names = FALSE, quote = FALSE)
  
  aucavgbyTrial = mean(DiDr.frame$auc)
  aucavgbyTrial.sd = sd(DiDr.frame$auc)
  
  #Store Scores and Labels for calculating and drawing AUC and ROC
  
  Result_by_Score_n_LabelFile = paste0("../Results/",Method,"_byTrial_",BipartiteNet,"_gamma_",gamma,"_delta_",delta,"_Score_n_Label.csv")
  cat(Result_by_Score_n_LabelFile,"\n")
  write.csv(res,Result_by_Score_n_LabelFile, row.names = FALSE, quote = FALSE)
  
  resultspred = prediction(res$Scores, res$Labels)
  auc.perf = performance(resultspred, measure = "auc")
  aucavgbyAll = auc.perf@y.values[[1]]
  
  cat("Method=",Method,'BipartiteNet=',BipartiteNet,'aucavgbyAll=',aucavgbyAll,'aucavgbyTrial=',aucavgbyTrial,"(+-",aucavgbyTrial.sd,")\n")
  
  
  end_time <- Sys.time()
  end_time - start_time
}
