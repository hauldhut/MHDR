# Methods = c("H1","H2","H3","MH12","MH13","MH23","MH123")
Method = "MH123"
delta = 0.5
gamma = 0.5#gamma
# for(delta in c(0.1,0.3, 0.7, 0.9)){#0.1,
for(gamma in c(0.1,0.3,0.7, 0.9)){#0.1,
  cat("PROCESSING FOR delta=", delta,"gamma=",gamma,"\n")
  start_time <- Sys.time()
  
  library(RandomWalkRestartMH)
  library(igraph)
  library(foreach)
  library(doParallel)
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
  
  if(Method == "MH123"){
    disease_MultiplexObject <- create.multiplex(list(DiSimNet1.g,DiSimNet2.g,DiSimNet3.g),Layers_Name = c("DiSimNet1","DiSimNet2","DiSimNet3"))  
    tau1 = 1
    tau2 = 1
    tau3 = 1
    tau <- c(tau1, tau2, tau3)
  }else if(Method == "MH12"){
    disease_MultiplexObject <- create.multiplex(list(DiSimNet1.g,DiSimNet2.g),Layers_Name = c("DiSimNet1","DiSimNet2"))
    tau1 = 1
    tau2 = 1
    tau <- c(tau1, tau2)
  }else if(Method == "MH13"){
    disease_MultiplexObject <- create.multiplex(list(DiSimNet1.g,DiSimNet3.g),Layers_Name = c("DiSimNet1","DiSimNet3"))
    tau1 = 1
    tau2 = 1
    tau <- c(tau1, tau2)
  }else if(Method == "MH23"){
    disease_MultiplexObject <- create.multiplex(list(DiSimNet2.g,DiSimNet3.g),Layers_Name = c("DiSimNet2","DiSimNet3"))
    tau1 = 1
    tau2 = 1
    tau <- c(tau1, tau2)
  }else if(Method == "H1"){
    disease_MultiplexObject <- create.multiplex(list(DiSimNet1.g),Layers_Name = c("DiSimNet1"))
    tau <- c(1)
  }else if(Method == "H2"){
    disease_MultiplexObject <- create.multiplex(list(DiSimNet2.g),Layers_Name = c("DiSimNet2"))
    tau <- c(1)
  }else{ 
    disease_MultiplexObject <- create.multiplex(list(DiSimNet3.g),Layers_Name = c("DiSimNet3"))
    tau <- c(1)
  }
  
  #DrugSimNet_PREDICT.txt
  #DrugSimNet_CHEM.txt
  DrugSimNet = "PREDICT"
  cat("PROCESSING FOR Method=", Method, "and DrugSimNet=", DrugSimNet,"\n")
  
  if(DrugSimNet == "PREDICT"){
    DrSimNet <- read.delim("../Data/DrugSimNet_PREDICT.txt",header = FALSE)
  }else{
    DrSimNet <- read.delim("../Data/DrugSimNet_CHEM.txt",header = FALSE)  
  }
  
  DrSimNet.frame <- data.frame(DrSimNet[[1]], DrSimNet[[3]])
  DrSimNet.weight = DrSimNet[[2]]
  
  DrSimNet.g <- graph.data.frame(d = DrSimNet.frame, directed = FALSE)
  E(DrSimNet.g)$weight <- DrSimNet.weight
  
  drug_MultiplexObject <- create.multiplex(list(DrSimNet.g),
                                              Layers_Name = c("DrSimNet"))
  
  #Add DiDrRelation
  DiDr.frame <- read.csv("../Data/Drug2Disease_PREDICT_BinaryInteraction.csv", header = TRUE)
  DiDr.frame <- DiDr.frame[which(DiDr.frame$diseaseid %in% disease_MultiplexObject$Pool_of_Nodes),]
  DiDr.frame <- DiDr.frame[which(DiDr.frame$drugid %in% drug_MultiplexObject$Pool_of_Nodes),]
  
  dim(DiDr.frame)
  
  #func
  do_something <- function(disease_MultiplexObject,drug_MultiplexObject,
                           DiDrRelation,SeedDisease, seeddrug, prd_diseases) {
    
    #Create multiplex-heterosgenous nw
    
    DiDrRelation_disease <- DiDrRelation[which(DiDrRelation$diseaseid %in% disease_MultiplexObject$Pool_of_Nodes),]
    
    #Create multiplex-heterosgenous nw
    disease_drug_Net <- create.multiplexHet(disease_MultiplexObject, drug_MultiplexObject, 
                                                DiDrRelation_disease)
    
    disease_drug_Net_TranMatrix <- compute.transition.matrix(disease_drug_Net, delta1=delta, delta2=delta)
    
    #compute 
    #tau <- c(1,1)
    Ranking_Results <- Random.Walk.Restart.MultiplexHet(disease_drug_Net_TranMatrix,
                                                                      disease_drug_Net,SeedDisease,
                                                                      seeddrug, r = gamma)
                                                                      
    
    #create labels for ranking results
    tf = Ranking_Results$RWRMH_Multiplex1
    tf$labels <- ifelse(tf$NodeNames==prd_disease, 1, 0)
    
    # calculating AUC
    resultspred = prediction(tf$Score, tf$labels)
    
    pauc.perf = performance(resultspred, measure = "auc")
    return(list(pauc.perf@y.values[[1]],data.frame(Scores=tf$Score, Labels=tf$labels)))
    # return(pauc.perf@y.values[[1]])
  }
  
  
  no_cores <- 15
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  #loop through
  res <- foreach(i=1:length(DiDr.frame$diseaseid), .combine = rbind) %dopar% {
    
    library(RandomWalkRestartMH)
    library(igraph)
    library(ROCR)
    
    prd_disease = DiDr.frame$diseaseid[[i]]
    seeddrug = DiDr.frame$drugid[[i]]
    
    disease_relation = DiDr.frame[which(DiDr.frame$drugid==seeddrug),]
    
    SeedDisease = disease_relation$diseaseid[-c(which(disease_relation$diseaseid==prd_disease))]
    
    DiDrRelation <- DiDrRelation <- DiDr.frame[-with(DiDr.frame, which(diseaseid %in% prd_disease & drugid %in% seeddrug)),]
    
    res <- do_something(disease_MultiplexObject,drug_MultiplexObject,
                      DiDrRelation,SeedDisease, seeddrug, prd_diseases)
  }
  
  dim(res)
  stopCluster(cl)
  
  #Store auc by Trial --> From this, auc by Disease can be obtain by aggregate function when summarizing AUCs
  DiDr.frame$auc <- unlist(res[,1])
  
  Result_byTrialFile = paste0("../Results/",Method,"_byTrial_",DrugSimNet,"_gamma_",gamma,"_delta_",delta,"_ROC.csv")
  cat(Result_byTrialFile,"\n")
  write.csv(DiDr.frame,Result_byTrialFile, row.names = FALSE, quote = FALSE)
  
  aucavgbyTrial = round(mean(DiDr.frame$auc),3)
  aucavgbyTrial.sd = round(sd(DiDr.frame$auc),3)
  
  
  res.final = NULL
  for(i in 1:nrow(res)){
    res.final = rbind(res.final, res[i,2][[1]])
  }
  
  Result_by_Score_n_LabelFile = paste0("../Results/",Method,"_byTrial_",DrugSimNet,"_gamma_",gamma,"_delta_",delta,"_Score_n_Label.csv")
  cat(Result_by_Score_n_LabelFile,"\n")
  write.csv(res.final,Result_by_Score_n_LabelFile, row.names = FALSE, quote = FALSE)
  # write.csv(data.frame(Scores = Scores, Labels = Labels),Result_by_Score_n_LabelFile, row.names = FALSE, quote = FALSE)
  
  resultspred = prediction(res.final$Scores, res.final$Labels)
  auc.perf = performance(resultspred, measure = "auc")
  aucavgbyAll = round(auc.perf@y.values[[1]],3)
  
  cat("Method=",Method,'DrugSimNet=',DrugSimNet,'aucavgbyAll=',aucavgbyAll,'aucavgbyTrial=',aucavgbyTrial,"(+-",aucavgbyTrial.sd,")\n")
  
  end_time <- Sys.time()
  timediff = end_time - start_time
  print(timediff)
}