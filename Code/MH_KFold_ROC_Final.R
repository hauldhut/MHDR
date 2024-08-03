Method = "MH123"#MH123/MH12/MH13/MH23/H1/H2/H3
gamma = 0.5

start_time <- Sys.time()

library(RandomWalkRestartMH)
library(igraph)

#need to install foreach and doParallel packages for this code to run
library(foreach)
library(doParallel)

setwd("~/Manuscripts/99HDR2/Code")

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

#add func for RWR on multiplex-heter nw
do_something <- function(disease_MultiplexObject,drug_MultiplexObject,
                         DiDrRelation,SeedDisease, seeddrug, prd_diseases) {
  
  #Create multiplex-heterosgenous nw
  
  DiDrRelation_disease <- DiDrRelation[which(DiDrRelation$diseaseid %in% disease_MultiplexObject$Pool_of_Nodes),]
  
  #Create multiplex-heterosgenous nw
  disease_drug_Net <- create.multiplexHet(disease_MultiplexObject, drug_MultiplexObject, 
                                              DiDrRelation_disease)
  
  disease_drug_Net_TranMatrix <- compute.transition.matrix(disease_drug_Net)
  
  #compute 
  # tau <- c(1,1)
  Ranking_Results <- Random.Walk.Restart.MultiplexHet(disease_drug_Net_TranMatrix,
                                                                    disease_drug_Net,SeedDisease,
                                                                    seeddrug, r = gamma)
  
  #create labels for ranking results
  tf = Ranking_Results$RWRMH_Multiplex1
  
  tf$labels <- ifelse(tf$NodeNames %in% prd_diseases, 1, 0)
  
  # calculating AUC
  resultspred = prediction(tf$Score, tf$labels)
  
  pauc.perf = performance(resultspred, measure = "auc")
  return(list(pauc.perf@y.values[[1]],data.frame(Scores=tf$Score, Labels=tf$labels)))
}

#count disease for each drug
sub_sum <- aggregate(diseaseid~drugid, data=DiDr.frame, FUN=function(x) c(count=length(x)))

#extract drug with only k or more diseases
k=10
sub_sum <- sub_sum[which(sub_sum$diseaseid>=k),]
sub_sum$drugid_no <- c(1:length(sub_sum$drugid))

#extract DiDr.frame with only disease from sub_sum
DiDr.frame1 <- DiDr.frame[which(DiDr.frame$drugid %in% sub_sum$drugid),]
rownames(DiDr.frame1) <- NULL #reset frame index

#func to assign k groups for each set of drug-disease, 
#as well as increment group no. for each group (k=3)
assign_group_no <- function(sub_sum,DiDr.frame1,k) {
  
  #set an empty data frame for a new DiDr.frame
  mylist.names <- c("drugid","diseaseid", "drugid_no","group_no")
  DiDr.frame2 <- sapply(mylist.names,function(x) NULL)
  
  
  for (j in 1:length(sub_sum$drugid)) {
    
    count = sub_sum$diseaseid[[j]]
    
    if(count<k) next
    
    set_no = floor(count/k)
    
    print(paste(k,count, set_no))
    
    group_vec = vector()
    for(gi  in 1:(k-1)){
      group_vec = c(group_vec,rep(gi,set_no))
    }
    group_vec = c(group_vec,rep(k,count-set_no*(k-1)))

    # group_vec <- c(rep(1,set_no), rep(2,set_no), rep(3,set_no), rep(4,set_no),rep(5,set_no), rep(6,set_no),rep(7,set_no), rep(8,set_no), rep(9,set_no),rep(10,(count-set_no*9)))
    
    subset <- DiDr.frame[which(DiDr.frame$drugid==sub_sum$drugid[[j]]),]
    subset$drugid_no <- rep(j,count)
    subset$group <- group_vec
    
    DiDr.frame2 <- rbind(DiDr.frame2,subset)
  }
  return(DiDr.frame2)
}

#assign each group of n/k disease-drug with an group id
out <- assign_group_no(sub_sum,DiDr.frame1,k)
out$group <- (out$drugid_no-1)*k+out$group
DiDr.frame2 <- out

#set an empty data frame for a new DiDr.frame
auc_results <- sapply(c("disease","group","auc"),function(x) NULL)

for (i in 1:max(DiDr.frame2$group)) {
  seeddrug = unique(DiDr.frame2$drugid[which(DiDr.frame2$group==i)])
  auc_results$drugid[i] <- seeddrug
  auc_results$group[i] <- i
}

#set up paralell processing (adjust the no_cores as per running system)
no_cores <- 6
cl <- makeCluster(no_cores)
registerDoParallel(cl)

#loop through
res <- foreach(i = 1:max(DiDr.frame2$group), .combine = rbind) %dopar% {
  
  library(RandomWalkRestartMH)
  library(igraph)
  library(ROCR)
  
  prd_diseases = DiDr.frame2$diseaseid[which(DiDr.frame2$group==i)]
  seeddrug = unique(DiDr.frame2$drugid[which(DiDr.frame2$group==i)])
  
  disease_relation = DiDr.frame2[which(DiDr.frame2$drugid==seeddrug),]
  SeedDisease = disease_relation$diseaseid[-c(which(disease_relation$diseaseid %in% prd_diseases))]
  
  # get bipartite graph without prd_diseases - disease linkages
  DiDrRelation <- DiDr.frame2[-with(DiDr.frame2, which(diseaseid %in% prd_diseases & drugid %in% seeddrug)),][1:2]
  
  res <- do_something(disease_MultiplexObject,drug_MultiplexObject,
                    DiDrRelation,SeedDisease, seeddrug, prd_diseases)
  res <- append(res, seeddrug)
}

dim(res)
stopCluster(cl)


df.res = data.frame(trial = c(1:nrow(res)),auc = unlist(res[,1]))
df.res.final = merge(df.res, DiDr.frame2[c("disease","group")], by.x="trial", by.y="group", all.x = TRUE)
df.res.final = unique(df.res.final)

Result_byTrialFile = paste0("../Results/",Method,"_byTrial_",DrugSimNet,"_ROC_KFold10.csv")
cat(Result_byTrialFile,"\n")
write.csv(df.res.final,Result_byTrialFile, row.names = FALSE, quote = FALSE)

aucavgbyTrial = round(mean(df.res.final$auc),3)
aucavgbyTrial.sd = round(sd(df.res.final$auc),3)


res.final = NULL
for(i in 1:nrow(res)){
  res.final = rbind(res.final, res[i,2][[1]])
}

Result_by_Score_n_LabelFile = paste0("../Results/",Method,"_byTrial_",DrugSimNet,"_Score_n_Label_KFold10.csv")
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
