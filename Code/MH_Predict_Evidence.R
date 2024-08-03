Method = "MH123"#MH123/MH12/MH13/MH23/H1/H2/H3

start_time <- Sys.time()

library(RandomWalkRestartMH)
library(igraph)
library(foreach)
library(doParallel)
library(rclinicaltrials)

setwd("~/Manuscripts/99HDR2/Code")

#DiseaseSimNet_OMIM.txt
#DiseaseSimNet_HPO.sif
#DiseaseSimNet_GeneNet.txt

enhancer1 <- read.delim("../Data/DiseaseSimNet_OMIM.txt",header = FALSE)
enhancer1.frame <- data.frame(enhancer1[[1]], enhancer1[[3]])
enhancer1.g <- graph.data.frame(d = enhancer1.frame, directed = FALSE)
enhancer1.weight = enhancer1[[2]]
E(enhancer1.g)$weight <- enhancer1.weight


enhancer2 <- read.delim("../Data/DiseaseSimNet_HPO.sif",header = FALSE)
enhancer2.frame <- data.frame(enhancer2[[1]], enhancer2[[3]])
enhancer2.g <- graph.data.frame(d = enhancer2.frame, directed = FALSE)
enhancer2.weight = enhancer2[[2]]
E(enhancer2.g)$weight <- enhancer2.weight

enhancer3 <- read.delim("../Data/DiseaseSimNet_GeneNet.txt",header = FALSE)
enhancer3.frame <- data.frame(enhancer3[[1]], enhancer3[[3]])
enhancer3.g <- graph.data.frame(d = enhancer3.frame, directed = FALSE)
enhancer3.weight = enhancer3[[2]]
E(enhancer3.g)$weight <- enhancer3.weight

if(Method == "MH123"){
  enhancer_MultiplexObject <- create.multiplex(list(enhancer1.g,enhancer2.g,enhancer3.g),Layers_Name = c("enhancer1","enhancer2","enhancer3"))  
  tau1 = 1
  tau2 = 1
  tau3 = 1
  tau <- c(tau1, tau2, tau3)
}else if(Method == "MH12"){
  enhancer_MultiplexObject <- create.multiplex(list(enhancer1.g,enhancer2.g),Layers_Name = c("enhancer1","enhancer2"))
  tau1 = 1
  tau2 = 1
  tau <- c(tau1, tau2)
}else if(Method == "MH13"){
  enhancer_MultiplexObject <- create.multiplex(list(enhancer1.g,enhancer3.g),Layers_Name = c("enhancer1","enhancer3"))
  tau1 = 1
  tau2 = 1
  tau <- c(tau1, tau2)
}else if(Method == "MH23"){
  enhancer_MultiplexObject <- create.multiplex(list(enhancer2.g,enhancer3.g),Layers_Name = c("enhancer2","enhancer3"))
  tau1 = 1
  tau2 = 1
  tau <- c(tau1, tau2)
}else if(Method == "H1"){
  enhancer_MultiplexObject <- create.multiplex(list(enhancer1.g),Layers_Name = c("enhancer"))
  tau <- c(1)
}else if(Method == "H2"){
  enhancer_MultiplexObject <- create.multiplex(list(enhancer2.g),Layers_Name = c("enhancer"))
  tau <- c(1)
}else{ 
  enhancer_MultiplexObject <- create.multiplex(list(enhancer3.g),Layers_Name = c("enhancer"))
  tau <- c(1)
}

#DrugSimNet_PREDICT.txt
#DrugSimNet_CHEM.txt
DrugSimNet = "PREDICT"
if(DrugSimNet == "PREDICT"){
  disease <- read.delim("../Data/DrugSimNet_PREDICT.txt",header = FALSE)
}else{
  disease <- read.delim("../Data/DrugSimNet_CHEM.txt",header = FALSE)  
}

disease.frame <- data.frame(disease[[1]], disease[[3]])
disease.weight = disease[[2]]

disease.g <- graph.data.frame(d = disease.frame, directed = FALSE)
E(disease.g)$weight <- disease.weight

disease_MultiplexObject <- create.multiplex(list(disease.g),
                                             Layers_Name = c("disease"))

#Add EDRelation
ED.frame <- read.csv("../Data/Drug2Disease_Name_PREDICT_ID.txt_BinaryInteraction.csv", header = TRUE)
ED.frame <- ED.frame[which(ED.frame$enhancer %in% enhancer_MultiplexObject$Pool_of_Nodes),]
ED.frame <- ED.frame[which(ED.frame$disease %in% disease_MultiplexObject$Pool_of_Nodes),]


#Create multiplex-heterosgenous nw
enhancer_Disease_Net <- create.multiplexHet(enhancer_MultiplexObject, disease_MultiplexObject, ED.frame)

enhancerHetTranMatrix <- compute.transition.matrix(enhancer_Disease_Net, lambda = 0.5)

#RANK FOR EACH DRUG
Drug_list = unique(ED.frame$disease)

r=0.5

library(hash)


no_cores <- 4
cl <- makeCluster(no_cores)
registerDoParallel(cl)

dri=1
topRanking = 10
res <- foreach(dri = 1:length(Drug_list), .combine = cbind) %dopar% {#length(Drug_list)
  library(RandomWalkRestartMH)
# res = hash()
# dc = 0
# for(dr in Drug_list){
#   dc=dc+1
#   cat(dc,"/",length(Drug_list),": Ranking for drug:", dr, "\n")
  
  dr = Drug_list[dri]
  
  seeddisease = dr
  SeedEnhancer = ED.frame[ED.frame$disease==dr,]$enhancer

  
  RWRH_enhancer_Disease_Results <- Random.Walk.Restart.MultiplexHet(enhancerHetTranMatrix,
                                                                    enhancer_Disease_Net,SeedEnhancer,
                                                                    seeddisease, r = r)
  top10_di = head(RWRH_enhancer_Disease_Results$RWRMH_Multiplex1$NodeNames,topRanking)
  
  df.ranking <-  data.frame(top10_di)
  names(df.ranking)[1] <- dr
  
  res <- df.ranking
  
  # df.ranking <-  data.frame(DrID = top10_di)
  # names(df.ranking)[names(df.ranking) == "DrID"] <- dr
}
stopCluster(cl)

res

#COLLECT EVIDENCE FOR EACH DRUG
Gene_Info <- read.delim("../Data/EntrezGeneInfo_New.txt",header = TRUE)
Disease_Info <- read.delim("../Data/Phenotype2Genes_Full.txt",header = TRUE)
Drug_Info <- read.delim("../Data/Drug2Targets_Full_Detail.txt", header = TRUE)
Pathway_Info <- read.delim("../Data/Pathway2Genes_New.txt", header = TRUE)
Complex_Info <- read.delim("../Data/ProteinComplex2Genes.txt", header = TRUE)

# library(stringr)
# library(stringi)
# Disease_Info$GeneList = stringr::str_trim(Disease_Info$GeneList)
# Disease_Info$GeneList = stringi::stri_trim(Disease_Info$GeneList)
Disease_Info$GeneList = gsub("[\r\n ]", "", Disease_Info$GeneList)
Drug_Info$GeneList = gsub("[\r\n ]", "", Drug_Info$GeneList)
Pathway_Info$GeneList = gsub("[\r\n ]", "", Pathway_Info$GeneList)
Complex_Info$GeneList = gsub("[\r\n ]", "", Complex_Info$GeneList)

#Level #0
#ClinicalTrials.gov
h.evid_dr_CT= hash()
df.evid_dr_di_CT = NULL

dc = 0
for(dr in Drug_list){
  # if(dr!="D00067") next
  dc=dc+1
  if(dc<321) next#221/190/256/279/320
  cat(dc,"/",length(Drug_list),": Collecting evidence (ClinicalTrials.gov) for drug:", dr,"\n")
  
  top10_di = res[[dr]]
  dr_GeneList = Drug_Info[Drug_Info$KEGGID==dr,]$GeneList
  
  if(nrow(Drug_Info[Drug_Info$KEGGID==dr,])==0) next
  dr_Name_Full = Drug_Info[Drug_Info$KEGGID==dr,]$Name
  # cat("DrugName:",dr_Name_Full,"\n")
  
  dr_NameVec = strsplit(dr_Name_Full,", ")[[1]]
  
  h.evid_di_CT = hash()
  for(di in top10_di){
    if(nrow(Disease_Info[Disease_Info$MIMID==di,])==0) next
    
    di_Name_Full = Disease_Info[Disease_Info$MIMID==di,]$Name
    # cat("DiseaseName:",di_Name_Full,"\n")
    
    di_Name = strsplit(di_Name_Full,";")[[1]][1]
    
    NCT_Vec = NULL
    for(dr_Name in dr_NameVec){    
      # print(di)
      
      pos = which(strsplit(dr_Name, "")[[1]]=="(")[1]
      # cat(dr_Name,pos,"\n")
      
      dr_Name_Compact = dr_Name
      if(!is.na(pos)){
        dr_Name_Compact = substr(dr_Name,1,pos-2)  
      }
      
      term = dr_Name_Compact
      cond = di_Name
      
      query = c(paste0("term=",term), paste0("cond=",cond))
      # print(query)
      
      # clinicaltrials_count(query = query)
      
      #To download detailed study information, including results, use clinicaltrials_download(). 
      #Downloading lots of results may take a long time and use a substantial amount of hard drive space. 
      #You can limit the number of studies downloaded with the count option. By default, the count is limited to 20.
      y <- clinicaltrials_download(query = query, count = 20, include_results = TRUE)
      nci_vec = unique(y$study_information$study_info$nct_id)
      
      NCT_Vec = unique(append(NCT_Vec,nci_vec))
    }
    cat("\t=>","NCI IDs",dr,di,paste0(NCT_Vec,collape=", "),"\n")
    if(length(NCT_Vec)>0){
      h.evid_di_CT[[di]] = NCT_Vec
      df.evid_dr_di_CT = rbind(df.evid_dr_di_CT, data.frame(KEGGID = dr, DrugName = dr_Name_Full, MIMID = di, DiseaseName = di_Name_Full, NCTID= paste0(NCT_Vec,collapse = ", ")))  
    }
  }
  
  if(length(h.evid_di_CT)>0){
    h.evid_dr_CT[[dr]] = h.evid_di_CT  
  }
  
}
dim(df.evid_dr_di_CT)
df.evid_dr_di_CT$MIMID = gsub("MIM","",df.evid_dr_di_CT$MIMID)

CT_evidFile = paste0("../Results/Evidence/",Method,"_df.evid_dr_di_CT.txt")
write.table(df.evid_dr_di_CT,CT_evidFile, sep = "\t")

#Level #1: Shared Genes
h.evid_dr_shareGene= hash()
df.evid_shareGene = NULL


dc = 0
for(dr in Drug_list){
  dc=dc+1
  
  cat(dc,"/",length(Drug_list),": Collecting evidence (shared Gene) for drug:", dr,"\n")
  
  top10_di = res[[dr]]
  dr_GeneList = Drug_Info[Drug_Info$KEGGID==dr,]$GeneList
  
  if(nrow(Drug_Info[Drug_Info$KEGGID==dr,])==0 || nchar(dr_GeneList)<1) next
  
  dr_GeneVec = strsplit(dr_GeneList,",")[[1]]
  h.evid_di = hash()
  for(di in top10_di){
    di_GeneList = Disease_Info[Disease_Info$MIMID==di,]$GeneList
    # print(di)
    if(nrow(Disease_Info[Disease_Info$MIMID==di,])==0 || nchar(di_GeneList)<1) next
    
    di_GeneVec = strsplit(di_GeneList,",")[[1]]
    
    if(length(dr_GeneVec)>0 && length(di_GeneVec)>0){
      # cat("\t",dr,"\t",di,"\n")
      # print(dr_GeneVec)
      # print(di_GeneVec)
      dr_di_GeneVec = intersect(dr_GeneVec, di_GeneVec)
      
      if(length(dr_di_GeneVec)>0){
        h.evid_di[[di]] = dr_di_GeneVec
        
        cat("\t",dr,"\t",di,"\t",dr_di_GeneVec,"\n")
        
        dr_di_GeneVec_Symbol = NULL
        for(geneid in dr_di_GeneVec){
          dr_di_GeneVec_Symbol = append(dr_di_GeneVec_Symbol, Gene_Info[Gene_Info$EntrezID==geneid,]$Symbol)   
        }
        
        SharedGene = paste(dr_di_GeneVec_Symbol,collapse = ", ")
        df.evid_shareGene <- rbind(df.evid_shareGene, data.frame(KEGGID = dr, MIMID=di, SharedGene=SharedGene))
      }  
    }
  }
  
  if(length(h.evid_di)>0){
    h.evid_dr_shareGene[[dr]] = h.evid_di  
  }
  
}
df.evid_shareGene
df.evid_shareGene.Final = merge(df.evid_shareGene, Disease_Info[c("MIMID","Name")], by="MIMID")
df.evid_shareGene.Final = merge(df.evid_shareGene.Final, Drug_Info[c("KEGGID","Name")], by="KEGGID")
colnames(df.evid_shareGene.Final) = c("KEGGID","MIMID","SharedGene","DiseaseName","DrugName")

df.evid_shareGene.Final$MIMID = gsub("MIM","",df.evid_shareGene.Final$MIMID)

head(df.evid_shareGene.Final)

shareGene_evidFile = paste0("../Results/Evidence/",Method,"_df.evid_shareGene.txt")
write.table(df.evid_shareGene.Final[c("KEGGID","DrugName","MIMID","DiseaseName","SharedGene")],shareGene_evidFile, sep = "\t")

#Level #2: Shared Pathways
h.evid_dr_sharePathway= hash()
df.evid_sharePathway = NULL

dc=0
for(dr in Drug_list){
  dc=dc+1
  cat(dc,"/",length(Drug_list),": Collecting evidence (shared Pathways) for drug:", dr, "\n")
  
  top10_di = res[[dr]]
  dr_GeneList = Drug_Info[Drug_Info$KEGGID==dr,]$GeneList
  
  if(nrow(Drug_Info[Drug_Info$KEGGID==dr,])==0 || nchar(dr_GeneList)<1) next
  
  dr_GeneVec = strsplit(dr_GeneList,",")[[1]]
  h.evid_di = hash()
  for(di in top10_di){
    if(nrow(df.evid_shareGene[df.evid_shareGene$Drug==dr,][df.evid_shareGene$Disease==di,])>0){
      cat("EXIST",dr,"-",di,"in df.evid_shareGene\n")
      break
    }
    
    di_GeneList = Disease_Info[Disease_Info$MIMID==di,]$GeneList
    
    if(nrow(Disease_Info[Disease_Info$MIMID==di,])==0 || nchar(di_GeneList)<1) next
    
    di_GeneVec = strsplit(di_GeneList,",")[[1]]
    
    if(length(dr_GeneVec)>0 && length(di_GeneVec)>0){
      SharedPathway = vector()
      for(pi in 1:nrow(Pathway_Info)){
        p = Pathway_Info[pi,]
        p_GeneVec = strsplit(p$GeneList,",")[[1]]
        
        if(length(intersect(dr_GeneVec, p_GeneVec))>0 && length(intersect(di_GeneVec, p_GeneVec))>0){
          pid = gsub("path:","",p$KEGID)
          SharedPathway = rbind(SharedPathway,pid)
        }
      }
      if(length(SharedPathway)>0){
        h.evid_di[[di]] = SharedPathway
        
        SharedPathway = paste(SharedPathway,collapse = ", ")
        cat("\t",dr,"\t",di,"\t",SharedPathway,"\n")
        
        df.evid_sharePathway <- rbind(df.evid_sharePathway, data.frame(KEGGID = dr, MIMID=di, SharedPathway=SharedPathway))
      }
      
    }
  }
  h.evid_dr_sharePathway[[dr]] = h.evid_di
}
df.evid_sharePathway

df.evid_sharePathway.Final = merge(df.evid_sharePathway, Disease_Info[c("MIMID","Name")], by="MIMID")
df.evid_sharePathway.Final = merge(df.evid_sharePathway.Final, Drug_Info[c("KEGGID","Name")], by="KEGGID")
colnames(df.evid_sharePathway.Final) = c("KEGGID","MIMID","SharedPathway","DiseaseName","DrugName")

df.evid_sharePathway.Final$MIMID = gsub("MIM","",df.evid_sharePathway.Final$MIMID)

head(df.evid_sharePathway.Final)

sharePathway_evidFile = paste0("../Results/Evidence/",Method,"_df.evid_sharePathway.txt")
write.table(df.evid_sharePathway.Final[c("KEGGID","DrugName","MIMID","DiseaseName","SharedPathway")],sharePathway_evidFile, sep = "\t")


#Level #3: Shared Complexes
h.evid_dr_shareComplex= hash()
df.evid_shareComplex = NULL

dc=0
for(dr in Drug_list){
  dc=dc+1
  cat(dc,"/",length(Drug_list),": Collecting evidence (shared Complex) for drug:", dr, "\n")
  
  top10_di = res[[dr]]
  dr_GeneList = Drug_Info[Drug_Info$KEGGID==dr,]$GeneList
  
  if(nrow(Drug_Info[Drug_Info$KEGGID==dr,])==0 || nchar(dr_GeneList)<1) next
  
  dr_GeneVec = strsplit(dr_GeneList,",")[[1]]
  h.evid_di = hash()
  for(di in top10_di){
    if(nrow(df.evid_shareGene[df.evid_shareGene$Drug==dr,][df.evid_shareGene$Disease==di,])>0){
      cat("EXIST",dr,"-",di,"in df.evid_shareGene\n")
      break
    }
    
    di_GeneList = Disease_Info[Disease_Info$MIMID==di,]$GeneList
    
    if(nrow(Disease_Info[Disease_Info$MIMID==di,])==0 || nchar(di_GeneList)<1) next
    
    di_GeneVec = strsplit(di_GeneList,",")[[1]]
    
    if(length(dr_GeneVec)>0 && length(di_GeneVec)>0){
      SharedComplex = vector()
      for(ci in 1:nrow(Complex_Info)){
        c = Complex_Info[ci,]
        c_GeneVec = strsplit(c$GeneList,",")[[1]]
        
        if(length(intersect(dr_GeneVec, c_GeneVec))>0 && length(intersect(di_GeneVec, c_GeneVec))>0){
          SharedComplex = rbind(SharedComplex,c$ComplexID)
          
        }
      }
      if(length(SharedComplex)>0){
        h.evid_di[[di]] = SharedComplex
        
        SharedComplex = paste(SharedComplex,collapse = ", ")
        cat("\t",dr,"\t",di,"\t",SharedComplex,"\n")
        
        df.evid_shareComplex <- rbind(df.evid_shareComplex, data.frame(KEGGID = dr, MIMID=di, SharedComplex=SharedComplex))
      }
      
    }
  }
  h.evid_dr_shareComplex[[dr]] = h.evid_di
}
df.evid_shareComplex

df.evid_shareComplex.Final = merge(df.evid_shareComplex, Disease_Info[c("MIMID","Name")], by="MIMID")
df.evid_shareComplex.Final = merge(df.evid_shareComplex.Final, Drug_Info[c("KEGGID","Name")], by="KEGGID")
colnames(df.evid_shareComplex.Final) = c("KEGGID","MIMID","SharedComplex","DiseaseName","DrugName")

df.evid_shareComplex.Final$MIMID = gsub("MIM","",df.evid_shareComplex.Final$MIMID)

head(df.evid_shareComplex.Final)

shareComplex_evidFile = paste0("../Results/Evidence/",Method,"_df.evid_shareComplex.txt")
write.table(df.evid_shareComplex.Final[c("KEGGID","DrugName","MIMID","DiseaseName","SharedComplex")],shareComplex_evidFile, sep = "\t")

### POST-PROCESS ###
pw <- read.delim(sharePathway_evidFile,header = TRUE)
pc <- read.delim(shareComplex_evidFile,header = TRUE)

dim(pw)
dim(pc)

df.evid_sharePathwayComplex <- merge(pw[c("KEGGID","DrugName","MIMID","DiseaseName","SharedPathway")], pc[c("KEGGID", "MIMID","SharedComplex")], by=c("KEGGID", "MIMID"))

dim(df.evid_sharePathwayComplex)

sharePathwayComplex_evidFile = paste0("../Results/Evidence/",Method,"_df.evid_sharePathwayComplex.txt")
write.table(df.evid_sharePathwayComplex[c("KEGGID","DrugName","MIMID","DiseaseName","SharedPathway","SharedComplex")],sharePathwayComplex_evidFile, sep = "\t")
