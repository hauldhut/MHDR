library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)

setwd("~/Manuscripts/99MHDR/Code")

# Function to read and visualize a GraphML file
visualize_graph <- function(graph, graphname, maxEdgeWeight) {
  # # Read the GraphML file
  # graph <- read_graph(file_path, format = "graphml")
  # 
  # # Check available node attributes
  # print(vertex_attr(graph))
  
  # Extract node types (ensure correct attribute name)
  node_types <- as.factor(V(graph)$type)  # Convert to factor to match scale mappings
  
  # Extract edge weights
  edge_weights <- E(graph)$weight
  
  # Define colors and shapes for each node type
  type_colors <- c("Drug" = "red", "Disease" = "blue", "Pathway" = "pink", "ProteinComplex" = "purple")
  # type_shapes <- c("Drug" = 21, "Disease" = 22, "Pathway" = 23, "ProteinComplex" = 24)  # Different shape encoding
  type_shapes <- c("Drug" = "circle", "Disease" = "square", "Pathway" = "triangle", "ProteinComplex" = "diamond")
  
  # Ensure only existing levels are used
  type_colors <- type_colors[levels(node_types)]
  type_shapes <- type_shapes[levels(node_types)]
  
  # Plot using ggraph
  p = ggraph(graph, layout = "fr") +
    geom_edge_link(aes(edge_width = edge_weights), color = "green", alpha = 0.8) +
    geom_node_point(aes(color = node_types, shape = node_types), size = 5) +
    geom_node_text(aes(label = name), repel = TRUE) +
    scale_color_manual(values = type_colors) +
    scale_shape_manual(values = type_shapes) +
    scale_edge_width(range = c(0.2, 0.2*maxEdgeWeight)) +  # Scale edge width based on weight
    labs(
      color = "Node Type",
      shape = "Node Type",
      edge_width = "Edge Weight"
    ) +
    theme_void()
  
  return(p)
}

#COLLECT EVIDENCE FOR EACH DRUG
Gene_Info <- read.delim("../Data/EntrezGeneInfo_New.txt",header = TRUE)
Disease_Info <- read.delim("../Data/Phenotype2Genes_Full.txt",header = TRUE)
Drug_Info <- read.delim("../Data/Drug2Targets_Full_Detail.txt", header = TRUE)
Pathway_Info <- read.delim("../Data/Pathway2Genes_New.txt", header = TRUE)
Complex_Info <- read.delim("../Data/ProteinComplex2Genes_Human_Latest.txt", header = TRUE)

# library(stringr)
# library(stringi)
# Disease_Info$GeneList = stringr::str_trim(Disease_Info$GeneList)
# Disease_Info$GeneList = stringi::stri_trim(Disease_Info$GeneList)
Disease_Info$GeneList = gsub("[\r\n ]", "", Disease_Info$GeneList)
Drug_Info$GeneList = gsub("[\r\n ]", "", Drug_Info$GeneList)
Pathway_Info$GeneList = gsub("[\r\n ]", "", Pathway_Info$GeneList)
Complex_Info$GeneList = gsub("[\r\n ]", "", Complex_Info$GeneList)

library(hash)
hDi = hash(keys = Disease_Info$MIMID, values = Disease_Info$Name)
hDr = hash(keys = Drug_Info$KEGGID, values = Drug_Info$Name)
hP = hash(keys = Pathway_Info$KEGID, values = Pathway_Info$Name)
hC = hash(keys = Complex_Info$ComplexID, values = Complex_Info$Name)


# share = "sharePathway"
# DrOI = "D06320"
# DiOI = "114480"

# share = "shareComplex"
# DrOI = "D06320"
# DiOI = "114480"

# 
share = "sharePathwayComplex"
DrOI = "D00067"
DiOI = "114480"

# DrOI = "D00400"
# DiOI = "145500"


Evid_Info <- read.delim(paste0("../Results/Evidence/MH123_df.evid_",share,".txt"),header = TRUE)


graphname = paste0(DrOI,"-",DiOI)

Evid_Info_Pair <- Evid_Info %>%
  filter(KEGGID == DrOI, MIMID == DiOI)

if (nrow(Evid_Info_Pair)<1){
  print("NO PAIR")
}

# Split the SharedPathway and SharedComplex into lists of pathways and complexes
if (share=="sharePathway"){
  SharedPathway <- unique(strsplit(Evid_Info_Pair[1,]$SharedPathway, ", ")[[1]])
}else if (share=="shareComplex"){
  SharedComplex <- unique(strsplit(Evid_Info_Pair[1,]$SharedComplex, ", ")[[1]])
}else if (share=="sharePathwayComplex"){
  SharedPathway <- unique(strsplit(Evid_Info_Pair[1,]$SharedPathway, ", ")[[1]])
  SharedComplex <- unique(strsplit(Evid_Info_Pair[1,]$SharedComplex, ", ")[[1]])
}

# Create an empty graph
g <- graph.empty(directed = FALSE)

# drug <- hDr[[Evid_Info_Pair[1,]$KEGGID]]
# disease <- hDi[[paste0("MIM",Evid_Info_Pair[1,]$MIMID)]]
hNodeID2Name = hash()

dr_id = DrOI
di_id = paste0("MIM",DiOI)
hNodeID2Name[[dr_id]] = hDr[[dr_id]]
hNodeID2Name[[di_id]] = hDi[[di_id]]


DrAssocProteins = unique(strsplit(Drug_Info[Drug_Info$KEGGID==dr_id,]$GeneList, ",")[[1]])
DiAssocGenes = unique(strsplit(Disease_Info[Disease_Info$MIMID==di_id,]$GeneList, ",")[[1]])

g <- add_vertices(g, 1, name = dr_id, type = "Drug")  
g <- add_vertices(g, 1, name = di_id, type = "Disease")  


if (share=="sharePathway"){
  for(pathway in SharedPathway){
    p_id = paste0("path:",pathway)
    hNodeID2Name[[p_id]] = strsplit(hP[[p_id]]," - ")[[1]][1]
    g <- add_vertices(g, 1, name = p_id, type = "Pathway")  
  }
}else if (share=="shareComplex"){
  for (complex in SharedComplex) {
    c_id = complex
    hNodeID2Name[[c_id]] = hC[[c_id]]
    g <- add_vertices(g, 1, name = c_id, type = "ProteinComplex")  
  }
}else if (share=="sharePathwayComplex"){
  for(pathway in SharedPathway){
    p_id = paste0("path:",pathway)
    hNodeID2Name[[p_id]] = strsplit(hP[[p_id]]," - ")[[1]][1]
    g <- add_vertices(g, 1, name = p_id, type = "Pathway")  
  }
  for (complex in SharedComplex) {
    c_id = complex
    hNodeID2Name[[c_id]] = hC[[c_id]]
    g <- add_vertices(g, 1, name = c_id, type = "ProteinComplex")  
  }
}

print(V(g))

# Add edges based on shared pathways and protein complexes

if (share=="sharePathway"){
  for (pathway in SharedPathway) {
    p_id = paste0("path:",pathway)
    
    pathway = paste0("path:",pathway)
    PInvolvedGenes = unique(strsplit(Pathway_Info[Pathway_Info$KEGID==pathway,]$GeneList, ",")[[1]])
    NuSharePvsDr = length(intersect(PInvolvedGenes,DrAssocProteins))
    NuSharePvsDi = length(intersect(PInvolvedGenes,DiAssocGenes))
    
    g <- g + edge(dr_id, p_id, type = "Drug-Pathway", weight = NuSharePvsDr)
    g <- g + edge(di_id, p_id, type = "Disease-Pathway", weight = NuSharePvsDi)
  }
}else if (share=="shareComplex"){
  for (complex in SharedComplex) {
    c_id = complex
    
    CInvolvedProteins = unique(strsplit(Complex_Info[Complex_Info$ComplexID==c_id,]$GeneList, ",")[[1]])
    NuShareCvsDr = length(intersect(CInvolvedProteins,DrAssocProteins))
    NuShareCvsDi = length(intersect(CInvolvedProteins,DiAssocGenes))
    
    g <- g + edge(dr_id, c_id, type = "Drug-ProteinComplex", weight = NuShareCvsDr)
    g <- g + edge(di_id, c_id, type = "Disease-ProteinComplex", weight = NuShareCvsDi)
  }
}else if (share=="sharePathwayComplex"){
  for (pathway in SharedPathway) {
    p_id = paste0("path:",pathway)
    
    pathway = paste0("path:",pathway)
    PInvolvedGenes = unique(strsplit(Pathway_Info[Pathway_Info$KEGID==pathway,]$GeneList, ",")[[1]])
    NuSharePvsDr = length(intersect(PInvolvedGenes,DrAssocProteins))
    NuSharePvsDi = length(intersect(PInvolvedGenes,DiAssocGenes))
    
    g <- g + edge(dr_id, p_id, type = "Drug-Pathway", weight = NuSharePvsDr)
    g <- g + edge(di_id, p_id, type = "Disease-Pathway", weight = NuSharePvsDi)
  }
  
  for (complex in SharedComplex) {
    c_id = complex
    
    CInvolvedProteins = unique(strsplit(Complex_Info[Complex_Info$ComplexID==c_id,]$GeneList, ",")[[1]])
    NuShareCvsDr = length(intersect(CInvolvedProteins,DrAssocProteins))
    NuShareCvsDi = length(intersect(CInvolvedProteins,DiAssocGenes))
    
    g <- g + edge(dr_id, c_id, type = "Drug-ProteinComplex", weight = NuShareCvsDr)
    g <- g + edge(di_id, c_id, type = "Disease-ProteinComplex", weight = NuShareCvsDi)
  }
}

maxEdgeWeight = max(E(g)$weight)
V(g)$name <- sapply(V(g)$name, function(x) hNodeID2Name[[x]])

# Example usage
p = visualize_graph(g, graphname, maxEdgeWeight)
saveRDS(p, paste0(graphname,"_Figure.rdata"))
print(p)
ggsave(paste0(graphname,"_Figure.pdf"), plot = p, width = 5, height = 5, dpi = 300)

# #Load to combine
# library(cowplot)
# p1 = readRDS("~/Manuscripts/99MHDR/Code/D00067-114480_Figure.rdata")
# p2 = readRDS("~/Manuscripts/99MHDR/Code/D00400-145500_Figure.rdata")
# 
# plot_grid(p1, p2, labels = c("(a)","(b)"), ncol = 2, nrow = 1)
# ggsave(paste0("Figure5.pdf"), width = 10, height = 5, dpi = 600)

