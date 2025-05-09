For comparison with TP-NRWRH, DDAGDL, and RGLDR

- **Compare with TP-NRWRH** 
 - **MH_KFold_ROC_Final.R**: Run K-Fold cross-validation for Multiplex-Heterogeneous Networks of Drugs and Diseases
 - **Summarize_AUROC_AUPRC_KFold.R**: Summarize K-Fold AUROC and AUPRC 

- **Compare with DDAGDL and RGLDR** 
 - **gat_embeddings_Final.py**: Generate representations of drugs and diseases using Graph attention networks (GAT)
 - **transformer_embeddings_Final.py**: Generate representations of drugs and diseases using Graph transformer networks
 - **metapath2vec_embeddings_Final.py**: Generate representations of drugs and diseases using metapath2vec-based models
 - **drug_disease_classifier_XGB.py**: Evaluate prediction performance of the embedding methods (GAT, Graph transformer, and metapath-based) in terms of AUROC and AUPRC using 10-Fold cross validation

- **Setup**
 - Use Python ver 3.10.16
 - **GNN_Pytorch.yml**: Install the Python environment
