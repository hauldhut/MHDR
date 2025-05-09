# Code
- **M_LOOCV_Params_Final.R**: Run Leave-One-Out cross-validation for Monoplex/Multiplex Disease Networks for restart probability (γ) in range [0.1, 0.9], and between-disease-disease-network jumping probability (𝛿) in range [0.1, 0.9].
  - *M1, M2, M3*: For Monoplex Disease Networks (i.e., Disease Similarity Networks)
  - *M12, M23, M13, M123*: For Multiplex Disease Networks
 
- **MH_LOOCV_Params_Final.R**: Run Leave-One-Out cross-validation for Heterogeneous/Multiplex-Heterogeneous Networks for restart probability (γ) in range [0.1, 0.9], and (b) between-disease-disease-network jumping probability (𝛿) in range [0.1, 0.9].
  - *H1, H2, H3*: For Heterogeneous Networks which connect a Drug Similarity Network with a Monoplex Disease Network
  - *MH12, MH23, MH13, MH123*: For Multiplex-Heterogeneous Networks which connect a Drug Similarity Network with a Multiplex Disease Network

- **Investigate different kLN, and similarity thresholds on MimMiner-based Disease Similarity Networks**:
  - **M_LOOCV_MimMiner_Final.R**: Run Leave-One-Out cross-validation for Monoplex/Multiplex Disease Networks
  - **MH_LOOCV_MimMiner_Final.R**: Run Leave-One-Out cross-validation for Heterogeneous/Multiplex-Heterogeneous Networks

- **Summary**:
  - **Summarize_AUROC_Params_Final.R**: To investigate the prediction performance in terms of AUROC by parameters
  - **Summarize_AUROC_All_Final.R**: To summarize and draw ROC curves
  - **Summarize_AUPRC_All_Final.R**: To summarize and draw PRC curves
  - **Create_Figures.R**:

- **Prediction and Evidence Collection**:
  - **MH_Predict_Evidence_Final.R**: To predict and select top 10 highly ranked diseases for each drug, then find evidence supporting the promissing drug-disease associations
  - **Visualize_DrugDisease_from_Evidence_Final.R**: Visualize promissing drug-disease associations

- **Folder Comparison**: For the comparison with TP-NRWRH, DDAGDL, and RGLDR 
  
