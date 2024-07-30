# This data contains all source code used in the manuscript

- M_LOOCV_ROC_Final.R: Run Leave-One-Out cross-validation for Monoplex/Multiplex Disease Networks
  - M1, M2, M3: For Monoplex Disease Networks (i.e., Disease Similarity Networks)
  - M12, M23, M13, M123: For Multiplex Disease Networks
 
- MH_LOOCV_ROC_Final.R: Run Leave-One-Out cross-validation for Heterogeneous/Multiplex-Heterogeneous Networks of Drugs and Diseases
  - H1, H2, H3: For Heterogeneous Networks which connect a Drug Similarity Network with a Monoplex Disease Network
  - MH12, MH23, MH13, MH123: For Multiplex-Heterogeneous Networks which connect a Drug Similarity Network with a Multiplex Disease Network
 
- MH_KFold_ROC_Final.R: Run K-Fold cross-validation for Heterogeneous/Multiplex-Heterogeneous Networks of Drugs and Diseases

- Summarize_AUC_Params.R: To investigate the prediction performance resulted from M_LOOCV_ROC_Final.R and MH_LOOCV_ROC_Final.R by parameters

- Summarize_AUC_ROC.R: To summarize and draw ROC curves

- HDR2_MH_Predict_Evidence.R: To predict and select top 10 highly ranked diseases for each drug, then find evidence supporting the promissing drug-disease associations 
