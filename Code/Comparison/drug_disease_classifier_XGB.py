import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score, f1_score, accuracy_score
from sklearn.metrics import roc_curve, precision_recall_curve
from xgboost import XGBClassifier
from tqdm import tqdm
import itertools
import matplotlib.pyplot as plt
import os
import sys
import random

# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)

# Custom class to redirect print output to both console and file
class Tee:
    def __init__(self, filename):
        self.file = open(filename, 'w')
        self.stdout = sys.stdout
    
    def write(self, message):
        self.file.write(message)
        self.stdout.write(message)
    
    def flush(self):
        self.file.flush()
        self.stdout.flush()
    
    def close(self):
        self.file.close()

# Step 1: Load embeddings and drug-disease pairs, filter invalid pairs
def load_data(embeddings_file, drug_disease_file):
    print("Loading embeddings...")
    embeddings_df = pd.read_csv(embeddings_file)
    
    valid_drugs = set(embeddings_df[embeddings_df['type'] == 'drug']['node_id'])
    valid_diseases = set(embeddings_df[embeddings_df['type'] == 'disease']['node_id'])
    
    drug_emb = embeddings_df[embeddings_df['type'] == 'drug'].set_index('node_id')
    disease_emb = embeddings_df[embeddings_df['type'] == 'disease'].set_index('node_id')
    
    emb_cols = [col for col in embeddings_df.columns if col.startswith('dim_')]
    embedding_size = len(emb_cols)
    
    drug_emb = drug_emb[emb_cols].to_numpy()
    disease_emb = disease_emb[emb_cols].to_numpy()
    
    drugs = list(valid_drugs)
    diseases = list(valid_diseases)
    
    print(f"Number of drugs with embeddings: {len(drugs)}")
    print(f"Number of diseases with embeddings: {len(diseases)}")
    
    print("Loading positive drug-disease pairs...")
    drug_disease_df = pd.read_csv(drug_disease_file)
    
    positive_pairs = []
    total_pairs = len(drug_disease_df)
    for _, row in drug_disease_df.iterrows():
        drug, disease = row['drug'], row['disease']
        if drug in valid_drugs and disease in valid_diseases:
            positive_pairs.append((drug, disease))
    
    positive_pairs = set(positive_pairs)
    skipped_pairs = total_pairs - len(positive_pairs)
    print(f"Number of positive pairs loaded: {len(positive_pairs)}")
    print(f"Number of pairs skipped (missing embeddings): {skipped_pairs}")
    
    if not positive_pairs:
        raise ValueError("No valid positive pairs found after filtering. Check embeddings_file and drug_disease_file.")
    
    return drugs, diseases, drug_emb, disease_emb, positive_pairs, embedding_size

# Step 2: Generate feature vectors and labels with all possible pairs
def generate_features_labels(drugs, diseases, drug_emb, disease_emb, positive_pairs, embedding_size):
    print("Generating feature vectors and labels...")
    features = []
    labels = []
    pairs = []
    
    drug_idx = {drug: i for i, drug in enumerate(drugs)}
    disease_idx = {disease: i for i, disease in enumerate(diseases)}
    
    # Generate all possible pairs and label them
    for drug, disease in tqdm(list(itertools.product(drugs, diseases)), desc="Generating pairs"):
        pair = (drug, disease)
        pairs.append(pair)
        # Label: 1 if positive, 0 if negative
        label = 1 if pair in positive_pairs else 0
        labels.append(label)
        
        drug_vec = drug_emb[drug_idx[drug]]
        disease_vec = disease_emb[disease_idx[disease]]
        feature_vec = np.concatenate([drug_vec, disease_vec])
        features.append(feature_vec)
    
    return np.array(features), np.array(labels), pairs

# Step 3: Train and evaluate XGBoost with 10-fold cross-validation
def evaluate_model(features, labels, base_name):
    print("Training and evaluating XGBoost model...")
    # Compute scale_pos_weight to handle class imbalance
    neg_count = len(labels) - sum(labels)
    pos_count = sum(labels)
    scale_pos_weight = neg_count / pos_count if pos_count > 0 else 1
    
    xgb = XGBClassifier(
        n_estimators=100,
        max_depth=5,
        learning_rate=0.1,
        random_state=42,
        n_jobs=-1,
        eval_metric='logloss',
        scale_pos_weight=scale_pos_weight  # Handle class imbalance
    )
    
    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    
    auroc_scores = []
    auprc_scores = []
    f1_scores = []
    accuracy_scores = []
    
    for fold, (train_idx, test_idx) in enumerate(skf.split(features, labels), 1):
        print(f"Processing fold {fold}/10...")
        X_train, X_test = features[train_idx], features[test_idx]
        y_train, y_test = labels[train_idx], labels[test_idx]
        
        xgb.fit(X_train, y_train)
        y_pred_proba = xgb.predict_proba(X_test)[:, 1]
        
        auroc = roc_auc_score(y_test, y_pred_proba)
        auprc = average_precision_score(y_test, y_pred_proba)
        y_pred = (y_pred_proba >= 0.5).astype(int)
        f1 = f1_score(y_test, y_pred)
        accuracy = accuracy_score(y_test, y_pred)
        
        auroc_scores.append(auroc)
        auprc_scores.append(auprc)
        f1_scores.append(f1)
        accuracy_scores.append(accuracy)
        
        print(f"Fold {fold} - AUROC: {auroc:.4f}, AUPRC: {auprc:.4f}, F1: {f1:.4f}, Accuracy: {accuracy:.4f}")
    
    auroc_mean, auroc_std = np.mean(auroc_scores), np.std(auroc_scores)
    auprc_mean, auprc_std = np.mean(auprc_scores), np.std(auprc_scores)
    f1_mean, f1_std = np.mean(f1_scores), np.std(f1_scores)
    accuracy_mean, accuracy_std = np.mean(accuracy_scores), np.std(accuracy_scores)
    
    print("\nFinal Results:")
    print(f"AUROC: {auroc_mean:.4f} ± {auroc_std:.4f}")
    print(f"AUPRC: {auprc_mean:.4f} ± {auprc_std:.4f}")
    print(f"F1-score: {f1_mean:.4f} ± {f1_std:.4f}")
    print(f"Accuracy: {accuracy_mean:.4f} ± {accuracy_std:.4f}")
    
    return auroc_mean, auroc_std, auprc_mean, auprc_std, f1_mean, f1_std, accuracy_mean, accuracy_std

# Step 4: Sample data for t-SNE visualization
def sample_for_tsne(features, labels, sample_size=10000):
    if len(features) <= sample_size:
        return features, labels
    indices = random.sample(range(len(features)), sample_size)
    return features[indices], labels[indices]

# Main function
cl = "XGB"
def main():
    embedding_size = 64
    epochs = 50
    emb_methods = ["gat"]#gat/transformer
    drug_simNets = ["DrugSimNet_PREDICT"]#["DrugSimNet_PREDICT", "DrugSimNet_CHEM"]
    disease_simNets = ["DiseaseSimNet_OHG"]#["DiseaseSimNet_OHG", "DiseaseSimNet_OMIM", "DiseaseSimNet_HPO", "DiseaseSimNet_GeneNet"]
    
    # Collect results for summary
    results = []
    
    
    for emb in emb_methods:
        for drug_sim in drug_simNets:
            for disease_sim in disease_simNets:
                # Construct file paths
                embeddings_file = f"./Results/{emb}_{drug_sim}_{disease_sim}_d_{embedding_size}_e_{epochs}.csv"
                drug_disease_file = "../Data/Drug2Disease_Name_PREDICT_ID.txt_BinaryInteraction.csv"
                
                base_name = "./Results/" + os.path.splitext(os.path.basename(embeddings_file))[0] + f"_NotBalanced_{cl}"
                
                print(f"\nProcessing pair:")
                print(f"Embeddings file: {embeddings_file}")
                print(f"Drug-disease file: {drug_disease_file}")
                
                output_file = f'{base_name}_output.txt'
                tee = Tee(output_file)
                sys.stdout = tee
                
                try:
                    drugs, diseases, drug_emb, disease_emb, positive_pairs, emb_size = load_data(
                        embeddings_file, drug_disease_file
                    )
                    
                    features, labels, pairs = generate_features_labels(
                        drugs, diseases, drug_emb, disease_emb, positive_pairs, emb_size
                    )
                    
                    # Sample data for t-SNE
                    features_tsne, labels_tsne = sample_for_tsne(features, labels)
                    from sklearn.manifold import TSNE
                    import seaborn as sns
                    tsne = TSNE(n_components=2, random_state=42)
                    embeddings_2d = tsne.fit_transform(features_tsne)
                    sns.scatterplot(x=embeddings_2d[:, 0], y=embeddings_2d[:, 1], hue=labels_tsne)
                    plt.savefig(f'{base_name}_tsne.png')
                    plt.close()
                    
                    auroc_mean, auroc_std, auprc_mean, auprc_std, f1_mean, f1_std, accuracy_mean, accuracy_std = evaluate_model(features, labels, base_name)
                    
                    # Collect results
                    results.append({
                        'emb': emb,
                        'drug_sim': drug_sim,
                        'disease_sim': disease_sim,
                        'auroc_mean': auroc_mean,
                        'auroc_std': auroc_std,
                        'auprc_mean': auprc_mean,
                        'auprc_std': auprc_std,
                        'f1_mean': f1_mean,
                        'f1_std': f1_std,
                        'accuracy_mean': accuracy_mean,
                        'accuracy_std': accuracy_std
                    })
                
                except Exception as e:
                    print(f"Error processing emb={emb}, drug_sim={drug_sim}, disease_sim={disease_sim}: {str(e)}")
                
                finally:
                    sys.stdout = tee.stdout
                    tee.close()
    
                # Save summary results to a CSV file
                summary_df = pd.DataFrame(results)
                summary_file = base_name + "_summary_metrics.csv"
                summary_df.to_csv(summary_file, index=False)
                print(f"\nSummary results saved to {summary_file}")

# Execute
if __name__ == "__main__":
    main()
