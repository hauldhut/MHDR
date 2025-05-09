import numpy as np
import networkx as nx
from gensim.models import Word2Vec
from gensim.models.callbacks import CallbackAny2Vec
import random
from collections import defaultdict
import pandas as pd
from tqdm import tqdm
import argparse
import multiprocessing
import os


emb_method = "metapath2vec"
# disease_sim = ""
# drug_sim = ""
# Set random seed for reproducibility
random.seed(42)
np.random.seed(42)

# Custom callback to show training progress per epoch
class EpochProgress(CallbackAny2Vec):
    def __init__(self, epochs):
        self.epochs = epochs
        self.epoch = 0
        self.pbar = None
    
    def on_epoch_begin(self, model):
        if self.pbar is None:
            self.pbar = tqdm(total=self.epochs, desc="Training Word2Vec epochs")
    
    def on_epoch_end(self, model):
        self.pbar.update(1)
        self.epoch += 1
        if self.epoch >= self.epochs:
            self.pbar.close()

# Step 1: Load the drug and disease similarity networks
def load_network(drug_sim_file="../Data/DrugSimNet_PREDICT.txt", 
                 disease_sim_file="../Data/DiseaseSimNet_OHG.txt"):
    
    print(drug_sim_file)
    print(disease_sim_file)

    
    G = nx.Graph()
    
    # Load drug similarity network
    print("Loading drug similarity network...")
    drug_sim_df = pd.read_csv(drug_sim_file, sep='\t', header=None, names=['drug1', 'sim', 'drug2'])
    drugs = set(drug_sim_df['drug1']).union(set(drug_sim_df['drug2']))
    
    # Add drug nodes with type
    for drug in tqdm(drugs, desc="Adding drug nodes"):
        G.add_node(drug, type="drug")
    
    # Add drug-drug edges with weights
    for _, row in tqdm(drug_sim_df.iterrows(), total=len(drug_sim_df), desc="Adding drug-drug edges"):
        G.add_edge(row['drug1'], row['drug2'], weight=float(row['sim']))
    
    # Load disease similarity network
    print("Loading disease similarity network...")
    disease_sim_df = pd.read_csv(disease_sim_file, sep='\t', header=None, names=['disease1', 'sim', 'disease2'])
    diseases = set(disease_sim_df['disease1']).union(set(disease_sim_df['disease2']))
    
    # Add disease nodes with type
    for disease in tqdm(diseases, desc="Adding disease nodes"):
        G.add_node(disease, type="disease")
    
    # Add disease-disease edges with weights
    for _, row in tqdm(disease_sim_df.iterrows(), total=len(disease_sim_df), desc="Adding disease-disease edges"):
        G.add_edge(row['disease1'], row['disease2'], weight=float(row['sim']))
    
    return G, list(drugs), list(diseases)

# Step 2: Generate metapath-guided random walks with weights
def metapath_random_walk(G, metapaths, num_walks=100, walk_length=10):
    walks = []
    
    # Create adjacency list by node type
    adj_by_type = defaultdict(lambda: defaultdict(list))
    for node in tqdm(G.nodes(), desc="Building adjacency list"):
        node_type = G.nodes[node]["type"]
        for neighbor in G.neighbors(node):
            neighbor_type = G.nodes[neighbor]["type"]
            adj_by_type[node_type][neighbor].append(node)
    
    # Generate walks
    total_iterations = num_walks * len(G.nodes()) * len(metapaths)
    with tqdm(total=total_iterations, desc="Generating random walks") as pbar:
        for _ in range(num_walks):
            for start_node in G.nodes():
                for metapath in metapaths:
                    current_node = start_node
                    walk = [current_node]
                    current_type = G.nodes[current_node]["type"]
                    
                    # Ensure the starting node matches the first metapath type
                    if current_type != metapath[0]:
                        pbar.update(1)
                        continue
                    
                    # Generate one walk
                    for step in range(walk_length - 1):
                        # Get the expected type for the next node
                        next_type_idx = (step + 1) % len(metapath)
                        next_type = metapath[next_type_idx]
                        
                        # Get neighbors of the current node
                        neighbors = list(G.neighbors(current_node))
                        if not neighbors:
                            break
                        
                        # Filter neighbors by the expected type
                        valid_neighbors = [
                            n for n in neighbors if G.nodes[n]["type"] == next_type
                        ]
                        if not valid_neighbors:
                            break
                        
                        # Get weights for valid neighbors
                        weights = []
                        for neighbor in valid_neighbors:
                            weight = G[current_node][neighbor].get('weight', 1.0)
                            weights.append(weight)
                        
                        # Normalize weights to probabilities
                        weights = np.array(weights)
                        if weights.sum() == 0:
                            weights = np.ones_like(weights)  # Fallback to uniform if all weights are 0
                        probabilities = weights / weights.sum()
                        
                        # Choose a neighbor based on weights
                        current_node = np.random.choice(valid_neighbors, p=probabilities)
                        walk.append(current_node)
                    
                    if len(walk) > 1:  # Only add non-trivial walks
                        walks.append(walk)
                    
                    pbar.update(1)
    
    return walks

# Step 3: Train metapath2vec model
def train_metapath2vec(walks, embedding_size=128, window=5, epochs=5):
    # Convert walks to strings
    str_walks = [[str(n) for n in walk] for walk in tqdm(walks, desc="Converting walks to strings")]
    
    # Train Word2Vec model with multiple workers
    print("Training Word2Vec model...")
    workers = min(16, multiprocessing.cpu_count())
    model = Word2Vec(
        str_walks,
        vector_size=embedding_size,
        window=window,
        min_count=0,
        sg=1,  # Skip-gram
        workers=workers,
        epochs=epochs,
        callbacks=[EpochProgress(epochs)]
    )
    return model

# Step 4: Extract and save embeddings
def save_embeddings(embeddings, G, output_file="embeddings.csv"):
    print("Saving embeddings...")
    data = []
    for node, emb in embeddings.items():
        node_type = G.nodes[node]["type"]
        row = [node, node_type] + emb.tolist()
        data.append(row)
    
    columns = ['node_id', 'type'] + [f'dim_{i+1}' for i in range(len(emb))]
    df = pd.DataFrame(data, columns=columns)
    df.to_csv(output_file, index=False)
    print(f"Embeddings saved to {output_file}")

# Step 5: Extract embeddings
def extract_embeddings(model, drugs, diseases):
    embeddings = {}
    for node in tqdm(drugs + diseases, desc="Extracting embeddings"):
        embeddings[node] = model.wv[str(node)]
    return embeddings

# Main function
def main(drug_sim_file, disease_sim_file, embedding_size=128, epochs=5, output_file="embeddings.csv"):
    # Load network from files
    G, drugs, diseases = load_network(drug_sim_file, disease_sim_file)
    
    # Define metapaths for similarity networks
    metapaths = [
        ["drug", "drug", "drug"],
        ["disease", "disease", "disease"]
    ]
    
    # Generate random walks
    walks = metapath_random_walk(G, metapaths)
    
    # Train metapath2vec
    model = train_metapath2vec(walks, embedding_size=embedding_size, epochs=epochs)
    
    # Extract embeddings
    embeddings = extract_embeddings(model, drugs, diseases)
    
    # Save embeddings
    
    save_embeddings(embeddings, G, output_file=output_file)
    
    return embeddings

# Execute and get embeddings
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Learn embeddings from drug and disease similarity networks")
    parser.add_argument('--drug_sim_file', type=str, default=f"../Data/DrugSimNet_PREDICT.txt", help='Drug sim net')
    parser.add_argument('--disease_sim_file', type=str, default=f"../Data/DiseaseSimNet_OHG.txt", help='Disease sim net')
    parser.add_argument('--embedding_size', type=int, default=64, help='Size of the embedding vectors')
    parser.add_argument('--epochs', type=int, default=50, help='Number of training epochs for Word2Vec')
    parser.add_argument('--output_file', type=str, default=f'./Results/{emb_method}', help='Output file for embeddings')
    args = parser.parse_args()
    
    drug_sim = os.path.splitext(os.path.basename(args.drug_sim_file))[0]
    disease_sim = os.path.splitext(os.path.basename(args.disease_sim_file))[0]

    embeddings = main(
        drug_sim_file=args.drug_sim_file,
        disease_sim_file=args.disease_sim_file,
        embedding_size=args.embedding_size,
        epochs=args.epochs,
        output_file=args.output_file + "_" + drug_sim + "_" + disease_sim + "_d_" + str(args.embedding_size) + "_e_" + str(args.epochs) + ".csv"
    )
    # Print embeddings for demonstration
    for node, emb in embeddings.items():
        print(f"Node: {node}, Embedding: {emb[:5]}...")  # Show first 5 dimensions
    
    # for node, emb in list(drug_embeddings.items())[:5]:
    #     print(f"Drug Node: {node}, Embedding: {emb[:5]}...")
    # for node, emb in list(disease_embeddings.items())[:5]:
    #     print(f"Disease Node: {node}, Embedding: {emb[:5]}...")
