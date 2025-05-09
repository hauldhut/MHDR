import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import TransformerConv
import networkx as nx
import os

emb_method = "transformer"

# Set random seed for reproducibility
torch.manual_seed(42)
np.random.seed(42)

# Step 1: Load the drug and disease similarity networks separately
def load_networks(drug_sim_file="../Data/DiseaseSimNet_OMIM.txt", 
                  disease_sim_file="../Data/HomomiRNAIntegrated.txt"):
    print(drug_sim_file)
    print(disease_sim_file)
    
    print("Loading drug similarity network...")
    drug_sim_df = pd.read_csv(drug_sim_file, sep='\t', header=None, names=['drug1', 'sim', 'drug2'])
    drugs = set(drug_sim_df['drug1']).union(set(drug_sim_df['drug2']))
    
    drug_G = nx.Graph()
    for drug in tqdm(drugs, desc="Adding drug nodes"):
        drug_G.add_node(drug, type="drug")
    
    for _, row in tqdm(drug_sim_df.iterrows(), total=len(drug_sim_df), desc="Adding drug-drug edges"):
        drug_G.add_edge(row['drug1'], row['drug2'], weight=float(row['sim']))
    
    print("Loading disease similarity network...")
    disease_sim_df = pd.read_csv(disease_sim_file, sep='\t', header=None, names=['disease1', 'sim', 'disease2'])
    diseases = set(disease_sim_df['disease1']).union(set(disease_sim_df['disease2']))
    
    disease_G = nx.Graph()
    for disease in tqdm(diseases, desc="Adding disease nodes"):
        disease_G.add_node(disease, type="disease")
    
    for _, row in tqdm(disease_sim_df.iterrows(), total=len(disease_sim_df), desc="Adding disease-disease edges"):
        disease_G.add_edge(row['disease1'], row['disease2'], weight=float(row['sim']))
    
    return drug_G, disease_G, list(drugs), list(diseases)

# Step 2: Convert networkx graph to PyTorch Geometric Data object
def graph_to_pyg_data(G, embedding_size=128):
    print(f"Graph has {len(G.nodes())} nodes and {len(G.edges())} edges")
    node_to_idx = {node: idx for idx, node in enumerate(G.nodes())}
    print("Node index mapping created")
    
    for node in G.nodes():
        if not isinstance(node, (str, int)):
            print(f"Invalid node ID: {node} (type: {type(node)})")
    
    edge_index = []
    edge_weight = []
    for u, v, data in G.edges(data=True):
        edge_index.append([node_to_idx[u], node_to_idx[v]])
        edge_index.append([node_to_idx[v], node_to_idx[u]])
        weight = data.get('weight', 1.0)
        if not isinstance(weight, (int, float)) or weight <= 0 or np.isnan(weight) or np.isinf(weight):
            print(f"Invalid weight found for edge ({u}, {v}): {weight}, using default 1.0")
            weight = 1.0
        edge_weight.append(weight)
        edge_weight.append(weight)
    print(f"Edge index and weights created: {len(edge_index)} edges")
    
    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    print("Edge index tensor created")
    edge_weight = torch.tensor(edge_weight, dtype=torch.float)
    print("Edge weight tensor created")
    
    num_nodes = len(G.nodes())
    x = torch.ones((num_nodes, embedding_size), dtype=torch.float)
    print("Node features tensor created")
    
    data = Data(x=x, edge_index=edge_index)
    print("PyTorch Geometric Data object created without edge weights")
    data.edge_weight = edge_weight
    print("Edge weights added to Data object")
    return data, node_to_idx

# Step 3: Define a custom WeightedTransformerConv layer that supports edge weights
class WeightedTransformerConv(TransformerConv):
    def __init__(self, in_channels, out_channels, heads=1, concat=True, dropout=0.0, edge_dim=1):
        super(WeightedTransformerConv, self).__init__(in_channels, out_channels, heads=heads, concat=concat, dropout=dropout, edge_dim=edge_dim)
    
    def forward(self, x, edge_index, edge_weight=None):
        # TransformerConv uses edge_attr for edge features
        if edge_weight is not None:
            edge_attr = edge_weight.view(-1, 1)  # Reshape to (num_edges, 1) for edge_dim=1
        else:
            edge_attr = None
        return super().forward(x, edge_index, edge_attr=edge_attr)

# Step 4: Define GraphTransformer model with WeightedTransformerConv
class GraphTransformer(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim, heads=4, dropout=0.5):
        super(GraphTransformer, self).__init__()
        self.conv1 = WeightedTransformerConv(input_dim, hidden_dim, heads=heads, concat=True, dropout=dropout, edge_dim=1)
        # After concat, the output dimension of conv1 is hidden_dim * heads
        self.conv2 = WeightedTransformerConv(hidden_dim * heads, output_dim, heads=1, concat=False, dropout=dropout, edge_dim=1)
        self.dropout = dropout
    
    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_weight
        x = self.conv1(x, edge_index, edge_weight=edge_weight)
        x = F.relu(x)
        x = F.dropout(x, p=self.dropout, training=self.training)
        x = self.conv2(x, edge_index, edge_weight=edge_weight)
        return x

# Step 5: Train GraphTransformer model with unsupervised loss (graph reconstruction)
def train_transformer(data, embedding_size=128, hidden_dim=64, epochs=20, lr=0.01, node_type="nodes"):
    device = torch.device('cpu')
    data = data.to(device)
    
    model = GraphTransformer(input_dim=data.x.shape[1], hidden_dim=hidden_dim, output_dim=embedding_size, heads=4).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    
    model.train()
    for epoch in tqdm(range(epochs), desc=f"Training GraphTransformer for {node_type}", ascii=False, smoothing=0.1):
        optimizer.zero_grad()
        embeddings = model(data)
        
        u, v = data.edge_index[0], data.edge_index[1]
        edge_scores = torch.sum(embeddings[u] * embeddings[v], dim=-1)
        pos_loss = -torch.log(torch.sigmoid(edge_scores) + 1e-15).mean()
        
        num_neg_samples = data.edge_index.shape[1]
        neg_u = torch.randint(0, data.num_nodes, (num_neg_samples,), device=device)
        neg_v = torch.randint(0, data.num_nodes, (num_neg_samples,), device=device)
        neg_mask = torch.ones(num_neg_samples, dtype=torch.bool, device=device)
        for i in range(num_neg_samples):
            u, v = neg_u[i], neg_v[i]
            edge_exists = (data.edge_index[0].eq(u) & data.edge_index[1].eq(v)).any().item()
            if edge_exists:
                neg_mask[i] = 0
        neg_u, neg_v = neg_u[neg_mask.bool()], neg_v[neg_mask.bool()]
        
        # Compute loss
        loss = pos_loss
        if len(neg_u) > 0:
            neg_scores = torch.sum(embeddings[neg_u] * embeddings[neg_v], dim=-1)
            neg_loss = -torch.log(1 - torch.sigmoid(neg_scores) + 1e-15).mean()
            loss = loss + neg_loss
        else:
            print(f"Epoch {epoch+1}: No negative samples after filtering, skipping neg_loss")
        
        loss.backward()
        optimizer.step()
        
        if (epoch + 1) % max(1, epochs // 5) == 0:
            print(f"Epoch {epoch+1}/{epochs}, Loss: {loss.item():.4f}")
    
    model.eval()
    with torch.no_grad():
        embeddings = model(data).cpu().numpy()
    return embeddings

# Step 6: Extract and process embeddings
def extract_embeddings(embeddings, nodes, node_to_idx, node_type):
    embeddings_dict = {}
    for node in tqdm(nodes, desc=f"Extracting embeddings for {node_type}"):
        idx = node_to_idx[node]
        embeddings_dict[node] = embeddings[idx]
    return embeddings_dict

# Step 7: Save embeddings
def save_embeddings(drug_embeddings, disease_embeddings, drug_G, disease_G, output_file="embeddings.csv"):
    print("Saving embeddings...")
    data = []
    
    for node, emb in drug_embeddings.items():
        node_type = drug_G.nodes[node]["type"]
        row = [node, node_type] + emb.tolist()
        data.append(row)
    
    for node, emb in disease_embeddings.items():
        node_type = disease_G.nodes[node]["type"]
        row = [node, node_type] + emb.tolist()
        data.append(row)
    
    embedding_size = len(next(iter(drug_embeddings.values())))
    columns = ['node_id', 'type'] + [f'dim_{i+1}' for i in range(embedding_size)]
    df = pd.DataFrame(data, columns=columns)
    df.to_csv(output_file, index=False)
    print(f"Embeddings saved to {output_file}")

# Main function
def main(drug_sim_file, disease_sim_file, embedding_size=128, epochs=200, output_file="embeddings.csv"):
    drug_G, disease_G, drugs, diseases = load_networks(drug_sim_file, disease_sim_file)
    
    print("Converting drug graph to PyTorch Geometric format...")
    drug_data, drug_node_to_idx = graph_to_pyg_data(drug_G, embedding_size)
    print("Converting disease graph to PyTorch Geometric format...")
    disease_data, disease_node_to_idx = graph_to_pyg_data(disease_G, embedding_size)
    
    drug_embeddings = train_transformer(drug_data, embedding_size=embedding_size, epochs=epochs, node_type="drugs")
    disease_embeddings = train_transformer(disease_data, embedding_size=embedding_size, epochs=epochs, node_type="diseases")
    
    drug_embeddings_dict = extract_embeddings(drug_embeddings, drugs, drug_node_to_idx, "drugs")
    disease_embeddings_dict = extract_embeddings(disease_embeddings, diseases, disease_node_to_idx, "diseases")
    
    save_embeddings(drug_embeddings_dict, disease_embeddings_dict, drug_G, disease_G, output_file=output_file)
    
    return drug_embeddings_dict, disease_embeddings_dict

# Execute and get embeddings
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Learn Transformer-based GNN embeddings from drug and disease similarity networks")
    parser.add_argument('--drug_sim_file', type=str, default=f"../Data/DrugSimNet_PREDICT.txt", help='Drug sim net')
    parser.add_argument('--disease_sim_file', type=str, default=f"../Data/DiseaseSimNet_OHG.txt", help='Disease sim net')
    parser.add_argument('--embedding_size', type=int, default=64, help='Size of the embedding vectors')
    parser.add_argument('--epochs', type=int, default=50, help='Number of training epochs for GAT')
    parser.add_argument('--output_file', type=str, default=f'./Results/{emb_method}', help='Output file prefix for embeddings')
    args = parser.parse_args()

    print(args.drug_sim_file)
    print(args.disease_sim_file)
    print(args.embedding_size)
    print(args.epochs)
    print(args.output_file)

    drug_sim = os.path.splitext(os.path.basename(args.drug_sim_file))[0]
    disease_sim = os.path.splitext(os.path.basename(args.disease_sim_file))[0]
    
    drug_embeddings, disease_embeddings = main(
        drug_sim_file=args.drug_sim_file,
        disease_sim_file=args.disease_sim_file,
        embedding_size=args.embedding_size,
        epochs=args.epochs,
        output_file=args.output_file + "_" + drug_sim + "_" + disease_sim + "_d_" + str(args.embedding_size) + "_e_" + str(args.epochs) + ".csv"
    )
    for node, emb in list(drug_embeddings.items())[:5]:
        print(f"Drug Node: {node}, Embedding: {emb[:5]}...")
    for node, emb in list(disease_embeddings.items())[:5]:
        print(f"Disease Node: {node}, Embedding: {emb[:5]}...")
