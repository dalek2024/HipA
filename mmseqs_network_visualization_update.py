import subprocess
import os
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from networkx.algorithms.community import greedy_modularity_communities

def run_mmseqs(fasta_file, output_dir, max_seqs):
    # Define MMseqs2 commands and parameters
    db_name = os.path.join(output_dir, "sequenceDB")
    result_name = os.path.join(output_dir, "searchResult")
    tmp_dir = os.path.join(output_dir, "tmp")
    tsv_file = os.path.join(output_dir, "result.tsv")
    filtered_tsv_file = os.path.join(output_dir, "filtered_result.tsv")

    # Create database from fasta file
    subprocess.run(["mmseqs", "createdb", fasta_file, db_name], check=True)

    # Run MMseqs2 search
    subprocess.run(["mmseqs", "search", db_name, db_name, result_name, tmp_dir, "--max-seqs", str(max_seqs)], check=True)

    # Convert search results to tsv format with additional columns for coverage
    subprocess.run([
        "mmseqs", "convertalis", db_name, db_name, result_name, tsv_file,
        "--format-output", "query,target,fident,qcov,tcov"
    ], check=True)

    print(f"Results written to {tsv_file}")

    # Filter based on coverage
    #network_data = pd.read_csv(tsv_file, sep='\t', header=None, names=['Protein1', 'Protein2', 'Fident','Qcov','Tcov'])
    #filtered_data = network_data[(network_data['Qcov'] >= coverage) & (network_data['Tcov'] >= coverage)]
    #filtered_data.to_csv(filtered_tsv_file, sep='\t', index=False, header=False)

    #print(f"Filtered results written to {filtered_tsv_file}")
    return tsv_file

def read_label_node(file):
    label_info = {}
    with open(file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            node = parts[0]
            label = parts[1] if len(parts) > 1 else ""
            color = parts[2] if len(parts) > 2 else 'red'
            size = int(parts[3]) if len(parts) > 3 else 300
            label_info[node] = {'label': label, 'color': color, 'size': size}
    return label_info

def create_network(tsv_file, threshold, output_dir, coverage, label_node_file=None):
    # Load the network data
    network_data = pd.read_csv(tsv_file, sep='\t', header=None, names=['Protein1', 'Protein2', 'RawScore', 'Qcov','Tcov'])

    # Read label node information if provided
    label_nodes = read_label_node(label_node_file) if label_node_file else {}

    # Create the graph
    G = nx.Graph()

    # Add edges to the graph based on the threshold and coverage
    for _, row in network_data.iterrows():
        if row['RawScore'] >= threshold and row['Qcov'] >= coverage and row['Tcov'] >= coverage and row['Protein1'] != row['Protein2']:  # Remove self-loops
            G.add_edge(row['Protein1'], row['Protein2'], weight=row['RawScore'])

    # Detect communities using the modularity-based algorithm
    communities = greedy_modularity_communities(G)

    # Map nodes to their communities
    community_map = {}
    community_members = {}
    for i, community in enumerate(communities):
        community_members[i] = list(community)
        for node in community:
            community_map[node] = i

    # Output the detected communities to TSV files and count members
    community_members_file = os.path.join(output_dir, f'community_members_threshold_{threshold}.tsv')
    community_counts_file = os.path.join(output_dir, f'community_counts_threshold_{threshold}.tsv')
    community_threshold_file = os.path.join(output_dir, f'community_threshold_{threshold}.txt')
    network_image_file = os.path.join(output_dir, f'protein_similarity_network_clusters_threshold_{threshold}.png')

    with open(community_members_file, 'w') as f:
        f.write('Community\tProtein\n')
        for i, community in community_members.items():
            for node in community:
                f.write(f"{i}\t{node}\n")

    with open(community_counts_file, 'w') as f:
        f.write('Community\tMemberCount\n')
        for i, community in community_members.items():
            f.write(f"{i}\t{len(community)}\n")

    with open(community_threshold_file, 'w') as f:
        f.write(f"Threshold used for community detection: {threshold}\n")

    print(f"Threshold used for community detection: {threshold}")

    # Define node color based on clusters
    node_color = [community_map[node] for node in G]

    # Set a constant node size
    node_size = 300  # Adjust this value to make nodes larger or smaller

    # Remove edges between clusters
    edges_to_remove = [(u, v) for u, v in G.edges() if community_map[u] != community_map[v]]
    G.remove_edges_from(edges_to_remove)

    # Define layout (spring layout, spread clusters)
    pos = nx.spring_layout(G, k=0.2, iterations=50)

    # Draw the network
    plt.figure(figsize=(20, 20))
    nx.draw_networkx_edges(G, pos, alpha=0.5, edge_color='grey')
    nodes = nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color=node_color, cmap=plt.cm.tab20, alpha=0.8)

    # Highlight and label specific nodes
    for node, info in label_nodes.items():
        if node in G:
            nx.draw_networkx_nodes(G, pos, nodelist=[node], node_size=info['size'], node_color=info['color'], label=info['label'])
            nx.draw_networkx_labels(G, pos, labels={node: info['label']}, font_size=40, font_color='black')

    # Add color bar
    plt.colorbar(nodes)

    # Remove axis
    plt.axis('off')

    # Save the figure
    plt.savefig(network_image_file, format='PNG')
    plt.show()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate network file using MMseqs2 and visualize it")
    parser.add_argument("fasta_file", help="Input FASTA file")
    parser.add_argument("output_dir", help="Output directory for MMseqs2 results")
    parser.add_argument("--threshold", type=float, nargs='+', required=True, help="Space-separated list of score thresholds for network edges")
    parser.add_argument("--max-seqs", type=int, default=300, help="Maximum number of sequences to return per query")
    parser.add_argument("--coverage", type=float, required=True, help="Minimum coverage threshold for filtering alignments")
    parser.add_argument("--label-node", type=str, required=True, help="TSV file with nodes to label")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    thresholds = args.threshold
    tsv_file = run_mmseqs(args.fasta_file, args.output_dir, args.max_seqs)

    for threshold in thresholds:
        create_network(tsv_file, threshold, args.output_dir, args.coverage, args.label_node)
