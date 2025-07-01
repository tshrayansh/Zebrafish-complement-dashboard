"""import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Load filtered data
df = pd.read_csv("data/filtered_expression.csv")

# Optional: Create plots folder
Path("plots").mkdir(exist_ok=True)

# === Plot 1: Gene counts per tissue ===
plt.figure(figsize=(10, 5))
df["Super Structure Name"].value_counts().sort_values().plot(kind="barh", color="skyblue")
plt.title("Complement Gene Expression Across Tissues")
plt.xlabel("Number of Annotations")
plt.tight_layout()
plt.savefig("plots/gene_expression_by_tissue.png")
plt.close()

# === Plot 2: Gene counts per developmental stage ===
plt.figure(figsize=(12, 5))
df["Start Stage"].value_counts().sort_values().plot(kind="barh", color="lightgreen")
plt.title("Complement Gene Expression by Developmental Stage")
plt.xlabel("Number of Annotations")
plt.tight_layout()
plt.savefig("plots/gene_expression_by_stage.png")
plt.close()

print("✅ Plots saved in 'plots/' folder.") """

# src/plot_expression.py

import pandas as pd
import matplotlib.pyplot as plt
import os

# Load filtered expression data
df = pd.read_csv('data/filtered_expression.csv')

# Make sure the plots/ directory exists
os.makedirs('plots/genes', exist_ok=True)

# Get list of unique genes
genes = df['Gene Symbol'].unique()

# Loop through each gene
for gene in genes:
    gene_df = df[df['Gene Symbol'] == gene]

    # --- Plot 1: Expression by tissue ---
    tissue_counts = gene_df['Sub Structure Name'].value_counts()
    if not tissue_counts.empty:
        plt.figure(figsize=(10, 4))
        tissue_counts.plot(kind='bar', color='steelblue')
        plt.title(f'Expression of {gene} by Tissue')
        plt.xlabel('Tissue')
        plt.ylabel('Count')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(f'plots/genes/{gene}_by_tissue.png')
        plt.close()

    # --- Plot 2: Expression by stage ---
    stage_counts = gene_df['Start Stage'].value_counts().sort_index()
    if not stage_counts.empty:
        plt.figure(figsize=(10, 4))
        stage_counts.plot(kind='bar', color='darkorange')
        plt.title(f'Expression of {gene} by Developmental Stage')
        plt.xlabel('Stage')
        plt.ylabel('Count')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(f'plots/genes/{gene}_by_stage.png')
        plt.close()

print(f"✅ Done! Plots saved in plots/genes/")
