"""import pandas as pd
import matplotlib.pyplot as plt
import os

# Load filtered data
df = pd.read_csv("data/filtered_expression.csv")

# Clean column if needed
if 'Sub Structure Name' not in df.columns:
    raise KeyError("❌ Column 'Sub Structure Name' not found in CSV")

# Group by gene symbol
unique_genes = df["Gene Symbol"].dropna().unique()

# Make plots directory
os.makedirs("plots/genes_by_organ", exist_ok=True)

for gene in unique_genes:
    gene_df = df[df["Gene Symbol"] == gene]

    # Count expression per organ
    organ_counts = gene_df["Sub Structure Name"].value_counts()

    # Skip empty plots
    if organ_counts.empty:
        continue

    # Plot
    plt.figure(figsize=(10, 4))
    organ_counts.plot(kind="bar", color="teal")
    plt.title(f"🧬 {gene} Expression by Organ")
    plt.ylabel("Annotation Count")
    plt.xlabel("Organ (Sub Structure)")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(f"plots/genes_by_organ/{gene}_organ_expression.png")
    plt.close()

print("✅ All plots saved to plots/genes_by_organ/") """
# src/plot_by_organ.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load the filtered dataset
df = pd.read_csv("data/filtered_expression.csv")

# ✅ FILTER: Only adult-stage entries
df = df[(df["Start Stage"] == "Adult") & (df["End Stage"] == "Adult")]

# ✅ Drop missing values
df = df.dropna(subset=["Gene Symbol", "Sub Structure Name"])

# ✅ Create a pivot table: rows = genes, columns = organs, values = expression counts
heatmap_data = df.pivot_table(
    index="Gene Symbol",
    columns="Sub Structure Name",
    aggfunc="size",
    fill_value=0
)

# ✅ Plot heatmap
plt.figure(figsize=(24, 14))
sns.heatmap(heatmap_data, cmap="YlGnBu", linewidths=0.4, linecolor='gray')

plt.title("🧬 Complement Gene Expression in Organs (Adult Stage)", fontsize=18)
plt.xlabel("Organ (Sub Structure Name)", fontsize=12)
plt.ylabel("Gene Symbol", fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.yticks(fontsize=8)

# ✅ Save plot
os.makedirs("plots/summary", exist_ok=True)
plt.tight_layout()
plt.savefig("plots/summary/all_genes_by_organ_adult_heatmap.png", dpi=300)
plt.close()
