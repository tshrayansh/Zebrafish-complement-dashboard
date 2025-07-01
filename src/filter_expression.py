import pandas as pd
from pathlib import Path

# === CONFIGURATION ===
DATA_PATH = "data/wildtype-expression_fish_2025.06.30.txt"
OUTPUT_PATH = "data/filtered_expression.csv"

# ZFIN data column headers
column_names = [
    "Gene ID", "Gene Symbol", "Fish Name",
    "Super Structure ID", "Super Structure Name",
    "Sub Structure ID", "Sub Structure Name",
    "Start Stage", "End Stage", "Assay",
    "Assay MMO ID", "Publication ID",
    "Probe ID", "Antibody ID", "Fish ID"
]

# === Complement system target genes ===
target_genes = {
    "c3a.1", "c3a.2", "c3a.3", "c3a.4", "c3a.5", "c3a.6",
    "c3b.1", "c3b.2", "c5", "c5ar1", "c6.1", "c7b",
    "c8a", "c8b", "c8g", "c8gl", "c9",
    "cfb", "cfbl", "cfd", "cfh", "cfhl1", "cfhl2", "cfhl3", "cfhl4", "cfi", "cfp"
}

# === Load data ===
print("📥 Loading expression data...")
df = pd.read_csv(DATA_PATH, sep="\t", header=None, names=column_names, skiprows=2)

# === Standardize gene symbol casing ===
df["Gene Symbol"] = df["Gene Symbol"].astype(str).str.strip().str.lower()

# === Filter expression data for target complement genes ===
filtered_df = df[df["Gene Symbol"].isin(target_genes)]

# === Save result ===
Path(OUTPUT_PATH).parent.mkdir(parents=True, exist_ok=True)
filtered_df.to_csv(OUTPUT_PATH, index=False)

# === Summary ===
print("✅ Extracted", len(filtered_df), "rows for", len(target_genes), "complement genes")
print("📁 Saved to", OUTPUT_PATH)
