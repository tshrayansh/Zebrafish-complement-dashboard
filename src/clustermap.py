import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --- Configuration ---
DATA_FILE = "data/filtered_expression.csv"
OUTPUT_DIR = "plots/summary"
OUTPUT_FILENAME = "genes_vs_super_structure_clustermap_overall.png"

# Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Data Loading and Initial Preprocessing ---
df = None
try:
    df = pd.read_csv(DATA_FILE)
    print(f"Successfully loaded data from {DATA_FILE}")

    # Clean up column names by stripping whitespace
    df.columns = df.columns.str.strip()

    # Validate essential columns for the heatmap axes
    required_columns_for_heatmap_axes = ["Gene Symbol", "Super Structure Name"]
    if not all(col in df.columns for col in required_columns_for_heatmap_axes):
        missing_cols = [col for col in required_columns_for_heatmap_axes if col not in df.columns]
        print(f"Error: Missing required columns for heatmap axes in the CSV file: {', '.join(missing_cols)}")
        print(f"Available columns: {df.columns.tolist()}")
        exit()

except FileNotFoundError:
    print(f"Error: The file '{DATA_FILE}' was not found.")
    print("Please ensure the 'data' directory and 'filtered_expression.csv' exist.")
    exit()
except Exception as e:
    print(f"An error occurred while loading or processing the CSV: {e}")
    exit()

if df is None:
    print("DataFrame could not be loaded or initialized. Exiting.")
    exit()

# --- Comprehensive Data Inspection (Before dropping NaNs) ---
print("\n--- Comprehensive Data Inspection (Before dropping NaNs) ---")
print(f"Total rows in original DataFrame: {len(df)}")
print("\nDataFrame Info (Non-Null Counts and Dtypes):")
df.info()

print("\n'Gene Symbol' Column Value Counts (including NaNs):")
print(df['Gene Symbol'].value_counts(dropna=False).head(10))
print(f"Number of NaN values in 'Gene Symbol': {df['Gene Symbol'].isnull().sum()}")

print("\n'Super Structure Name' Column Value Counts (including NaNs):")
print(df['Super Structure Name'].value_counts(dropna=False).head(10))
print(f"Number of NaN values in 'Super Structure Name': {df['Super Structure Name'].isnull().sum()}")

print("\n'Sub Structure Name' Column Value Counts (including NaNs):")
print(df['Sub Structure Name'].value_counts(dropna=False).head(10))
print(f"Number of NaN values in 'Sub Structure Name': {df['Sub Structure Name'].isnull().sum()}")

print("\nFirst 10 rows of the original DataFrame:")
print(df.head(10))
print("---------------------------------------\n")


# --- Data Preparation for Clustermap ---
# Drop rows only where 'Gene Symbol' or 'Super Structure Name' are missing.
df_for_heatmap = df.dropna(subset=required_columns_for_heatmap_axes)

# --- Debugging Prints (After dropping NaNs) ---
print("\n--- Debugging Information (After dropping NaNs for Clustermap) ---")
print(f"Rows remaining after dropping NaNs in 'Gene Symbol'/'Super Structure Name': {len(df_for_heatmap)}")
print(f"Unique genes for clustermap: {df_for_heatmap['Gene Symbol'].nunique()}")
print(f"Unique super structures for clustermap: {df_for_heatmap['Super Structure Name'].nunique()}")
print("-------------------------------------------\n")


# Create a pivot table for the clustermap data:
# Index: Gene Symbol
# Columns: Super Structure Name
# Values: Count of observations
heatmap_data = df_for_heatmap.pivot_table(
    index="Gene Symbol",
    columns="Super Structure Name",
    aggfunc="size",
    fill_value=0
)

# Check if heatmap_data is empty
if heatmap_data.empty:
    print("No data found for clustermap after removing rows with missing Gene Symbol or Super Structure Name.")
    print("Please inspect and clean your 'filtered_expression.csv' file.")
    exit()

# --- Clustermap Plotting ---
# Clustermap automatically handles figure sizing and layout, but we can
# pass `figsize` to control the overall size of the plot including dendrograms.
# Adjust figsize based on the number of genes and super structures for better readability.
# A heuristic: 0.5 inches per gene row, 0.3 inches per super structure column.
fig_height_base = heatmap_data.shape[0] * 0.5
fig_width_base = heatmap_data.shape[1] * 0.3

# Ensure a reasonable minimum size
fig_height = max(10, fig_height_base)
fig_width = max(12, fig_width_base)


# Create the clustermap
# `row_cluster=True` and `col_cluster=True` (default) perform clustering.
# `cmap` sets the color scheme.
# `annot=True` displays the values on the heatmap cells.
# `fmt="d"` formats annotations as integers.
# `cbar_kws` customizes the color bar label.
# `figsize` controls the overall figure size.
g = sns.clustermap(
    heatmap_data,
    cmap="YlGnBu",
    linewidths=0.4,
    linecolor='gray',
    annot=True,
    fmt="d",
    cbar_kws={"label": "Number of Expression Annotations"},
    figsize=(fig_width, fig_height),
    # Optional: You can disable clustering for rows or columns if desired
    # row_cluster=False,
    # col_cluster=False,
    # Optional: Standardize data before clustering/plotting (e.g., by row or column)
    # z_score=0, # Z-score normalization across columns (organs)
    # z_score=1, # Z-score normalization across rows (genes)
    # standard_scale=0, # Scale each column to range [0, 1]
    # standard_scale=1, # Scale each row to range [0, 1]
)

# Adjust title and labels. Clustermap handles x/y tick labels automatically.
# We set the main title of the entire figure.
g.ax_row_dendrogram.set_title("Overall Gene Expression Clustermap (All Zebrafish Stages)", fontsize=16)
g.ax_col_dendrogram.set_title("Clustering by Super Structure", fontsize=12)
g.ax_heatmap.set_ylabel("Gene Symbol", fontsize=12)
g.ax_heatmap.set_xlabel("Super Structure Name", fontsize=12)

# Rotate x-axis labels for better readability
# Access the heatmap axes and rotate the tick labels
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=10)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)


# Save the plot
output_path = os.path.join(OUTPUT_DIR, OUTPUT_FILENAME)
plt.savefig(output_path, dpi=300, bbox_inches='tight') # Use bbox_inches='tight' for better saving with dendrograms
plt.close() # Close the plot to free up memory

print(f"✅ Clustermap saved successfully to {output_path}")