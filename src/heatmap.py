import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --- Configuration ---
DATA_FILE = "data/filtered_expression.csv"
OUTPUT_DIR = "plots/summary"
OUTPUT_FILENAME = "all_genes_by_super_sub_organ_overall_heatmap.png" # Filename from your last output

# Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Data Loading and Initial Preprocessing ---
df = None
try:
    df = pd.read_csv(DATA_FILE)
    print(f"Successfully loaded data from {DATA_FILE}")

    # Clean up column names by stripping whitespace
    df.columns = df.columns.str.strip()

    # Validate essential columns for the heatmap.
    # We will now use 'Gene Symbol' and 'Super Structure Name' for the heatmap axes.
    # 'Sub Structure Name' is excluded from this check for dropping NaNs,
    # because it's mostly empty in your data.
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
print(df['Gene Symbol'].value_counts(dropna=False).head(10)) # Show top 10 for brevity
print(f"Number of NaN values in 'Gene Symbol': {df['Gene Symbol'].isnull().sum()}")

print("\n'Super Structure Name' Column Value Counts (including NaNs):")
print(df['Super Structure Name'].value_counts(dropna=False).head(10)) # Show top 10
print(f"Number of NaN values in 'Super Structure Name': {df['Super Structure Name'].isnull().sum()}")

print("\n'Sub Structure Name' Column Value Counts (including NaNs):")
print(df['Sub Structure Name'].value_counts(dropna=False).head(10)) # Show top 10
print(f"Number of NaN values in 'Sub Structure Name': {df['Sub Structure Name'].isnull().sum()}")

print("\nFirst 10 rows of the original DataFrame:")
print(df.head(10))
print("---------------------------------------\n")


# --- Data Preparation for Heatmap ---
# Drop rows only where 'Gene Symbol' or 'Super Structure Name' are missing.
# We are intentionally NOT dropping based on 'Sub Structure Name' here,
# as it has too many missing values to be a primary axis.
df_for_heatmap = df.dropna(subset=required_columns_for_heatmap_axes)

# --- Debugging Prints (After dropping NaNs) ---
print("\n--- Debugging Information (After dropping NaNs for Heatmap) ---")
print(f"Rows remaining after dropping NaNs in 'Gene Symbol'/'Super Structure Name': {len(df_for_heatmap)}")
print(f"Unique genes for heatmap: {df_for_heatmap['Gene Symbol'].nunique()}")
print(f"Unique super structures for heatmap: {df_for_heatmap['Super Structure Name'].nunique()}")
print("-------------------------------------------\n")


# Create a pivot table for the heatmap data:
# Index: Gene Symbol
# Columns: Super Structure Name (only, as Sub Structure Name is too sparse)
# Values: Count of observations
heatmap_data = df_for_heatmap.pivot_table(
    index="Gene Symbol",
    columns="Super Structure Name", # Now only using Super Structure Name
    aggfunc="size",
    fill_value=0
)

# Check if heatmap_data is empty
if heatmap_data.empty:
    print("No data found for heatmap after removing rows with missing Gene Symbol or Super Structure Name.")
    print("Please inspect and clean your 'filtered_expression.csv' file.")
    exit()

# --- Heatmap Plotting ---
# Adjust figure size dynamically based on the number of columns (super structures)
width_per_column = 0.8
fig_width = max(15, heatmap_data.shape[1] * width_per_column) # Min width 15
fig_height = max(10, heatmap_data.shape[0] * 0.4) # Adjust height based on number of genes

plt.figure(figsize=(fig_width, fig_height))

sns.heatmap(
    heatmap_data,
    cmap="YlGnBu",
    linewidths=0.4,
    linecolor='gray',
    annot=True,
    fmt="d",
    cbar_kws={"label": "Number of Expression Annotations"}
)

plt.title("Overall Gene Expression Across Super Structures (All Zebrafish Stages)", fontsize=18)
plt.xlabel("Super Structure Name", fontsize=12)
plt.ylabel("Gene Symbol", fontsize=12)

plt.xticks(rotation=45, ha='right', fontsize=10) # Adjust rotation and fontsize for single level
plt.yticks(fontsize=8)

plt.tight_layout()

output_path = os.path.join(OUTPUT_DIR, OUTPUT_FILENAME)
plt.savefig(output_path, dpi=300)
plt.close()

print(f"✅ Heatmap saved successfully to {output_path}")
