import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np # Still useful for potential future needs, but jitter removed for now

# --- Configuration ---
DATA_FILE = "data/filtered_expression.csv"
OUTPUT_DIR = "plots/summary"
OUTPUT_FILENAME = "genes_vs_super_structure_bubble_heatmap_grid_aligned.png" # Updated filename

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
        print(f"Error: Missing required columns for bubble heatmap axes in the CSV file: {', '.join(missing_cols)}")
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


# --- Data Preparation for Bubble Heatmap ---
# Drop rows only where 'Gene Symbol' or 'Super Structure Name' are missing.
df_for_plotting = df.dropna(subset=required_columns_for_heatmap_axes)

# --- Debugging Prints (After dropping NaNs) ---
print("\n--- Debugging Information (After dropping NaNs for Bubble Heatmap) ---")
print(f"Rows remaining after dropping NaNs in 'Gene Symbol'/'Super Structure Name': {len(df_for_plotting)}")
print(f"Unique genes for bubble heatmap: {df_for_plotting['Gene Symbol'].nunique()}")
print(f"Unique super structures for bubble heatmap: {df_for_plotting['Super Structure Name'].nunique()}")
print("-------------------------------------------\n")


# Create a pivot table to get counts for each Gene Symbol-Super Structure Name pair
counts_pivot = df_for_plotting.pivot_table(
    index="Gene Symbol",
    columns="Super Structure Name",
    aggfunc="size",
    fill_value=0
)

# Check if counts_pivot is empty
if counts_pivot.empty:
    print("No data found for bubble heatmap after removing rows with missing Gene Symbol or Super Structure Name.")
    print("Please inspect and clean your 'filtered_expression.csv' file.")
    exit()

# Melt the pivot table back into a long format suitable for scatterplot
bubble_data = counts_pivot.reset_index().melt(
    id_vars="Gene Symbol",
    var_name="Super Structure Name",
    value_name="Expression Count"
)

# Filter out rows where Expression Count is 0, as we don't need to plot empty bubbles
# This is crucial, as 0-count bubbles would still occupy space on the grid.
bubble_data = bubble_data[bubble_data["Expression Count"] > 0]

# Check if bubble_data is empty after filtering for counts > 0
if bubble_data.empty:
    print("No actual expression data (counts > 0) found to plot for the bubble heatmap.")
    print("This might mean all gene-super structure pairs have 0 expression counts.")
    exit()

# --- Bubble Heatmap Plotting ---
# Adjust figure size dynamically based on the number of unique genes and super structures.
num_genes = bubble_data['Gene Symbol'].nunique()
num_super_structures = bubble_data['Super Structure Name'].nunique()

# Base size per item - adjusted for more space and to prevent overlap
# We need more space per column for distinct bubbles without jitter.
base_width_per_col = 1.0 # Adjusted for direct alignment, may need tuning
base_height_per_row = 0.6

fig_width = max(20, num_super_structures * base_width_per_col)
fig_height = max(15, num_genes * base_height_per_row)

plt.figure(figsize=(fig_width, fig_height))

# Create the scatter plot (bubble heatmap)
# Crucially, we are now using the original categorical columns for x and y.
# This makes the bubbles align to a grid.
sns.scatterplot(
    data=bubble_data,
    x="Super Structure Name", # Original categorical X-axis for grid alignment
    y="Gene Symbol",          # Original categorical Y-axis for grid alignment
    size="Expression Count",
    sizes=(150, 2000), # Adjust min and max bubble sizes carefully to prevent overlap
                               # Increased min size slightly to ensure visibility of smallest bubbles
    hue="Expression Count", # Color by expression count for visual emphasis
    palette="Reds", # Red palette: light to dark for increasing values
    legend="full", # Show the full legend for size and color
    alpha=0.8, # Increased transparency slightly to help with potential overlaps
    edgecolor="black", # Border around bubbles
    linewidth=0.5 # Width of the border
)

plt.title("Gene Expression Bubble Heatmap (Overall, Gene vs. Super Structure)", fontsize=18)
plt.xlabel("Super Structure Name", fontsize=12)
plt.ylabel("Gene Symbol", fontsize=12)

# Rotate x-axis labels for better readability
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.yticks(fontsize=8)

# Place the legend outside the plot area if it overlaps
plt.legend(title="Expression Count", bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)

plt.tight_layout()

# Save the plot
output_path = os.path.join(OUTPUT_DIR, OUTPUT_FILENAME)
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"✅ Bubble Heatmap saved successfully to {output_path}")