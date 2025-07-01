import pandas as pd

df = pd.read_csv("data/filtered_expression.csv")
print(df.columns)  # Show available columns

# We'll pick these essential fields (adjust names as needed):
# 'Gene Symbol', 'Start Stage', 'End Stage', 'Assay', 'Super Structure Name'

df['Stage'] = (df['Start Stage'] + df['End Stage']) / 2  # Midpoint stage
df = df[['Gene Symbol', 'Stage', 'Super Structure Name', 'Assay']].dropna()
df = df.rename(columns={
    'Gene Symbol': 'Gene',
    'Super Structure Name': 'Anatomy'
})
df.to_csv("data/clean_expression.csv", index=False)
print("✅ Clean data saved to data/clean_expression.csv")
