import requests
from bs4 import BeautifulSoup
import pandas as pd

def get_expression_data(gene_symbol):
    url = f"https://zfin.org/action/marker/{gene_symbol}/expression"
    response = requests.get(url)
    soup = BeautifulSoup(response.content, "html.parser")
    table = soup.find("table")
    
    data = []
    if table:
        rows = table.find_all("tr")[1:]
        for row in rows:
            cols = row.find_all("td")
            data.append([col.text.strip() for col in cols])
    
    df = pd.DataFrame(data, columns=["Stage", "Anatomical Structure", "Assay", "Pattern", "Reference"])
    df.to_csv(f"data/{gene_symbol}_expression.csv", index=False)
    return df

# Test with one gene
get_expression_data("c3a.1")
