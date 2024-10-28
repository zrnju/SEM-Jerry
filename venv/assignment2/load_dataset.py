import urllib
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Step 1: Load the dataset
urllib.request.urlretrieve('https://datasets.cellxgene.cziscience.com/5d871206-9489-4d9f-8106-94305ccb1c3a.h5ad', 'dataset.h5ad')
adata = sc.read_h5ad('dataset.h5ad')
print(adata)

# Step 2: Load genes data
urllib.request.urlretrieve("https://raw.githubusercontent.com/D3sert650/SEM_Assignment_1/refs/heads/main/data/differential_expression_AD_Normal.csv", "differential_expression_AD_Normal.csv")
genes = pd.read_csv("differential_expression_AD_Normal.csv", comment="#")
genes.head()

# Mapping genes to features
lookup = dict(zip(adata.var.index, adata.var['feature_name']))
reverse_lookup = {v: k for k, v in lookup.items()}

# Import genes of interest
genes_of_interest = [
    "SLC26A3", "RASGEF1B", "RP11-701H24.9", "LINGO1", "PDE4DIP",
    "AC159540.1", "RP11-289H16.1", "RP11-219A15.1", "LINC01609",
    "PHYHIP", "RP11-745L13.2"
]

# Map genes to their indices
genes_of_interest = [reverse_lookup[gene] for gene in genes_of_interest if gene in reverse_lookup]
print("Genes mapped to indices:", genes_of_interest)

# Export current dataset
print(adata.obs)
adata.obs.to_csv('obs_data.csv',index=False)


# Step 3: Experiment and find the differences between the two cohorts in disease category (purple Alzheimer subject, green normal subject).

# Split into normal and AD
obs = adata.obs
normal = obs.loc[obs.disease == "normal"]
ad = obs.loc[obs.disease == "Alzheimer disease"]
print(normal)
print(ad)

# Comparison
normal.reset_index(inplace=True)
ad.reset_index(inplace=True)
comparison = normal.describe().compare(ad.describe())
comparison_plot = comparison.T
comparison_plot.plot(kind='bar', figsize=(10, 6))
plt.title('Comparison of Descriptive Statistics Between Normal and AD')
plt.ylabel('Values')
plt.xlabel('Statistics')
plt.xticks(rotation=45)
plt.legend(title='Group')
plt.tight_layout()
plt.show()