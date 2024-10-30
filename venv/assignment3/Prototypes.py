import urllib
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Import necessary libraries
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
#print(adata.obs)
#adata.obs.to_csv('obs_data.csv',index=False)


# Step 3: Experiment and find the differences between the two cohorts in disease category (purple Alzheimer subject, green normal subject).

# Split into normal and AD
obs = adata.obs
normal = obs.loc[obs.disease == "normal"]
ad = obs.loc[obs.disease == "Alzheimer disease"]
#print(normal)
#print(ad)

#subjects of gender cohorts
plt.figure(figsize=(10, 6))
for i, gene in enumerate(genes_of_interest, 1):
    plt.subplot(2, 3, i)
    sns.boxplot(
        x=obs['sex'],
        y=adata[:, gene].X.toarray().flatten(),
        palette={"male": "skyblue", "female": "salmon"}
    )
    plt.title(f"Expression of {lookup[gene]}")
    plt.xlabel('Gender')
    plt.ylabel('Expression')

plt.tight_layout()
plt.show()


# Set up the figure size and layout
fig, axs = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('Gene Expression Levels by Gender (Male vs. Female)', fontsize=16)

# Plot each gene as a violin plot
for i, gene in enumerate(genes_of_interest):
    row, col = divmod(i, 3)
    sns.violinplot(x=obs['sex'], y=adata[:, gene].X.toarray().flatten(), ax=axs[row, col],
                   palette={"female": "lightcoral", "male": "skyblue"})
    axs[row, col].set_title(f'Expression of {gene}')
    axs[row, col].set_xlabel('Gender')
    axs[row, col].set_ylabel('Expression')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout to make room for the title
plt.show()