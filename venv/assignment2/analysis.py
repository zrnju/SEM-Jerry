# Checking available fields in adata.obs to identify disease category information
print(adata.obs.columns)

# Plotting cohort distribution if the field is named 'disease_category'
plt.figure(figsize=(6, 4))
sns.countplot(data=adata.obs, x='disease_category', palette={'Alzheimer': 'purple', 'Normal': 'green'})
plt.title("Cell Count per Disease Category")
plt.show()