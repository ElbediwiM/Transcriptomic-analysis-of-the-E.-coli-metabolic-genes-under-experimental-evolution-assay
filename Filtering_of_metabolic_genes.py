#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd

# Step 1: Load the CSV file
csv_file_path = 'e:/metabolic_tr/transcriptomic/A.baumannii/file.csv'  


# Step 2: Identify duplicates in the 'genes' column 
genes_df = pd.read_csv('e:/metabolic_tr/transcriptomic/ E.coli /metabolic_genes.csv')
genes_list = genes_df['gene'].tolist() 

# Step 3: Filter the DataFrame

filtered_df = df[df['gene'].isin(genes_list)]  

# Step 4: Save the filtered DataFrame to a new CSV
output_file_path = 'e:/metabolic_tr/transcriptomic/ E.coli /Filtered data.csv'  
filtered_df.to_csv(output_file_path, index=False)

print(f"Filtered data saved to {output_file_path}")

