# phyton-script-for-microarray-ID-s
<<<<<<< HEAD

<br>
gse_id = "GSE42977"   # <-- Change to any GEO Series ID
gse = GEOparse.get_GEO(geo=gse_id, destdir=".")

print(f"Dataset: {gse_id}")
print(f"Title: {gse.metadata['title'][0]}")
#print(f"Organism: {gse.metadata['organism_ch1'][0]}")
print(f"Samples: {len(gse.gsms)}")

# Extract expression data# Step 2: Build expression matrix

dfs = []
for gsm_name, gsm in gse.gsms.items():
if gsm.table is not None and "VALUE" in gsm.table.columns:
df = gsm.table[["ID_REF", "VALUE"]].copy()
df.rename(columns={"VALUE": gsm_name}, inplace=True)
dfs.append(df)

expr_df = dfs[0]
for df in dfs[1:]:
expr_df = expr_df.merge(df, on="ID_REF")

expr_df.set_index("ID_REF", inplace=True)

print("Expression Matrix Shape:", expr_df.shape)

# Save raw expression matrix

expr_df.to_csv("expression.csv")
print("✅ Saved: expression.csv")

# Step 3: Differential Expression Analysis

n = len(expr_df.columns)
control = expr_df.iloc[:, :n//2]
treated = expr_df.iloc[:, n//2:]

pvals, logFCs, t_stats, ave_expr_vals, = [], [], [], [],
for gene_id in expr_df.index: # Renamed gene to gene_id for clarity
stat, p = ttest_ind(control.loc[gene_id], treated.loc[gene_id], nan_policy="omit")
logFC = treated.loc[gene_id].mean() - control.loc[gene_id].mean()
pvals.append(p)
logFCs.append(logFC)
t_stats.append(stat)
ave_expr_vals.append(expr_df.loc[gene_id].mean()) # Correct calculation of average expression

deg_df = pd.DataFrame({
"Gene": expr_df.index,
"logFC": logFCs,
"pval": pvals,
"significance": ["significant" if p_val < 0.05 else "non-significant" for p_val in pvals],
"description": ["up-regulated" if fc > 0 else "down-regulated" for fc in logFCs],
"t": t_stats, # Use the collected t-statistics
"AveExpr" : ave_expr_vals # Use the collected average expression values
})
deg_df["-log10p"] = -np.log10(deg_df["pval"])

# Save DEG results

deg_df.to_csv("DEG.csv", index=False)
print("✅ Saved: DEG.csv")
platform = list(gse.gpls.values())[0]
annot = platform.table

print(annot.columns)
annot = platform.table[['ID', 'Symbol']].copy() # Get both ID and Symbol
annot.rename(columns={'ID': 'Probe', 'Symbol': 'Gene_Symbol'}, inplace=True)

# The 'ID' column for this dataset contains string identifiers (e.g., 'AATK_E63_R'),

# so it should not be converted to numeric.

# annot['Probe'] = pd.to_numeric(annot['Probe'], errors='coerce')

# annot = annot.dropna(subset=['Probe'])

# annot['Probe'] = annot['Probe'].astype(int)

print(annot.head())

# The 'Gene' column in deg_df is currently the 'ID_REF' which can be mixed type (string/numeric).

# Rename the original 'Gene' column (which contains probe IDs) to 'Probe_ID'.

if 'Gene' in deg_df.columns and 'Probe_ID' not in deg_df.columns:
deg_df.rename(columns={'Gene': 'Probe_ID'}, inplace=True)

# Ensure 'Probe_ID' is string type for merging

deg_df['Probe_ID'] = deg_df['Probe_ID'].astype(str)
print("✅ Ensured 'Probe_ID' column in deg_df is string type for merging.")
deg_df = deg_df.merge(annot, left_on='Probe_ID', right_on='Probe', how='left')

# Rename 'Gene_Symbol' (which contains the gene names from annotation) to 'Gene'.

# This makes the 'Gene' column represent the gene symbol.

if 'Gene_Symbol' in deg_df.columns:
deg_df.rename(columns={'Gene_Symbol': 'Gene'}, inplace=True)

# Drop any redundant 'Probe' column that was introduced from the 'annot' DataFrame during the merge.

deg_df.drop(columns=['Probe'], errors='ignore', inplace=True)

# Drop any temporary columns if they were created in a previous attempt

deg_df.drop(columns=['Gene_Symbol_temp'], errors='ignore', inplace=True)

# Now, deg_df should contain 'Gene' (gene symbols) and 'Probe_ID' (probe identifiers),

# along with all other differential expression results.

deg_df = deg_df.sort_values('pval')
deg_df = deg_df.drop_duplicates(subset='Probe_ID')
deg_final = deg_df[['Gene', 'Probe_ID', 'logFC', 'pval', '-log10p','significance', 'description', "t","AveExpr"]]
deg_final.to_csv("DEG_gene.csv", index=False)

# Step 4: Filter significant DEGs

# Thresholds: p < 0.05 and |logFC| > 1

sig_df = deg_df[(deg_df["pval"] < 0.05) & (abs(deg_df["logFC"]) >1)]

sig_df.to_csv("significant.csv", index=False)
print(f"✅ Saved: significant.csv ({sig_df.shape[0]} significant genes)")

# Step 5: Volcano plot

plt.figure(figsize=(8,6))
plt.scatter(deg_df["logFC"], deg_df["-log10p"], alpha=0.5, label="All genes")
plt.scatter(sig_df["logFC"], sig_df["-log10p"], color="red", alpha=0.7, label="Significant")
plt.axhline(-np.log10(0.05), color="blue", linestyle="--")
plt.axvline(-1, color="green", linestyle="--")
plt.axvline(1, color="green", linestyle="--")
plt.xlabel("Log Fold Change")
plt.ylabel("-log10(p-value)")
plt.title(f"Volcano Plot - {gse_id}")
plt.legend()
plt.savefig("volcano_plot.png", dpi=300)
plt.show()
=======
>>>>>>> cef2397ec78edd9b3862a9a476694c34c60bc0e0
