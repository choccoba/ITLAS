########################################################################
# V18-C5: INDIVIDUAL GENE MARKERS — Tissue-Separated, Liver & Blood
# ============================================================
# Step (3) of 7-step pipeline
# Prereq: C3 (proportions), C4 (pathways) completed
# Input:  c5_genes list from C4, adata_liver, adata_blood
# Output: version18-analysis/C5_genes/
########################################################################

# %% [markdown]
# # V18-C5: Individual Gene Markers (Tissue-Separated)
# **Step 3/7**: Analyze individual genes from significant pathways, Liver AND Blood simultaneously.

# %%
# ============================================================
#  CELL 21: SETUP C5
# ============================================================
C5_SAVE = '/content/drive/MyDrive/ITLAS/results/version18-analysis/C5_genes/'
os.makedirs(C5_SAVE, exist_ok=True)

print("=" * 70)
print("  C5: INDIVIDUAL GENE MARKERS — Tissue-Separated")
print("=" * 70)

# Load selected genes from C4 (or use c5_genes from memory)
try:
    c5_gene_df = pd.read_csv(f'{C4_SAVE}C4_selected_genes_for_C5.csv')
    c5_genes = sorted(c5_gene_df['gene'].tolist())
except:
    print("  ⚠️ Using c5_genes from memory")

print(f"  Target genes: {len(c5_genes)}")

# Gene → pathway mapping
gene_to_pathway = {}
for gs_name, gs_genes in gene_sets.items():
    for g in gs_genes:
        if g in gene_to_pathway:
            gene_to_pathway[g] += f', {gs_name}'
        else:
            gene_to_pathway[g] = gs_name

# %%
# ============================================================
#  CELL 22: GENE SCREENING ENGINE (Tissue-Separated)
# ============================================================
print("=" * 70)
print("  GENE SCREENING ENGINE")
print("=" * 70)

def screen_gene_tissue(ad, tissue_name, gene, lin_data_name):
    """Screen one gene × one lineage across all comparisons in one tissue."""
    results = []
    
    lin_mask = ad.obs[lineage_col] == lin_data_name
    sub = ad[lin_mask]
    
    if gene not in sub.var_names:
        return results
    
    # Pre-extract expression
    gene_expr = sub[:, gene].X
    if hasattr(gene_expr, 'toarray'):
        gene_expr = gene_expr.toarray().flatten()
    else:
        gene_expr = np.array(gene_expr).flatten()
    
    for s1, s2, comp_name in COMPARISONS:
        def get_donor_means(stage):
            stage_mask = sub.obs[stage_col].values == stage
            donors = sub.obs.loc[stage_mask, donor_col].unique()
            means = []
            for d in donors:
                d_mask = (sub.obs[stage_col].values == stage) & (sub.obs[donor_col].values == d)
                if d_mask.sum() >= 5:
                    means.append(float(np.mean(gene_expr[d_mask])))
            return means
        
        vals1 = get_donor_means(s1)
        vals2 = get_donor_means(s2)
        
        if len(vals1) < 2 or len(vals2) < 2:
            continue
        
        mean1, mean2 = np.mean(vals1), np.mean(vals2)
        pct_change = ((mean2 - mean1) / abs(mean1) * 100) if abs(mean1) > 1e-10 else (99999 if mean2 > 0 else 0)
        
        try:
            stat, pval = mannwhitneyu(vals1, vals2, alternative='two-sided')
        except:
            pval = 1.0
        
        n_pairs = len(vals1) * len(vals2)
        if mean2 >= mean1:
            consist = sum(1 for v2 in vals2 for v1 in vals1 if v2 >= v1)
        else:
            consist = sum(1 for v2 in vals2 for v1 in vals1 if v2 < v1)
        
        direction = '↑' if mean2 > mean1 else ('↓' if mean2 < mean1 else '=')
        sig = '***' if pval<0.001 else ('**' if pval<0.01 else ('*' if pval<0.05 else ('†' if pval<0.10 else 'NS')))
        
        results.append({
            'tissue': tissue_name,
            'gene': gene,
            'lineage': lin_std,
            'pathway': gene_to_pathway.get(gene, 'other'),
            'comparison': comp_name,
            'stage1': s1, 'stage2': s2,
            'mean_s1': round(mean1, 6),
            'mean_s2': round(mean2, 6),
            'pct_change': round(pct_change, 1),
            'direction': direction,
            'p_value': round(pval, 6),
            'sig': sig,
            'consistency': f"{consist}/{n_pairs}",
            'consist_pct': round(consist/n_pairs*100, 1) if n_pairs > 0 else 0,
            'n_s1': len(vals1), 'n_s2': len(vals2),
        })
    
    return results

print(f"  Engine ready: {len(c5_genes)} genes × {len(lineage_map)} lineages × {len(COMPARISONS)} comparisons × 2 tissues")
print(f"  Expected tests: ~{len(c5_genes)*len(lineage_map)*len(COMPARISONS)*2:,}")

# %%
# ============================================================
#  CELL 23: RUN LIVER SCREENING
# ============================================================
print("=" * 70)
print("  🔴 LIVER GENE SCREENING")
print("=" * 70)

results_liver = []
t0 = time.time()
total = len(c5_genes)

for gi, gene in enumerate(c5_genes):
    if (gi+1) % 20 == 0 or gi == 0:
        elapsed = time.time() - t0
        rate = (gi+1)/elapsed if elapsed > 0 else 1
        eta = (total-gi-1)/rate/60 if rate > 0 else 0
        print(f"  [{gi+1:3d}/{total}] {gene:12s} | {elapsed:.0f}s | ETA: {eta:.1f}min")
    
    for lin_std, lin_data in lineage_map.items():
        results_liver.extend(screen_gene_tissue(adata_liver, 'Liver', gene, lin_data))

df_gene_liver = pd.DataFrame(results_liver)
elapsed = time.time() - t0
print(f"\n✅ Liver: {len(df_gene_liver):,} tests, {(df_gene_liver['p_value']<0.05).sum()} sig ({elapsed:.0f}s)")
df_gene_liver.to_csv(f'{C5_SAVE}C5_genes_liver.csv', index=False)

# %%
# ============================================================
#  CELL 24: RUN BLOOD SCREENING
# ============================================================
print("=" * 70)
print("  🔵 BLOOD GENE SCREENING")
print("=" * 70)

results_blood = []
t0 = time.time()

for gi, gene in enumerate(c5_genes):
    if (gi+1) % 20 == 0 or gi == 0:
        elapsed = time.time() - t0
        rate = (gi+1)/elapsed if elapsed > 0 else 1
        eta = (total-gi-1)/rate/60 if rate > 0 else 0
        print(f"  [{gi+1:3d}/{total}] {gene:12s} | {elapsed:.0f}s | ETA: {eta:.1f}min")
    
    for lin_std, lin_data in lineage_map.items():
        results_blood.extend(screen_gene_tissue(adata_blood, 'Blood', gene, lin_data))

df_gene_blood = pd.DataFrame(results_blood)
elapsed = time.time() - t0
print(f"\n✅ Blood: {len(df_gene_blood):,} tests, {(df_gene_blood['p_value']<0.05).sum()} sig ({elapsed:.0f}s)")
df_gene_blood.to_csv(f'{C5_SAVE}C5_genes_blood.csv', index=False)

# %%
# ============================================================
#  CELL 25: TISSUE DISCREPANCY ANALYSIS — NL→IT
# ============================================================
print("=" * 70)
print("  🔴🔵 TISSUE DISCREPANCY: NL→IT Gene Changes")
print("=" * 70)

# Merge Liver and Blood for NL→IT
l_it = df_gene_liver[df_gene_liver['comparison']=='NL→IT'][['gene','lineage','pct_change','p_value','sig','direction']].copy()
l_it.columns = ['gene','lineage','L_pct','L_p','L_sig','L_dir']

b_it = df_gene_blood[df_gene_blood['comparison']=='NL→IT'][['gene','lineage','pct_change','p_value','sig','direction']].copy()
b_it.columns = ['gene','lineage','B_pct','B_p','B_sig','B_dir']

disc = l_it.merge(b_it, on=['gene','lineage'], how='outer')
disc['dir_same'] = disc['L_dir'] == disc['B_dir']
disc['both_sig'] = (disc['L_p'] < 0.05) & (disc['B_p'] < 0.05)
disc['liver_only'] = (disc['L_p'] < 0.05) & (disc['B_p'] >= 0.05)
disc['blood_only'] = (disc['L_p'] >= 0.05) & (disc['B_p'] < 0.05)

disc.to_csv(f'{C5_SAVE}C5_tissue_discrepancy_NL_IT.csv', index=False)

# Summary
print(f"\n  Total gene-lineage pairs: {len(disc)}")
print(f"  Direction concordant: {disc['dir_same'].sum()} ({disc['dir_same'].mean()*100:.1f}%)")
print(f"  Both significant: {disc['both_sig'].sum()}")
print(f"  Liver-only sig: {disc['liver_only'].sum()}")
print(f"  Blood-only sig: {disc['blood_only'].sum()}")

# Direction reversals with significance
rev = disc[~disc['dir_same'] & ((disc['L_p']<0.05)|(disc['B_p']<0.05))]
if len(rev) > 0:
    print(f"\n  🚨 DIRECTION REVERSALS ({len(rev)}):")
    for _, row in rev.sort_values('L_p').head(20).iterrows():
        pw = gene_to_pathway.get(row['gene'], '')
        print(f"    {row['gene']:12s} {row['lineage']:8s} | 🔴{row['L_dir']}{abs(row['L_pct']):.0f}%({row['L_sig']}) vs 🔵{row['B_dir']}{abs(row['B_pct']):.0f}%({row['B_sig']}) [{pw}]")

# %%
# ============================================================
#  CELL 26: NL→IT SIGNIFICANT GENES — Liver vs Blood SIMULTANEOUS
# ============================================================
print("=" * 70)
print("  🔴🔵 NL → IT: SIGNIFICANT GENES (Simultaneous)")
print("=" * 70)

# Prepare discrepancy dataframe (disc)
l_dat = df_gene_liver[df_gene_liver['comparison'] == 'NL→IT'].copy()
b_dat = df_gene_blood[df_gene_blood['comparison'] == 'NL→IT'].copy()

disc = pd.merge(l_dat, b_dat, on=['lineage', 'gene'], suffixes=('_L', '_B'))
disc = disc.rename(columns={
    'p_value_L': 'L_p', 'p_value_B': 'B_p',
    'direction_L': 'L_dir', 'direction_B': 'B_dir',
    'pct_change_L': 'L_pct', 'pct_change_B': 'B_pct'
})
disc['dir_same'] = (disc['L_dir'] == disc['B_dir'])

for lin in LINEAGES:
    print(f"\n  ── {lin} ──")
    lin_disc = disc[(disc['lineage']==lin) & ((disc['L_p']<0.05)|(disc['B_p']<0.05))]
    lin_disc = lin_disc.sort_values('L_p')
    
    if len(lin_disc) == 0:
        print(f"    (no significant genes)")
        continue
    
    for _, row in lin_disc.iterrows():
        l_m = '★' if row['L_p']<0.05 else ('†' if row['L_p']<0.10 else ' ')
        b_m = '★' if row['B_p']<0.05 else ('†' if row['B_p']<0.10 else ' ')
        match = '✅' if row['dir_same'] else '⚠️↔'
        pw = gene_to_pathway.get(row['gene'], '')[:20]
        
        print(f"    {row['gene']:12s} | 🔴{l_m}{row['L_dir']}{abs(row['L_pct']):7.1f}% p={row['L_p']:.4f} | 🔵{b_m}{row['B_dir']}{abs(row['B_pct']):7.1f}% p={row['B_p']:.4f} | {match} [{pw}]")


# %%
# ============================================================
#  CELL 27: SELECT SIGNIFICANT GENES → Pass to C6/C7
# ============================================================
print("=" * 70)
print("  🎯 SIGNIFICANT GENE SELECTION FOR C6/C7")
print("=" * 70)

# All comparisons, collect significant genes per lineage per tissue
sig_gene_summary = []

for comp_name in [c[2] for c in COMPARISONS]:
    for lin in LINEAGES_STD:
        l_sig = df_gene_liver[(df_gene_liver['comparison']==comp_name) & (df_gene_liver['lineage']==lin) & (df_gene_liver['p_value']<0.05)]
        b_sig = df_gene_blood[(df_gene_blood['comparison']==comp_name) & (df_gene_blood['lineage']==lin) & (df_gene_blood['p_value']<0.05)]
        
        for _, row in l_sig.iterrows():
            sig_gene_summary.append({**row.to_dict(), 'sig_tissue': 'Liver'})
        for _, row in b_sig.iterrows():
            sig_gene_summary.append({**row.to_dict(), 'sig_tissue': 'Blood'})

df_sig_all = pd.DataFrame(sig_gene_summary)
df_sig_all.to_csv(f'{C5_SAVE}C5_all_significant_genes.csv', index=False)

# Top genes by frequency of significance across lineages/comparisons
gene_sig_counts = df_sig_all.groupby('gene').agg(
    n_sig=('comparison', 'count'),
    n_lineages=('lineage', 'nunique'),
    n_comparisons=('comparison', 'nunique'),
    tissues=('sig_tissue', lambda x: ','.join(sorted(set(x)))),
).sort_values('n_sig', ascending=False)

print(f"\n  Top 20 Most Frequently Significant Genes:")
print(f"  {'Gene':12s} {'#Sig':>5s} {'#Lin':>5s} {'#Comp':>6s} {'Tissues'}")
for gene, row in gene_sig_counts.head(20).iterrows():
    pw = gene_to_pathway.get(gene, '')[:25]
    print(f"  {gene:12s} {row['n_sig']:5d} {row['n_lineages']:5d} {row['n_comparisons']:6d} {row['tissues']:15s} [{pw}]")

gene_sig_counts.to_csv(f'{C5_SAVE}C5_gene_significance_ranking.csv')

# %%
# ============================================================
#  CELL 28: C5 COMPLETE
# ============================================================
print("=" * 70)
print("  ✅ C5 COMPLETE — Individual Gene Markers")
print("=" * 70)

print(f"""
  Files saved to: {C5_SAVE}
  ├── C5_genes_liver.csv            (all gene results, Liver)
  ├── C5_genes_blood.csv            (all gene results, Blood)
  ├── C5_tissue_discrepancy_NL_IT.csv
  ├── C5_all_significant_genes.csv
  └── C5_gene_significance_ranking.csv

  ▶ NEXT: C6 — Tissue Discrepancy Master Compilation
  ▶ THEN: C7 — Inter-relationship Pattern Discovery
""")
