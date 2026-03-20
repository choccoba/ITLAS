########################################################################
# V18-C4: 26 PATHWAY SCORING — Tissue-Separated, Liver & Blood 동시
# ============================================================
# Step (2) of 7-step pipeline
# Prereq: C3 completed (proportions)
# Input:  adata_liver, adata_blood, lineage_map (from C3)
# Output: version18-analysis/C4_pathway/
########################################################################

# %% [markdown]
# # V18-C4: 26 Pathway × 6 Lineage Scoring (Tissue-Separated)
# **Step 2/7**: Score 26 gene pathways per lineage, compare across disease spectrum, Liver AND Blood simultaneously.

# %%
# ============================================================
#  CELL 14: DEFINE 26 GENE SETS
# ============================================================
print("=" * 70)
print("  C4: 26 PATHWAY GENE SETS")
print("=" * 70)

C4_SAVE = '/content/drive/MyDrive/ITLAS/results/version18-analysis/C4_pathway/'
os.makedirs(C4_SAVE, exist_ok=True)

gene_sets = {
    'inflammasome':       ['NLRP3','AIM2','NLRC4','MEFV','PYCARD','CASP1'],
    'cytotoxicity':       ['GZMA','GZMB','GZMK','GNLY','PRF1','FASLG','NKG7'],
    'checkpoint':         ['PDCD1','LAG3','HAVCR2','TIGIT','CTLA4','CD160'],
    'exhaustion':         ['TOX','TOX2','ENTPD1','EOMES','BATF','LAYN'],
    'nk_function':        ['NCAM1','NCR1','NCR3','KLRK1','KLRD1','FCGR3A'],
    'nk_il15_dual':       ['GZMB','GNLY','PRF1','NKG7','KLRD1','FCGR3A','FCER1G','TYROBP'],
    'il15_mtor':          ['MTOR','RPTOR','RICTOR','EIF4EBP1','JAK1','JAK3','STAT5A','STAT5B','IL2RG','IL15RA'],
    'immune_evasion':     ['CD274','PDCD1LG2','LGALS9','VTCN1','HAVCR2','TGFB1','IDO1','TNFAIP3'],
    'treg':               ['FOXP3','IL2RA','CTLA4','TIGIT','IKZF2'],
    'naive_t':            ['TCF7','LEF1','SELL','CCR7','IL7R'],
    'memory_t':           ['GZMK','CD44','ID3','PRDM1','BCL6','BATF','CD27','CXCL13'],
    'tf_programs':        ['TBX21','GATA3','RORC','EOMES','BCL6','MAF','IRF4','BATF','PRDM1'],
    'tissue_resident':    ['ZNF683','ITGAE','CXCR6','PRDM1','CD69'],
    'stemness':           ['TCF7','LEF1','MYB','SELL','SLAMF6','BCL6','ID3'],
    'glycolysis':         ['HK1','HK2','PFKM','PKM','LDHA','ENO1','GAPDH','SLC2A1'],
    'oxphos':             ['NDUFA1','NDUFB1','SDHB','UQCRC1','COX4I1','ATP5F1A','ATP5F1B'],
    'mito_dysfunction':   ['MT-CO1','MT-ND1','MT-ND2','MT-CYB','MT-ATP6','PINK1','PRKN','TFAM'],
    'metabolic_recovery': ['TFAM','ESRRA','NRF1','SIRT3','PPARGC1A','CPT1A'],
    'apoptosis':          ['FAS','FASLG','BAK1','CASP3','CASP8','BCL2'],
    'senescence':         ['CDKN1A','CDKN2A','TP53','GLB1','IGFBP7','SERPINE1','RB1'],
    'cell_cycle':         ['CDK2','CDK4','CCND1','RB1','MKI67','TOP2A','PCNA','MCM2','CDK1','BIRC5'],
    'proliferation':      ['MKI67','TOP2A','PCNA','MCM2','CDK1'],
    'epigenetics':        ['DNMT1','DNMT3A','DNMT3B','TET2','EZH2','KDM6A','HDAC1','SIRT1'],
    'fibrosis':           ['TGFB1','TGFB2','TGFBR1','TGFBR2','COL1A1','COL3A1','ACTA2',
                           'FN1','TIMP1','MMP2'],
    'cancer_associated':  ['ATM','VIM','CDH1','SNAI1','ZEB1','MYC','CCND1','CDK4','CDKN2A'],
    'angiogenesis':       ['VEGFA','VEGFB','ANGPT2','PECAM1','FLT1','KDR','NRP1','HIF1A','EPAS1'],
}

# Check availability
gene_names = set(adata.var_names)
print(f"\n  Gene set availability:")
total_available = 0
total_defined = 0
for gs_name, gs_genes in gene_sets.items():
    avail = [g for g in gs_genes if g in gene_names]
    total_available += len(avail)
    total_defined += len(gs_genes)
    status = '✅' if len(avail) == len(gs_genes) else '⚠️'
    missing = [g for g in gs_genes if g not in gene_names]
    miss_str = f" missing: {missing}" if missing else ""
    print(f"  {status} {gs_name:20s} {len(avail):2d}/{len(gs_genes):2d}{miss_str}")

print(f"\n  Total: {total_available}/{total_defined} genes available")

# %%
# ============================================================
#  CELL 15: COMPUTE PATHWAY SCORES PER CELL (Mean Expression)
# ============================================================
print("=" * 70)
print("  COMPUTING PATHWAY SCORES (Mean Expression Proxy)")
print("  Method: mean(normalized expression) per gene set per cell")
print("  Validated: ρ = 0.922 concordance with AUCell (v17)")
print("=" * 70)

t0 = time.time()

# --- Split AnnData by Tissue (Fix for NameError) ---
print("  Splitting adata into Liver and Blood subsets...")
adata_liver = adata[adata.obs[tissue_col] == 'Liver'].copy()
adata_blood = adata[adata.obs[tissue_col] == 'Blood'].copy()
print(f"  adata_liver: {adata_liver.shape}")
print(f"  adata_blood: {adata_blood.shape}")

def compute_pathway_scores(ad, gene_sets_dict):
    """Compute mean expression-based pathway scores for each cell."""
    gene_names_ds = set(ad.var_names)
    for gs_name, gs_genes in gene_sets_dict.items():
        avail = [g for g in gs_genes if g in gene_names_ds]
        if len(avail) >= 2:
            gene_indices = [list(ad.var_names).index(g) for g in avail]
            # Extract expression data (handle sparse matrix)
            expr = ad.X[:, gene_indices]
            if hasattr(expr, 'toarray'):
                expr = expr.toarray()
            elif hasattr(expr, 'A'): # matrix subclass
                expr = expr.A
            
            # Compute mean per cell
            ad.obs[f'pw_{gs_name}'] = np.mean(expr, axis=1)
        else:
            ad.obs[f'pw_{gs_name}'] = 0.0
    return ad

adata_liver = compute_pathway_scores(adata_liver, gene_sets)
adata_blood = compute_pathway_scores(adata_blood, gene_sets)

print(f"  Computed {len(gene_sets)} pathway scores in {time.time()-t0:.1f}s")

# Verify
pw_cols = [c for c in adata_liver.obs.columns if c.startswith('pw_')]
print(f"  Pathway columns: {len(pw_cols)}")

# %%
# ============================================================
#  CELL 16: DONOR-LEVEL PATHWAY AGGREGATION & MWU TESTS
# ============================================================
print("=" * 70)
print("  DONOR-LEVEL PATHWAY SCORING (Tissue-Separated)")
print("=" * 70)

# Define Comparisons
COMPARISONS = [
    ('NL', 'IT', 'NL→IT'),
    ('NL', 'IA', 'NL→IA'),
    ('NL', 'AR', 'NL→AR'),
    ('NL', 'CR', 'NL→CR'),
    ('IT', 'IA', 'IT→IA'),
    ('IA', 'AR', 'IA→AR'),
    ('CR', 'AR', 'CR→AR')
]

def screen_pathways_tissue(ad, tissue_name, lineage_map, gene_sets, comparisons):
    """Screen all pathways × lineages × comparisons for one tissue."""
    results = []
    pw_names = list(gene_sets.keys())
    
    for lin_std, lin_data in lineage_map.items():
        lin_mask = ad.obs[lineage_col] == lin_data
        sub = ad[lin_mask]
        
        for pw_name in pw_names:
            pw_col = f'pw_{pw_name}'
            if pw_col not in sub.obs.columns:
                continue
            
            pw_vals = sub.obs[pw_col].values
            
            for s1, s2, comp_name in comparisons:
                # Donor-level aggregation
                def get_donor_means(stage):
                    stage_mask = sub.obs[stage_col].values == stage
                    donors = sub.obs.loc[stage_mask, donor_col].unique()
                    means = []
                    for d in donors:
                        d_mask = (sub.obs[stage_col].values == stage) & (sub.obs[donor_col].values == d)
                        if d_mask.sum() >= 5:
                            means.append(float(np.mean(pw_vals[d_mask])))
                    return means
                
                vals1 = get_donor_means(s1)
                vals2 = get_donor_means(s2)
                
                if len(vals1) < 2 or len(vals2) < 2:
                    continue
                
                mean1, mean2 = np.mean(vals1), np.mean(vals2)
                pct_change = ((mean2 - mean1) / abs(mean1) * 100) if abs(mean1) > 1e-10 else 0
                
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
                    'lineage': lin_std,
                    'pathway': pw_name,
                    'comparison': comp_name,
                    'stage1': s1, 'stage2': s2,
                    'mean_s1': round(mean1, 6),
                    'mean_s2': round(mean2, 6),
                    'pct_change': round(pct_change, 1),
                    'direction': direction,
                    'p_value': round(pval, 6),
                    'sig': sig,
                    'consistency': f"{consist}/{n_pairs}",
                    'n_s1': len(vals1),
                    'n_s2': len(vals2),
                })
    
    return pd.DataFrame(results)

t0 = time.time()
df_pw_liver = screen_pathways_tissue(adata_liver, 'Liver', lineage_map, gene_sets, COMPARISONS)
print(f"  Liver: {len(df_pw_liver)} tests, {(df_pw_liver['p_value']<0.05).sum()} sig ({time.time()-t0:.0f}s)")

t0 = time.time()
df_pw_blood = screen_pathways_tissue(adata_blood, 'Blood', lineage_map, gene_sets, COMPARISONS)
print(f"  Blood: {len(df_pw_blood)} tests, {(df_pw_blood['p_value']<0.05).sum()} sig ({time.time()-t0:.0f}s)")

df_pw_liver.to_csv(f'{C4_SAVE}C4_pathway_liver.csv', index=False)
df_pw_blood.to_csv(f'{C4_SAVE}C4_pathway_blood.csv', index=False)
print(f"\n💾 Saved: C4_pathway_liver.csv, C4_pathway_blood.csv")

# %%
# ============================================================
#  CELL 17: NL→IT PATHWAY CHANGES — Liver vs Blood SIMULTANEOUS
# ============================================================
print("=" * 70)
print("  🔴🔵 NL → IT: PATHWAY CHANGES (Liver vs Blood)")
print("=" * 70)

for lin in LINEAGES:  # Corrected from LINEAGES_STD to LINEAGES
    print(f"\n  ── {lin} ──")
    
    l_data = df_pw_liver[(df_pw_liver['comparison'] == 'NL→IT') & (df_pw_liver['lineage'] == lin)]
    b_data = df_pw_blood[(df_pw_blood['comparison'] == 'NL→IT') & (df_pw_blood['lineage'] == lin)]
    
    if len(l_data) == 0:
        print(f"    (no data)")
        continue
    
    merged = l_data[['pathway', 'pct_change', 'p_value', 'sig', 'direction']].merge(
        b_data[['pathway', 'pct_change', 'p_value', 'sig', 'direction']],
        on='pathway', suffixes=('_L', '_B')
    )
    
    # Show significant in either tissue first
    merged['any_sig'] = (merged['p_value_L'] < 0.05) | (merged['p_value_B'] < 0.05)
    merged = merged.sort_values(['any_sig', 'p_value_L'], ascending=[False, True])
    
    for _, row in merged.iterrows():
        if not row['any_sig']:
            continue
        
        l_mark = '★' if row['p_value_L'] < 0.05 else ('†' if row['p_value_L'] < 0.10 else ' ')
        b_mark = '★' if row['p_value_B'] < 0.05 else ('†' if row['p_value_B'] < 0.10 else ' ')
        match = '✅' if row['direction_L'] == row['direction_B'] else '⚠️↔'
        
        print(f"    {row['pathway']:22s} | 🔴{l_mark}{row['direction_L']}{abs(row['pct_change_L']):6.1f}% p={row['p_value_L']:.4f} | 🔵{b_mark}{row['direction_B']}{abs(row['pct_change_B']):6.1f}% p={row['p_value_B']:.4f} | {match}")

# %%
# ============================================================
#  CELL 18: ALL COMPARISONS — Significant Pathways Summary
# ============================================================
print("=" * 70)
print("  ALL COMPARISONS: Significant Pathways Count")
print("=" * 70)

print(f"\n  {'Lineage':<10s} {'Comparison':<10s} {'🔴Liver':>8s} {'🔵Blood':>8s} {'Both':>6s}")
print(f"  {'─'*48}")

for comp_name in [c[2] for c in COMPARISONS]:
    for lin in LINEAGES:  # Corrected from LINEAGES_STD to LINEAGES
        l_n = len(df_pw_liver[(df_pw_liver['comparison']==comp_name) & (df_pw_liver['lineage']==lin) & (df_pw_liver['p_value']<0.05)])
        b_n = len(df_pw_blood[(df_pw_blood['comparison']==comp_name) & (df_pw_blood['lineage']==lin) & (df_pw_blood['p_value']<0.05)])
        
        if l_n > 0 or b_n > 0:
            both = '✅' if l_n > 0 and b_n > 0 else ''
            print(f"  {lin:<10s} {comp_name:<10s} {l_n:>8d} {b_n:>8d} {both:>6s}")

# %%
# ============================================================
#  CELL 19: SELECT SIGNIFICANT PATHWAYS → Pass to C5
# ============================================================
print("=" * 70)
print("  🎯 SIGNIFICANT PATHWAY SELECTION FOR C5")
print("=" * 70)

# For each lineage, collect pathways significant in NL→IT (either tissue)
c5_selections = {}
for lin in LINEAGES:
    l_sig = set(df_pw_liver[(df_pw_liver['comparison']=='NL→IT') & (df_pw_liver['lineage']==lin) & (df_pw_liver['p_value']<0.05)]['pathway'])
    b_sig = set(df_pw_blood[(df_pw_blood['comparison']=='NL→IT') & (df_pw_blood['lineage']==lin) & (df_pw_blood['p_value']<0.05)]['pathway'])
    
    # Also include IT→IA and NL→CR significant
    l_sig2 = set(df_pw_liver[(df_pw_liver['comparison'].isin(['IT→IA','NL→CR'])) & (df_pw_liver['lineage']==lin) & (df_pw_liver['p_value']<0.05)]['pathway'])
    b_sig2 = set(df_pw_blood[(df_pw_blood['comparison'].isin(['IT→IA','NL→CR'])) & (df_pw_blood['lineage']==lin) & (df_pw_blood['p_value']<0.05)]['pathway'])
    
    all_sig = sorted(l_sig | b_sig | l_sig2 | b_sig2)
    c5_selections[lin] = all_sig
    
    print(f"\n  {lin}: {len(all_sig)} pathways")
    for pw in all_sig:
        in_l = '🔴' if pw in l_sig else '  '
        in_b = '🔵' if pw in b_sig else '  '
        print(f"    {in_l}{in_b} {pw}")

# Collect all genes from selected pathways
c5_genes = set()
for lin, pws in c5_selections.items():
    for pw in pws:
        if pw in gene_sets:
            c5_genes.update(gene_sets[pw])

c5_genes = sorted([g for g in c5_genes if g in gene_names])
print(f"\n  ▶ Total unique genes for C5: {len(c5_genes)}")

# Save
selection_rows = []
for lin, pws in c5_selections.items():
    for pw in pws:
        selection_rows.append({'lineage': lin, 'pathway': pw})
pd.DataFrame(selection_rows).to_csv(f'{C4_SAVE}C4_selected_pathways_for_C5.csv', index=False)
pd.DataFrame({'gene': c5_genes}).to_csv(f'{C4_SAVE}C4_selected_genes_for_C5.csv', index=False)

print(f"\n💾 Saved: C4_selected_pathways_for_C5.csv, C4_selected_genes_for_C5.csv")

# %%
# ============================================================
#  CELL 20: C4 SUMMARY
# ============================================================
print("=" * 70)
print("  ✅ C4 COMPLETE — Pathway Scoring Foundation")
print("=" * 70)

print(f"""
  Files saved to: {C4_SAVE}
  ├── C4_pathway_liver.csv          (all results, Liver)
  ├── C4_pathway_blood.csv          (all results, Blood)
  ├── C4_selected_pathways_for_C5.csv
  └── C4_selected_genes_for_C5.csv  ({len(c5_genes)} genes)

  ▶ NEXT: C5 — Individual Gene Markers for selected pathways
""")
