########################################################################
# V18-C6/C7/C8: TISSUE DISCREPANCY → PATTERNS → DISEASE CHARACTERISTICS
# ============================================================
# Steps (4), (5), (6) of 7-step pipeline
# Prereq: C3-C5 completed
# Output: version18-analysis/C6_discrepancy/, C7_patterns/, C8_characteristics/
########################################################################

# %% [markdown]
# # V18-C6: Tissue Discrepancy Master Compilation (Step 4/7)
# # V18-C7: Inter-relationship Pattern Discovery (Step 5/7)
# # V18-C8: Disease Spectrum-Specific Characteristics (Step 6/7)

# %%
# ============================================================
#  C6: TISSUE DISCREPANCY MASTER COMPILATION
#  = User Step (4)
# ============================================================
print("=" * 70)
print("  C6: MASTER TISSUE DISCREPANCY COMPILATION")
print("=" * 70)

C6_SAVE = '/content/drive/MyDrive/ITLAS/results/version18-analysis/C6_discrepancy/'
os.makedirs(C6_SAVE, exist_ok=True)

# Load C3-C5 results
# df_mwu_liver/blood (from C3), df_pw_liver/blood (from C4), df_gene_liver/blood (from C5)

# ── C6-1: Proportion discrepancy (Skipped as df_simul is missing) ──
print("\n  ── C6-1: Proportion Discrepancy (Skipped) ──")
prop_disc = pd.DataFrame() # Empty placeholder

# ── C6-2: Pathway discrepancy ──
print("  ── C6-2: Pathway Discrepancy ──")
pw_disc_rows = []
for comp_name in [c[2] for c in COMPARISONS]:
    for lin in LINEAGES: # Corrected from LINEAGES_STD
        l_pws = df_pw_liver[(df_pw_liver['comparison']==comp_name) & (df_pw_liver['lineage']==lin)]
        b_pws = df_pw_blood[(df_pw_blood['comparison']==comp_name) & (df_pw_blood['lineage']==lin)]
        
        merged = l_pws[['pathway','pct_change','p_value','direction']].merge(
            b_pws[['pathway','pct_change','p_value','direction']],
            on='pathway', suffixes=('_L','_B'), how='outer'
        )
        
        for _, row in merged.iterrows():
            pw_disc_rows.append({
                'level': 'pathway',
                'item': f"{lin}/{row['pathway']}",
                'lineage': lin,
                'comparison': comp_name,
                'L_pct': row.get('pct_change_L', np.nan),
                'B_pct': row.get('pct_change_B', np.nan),
                'L_p': row.get('p_value_L', np.nan),
                'B_p': row.get('p_value_B', np.nan),
                'L_dir': row.get('direction_L', ''),
                'B_dir': row.get('direction_B', ''),
                'dir_same': row.get('direction_L','') == row.get('direction_B',''),
                'both_sig': (row.get('p_value_L',1) < 0.05) and (row.get('p_value_B',1) < 0.05),
                'liver_only': (row.get('p_value_L',1) < 0.05) and (row.get('p_value_B',1) >= 0.05),
                'blood_only': (row.get('p_value_L',1) >= 0.05) and (row.get('p_value_B',1) < 0.05),
            })

df_pw_disc = pd.DataFrame(pw_disc_rows)

# ── C6-3: Gene discrepancy ── (already computed in C5 for NL→IT; extend to all comparisons)
print("  ── C6-3: Gene Discrepancy (all comparisons) ──")
gene_disc_rows = []
for comp_name in [c[2] for c in COMPARISONS]:
    l_genes = df_gene_liver[df_gene_liver['comparison']==comp_name][['gene','lineage','pct_change','p_value','direction']].copy()
    l_genes.columns = ['gene','lineage','L_pct','L_p','L_dir']
    
    b_genes = df_gene_blood[df_gene_blood['comparison']==comp_name][['gene','lineage','pct_change','p_value','direction']].copy()
    b_genes.columns = ['gene','lineage','B_pct','B_p','B_dir']
    
    merged = l_genes.merge(b_genes, on=['gene','lineage'], how='outer')
    merged['comparison'] = comp_name
    merged['level'] = 'gene'
    merged['item'] = merged['gene'] + '/' + merged['lineage']
    merged['dir_same'] = merged['L_dir'] == merged['B_dir']
    merged['both_sig'] = (merged['L_p'] < 0.05) & (merged['B_p'] < 0.05)
    merged['liver_only'] = (merged['L_p'] < 0.05) & (merged['B_p'] >= 0.05)
    merged['blood_only'] = (merged['L_p'] >= 0.05) & (merged['B_p'] < 0.05)
    
    gene_disc_rows.append(merged)

df_gene_disc = pd.concat(gene_disc_rows, ignore_index=True)

# ── C6-4: Master summary ──
print("\n  ── MASTER DISCREPANCY SUMMARY ──")

for level, df in [('Proportion', prop_disc), # Removed condition
                   ('Pathway', df_pw_disc),
                   ('Gene', df_gene_disc)]:
    if len(df) == 0:
        continue
    
    total = len(df)
    dir_same = df['dir_same'].sum() if 'dir_same' in df else 0
    both = df['both_sig'].sum() if 'both_sig' in df else 0
    l_only = df['liver_only'].sum() if 'liver_only' in df else 0
    b_only = df['blood_only'].sum() if 'blood_only' in df else 0
    
    print(f"\n  {level} Level:")
    print(f"    Total pairs: {total}")
    print(f"    Direction concordant: {dir_same} ({dir_same/total*100:.0f}%)")
    print(f"    Both significant: {both}")
    print(f"    Liver-only sig: {l_only}")
    print(f"    Blood-only sig: {b_only}")

# Save master file
df_gene_disc.to_csv(f'{C6_SAVE}C6_gene_discrepancy_master.csv', index=False)
df_pw_disc.to_csv(f'{C6_SAVE}C6_pathway_discrepancy_master.csv', index=False)
print(f"\n💾 Saved: C6_gene_discrepancy_master.csv, C6_pathway_discrepancy_master.csv")

# %%
# ============================================================
#  C7: INTER-RELATIONSHIP PATTERN DISCOVERY
#  = User Step (5)
# ============================================================
print("=" * 70)
print("  C7: INTER-RELATIONSHIP PATTERN DISCOVERY")
print("=" * 70)

C7_SAVE = '/content/drive/MyDrive/ITLAS/results/version18-analysis/C7_patterns/'
os.makedirs(C7_SAVE, exist_ok=True)

# ── C7-1: Select top genes (most frequently significant) ──
print("  ── C7-1: Top Gene Selection ──")

# Combine liver and blood significant genes
all_sig = pd.concat([
    df_gene_liver[df_gene_liver['p_value']<0.05].assign(sig_tissue='Liver'),
    df_gene_blood[df_gene_blood['p_value']<0.05].assign(sig_tissue='Blood'),
])

gene_freq = all_sig.groupby('gene').agg(
    n_sig=('comparison', 'count'),
    n_lineages=('lineage', 'nunique'),
    n_comparisons=('comparison', 'nunique'),
    lineages=('lineage', lambda x: ','.join(sorted(set(x)))),
    tissues=('sig_tissue', lambda x: ','.join(sorted(set(x)))),
    comps=('comparison', lambda x: ','.join(sorted(set(x)))),
).sort_values('n_sig', ascending=False)

# Select top 10-15 genes
TOP_N = min(15, len(gene_freq))
top_genes = gene_freq.head(TOP_N).index.tolist()

print(f"\n  Top {TOP_N} genes selected for pattern analysis:")
for i, gene in enumerate(top_genes):
    row = gene_freq.loc[gene]
    pw = gene_to_pathway.get(gene, '')[:30]
    print(f"  {i+1:2d}. {gene:12s} {row['n_sig']:3d}sig in {row['n_lineages']}lin/{row['n_comparisons']}comp | {row['tissues']:15s} [{pw}]")

pd.DataFrame({'rank': range(1, TOP_N+1), 'gene': top_genes}).to_csv(f'{C7_SAVE}C7_top_genes.csv', index=False)

# ── C7-2: Gene-to-Gene Correlation (within lineage) ──
print("\n  ── C7-2: Gene-to-Gene Correlations ──")
from scipy.stats import spearmanr

def compute_gene_gene_corr(ad, tissue_name, genes, lin_data):
    """Donor-level Spearman correlation between gene pairs."""
    lin_mask = ad.obs[lineage_col] == lin_data
    sub = ad[lin_mask]
    
    # Donor-level means for each gene
    donors = sub.obs[donor_col].unique()
    donor_means = {}
    
    for gene in genes:
        if gene not in sub.var_names:
            continue
        expr = sub[:, gene].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray().flatten()
        else:
            expr = np.array(expr).flatten()
        
        d_means = []
        for d in donors:
            d_mask = sub.obs[donor_col].values == d
            if d_mask.sum() >= 5:
                d_means.append(float(np.mean(expr[d_mask])))
            else:
                d_means.append(np.nan)
        donor_means[gene] = d_means
    
    # Pairwise correlations
    results = []
    gene_list = [g for g in genes if g in donor_means]
    
    for i in range(len(gene_list)):
        for j in range(i+1, len(gene_list)):
            g1, g2 = gene_list[i], gene_list[j]
            v1 = np.array(donor_means[g1])
            v2 = np.array(donor_means[g2])
            
            # Remove NaN
            valid = ~np.isnan(v1) & ~np.isnan(v2)
            if valid.sum() < 5:
                continue
            
            rho, pval = spearmanr(v1[valid], v2[valid])
            results.append({
                'tissue': tissue_name,
                'lineage': lin_std,
                'gene1': g1, 'gene2': g2,
                'spearman_rho': round(rho, 4),
                'p_value': round(pval, 6),
                'n_donors': int(valid.sum()),
                'sig': '*' if pval < 0.05 else 'NS',
            })
    
    return results

corr_results = []
for lin_std, lin_data in lineage_map.items():
    for tissue_name, ad in [('Liver', adata_liver), ('Blood', adata_blood)]:
        corr_results.extend(compute_gene_gene_corr(ad, tissue_name, top_genes, lin_data))

df_corr = pd.DataFrame(corr_results)
if len(df_corr) > 0:
    sig_corr = df_corr[df_corr['p_value'] < 0.05]
    print(f"\n  Total correlations: {len(df_corr)}, Significant: {len(sig_corr)}")
    
    # Show strongest correlations
    print(f"\n  Top 20 strongest significant correlations:")
    for _, row in sig_corr.sort_values('spearman_rho', ascending=False).head(20).iterrows():
        print(f"    {row['gene1']:10s}↔{row['gene2']:10s} ρ={row['spearman_rho']:+.3f} p={row['p_value']:.4f} | {row['tissue']}/{row['lineage']}")
    
    df_corr.to_csv(f'{C7_SAVE}C7_gene_gene_correlations.csv', index=False)

# ── C7-3: Cell-to-Cell Patterns (same gene, multiple lineages) ──
print("\n  ── C7-3: Cell-to-Cell Patterns (Cross-Lineage) ──")

for gene in top_genes[:10]:
    l_data = df_gene_liver[(df_gene_liver['gene']==gene) & (df_gene_liver['comparison']=='NL→IT')]
    b_data = df_gene_blood[(df_gene_blood['gene']==gene) & (df_gene_blood['comparison']=='NL→IT')]
    
    l_sig = l_data[l_data['p_value'] < 0.05]
    b_sig = b_data[b_data['p_value'] < 0.05]
    
    if len(l_sig) >= 2 or len(b_sig) >= 2:  # Multi-lineage
        print(f"\n  {gene} (NL→IT):")
        for lin in LINEAGES:
            lr = l_data[l_data['lineage']==lin]
            br = b_data[b_data['lineage']==lin]
            if len(lr) > 0 and len(br) > 0:
                lr, br = lr.iloc[0], br.iloc[0]
                l_m = '★' if lr['p_value']<0.05 else ' '
                b_m = '★' if br['p_value']<0.05 else ' '
                print(f"    {lin:10s} 🔴{l_m}{lr['direction']}{abs(lr['pct_change']):6.1f}% 🔵{b_m}{br['direction']}{abs(br['pct_change']):6.1f}%")

# %%
# ============================================================
#  C8: DISEASE SPECTRUM-SPECIFIC CHARACTERISTICS
#  = User Step (6)
# ============================================================
print("=" * 70)
print("  C8: DISEASE SPECTRUM-SPECIFIC CHARACTERISTICS")
print("=" * 70)

C8_SAVE = '/content/drive/MyDrive/ITLAS/results/version18-analysis/C8_characteristics/'
os.makedirs(C8_SAVE, exist_ok=True)

# ── IT-Specific: NL→IT sig, NL→IA NS ──
print("\n  ── IT-SPECIFIC GENES (NL→IT sig, NL→IA NS) ──")
for tissue_name, df in [('Liver', df_gene_liver), ('Blood', df_gene_blood)]:
    nl_it = df[df['comparison']=='NL→IT'].set_index(['gene','lineage'])
    nl_ia = df[df['comparison']=='NL→IA'].set_index(['gene','lineage'])
    
    common = nl_it.index.intersection(nl_ia.index)
    it_specific = []
    for idx in common:
        if nl_it.loc[idx, 'p_value'] < 0.05 and nl_ia.loc[idx, 'p_value'] > 0.10:
            it_specific.append({
                'gene': idx[0], 'lineage': idx[1], 'tissue': tissue_name,
                'IT_pct': nl_it.loc[idx, 'pct_change'],
                'IT_p': nl_it.loc[idx, 'p_value'],
                'IA_p': nl_ia.loc[idx, 'p_value'],
            })
    
    df_it_spec = pd.DataFrame(it_specific)
    if len(df_it_spec) > 0:
        print(f"\n  {tissue_name}: {len(df_it_spec)} IT-specific genes")
        for _, row in df_it_spec.sort_values('IT_p').head(10).iterrows():
            print(f"    {row['gene']:12s} {row['lineage']:8s} IT:{row['IT_pct']:+.1f}% p={row['IT_p']:.4f} | IA:p={row['IA_p']:.3f}(NS)")

# ── IA Transition: IT→IA sig ──
print("\n  ── IA TRANSITION GENES (IT→IA sig) ──")
for tissue_name, df in [('Liver', df_gene_liver), ('Blood', df_gene_blood)]:
    it_ia_sig = df[(df['comparison']=='IT→IA') & (df['p_value']<0.05)]
    if len(it_ia_sig) > 0:
        print(f"\n  {tissue_name}: {len(it_ia_sig)} IT→IA transition genes")
        for _, row in it_ia_sig.sort_values('p_value').head(10).iterrows():
            print(f"    {row['gene']:12s} {row['lineage']:8s} {row['direction']}{abs(row['pct_change']):.1f}% p={row['p_value']:.4f}")

# ── CR Scar: NL→CR sig, NL→IA NS ──
print("\n  ── CR SCAR GENES (NL→CR sig, NL→IA NS) ──")
for tissue_name, df in [('Liver', df_gene_liver), ('Blood', df_gene_blood)]:
    nl_cr = df[df['comparison']=='NL→CR'].set_index(['gene','lineage'])
    nl_ia = df[df['comparison']=='NL→IA'].set_index(['gene','lineage'])
    
    common = nl_cr.index.intersection(nl_ia.index)
    cr_scar = []
    for idx in common:
        if nl_cr.loc[idx, 'p_value'] < 0.05 and nl_ia.loc[idx, 'p_value'] > 0.10:
            cr_scar.append({
                'gene': idx[0], 'lineage': idx[1], 'tissue': tissue_name,
                'CR_pct': nl_cr.loc[idx, 'pct_change'],
                'CR_p': nl_cr.loc[idx, 'p_value'],
                'IA_p': nl_ia.loc[idx, 'p_value'],
            })
    
    df_cr_scar = pd.DataFrame(cr_scar)
    if len(df_cr_scar) > 0:
        print(f"\n  {tissue_name}: {len(df_cr_scar)} CR scar genes")
        for _, row in df_cr_scar.sort_values('CR_p').head(10).iterrows():
            print(f"    {row['gene']:12s} {row['lineage']:8s} CR:{row['CR_pct']:+.1f}% p={row['CR_p']:.4f} | IA:p={row['IA_p']:.3f}(NS)")

# ── IA vs AR: What chronic fails but acute achieves ──
print("\n  ── IA vs AR DISCRIMINATORS (chronic vs acute resolution) ──")
for tissue_name, df in [('Liver', df_gene_liver), ('Blood', df_gene_blood)]:
    ia_ar_sig = df[(df['comparison']=='IA→AR') & (df['p_value']<0.05)]
    if len(ia_ar_sig) > 0:
        print(f"\n  {tissue_name}: {len(ia_ar_sig)} IA-AR discriminators")
        for _, row in ia_ar_sig.sort_values('p_value').head(10).iterrows():
            # AR direction: ↑AR means AR higher, which could be recovery
            ar_label = '↑AR' if row['direction'] == '↑' else '↓AR'
            print(f"    {row['gene']:12s} {row['lineage']:8s} {ar_label}{abs(row['pct_change']):.1f}% p={row['p_value']:.4f}")

# %%
# ============================================================
#  C6-C8 COMPLETE
# ============================================================
print("=" * 70)
print("  ✅ C6-C8 COMPLETE")
print("=" * 70)

print(f"""
  C6 Files: {C6_SAVE}
  C7 Files: {C7_SAVE}
  C8 Files: {C8_SAVE}

  ▶ NEXT: C9 — Manuscript Preparation
  ▶ Order: Results → Discussion → Figures/Tables → M&M → 
           Introduction → References → Title/Abstract/Key Points
""")

# %%
# ============================================================
#  SAVE C8 RESULTS TO CSV
# ============================================================
print("=" * 70)
print("  💾 SAVING C8 CHARACTERISTICS TO FILES")
print("=" * 70)

# Containers for aggregated results
c8_it_specific = []
c8_ia_transition = []
c8_cr_scar = []
c8_ia_ar = []

# Re-extract and accumulate data from both tissues
for tissue_name, df in [('Liver', df_gene_liver), ('Blood', df_gene_blood)]:
    # 1. IT-Specific (NL->IT sig, NL->IA NS)
    nl_it = df[df['comparison']=='NL→IT'].set_index(['gene','lineage'])
    nl_ia = df[df['comparison']=='NL→IA'].set_index(['gene','lineage'])
    common_it = nl_it.index.intersection(nl_ia.index)
    
    for idx in common_it:
        if nl_it.loc[idx, 'p_value'] < 0.05 and nl_ia.loc[idx, 'p_value'] > 0.10:
            c8_it_specific.append({
                'tissue': tissue_name, 'lineage': idx[1], 'gene': idx[0],
                'pathway': nl_it.loc[idx, 'pathway'],
                'IT_pct': nl_it.loc[idx, 'pct_change'], 'IT_p': nl_it.loc[idx, 'p_value'],
                'IA_p': nl_ia.loc[idx, 'p_value']
            })

    # 2. IA Transition (IT->IA sig)
    it_ia = df[(df['comparison']=='IT→IA') & (df['p_value']<0.05)]
    for _, row in it_ia.iterrows():
        c8_ia_transition.append(row.to_dict())

    # 3. CR Scar (NL->CR sig, NL->IA NS)
    nl_cr = df[df['comparison']=='NL→CR'].set_index(['gene','lineage'])
    common_cr = nl_cr.index.intersection(nl_ia.index)
    
    for idx in common_cr:
        if nl_cr.loc[idx, 'p_value'] < 0.05 and nl_ia.loc[idx, 'p_value'] > 0.10:
            c8_cr_scar.append({
                'tissue': tissue_name, 'lineage': idx[1], 'gene': idx[0],
                'pathway': nl_cr.loc[idx, 'pathway'],
                'CR_pct': nl_cr.loc[idx, 'pct_change'], 'CR_p': nl_cr.loc[idx, 'p_value'],
                'IA_p': nl_ia.loc[idx, 'p_value']
            })

    # 4. IA vs AR Discriminators (IA->AR sig)
    ia_ar = df[(df['comparison']=='IA→AR') & (df['p_value']<0.05)]
    for _, row in ia_ar.iterrows():
        c8_ia_ar.append(row.to_dict())

# Convert to DataFrames and Save
df_c8_it = pd.DataFrame(c8_it_specific)
df_c8_ia_trans = pd.DataFrame(c8_ia_transition)
df_c8_cr = pd.DataFrame(c8_cr_scar)
df_c8_ar = pd.DataFrame(c8_ia_ar)

df_c8_it.to_csv(f'{C8_SAVE}C8_IT_specific_genes.csv', index=False)
df_c8_ia_trans.to_csv(f'{C8_SAVE}C8_IA_transition_genes.csv', index=False)
df_c8_cr.to_csv(f'{C8_SAVE}C8_CR_scar_genes.csv', index=False)
df_c8_ar.to_csv(f'{C8_SAVE}C8_IA_AR_discriminators.csv', index=False)

# Verify
print(f"✅ Saved files to {C8_SAVE}:")
for f in sorted(os.listdir(C8_SAVE)):
    print(f"  - {f}")

print(f"\nSummary of saved features:")
print(f"  IT-Specific: {len(df_c8_it)} rows")
print(f"  IA-Transition: {len(df_c8_ia_trans)} rows")
print(f"  CR-Scar: {len(df_c8_cr)} rows")
print(f"  IA-AR Discriminators: {len(df_c8_ar)} rows")

# %%
# ============================================================
#  C6-C8 COMPLETE
# ============================================================
print("=" * 70)
print("  ✅ C6-C8 COMPLETE")
print("=" * 70)

print(f"""
  C6 Files: {C6_SAVE}
  C7 Files: {C7_SAVE}
  C8 Files: {C8_SAVE}

  ▶ NEXT: C9 — Manuscript Preparation
  ▶ Order: Results → Discussion → Figures/Tables → M&M → 
           Introduction → References → Title/Abstract/Key Points
""")


