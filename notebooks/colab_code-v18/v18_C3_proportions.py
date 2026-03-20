########################################################################
# V18-C3: CELL LINEAGE PROPORTIONS — Tissue-Separated
# ============================================================
# Step (1) of 7-step pipeline
# Dataset: GSE182159 (243K cells, 23 donors, 5 groups)
# Output: version18-analysis/C3_proportions/
#
# ⚠️ RULES:
#   ❌ NO combined (Liver+Blood) data ever
#   ❌ NO line graphs → dot/box plot only
#   ❌ NEVER "five disease stages"
#   ✅ Liver AND Blood SIMULTANEOUS comparison
#   ✅ Liver = red ●, Blood = blue ▲
#   ✅ P-values always shown
#   ✅ Bottom-up: data speaks first
#
# Disease terminology:
#   "disease spectrum/groups" = NL, IT, IA, AR, CR (all 5)
#   "disease stages" = IT→IA→CR ONLY (chronic HBV)
#   AR = independent comparator (adult acute self-limited HBV)
#   NL = healthy normal control
#
# Comparisons (6 primary):
#   NL→IT, NL→IA, NL→AR, NL→CR, IA-AR(+NL), CR-AR(+NL)
#   Plus: IT→IA, IT→CR (transition/endpoint)
########################################################################

# %% [markdown]
# # V18-C3: Cell Lineage Proportions (Tissue-Separated)
# **Step 1/7**: Identify which cell lineages show significant proportion changes across disease spectrum, separately in Liver and Blood.

# %%
# ============================================================
#  CELL 1: SETUP & DATA LOADING
# ============================================================
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['figure.dpi'] = 150
matplotlib.rcParams['font.family'] = 'Arial'
import warnings
warnings.filterwarnings('ignore')
import os, time

# ── Paths ──
DATA_PATH = '/content/drive/MyDrive/ITLAS/data/processed/GSE182159_gut2021_annotated.h5ad'
SAVE_DIR  = '/content/drive/MyDrive/ITLAS/results/version18-analysis/C3_proportions/'
os.makedirs(SAVE_DIR, exist_ok=True)

print("=" * 70)
print("  V18-C3: CELL LINEAGE PROPORTIONS — Tissue-Separated")
print("=" * 70)

t0 = time.time()
adata = sc.read_h5ad(DATA_PATH)
print(f"✅ Loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes ({time.time()-t0:.1f}s)")
print(f"\nColumns: {list(adata.obs.columns)}")

# %%
# ============================================================
#  CELL 2: IDENTIFY KEY COLUMNS
# ============================================================
# ── Auto-detect columns ──
print("=" * 70)
print("  COLUMN IDENTIFICATION")
print("=" * 70)

# Tissue column
tissue_col = None
for col in ['tissue', 'Tissue', 'source', 'Source', 'compartment', 'Compartment', 'sample_type']:
    if col in adata.obs.columns:
        tissue_col = col
        break
if tissue_col:
    print(f"✅ Tissue: '{tissue_col}' → {adata.obs[tissue_col].value_counts().to_dict()}")
else:
    print("⚠️ No tissue column found! Check columns above and set manually:")
    print("   tissue_col = 'YOUR_COLUMN_NAME'")

# Lineage column
lineage_col = None
for col in ['lineage', 'Lineage', 'cell_type', 'celltype', 'major_cluster', 'lineage_v2']:
    if col in adata.obs.columns:
        lineage_col = col
        break
if lineage_col:
    print(f"✅ Lineage: '{lineage_col}' → {sorted(adata.obs[lineage_col].unique())}")

# Stage and donor
stage_col = None
for col in ['Stage', 'stage', 'disease_stage', 'condition', 'group']:
    if col in adata.obs.columns:
        stage_col = col
        break
donor_col = None
for col in ['donor', 'Donor', 'patient', 'sample_id', 'donor_id']:
    if col in adata.obs.columns:
        donor_col = col
        break

print(f"✅ Stage: '{stage_col}' → {adata.obs[stage_col].value_counts().to_dict()}")
print(f"✅ Donor: '{donor_col}' → {adata.obs[donor_col].nunique()} unique")

# %%
# ============================================================
#  CELL 3: TISSUE SEPARATION & VALIDATION
# ============================================================
print("=" * 70)
print("  TISSUE SEPARATION")
print("=" * 70)

# ⚠️ ADJUST THESE IF NEEDED based on Cell 2 output
LIVER_VAL = 'Liver'   # e.g., 'liver', 'Intrahepatic', 'Liver'
BLOOD_VAL = 'Blood'   # e.g., 'blood', 'PBMC', 'Peripheral', 'Blood'

adata_liver = adata[adata.obs[tissue_col] == LIVER_VAL].copy()
adata_blood = adata[adata.obs[tissue_col] == BLOOD_VAL].copy()

print(f"  Liver: {adata_liver.shape[0]:,} cells ({adata_liver.shape[0]/adata.shape[0]*100:.1f}%)")
print(f"  Blood: {adata_blood.shape[0]:,} cells ({adata_blood.shape[0]/adata.shape[0]*100:.1f}%)")

# Donor × Stage × Tissue breakdown
print(f"\n{'─'*70}")
print(f"  {'Group':<6} {'Liver donors':>14} {'Liver cells':>13} {'Blood donors':>14} {'Blood cells':>13} {'L:B ratio':>10}")
print(f"{'─'*70}")

for stage in ['NL', 'IT', 'IA', 'AR', 'CR']:
    l_mask = adata_liver.obs[stage_col] == stage
    b_mask = adata_blood.obs[stage_col] == stage
    l_donors = adata_liver.obs.loc[l_mask, donor_col].nunique()
    b_donors = adata_blood.obs.loc[b_mask, donor_col].nunique()
    l_cells = l_mask.sum()
    b_cells = b_mask.sum()
    ratio = f"1:{b_cells/l_cells:.1f}" if l_cells > 0 else "N/A"
    print(f"  {stage:<6} {l_donors:>14} {l_cells:>13,} {b_donors:>14} {b_cells:>13,} {ratio:>10}")

# %%
# ============================================================
#  CELL 4: LINEAGE MAPPING
# ============================================================
print("=" * 70)
print("  LINEAGE MAPPING")
print("=" * 70)

# Standard 6 lineages for analysis (exclude gdT — too few cells)
LINEAGES_STD = ['Myeloid', 'CD4T', 'CD8T', 'NK', 'B', 'PlasmaB']

lineage_values = sorted(adata.obs[lineage_col].unique())
print(f"Dataset lineages: {lineage_values}")

# Auto-map (adjust manually if needed)
lineage_map = {}
for lv in lineage_values:
    lvl = lv.lower().replace('_', '').replace(' ', '')
    if 'cd4' in lvl:
        lineage_map['CD4T'] = lv
    elif 'cd8' in lvl:
        lineage_map['CD8T'] = lv
    elif 'myeloid' in lvl or 'macro' in lvl or 'mono' in lvl or 'dc' in lvl:
        lineage_map['Myeloid'] = lv
    elif lvl == 'nk':
        lineage_map['NK'] = lv
    elif 'plasma' in lvl:
        lineage_map['PlasmaB'] = lv
    elif lvl in ['b', 'bcell']:
        lineage_map['B'] = lv

print(f"\nMapping:")
for std, ds in lineage_map.items():
    # Count cells
    n = (adata.obs[lineage_col] == ds).sum()
    print(f"  {std:10s} → {ds:15s} ({n:,} cells)")

unmapped = [lv for lv in lineage_values if lv not in lineage_map.values()]
if unmapped:
    print(f"\n⚠️ Excluded: {unmapped}")

# %%
# ============================================================
#  CELL 5: COMPUTE DONOR-LEVEL PROPORTIONS (TISSUE-SEPARATED)
# ============================================================
print("=" * 70)
print("  DONOR-LEVEL PROPORTIONS (TISSUE-SEPARATED)")
print("=" * 70)

def compute_proportions_by_tissue(adata_tissue, tissue_name):
    """Compute per-donor lineage proportions within a specific tissue."""
    rows = []
    donors = adata_tissue.obs[donor_col].unique()
    
    for d in donors:
        d_mask = adata_tissue.obs[donor_col] == d
        d_stage = adata_tissue.obs.loc[d_mask, stage_col].iloc[0]
        d_total = d_mask.sum()
        
        for lin_std, lin_data in lineage_map.items():
            lin_mask = d_mask & (adata_tissue.obs[lineage_col] == lin_data)
            count = lin_mask.sum()
            pct = (count / d_total * 100) if d_total > 0 else 0
            
            rows.append({
                'donor': d,
                'stage': d_stage,
                'tissue': tissue_name,
                'lineage': lin_std,
                'count': count,
                'total': d_total,
                'pct': round(pct, 2),
            })
    
    return pd.DataFrame(rows)

df_prop_liver = compute_proportions_by_tissue(adata_liver, 'Liver')
df_prop_blood = compute_proportions_by_tissue(adata_blood, 'Blood')

print(f"  Liver proportions: {len(df_prop_liver)} rows ({df_prop_liver['donor'].nunique()} donors)")
print(f"  Blood proportions: {len(df_prop_blood)} rows ({df_prop_blood['donor'].nunique()} donors)")

# Save
df_prop_liver.to_csv(f'{SAVE_DIR}C3_proportions_liver.csv', index=False)
df_prop_blood.to_csv(f'{SAVE_DIR}C3_proportions_blood.csv', index=False)
print(f"\n💾 Saved: C3_proportions_liver.csv, C3_proportions_blood.csv")

# %%
# ============================================================
#  CELL 6: SUMMARY TABLE — Liver vs Blood (Mean ± SD)
# ============================================================
print("=" * 70)
print("  PROPORTION SUMMARY: Liver vs Blood (Mean ± SD)")
print("=" * 70)

def make_summary(df_prop):
    return df_prop.groupby(['stage', 'lineage'])['pct'].agg(['mean', 'std']).round(1)

liver_summary = make_summary(df_prop_liver)
blood_summary = make_summary(df_prop_blood)

STAGES = ['NL', 'IT', 'IA', 'AR', 'CR']

for lin in LINEAGES_STD:
    print(f"\n  ── {lin} ──")
    print(f"  {'':8s}", end='')
    for stage in STAGES:
        print(f"  {stage:>12s}", end='')
    print()
    
    # Liver row
    print(f"  {'Liver':8s}", end='')
    for stage in STAGES:
        try:
            m, s = liver_summary.loc[(stage, lin), ['mean', 'std']]
            print(f"  {m:5.1f}±{s:4.1f}", end='')
        except:
            print(f"  {'N/A':>12s}", end='')
    print()
    
    # Blood row
    print(f"  {'Blood':8s}", end='')
    for stage in STAGES:
        try:
            m, s = blood_summary.loc[(stage, lin), ['mean', 'std']]
            print(f"  {m:5.1f}±{s:4.1f}", end='')
        except:
            print(f"  {'N/A':>12s}", end='')
    print()

# %%
# ============================================================
#  CELL 7: MANN-WHITNEY U TESTS — Liver AND Blood (SIMULTANEOUS)
# ============================================================
print("=" * 70)
print("  STATISTICAL TESTING: All Comparisons × Liver/Blood")
print("=" * 70)

COMPARISONS = [
    ('NL', 'IT', 'NL→IT'),
    ('NL', 'IA', 'NL→IA'),
    ('NL', 'AR', 'NL→AR'),
    ('NL', 'CR', 'NL→CR'),
    ('IT', 'IA', 'IT→IA'),
    ('IT', 'CR', 'IT→CR'),
    ('IA', 'AR', 'IA→AR'),
    ('CR', 'AR', 'CR→AR'),
]

def run_mwu_proportions(df_prop, tissue_name):
    """Run MWU tests for all lineage × comparison combinations."""
    results = []
    for s1, s2, comp_name in COMPARISONS:
        for lin in LINEAGES_STD:
            vals1 = df_prop[(df_prop['stage'] == s1) & (df_prop['lineage'] == lin)]['pct'].values
            vals2 = df_prop[(df_prop['stage'] == s2) & (df_prop['lineage'] == lin)]['pct'].values
            
            if len(vals1) < 2 or len(vals2) < 2:
                continue
            
            mean1, mean2 = np.mean(vals1), np.mean(vals2)
            diff_ppt = mean2 - mean1  # percentage point difference
            pct_change = ((mean2 - mean1) / abs(mean1) * 100) if mean1 != 0 else 99999
            
            try:
                stat, pval = mannwhitneyu(vals1, vals2, alternative='two-sided')
            except:
                pval = 1.0
            
            # Consistency
            n_pairs = len(vals1) * len(vals2)
            if mean2 >= mean1:
                consist = sum(1 for v2 in vals2 for v1 in vals1 if v2 >= v1)
            else:
                consist = sum(1 for v2 in vals2 for v1 in vals1 if v2 < v1)
            
            direction = '↑' if mean2 > mean1 else ('↓' if mean2 < mean1 else '=')
            sig = '***' if pval<0.001 else ('**' if pval<0.01 else ('*' if pval<0.05 else ('†' if pval<0.10 else 'NS')))
            
            results.append({
                'tissue': tissue_name,
                'lineage': lin,
                'comparison': comp_name,
                'stage1': s1, 'stage2': s2,
                'mean_s1': round(mean1, 2),
                'mean_s2': round(mean2, 2),
                'diff_ppt': round(diff_ppt, 1),
                'pct_change': round(pct_change, 1),
                'direction': direction,
                'p_value': round(pval, 6),
                'sig': sig,
                'consistency': f"{consist}/{n_pairs}",
                'n_s1': len(vals1),
                'n_s2': len(vals2),
            })
    
    return pd.DataFrame(results)

df_mwu_liver = run_mwu_proportions(df_prop_liver, 'Liver')
df_mwu_blood = run_mwu_proportions(df_prop_blood, 'Blood')

print(f"  Liver: {len(df_mwu_liver)} tests, {(df_mwu_liver['p_value'] < 0.05).sum()} significant")
print(f"  Blood: {len(df_mwu_blood)} tests, {(df_mwu_blood['p_value'] < 0.05).sum()} significant")

# Save
df_mwu_liver.to_csv(f'{SAVE_DIR}C3_mwu_liver.csv', index=False)
df_mwu_blood.to_csv(f'{SAVE_DIR}C3_mwu_blood.csv', index=False)

# %%
# ============================================================
#  CELL 8: SIMULTANEOUS COMPARISON TABLE — Liver vs Blood
# ============================================================
print("=" * 70)
print("  🔴🔵 SIMULTANEOUS COMPARISON: Liver vs Blood")
print("=" * 70)

# Merge liver and blood results
liver_cols = df_mwu_liver[['lineage', 'comparison', 'mean_s1', 'mean_s2', 'diff_ppt', 'pct_change', 'p_value', 'sig', 'direction']].copy()
liver_cols.columns = ['lineage', 'comparison', 'L_mean1', 'L_mean2', 'L_diff', 'L_pct', 'L_p', 'L_sig', 'L_dir']

blood_cols = df_mwu_blood[['lineage', 'comparison', 'mean_s1', 'mean_s2', 'diff_ppt', 'pct_change', 'p_value', 'sig', 'direction']].copy()
blood_cols.columns = ['lineage', 'comparison', 'B_mean1', 'B_mean2', 'B_diff', 'B_pct', 'B_p', 'B_sig', 'B_dir']

df_simul = liver_cols.merge(blood_cols, on=['lineage', 'comparison'], how='outer')

# Concordance analysis
df_simul['dir_same'] = df_simul['L_dir'] == df_simul['B_dir']
df_simul['both_sig'] = (df_simul['L_p'] < 0.05) & (df_simul['B_p'] < 0.05)
df_simul['liver_only'] = (df_simul['L_p'] < 0.05) & (df_simul['B_p'] >= 0.05)
df_simul['blood_only'] = (df_simul['L_p'] >= 0.05) & (df_simul['B_p'] < 0.05)
df_simul['neither'] = (df_simul['L_p'] >= 0.05) & (df_simul['B_p'] >= 0.05)

# Print simultaneous table for each comparison
for comp_name in [c[2] for c in COMPARISONS]:
    comp_data = df_simul[df_simul['comparison'] == comp_name].sort_values('lineage')
    
    # Only show if at least one significant in either tissue
    any_sig = comp_data[(comp_data['L_p'] < 0.05) | (comp_data['B_p'] < 0.05)]
    
    print(f"\n  ── {comp_name} ──")
    if len(any_sig) == 0:
        print(f"    (No significant findings in either tissue)")
        continue
    
    print(f"  {'Lineage':<10s} {'🔴Liver':>28s}  {'🔵Blood':>28s}  {'Match':>6s}")
    print(f"  {'':─<10s} {'':─>28s}  {'':─>28s}  {'':─>6s}")
    
    for _, row in comp_data.iterrows():
        l_str = f"{row['L_dir']}{abs(row['L_diff']):.1f}pp ({row['L_sig']})"
        b_str = f"{row['B_dir']}{abs(row['B_diff']):.1f}pp ({row['B_sig']})"
        match = "✅" if row['dir_same'] else "⚠️↔"
        
        # Highlight significant
        l_mark = "★" if row['L_p'] < 0.05 else " "
        b_mark = "★" if row['B_p'] < 0.05 else " "
        
        print(f"  {row['lineage']:<10s} {l_mark}{l_str:>27s}  {b_mark}{b_str:>27s}  {match:>6s}")

# Save simultaneous comparison
df_simul.to_csv(f'{SAVE_DIR}C3_simultaneous_comparison.csv', index=False)
print(f"\n💾 Saved: C3_simultaneous_comparison.csv")

# %%
# ============================================================
#  CELL 9: FIGURE — Dot Plots (Liver vs Blood, Side by Side)
# ============================================================
print("=" * 70)
print("  FIGURE: Proportion Dot Plots")
print("=" * 70)

COLORS = {'NL': '#2ecc71', 'IT': '#e74c3c', 'IA': '#e67e22', 'AR': '#3498db', 'CR': '#9b59b6'}
STAGE_ORDER = ['NL', 'IT', 'IA', 'AR', 'CR']

fig, axes = plt.subplots(len(LINEAGES_STD), 2, figsize=(14, 3.5 * len(LINEAGES_STD)), 
                         sharex=False, sharey='row')

for i, lin in enumerate(LINEAGES_STD):
    for j, (tissue_name, df_prop) in enumerate([('Liver', df_prop_liver), ('Blood', df_prop_blood)]):
        ax = axes[i, j]
        lin_data = df_prop[df_prop['lineage'] == lin]
        
        for k, stage in enumerate(STAGE_ORDER):
            stage_vals = lin_data[lin_data['stage'] == stage]['pct'].values
            # Jitter
            x = np.random.normal(k, 0.08, len(stage_vals))
            marker = 'o' if tissue_name == 'Liver' else '^'
            ax.scatter(x, stage_vals, c=COLORS[stage], marker=marker, s=60, alpha=0.8, 
                      edgecolors='white', linewidth=0.5, zorder=3)
            # Mean bar
            if len(stage_vals) > 0:
                ax.hlines(np.mean(stage_vals), k-0.2, k+0.2, colors='black', linewidth=2, zorder=4)
        
        ax.set_xticks(range(len(STAGE_ORDER)))
        ax.set_xticklabels(STAGE_ORDER)
        ax.set_ylabel(f'{lin} (%)' if j == 0 else '')
        
        # Title with tissue indicator
        color = '#c0392b' if tissue_name == 'Liver' else '#2980b9'
        ax.set_title(f"{'🔴 ' if tissue_name=='Liver' else '🔵 '}{tissue_name}", 
                     color=color, fontweight='bold', fontsize=11)
        ax.grid(axis='y', alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Add significance stars
        df_mwu = df_mwu_liver if tissue_name == 'Liver' else df_mwu_blood
        sig_tests = df_mwu[(df_mwu['lineage'] == lin) & (df_mwu['comparison'] == 'NL→IT') & (df_mwu['p_value'] < 0.05)]
        if len(sig_tests) > 0:
            p = sig_tests.iloc[0]['p_value']
            star = '***' if p < 0.001 else ('**' if p < 0.01 else '*')
            y_max = ax.get_ylim()[1]
            ax.annotate(star, xy=(0.5, y_max * 0.95), fontsize=12, ha='center', color='red')

plt.suptitle('C3: Lineage Proportions — Liver (●) vs Blood (▲)', 
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(f'{SAVE_DIR}C3_proportion_dotplot_liver_blood.png', dpi=300, bbox_inches='tight')
plt.show()
print(f"💾 Saved: C3_proportion_dotplot_liver_blood.png")

# %%
# ============================================================
#  CELL 10: SIGNIFICANT LINEAGE SELECTION → Pass to C4
# ============================================================
print("=" * 70)
print("  🎯 SIGNIFICANT LINEAGE SELECTION FOR C4")
print("=" * 70)

# Criteria: p<0.05 in at least one comparison, in either tissue
sig_liver = df_mwu_liver[df_mwu_liver['p_value'] < 0.05]['lineage'].unique()
sig_blood = df_mwu_blood[df_mwu_blood['p_value'] < 0.05]['lineage'].unique()
sig_either = sorted(set(sig_liver) | set(sig_blood))
sig_both   = sorted(set(sig_liver) & set(sig_blood))

print(f"  Liver significant lineages: {sorted(sig_liver)}")
print(f"  Blood significant lineages: {sorted(sig_blood)}")
print(f"  Either tissue: {sig_either}")
print(f"  Both tissues: {sig_both}")

# Detailed breakdown
print(f"\n  {'Lineage':<10s} {'Liver sig':>10s} {'Blood sig':>10s} {'Selected':>10s}")
print(f"  {'─'*45}")
for lin in LINEAGES_STD:
    l_n = len(df_mwu_liver[(df_mwu_liver['lineage'] == lin) & (df_mwu_liver['p_value'] < 0.05)])
    b_n = len(df_mwu_blood[(df_mwu_blood['lineage'] == lin) & (df_mwu_blood['p_value'] < 0.05)])
    selected = '✅' if lin in sig_either else '❌'
    print(f"  {lin:<10s} {l_n:>10d} {b_n:>10d} {selected:>10s}")

# But also include trend lineages (p<0.10) for completeness
trend_liver = df_mwu_liver[df_mwu_liver['p_value'] < 0.10]['lineage'].unique()
trend_blood = df_mwu_blood[df_mwu_blood['p_value'] < 0.10]['lineage'].unique()
trend_either = sorted(set(trend_liver) | set(trend_blood))

print(f"\n  Including trends (p<0.10): {trend_either}")

# DECISION: Pass ALL 6 lineages to C4 (even non-significant may show pathway changes)
# Rationale: C1(v17) showed no significant NL→IT proportions, but C2 showed massive pathway changes
C4_LINEAGES = LINEAGES_STD.copy()
print(f"\n  ▶ DECISION: All 6 lineages passed to C4 (proportion ≠ molecular change)")
print(f"    Rationale: v17 showed no NL→IT proportion sig, but massive pathway changes")

# Save selection
pd.DataFrame({'lineage': C4_LINEAGES, 'selected_for_C4': True}).to_csv(
    f'{SAVE_DIR}C3_selected_lineages_for_C4.csv', index=False)

# %%
# ============================================================
#  CELL 11: KEY OBSERVATIONS (Bottom-Up, No Preconception)
# ============================================================
print("=" * 70)
print("  C3 KEY OBSERVATIONS — What the Data Says")
print("=" * 70)

for comp_name in ['NL→IT', 'NL→IA', 'IT→IA', 'NL→CR', 'IA→AR', 'CR→AR']:
    print(f"\n  📊 {comp_name}:")
    
    for lin in LINEAGES_STD:
        l_row = df_mwu_liver[(df_mwu_liver['comparison'] == comp_name) & (df_mwu_liver['lineage'] == lin)]
        b_row = df_mwu_blood[(df_mwu_blood['comparison'] == comp_name) & (df_mwu_blood['lineage'] == lin)]
        
        if len(l_row) == 0 or len(b_row) == 0:
            continue
        
        lr, br = l_row.iloc[0], b_row.iloc[0]
        
        # Direction concordance
        match = "✅" if lr['direction'] == br['direction'] else "⚠️↔"
        
        # Significance markers
        l_mark = f"p={lr['p_value']:.4f}" + (' ★' if lr['p_value'] < 0.05 else (' †' if lr['p_value'] < 0.10 else ''))
        b_mark = f"p={br['p_value']:.4f}" + (' ★' if br['p_value'] < 0.05 else (' †' if br['p_value'] < 0.10 else ''))
        
        print(f"    {lin:10s} | Liver:{lr['direction']}{abs(lr['diff_ppt']):5.1f}pp {l_mark} | Blood:{br['direction']}{abs(br['diff_ppt']):5.1f}pp {b_mark} | {match}")

# %%
# ============================================================
#  CELL 12: TISSUE DISCREPANCY SUMMARY
# ============================================================
print("=" * 70)
print("  TISSUE DISCREPANCY: Direction Reversals & Significance Divergence")
print("=" * 70)

# Direction reversals (opposite in Liver vs Blood)
reversals = df_simul[~df_simul['dir_same'] & ((df_simul['L_p'] < 0.10) | (df_simul['B_p'] < 0.10))]
if len(reversals) > 0:
    print(f"\n  🚨 DIRECTION REVERSALS (with trend significance):")
    for _, row in reversals.iterrows():
        print(f"    {row['lineage']:10s} {row['comparison']:8s} | Liver:{row['L_dir']}{abs(row['L_diff']):.1f}pp({row['L_sig']}) vs Blood:{row['B_dir']}{abs(row['B_diff']):.1f}pp({row['B_sig']})")
else:
    print(f"\n  ✅ No direction reversals with trend significance")

# Liver-only significant
liver_only = df_simul[df_simul['liver_only']]
if len(liver_only) > 0:
    print(f"\n  🔴 LIVER-ONLY SIGNIFICANT ({len(liver_only)}):")
    for _, row in liver_only.iterrows():
        print(f"    {row['lineage']:10s} {row['comparison']:8s} | Liver:{row['L_dir']}{abs(row['L_diff']):.1f}pp({row['L_sig']}) vs Blood:({row['B_sig']})")

# Blood-only significant
blood_only = df_simul[df_simul['blood_only']]
if len(blood_only) > 0:
    print(f"\n  🔵 BLOOD-ONLY SIGNIFICANT ({len(blood_only)}):")
    for _, row in blood_only.iterrows():
        print(f"    {row['lineage']:10s} {row['comparison']:8s} | Liver:({row['L_sig']}) vs Blood:{row['B_dir']}{abs(row['B_diff']):.1f}pp({row['B_sig']})")

# Both significant
both_sig = df_simul[df_simul['both_sig']]
if len(both_sig) > 0:
    print(f"\n  🟢 BOTH TISSUES SIGNIFICANT ({len(both_sig)}):")
    for _, row in both_sig.iterrows():
        print(f"    {row['lineage']:10s} {row['comparison']:8s} | Liver:{row['L_dir']}{abs(row['L_diff']):.1f}pp({row['L_sig']}) & Blood:{row['B_dir']}{abs(row['B_diff']):.1f}pp({row['B_sig']})")

# %%
# ============================================================
#  CELL 13: FINAL SUMMARY & FILES
# ============================================================
print("=" * 70)
print("  ✅ C3 COMPLETE — Proportion Foundation")
print("=" * 70)

print(f"""
  Files saved to: {SAVE_DIR}
  ├── C3_proportions_liver.csv     (donor-level proportions, Liver)
  ├── C3_proportions_blood.csv     (donor-level proportions, Blood)
  ├── C3_mwu_liver.csv             (MWU test results, Liver)
  ├── C3_mwu_blood.csv             (MWU test results, Blood)
  ├── C3_simultaneous_comparison.csv (Liver vs Blood side-by-side)
  ├── C3_proportion_dotplot_liver_blood.png
  └── C3_selected_lineages_for_C4.csv

  ▶ NEXT: C4 — 26 Pathway Scoring for all 6 lineages
  ▶ Pass: C4_LINEAGES = {C4_LINEAGES}
""")
