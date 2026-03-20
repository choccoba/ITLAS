#!/usr/bin/env python3
"""
==========================================================================
V18 C3: INDIVIDUAL GENE EXPRESSION ANALYSIS
==========================================================================
169 genes × 6 lineages × Liver/Blood (tissue-separated)
Donor-level aggregation + Mann-Whitney U tests
GPU-accelerated (A100) with CuPy sparse matrix operations

DESIGN PRINCIPLES:
  - Bottom-up only: data shows what it shows
  - No combined Liver+Blood analysis
  - All comparisons include NL baseline
  - Dot/box plots only (no line graphs)
  - Liver=red circle, Blood=blue triangle
  - P-values always displayed

OUTPUT: /content/drive/MyDrive/ITLAS/results/version18-analysis/C3_gene_expression/
==========================================================================
"""

# ============================================================
# CELL 0: GPU SETUP & VERIFICATION
# ============================================================
import subprocess
import sys

print("=" * 70)
print("  V18 C3: GPU ENVIRONMENT SETUP")
print("=" * 70)

# GPU detection
try:
    result = subprocess.run(['nvidia-smi', '--query-gpu=name,memory.total,compute_cap',
                            '--format=csv,noheader'], capture_output=True, text=True)
    gpu_info = result.stdout.strip()
    print(f"✅ GPU detected: {gpu_info}")
except:
    print("⚠️ No GPU detected — will use CPU fallback")

# Install CuPy for GPU acceleration
try:
    import cupy as cp
    print(f"✅ CuPy {cp.__version__} ready")
    print(f"   GPU memory: {cp.cuda.Device(0).mem_info[1] / 1e9:.1f} GB total")
except ImportError:
    print("📦 Installing CuPy...")
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'cupy-cuda12x', '-q'])
    import cupy as cp
    print(f"✅ CuPy {cp.__version__} installed")

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import warnings, os, time, gc
warnings.filterwarnings('ignore')

print("✅ All libraries loaded")
print("=" * 70)


# ============================================================
# CELL 1: DATA LOADING
# ============================================================
print("\n" + "=" * 70)
print("  CELL 1: DATA LOADING")
print("=" * 70)

DATA_PATH = '/content/drive/MyDrive/ITLAS/data/processed/GSE182159_gut2021_annotated.h5ad'
RESULTS_DIR = '/content/drive/MyDrive/ITLAS/results/version18-analysis/C3_gene_expression'
os.makedirs(RESULTS_DIR, exist_ok=True)

t0 = time.time()
adata = sc.read_h5ad(DATA_PATH)
print(f"✅ Loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes in {time.time()-t0:.1f}s")

# --- Auto-detect column names (robust to naming variations) ---
print("  Available columns:", list(adata.obs.columns))

# Tissue column
tissue_col = None
for c in ['tissue', 'Tissue', 'tissue_type']:
    if c in adata.obs.columns:
        tissue_col = c
        break
assert tissue_col is not None, f"❌ No tissue column found! Available: {list(adata.obs.columns)}"

# Stage column
stage_col = None
for c in ['Stage', 'stage', 'disease_stage', 'disease_group', 'group', 'condition']:
    if c in adata.obs.columns:
        stage_col = c
        break
assert stage_col is not None, f"❌ No stage column found! Available: {list(adata.obs.columns)}"

# Donor column (critical for donor-level aggregation)
donor_col = None
for c in ['donor', 'sample_id', 'donor_id', 'orig.ident', 'patient', 'subject', 'Patient', 'Donor', 'sample']:
    if c in adata.obs.columns:
        donor_col = c
        break
assert donor_col is not None, f"❌ No donor column found! Available: {list(adata.obs.columns)}"

# Lineage column
lineage_col = None
for c in ['lineage', 'Lineage', 'cell_lineage', 'celltype_major']:
    if c in adata.obs.columns:
        lineage_col = c
        break
assert lineage_col is not None, f"❌ No lineage column found! Available: {list(adata.obs.columns)}"

# Subcluster column
subcluster_col = None
for c in ['gut2021_subcluster_v2', 'subcluster', 'celltype', 'cell_subtype']:
    if c in adata.obs.columns:
        subcluster_col = c
        break
if subcluster_col is None:
    print("⚠️ No subcluster column found — will use lineage-level only")

print(f"  Columns: tissue={tissue_col}, stage={stage_col}, donor={donor_col}, lineage={lineage_col}")
print(f"  Tissues: {adata.obs[tissue_col].value_counts().to_dict()}")
print(f"  Stages: {adata.obs[stage_col].value_counts().to_dict()}")
print(f"  Lineages: {sorted(adata.obs[lineage_col].unique())}")
print(f"  Donors: {adata.obs[donor_col].nunique()} unique")


# ============================================================
# CELL 2: DEFINE 169 GENES
# ============================================================
print("\n" + "=" * 70)
print("  CELL 2: DEFINING 169 TARGET GENES")
print("=" * 70)

# C3 original 152 genes (from v17 analysis)
C3_GENES = [
    # Inflammasome & innate sensing
    'AIM2', 'NLRC4', 'MEFV', 'CASP1', 'NLRP3', 'IL1B', 'GSDMD', 'PYCARD',
    'CASP4', 'CASP5', 'NAIP',
    'TLR2', 'TLR4', 'TLR7', 'TLR8', 'TLR9', 'MYD88', 'TICAM1',
    'DDX58', 'IFIH1', 'MAVS', 'CGAS', 'STING1',
    # Type I IFN
    'IFNA1', 'IFNB1', 'IFNAR1', 'IFNAR2', 'STAT1', 'STAT2', 'IRF3', 'IRF7',
    'IRF9', 'MX1', 'MX2', 'OAS1', 'ISG15', 'IFIT1',
    # Tolerance / suppression
    'LGALS9', 'TGFB1', 'IL10', 'IDO1', 'TNFAIP3', 'HAVCR2', 'CD274',
    'PDCD1LG2', 'VSIR', 'SIGLEC10', 'LILRB1', 'LILRB2',
    # Treg
    'FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'TNFRSF18', 'TIGIT',
    'ENTPD1', 'NT5E', 'LRRC32', 'BACH2', 'PRDM1', 'IL7R',
    # T cell exhaustion
    'PDCD1', 'LAG3', 'TOX', 'EOMES', 'BATF', 'NFATC1',
    'TBX21', 'TCF7', 'SLAMF6', 'CXCR5', 'GZMB',
    # Terminal/progenitor exhaustion
    'ENTPD1', 'HAVCR2', 'CX3CR1', 'PRF1', 'GNLY', 'FGFBP2',
    'XCL1', 'XCL2', 'SELL',
    # Cytotoxicity & NK
    'GZMA', 'GZMK', 'GZMH', 'GZMM', 'NKG7', 'KLRK1', 'KLRD1',
    'NCR1', 'NCR3', 'FCGR3A', 'CD160', 'KIR2DL4',
    # mTOR / metabolism
    'MTOR', 'RPTOR', 'RICTOR', 'RPS6KB1', 'EIF4EBP1', 'AKT1',
    'HIF1A', 'LDHA', 'PKM', 'SLC2A1', 'PFKFB3',
    'NDUFS1', 'COX5A', 'ATP5F1A', 'UQCRC1', 'SDHB',
    'CPT1A', 'ACADVL', 'HADHA', 'PPARGC1A',
    # JAK-STAT
    'JAK1', 'JAK2', 'JAK3', 'TYK2', 'STAT3', 'STAT4', 'STAT5A',
    'STAT5B', 'STAT6', 'SOCS1', 'SOCS3', 'CISH', 'PIAS1',
    # B cell
    'CD19', 'MS4A1', 'CD79A', 'CD79B', 'PAX5', 'BCL6',
    'IRF4', 'XBP1', 'SDC1', 'TNFRSF17', 'MZB1',
    # Apoptosis
    'BCL2', 'MCL1', 'BCL2L1', 'BIRC3', 'CFLAR',
    'BAX', 'BAK1', 'BID', 'BBC3', 'CASP3', 'CASP8', 'FAS', 'FASLG',
    # DNA damage / epigenetic
    'TP53', 'ATM', 'ATR', 'BRCA1', 'CHEK1', 'CHEK2',
    'DNMT1', 'DNMT3A', 'TET2', 'HDAC1', 'EZH2', 'KDM6A',
    # Chemotaxis
    'CCR7', 'CXCR3', 'CXCR4', 'CXCR6', 'CCR2', 'CCR5',
    'CCL3', 'CCL4', 'CCL5', 'CXCL10', 'CXCL13', 'CX3CR1',
    # Antigen presentation
    'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRA', 'HLA-DRB1', 'HLA-DPA1',
    'HLA-DPB1', 'B2M', 'TAP1', 'TAP2', 'CIITA', 'CD74',
    # Cancer-associated
    'TERT', 'MYC', 'VEGFA', 'HGF', 'MET', 'CTNNB1',
    'APC', 'AXIN1', 'TP53', 'RB1', 'CDKN2A', 'MDM2',
]

# C3b supplementary 17 genes
C3B_GENES = [
    'AICDA', 'JCHAIN', 'IL1RN', 'CD27', 'TOX2',
    'RICTOR',  # already in C3 but ensure presence
    'AIM2',    # already in C3 - B cell context
    'MEFV',    # already in C3
    'NLRC4',   # already in C3
    'CASP1',   # already in C3
    'LGALS9',  # already in C3
    'TGFB1',   # already in C3
    'MTOR',    # already in C3
    'JAK1',    # already in C3
    'RPTOR',   # already in C3
    'PYCARD',  # already in C3
    'PRDM1',   # already in C3
]

# Deduplicate and create final list
ALL_GENES = sorted(set(C3_GENES + C3B_GENES))
print(f"Total unique target genes: {len(ALL_GENES)}")

# Check which genes exist in the dataset
available_genes = [g for g in ALL_GENES if g in adata.var_names]
missing_genes = [g for g in ALL_GENES if g not in adata.var_names]
print(f"✅ Available in dataset: {len(available_genes)}")
if missing_genes:
    print(f"⚠️ Missing ({len(missing_genes)}): {missing_genes}")

# Save gene list
gene_df = pd.DataFrame({
    'gene': available_genes,
    'in_C3': [g in C3_GENES for g in available_genes],
    'in_C3b': [g in C3B_GENES for g in available_genes],
    'source': ['C3+C3b' if (g in C3_GENES and g in C3B_GENES) else 'C3' if g in C3_GENES else 'C3b'
               for g in available_genes]
})
gene_df.to_csv(f'{RESULTS_DIR}/C3_gene_list_{len(available_genes)}genes.csv', index=False)
print(f"✅ Gene list saved")


# ============================================================
# CELL 3: GPU-ACCELERATED GENE EXPRESSION EXTRACTION
# ============================================================
print("\n" + "=" * 70)
print("  CELL 3: GPU-ACCELERATED GENE EXPRESSION EXTRACTION")
print("=" * 70)

# Define analysis groups
LINEAGES = ['Myeloid', 'CD4_T', 'CD8_T', 'NK', 'B', 'PlasmaB']
# Map potential alternative names
lineage_map = {}
for lin in adata.obs[lineage_col].unique():
    if 'CD4' in lin: lineage_map[lin] = 'CD4_T'
    elif 'CD8' in lin: lineage_map[lin] = 'CD8_T'
    elif lin.lower() in ['myeloid', 'nk', 'b', 'plasmab', 'gdt']:
        lineage_map[lin] = lin
    else:
        lineage_map[lin] = lin

STAGES = ['NL', 'IT', 'IA', 'AR', 'CR']
TISSUES = ['Liver', 'Blood']

def extract_gene_expression_gpu(adata, genes, tissue_val, batch_size=20000):
    """
    GPU-accelerated gene expression extraction.
    Uses CuPy for sparse matrix slicing and mean computation.
    Returns donor-level mean expression for each gene × lineage.
    """
    # Filter by tissue
    tissue_mask = adata.obs[tissue_col] == tissue_val
    adata_tissue = adata[tissue_mask]
    n_cells = adata_tissue.shape[0]
    print(f"\n  [{tissue_val}] {n_cells:,} cells")

    # Get gene indices
    gene_indices = []
    valid_genes = []
    for g in genes:
        if g in adata_tissue.var_names:
            gene_indices.append(list(adata_tissue.var_names).index(g))
            valid_genes.append(g)

    print(f"  Valid genes: {len(valid_genes)}/{len(genes)}")

    # Extract expression matrix (cells × genes) — sparse
    t0 = time.time()
    X = adata_tissue.X  # sparse matrix

    # Convert to CSC for efficient column slicing
    if sp.issparse(X):
        X_csc = X.tocsc()
    else:
        X_csc = sp.csc_matrix(X)

    # Slice only target genes → dense array (cells × n_genes)
    # Process in batches to avoid memory overflow
    n_genes = len(gene_indices)
    print(f"  Extracting {n_genes} genes from {n_cells:,} cells...")

    # GPU-accelerated: transfer gene columns in batches
    expr_matrix = np.zeros((n_cells, n_genes), dtype=np.float32)

    GENE_BATCH = 50  # Process 50 genes at a time
    for i in range(0, n_genes, GENE_BATCH):
        batch_idx = gene_indices[i:i+GENE_BATCH]
        batch_data = X_csc[:, batch_idx].toarray().astype(np.float32)
        expr_matrix[:, i:i+len(batch_idx)] = batch_data

    print(f"  ✅ Expression extracted in {time.time()-t0:.1f}s")

    # GPU-accelerated donor-level aggregation
    t1 = time.time()
    obs = adata_tissue.obs
    results = []

    for lineage in LINEAGES:
        # Get lineage mask
        lin_mask = obs[lineage_col].map(lineage_map).values == lineage
        if lin_mask.sum() == 0:
            continue

        lin_expr = expr_matrix[lin_mask]  # (cells_in_lineage × n_genes)
        lin_obs = obs[lin_mask]

        # Transfer to GPU for fast aggregation
        lin_expr_gpu = cp.asarray(lin_expr)

        for donor_id in lin_obs[donor_col].unique():
            donor_mask = (lin_obs[donor_col].values == donor_id)
            n_donor_cells = donor_mask.sum()
            if n_donor_cells == 0:
                continue

            stage = lin_obs[lin_obs[donor_col] == donor_id][stage_col].iloc[0]

            # GPU mean computation
            donor_expr_gpu = lin_expr_gpu[cp.asarray(donor_mask)]
            donor_means = cp.asnumpy(cp.mean(donor_expr_gpu, axis=0))

            row = {
                'tissue': tissue_val,
                'lineage': lineage,
                'donor': donor_id,
                'Stage': stage,
                'n_cells': int(n_donor_cells),
            }
            for j, gene in enumerate(valid_genes):
                row[gene] = float(donor_means[j])

            results.append(row)

        # Free GPU memory after each lineage
        del lin_expr_gpu
        cp.get_default_memory_pool().free_all_blocks()

    print(f"  ✅ Donor aggregation in {time.time()-t1:.1f}s ({len(results)} rows)")
    return pd.DataFrame(results), valid_genes


# Run extraction
print("\n" + "=" * 50)
print("  EXTRACTING: LIVER")
print("=" * 50)
liver_df, valid_genes = extract_gene_expression_gpu(adata, available_genes, 'Liver')

print("\n" + "=" * 50)
print("  EXTRACTING: BLOOD")
print("=" * 50)
blood_df, _ = extract_gene_expression_gpu(adata, available_genes, 'Blood')

# Save raw donor-level data
liver_df.to_csv(f'{RESULTS_DIR}/C3_liver_donor_gene_expression.csv.gz',
                index=False, compression='gzip')
blood_df.to_csv(f'{RESULTS_DIR}/C3_blood_donor_gene_expression.csv.gz',
                index=False, compression='gzip')
print(f"\n✅ Saved: Liver ({liver_df.shape}), Blood ({blood_df.shape})")
print(f"✅ Valid genes for analysis: {len(valid_genes)}")


# ============================================================
# CELL 4: STATISTICAL TESTING — Mann-Whitney U (All comparisons)
# ============================================================
print("\n" + "=" * 70)
print("  CELL 4: STATISTICAL TESTING")
print("=" * 70)

def run_statistics(df, tissue_name, genes, comparisons=None):
    """
    Run Mann-Whitney U tests for all gene × lineage × comparison combinations.
    Default comparisons: NL vs IT, NL vs IA, NL vs AR, NL vs CR, IT vs IA, IA vs AR, CR vs AR
    """
    if comparisons is None:
        comparisons = [
            ('NL', 'IT'), ('NL', 'IA'), ('NL', 'AR'), ('NL', 'CR'),
            ('IT', 'IA'), ('IA', 'AR'), ('CR', 'AR'),
        ]

    results = []
    total = len(LINEAGES) * len(genes) * len(comparisons)
    done = 0

    for lineage in LINEAGES:
        lin_df = df[df['lineage'] == lineage]
        if len(lin_df) == 0:
            continue

        for gene in genes:
            if gene not in lin_df.columns:
                continue

            for grp1, grp2 in comparisons:
                vals1 = lin_df[lin_df['Stage'] == grp1][gene].dropna().values
                vals2 = lin_df[lin_df['Stage'] == grp2][gene].dropna().values

                n1, n2 = len(vals1), len(vals2)
                if n1 < 2 or n2 < 2:
                    continue

                mean1, mean2 = np.mean(vals1), np.mean(vals2)

                # Percentage change
                if mean1 > 0:
                    pct_change = ((mean2 - mean1) / mean1) * 100
                elif mean2 > 0:
                    pct_change = 999.0  # 0→positive
                else:
                    pct_change = 0.0

                direction = '↑' if mean2 > mean1 else '↓' if mean2 < mean1 else '='

                # Mann-Whitney U test
                try:
                    U, p = mannwhitneyu(vals1, vals2, alternative='two-sided')
                except:
                    p = 1.0
                    U = 0

                # Consistency: count donor pairs where direction holds
                n_pairs = n1 * n2
                concordant = 0
                for v1 in vals1:
                    for v2 in vals2:
                        if direction == '↑' and v2 > v1:
                            concordant += 1
                        elif direction == '↓' and v2 < v1:
                            concordant += 1
                        elif direction == '=' and v2 == v1:
                            concordant += 1
                consistency = f"{concordant}/{n_pairs}"

                # Significance annotation
                if p < 0.01:
                    sig = '**'
                elif p < 0.05:
                    sig = '*'
                elif p < 0.10:
                    sig = '†'
                else:
                    sig = 'NS'

                results.append({
                    'tissue': tissue_name,
                    'lineage': lineage,
                    'gene': gene,
                    'comparison': f'{grp1}→{grp2}',
                    'grp1': grp1,
                    'grp2': grp2,
                    'n1': n1,
                    'n2': n2,
                    'mean_grp1': round(mean1, 6),
                    'mean_grp2': round(mean2, 6),
                    'pct_change': round(pct_change, 1),
                    'direction': direction,
                    'consistency': consistency,
                    'U_stat': U,
                    'p_value': round(p, 6),
                    'sig': sig,
                })

            done += len(comparisons)
            if done % 500 == 0:
                print(f"  [{tissue_name}] {done}/{total} tests ({done*100//total}%)")

    return pd.DataFrame(results)


t0 = time.time()

print(f"\n--- Liver Statistics ---")
liver_stats = run_statistics(liver_df, 'Liver', valid_genes)

print(f"\n--- Blood Statistics ---")
blood_stats = run_statistics(blood_df, 'Blood', valid_genes)

# Combine
all_stats = pd.concat([liver_stats, blood_stats], ignore_index=True)

# Save
all_stats.to_csv(f'{RESULTS_DIR}/C3_all_statistics.csv.gz', index=False, compression='gzip')

# Summary
n_sig_liver = (liver_stats['p_value'] < 0.05).sum()
n_sig_blood = (blood_stats['p_value'] < 0.05).sum()
n_trend_liver = ((liver_stats['p_value'] >= 0.05) & (liver_stats['p_value'] < 0.10)).sum()
n_trend_blood = ((blood_stats['p_value'] >= 0.05) & (blood_stats['p_value'] < 0.10)).sum()

print(f"\n✅ Statistics complete in {time.time()-t0:.1f}s")
print(f"  Liver: {n_sig_liver} significant (p<0.05), {n_trend_liver} trend (p<0.10), total={len(liver_stats)}")
print(f"  Blood: {n_sig_blood} significant (p<0.05), {n_trend_blood} trend (p<0.10), total={len(blood_stats)}")


# ============================================================
# CELL 5: TOP FINDINGS — NL→IT (IT-specific candidates)
# ============================================================
print("\n" + "=" * 70)
print("  CELL 5: TOP NL→IT FINDINGS (IT-SPECIFIC CANDIDATES)")
print("=" * 70)

it_stats = all_stats[all_stats['comparison'] == 'NL→IT'].copy()
it_sig = it_stats[it_stats['p_value'] < 0.05].sort_values('p_value')

print(f"\nSignificant NL→IT changes (p<0.05): {len(it_sig)}")
print(f"\n{'Tissue':<8} {'Lineage':<10} {'Gene':<12} {'Change':>10} {'Dir':>4} {'Consistency':>12} {'p-value':>10}")
print("-" * 70)
for _, row in it_sig.head(40).iterrows():
    print(f"{row['tissue']:<8} {row['lineage']:<10} {row['gene']:<12} "
          f"{row['pct_change']:>8.1f}% {row['direction']:>4} "
          f"{row['consistency']:>12} {row['p_value']:>10.4f}{row['sig']}")

# IT trends
it_trend = it_stats[(it_stats['p_value'] >= 0.05) & (it_stats['p_value'] < 0.10)].sort_values('p_value')
print(f"\n\nNL→IT Trends (0.05 ≤ p < 0.10): {len(it_trend)}")
for _, row in it_trend.head(20).iterrows():
    print(f"{row['tissue']:<8} {row['lineage']:<10} {row['gene']:<12} "
          f"{row['pct_change']:>8.1f}% {row['direction']:>4} "
          f"{row['consistency']:>12} {row['p_value']:>10.4f}{row['sig']}")

# Save IT-specific results
it_stats.to_csv(f'{RESULTS_DIR}/C3_NL_vs_IT_all.csv', index=False)
it_sig.to_csv(f'{RESULTS_DIR}/C3_NL_vs_IT_significant.csv', index=False)


# ============================================================
# CELL 6: CROSS-TISSUE DISCREPANCY ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("  CELL 6: LIVER vs BLOOD DISCREPANCY (NL→IT)")
print("=" * 70)

# Pivot: for each gene × lineage, compare Liver vs Blood NL→IT
liver_it = it_stats[it_stats['tissue'] == 'Liver'][['lineage', 'gene', 'pct_change', 'p_value', 'sig']].copy()
liver_it.columns = ['lineage', 'gene', 'liver_pct', 'liver_p', 'liver_sig']

blood_it = it_stats[it_stats['tissue'] == 'Blood'][['lineage', 'gene', 'pct_change', 'p_value', 'sig']].copy()
blood_it.columns = ['lineage', 'gene', 'blood_pct', 'blood_p', 'blood_sig']

cross = liver_it.merge(blood_it, on=['lineage', 'gene'], how='outer')
cross['discrepancy'] = abs(cross['liver_pct'].fillna(0) - cross['blood_pct'].fillna(0))
cross['direction_match'] = (
    (cross['liver_pct'].fillna(0) > 0) == (cross['blood_pct'].fillna(0) > 0)
)

# Flag high discrepancies
cross = cross.sort_values('discrepancy', ascending=False)
cross.to_csv(f'{RESULTS_DIR}/C3_liver_blood_discrepancy_NL_IT.csv', index=False)

print(f"\n--- Top 20 Liver vs Blood Discrepancies (NL→IT) ---")
print(f"{'Lineage':<10} {'Gene':<12} {'Liver%':>10} {'Blood%':>10} {'Disc':>8} {'Dir Match':>10}")
print("-" * 65)
for _, row in cross.head(20).iterrows():
    lp = f"{row['liver_pct']:.1f}%{row['liver_sig']}" if pd.notna(row['liver_pct']) else "N/A"
    bp = f"{row['blood_pct']:.1f}%{row['blood_sig']}" if pd.notna(row['blood_pct']) else "N/A"
    match = "✅" if row['direction_match'] else "⚠️"
    print(f"{row['lineage']:<10} {row['gene']:<12} {lp:>10} {bp:>10} {row['discrepancy']:>8.1f} {match:>10}")


# ============================================================
# CELL 7: PATTERN CLASSIFICATION (All comparisons)
# ============================================================
print("\n" + "=" * 70)
print("  CELL 7: PATTERN CLASSIFICATION")
print("=" * 70)

def classify_pattern(gene_lineage_stats):
    """
    Classify gene-lineage combination into disease spectrum patterns.
    Based on NL→X significance patterns:
      - IT-specific: NL→IT sig, NL→IA NS
      - Chronic_persistent: NL→IT sig AND NL→IA sig (same direction)
      - IA-emergent: NL→IA sig, NL→IT NS
      - AR-advantage: NL→AR sig, NL→IT NS, NL→IA NS
      - CR-scar: NL→CR sig (unique to CR)
    """
    nl_it = gene_lineage_stats.get('NL→IT', {})
    nl_ia = gene_lineage_stats.get('NL→IA', {})
    nl_ar = gene_lineage_stats.get('NL→AR', {})
    nl_cr = gene_lineage_stats.get('NL→CR', {})

    it_sig = nl_it.get('p', 1) < 0.05
    ia_sig = nl_ia.get('p', 1) < 0.05
    ar_sig = nl_ar.get('p', 1) < 0.05
    cr_sig = nl_cr.get('p', 1) < 0.05

    it_trend = nl_it.get('p', 1) < 0.10
    ia_trend = nl_ia.get('p', 1) < 0.10

    if it_sig and not ia_sig:
        return 'IT-specific'
    elif it_sig and ia_sig:
        # Check if same direction
        it_dir = nl_it.get('dir', '=')
        ia_dir = nl_ia.get('dir', '=')
        if it_dir == ia_dir:
            return 'Chronic_persistent'
        else:
            return 'IT-specific (reversed in IA)'
    elif not it_sig and ia_sig:
        return 'IA-emergent'
    elif not it_sig and not ia_sig and ar_sig:
        return 'AR-advantage'
    elif not it_sig and not ia_sig and not ar_sig and cr_sig:
        return 'CR-scar'
    elif it_trend and not ia_trend:
        return 'IT-trend'
    else:
        return 'NS'


# Classify for each tissue
pattern_results = []

for tissue in TISSUES:
    tissue_stats = all_stats[all_stats['tissue'] == tissue]

    for lineage in LINEAGES:
        lin_stats = tissue_stats[tissue_stats['lineage'] == lineage]

        for gene in valid_genes:
            gene_stats = lin_stats[lin_stats['gene'] == gene]

            comp_dict = {}
            for _, row in gene_stats.iterrows():
                comp_dict[row['comparison']] = {
                    'p': row['p_value'],
                    'pct': row['pct_change'],
                    'dir': row['direction'],
                    'sig': row['sig'],
                    'consistency': row['consistency'],
                }

            pattern = classify_pattern(comp_dict)

            nl_it = comp_dict.get('NL→IT', {})
            nl_ia = comp_dict.get('NL→IA', {})
            nl_cr = comp_dict.get('NL→CR', {})

            pattern_results.append({
                'tissue': tissue,
                'lineage': lineage,
                'gene': gene,
                'pattern': pattern,
                'NL_IT_pct': nl_it.get('pct', None),
                'NL_IT_p': nl_it.get('p', None),
                'NL_IT_sig': nl_it.get('sig', 'N/A'),
                'NL_IA_pct': nl_ia.get('pct', None),
                'NL_IA_p': nl_ia.get('p', None),
                'NL_IA_sig': nl_ia.get('sig', 'N/A'),
                'NL_CR_pct': nl_cr.get('pct', None),
                'NL_CR_p': nl_cr.get('p', None),
                'NL_CR_sig': nl_cr.get('sig', 'N/A'),
            })

pattern_df = pd.DataFrame(pattern_results)
pattern_df.to_csv(f'{RESULTS_DIR}/C3_pattern_classification.csv', index=False)

# Summary
print("\n--- Pattern Distribution ---")
for tissue in TISSUES:
    t_df = pattern_df[pattern_df['tissue'] == tissue]
    counts = t_df['pattern'].value_counts()
    print(f"\n  [{tissue}]")
    for pat, cnt in counts.items():
        if pat != 'NS':
            print(f"    {pat:<30} {cnt:>5}")

# IT-specific genes
it_specific = pattern_df[pattern_df['pattern'] == 'IT-specific']
print(f"\n\n--- IT-SPECIFIC Genes (tissue-separated) ---")
for tissue in TISSUES:
    t_df = it_specific[it_specific['tissue'] == tissue].sort_values('NL_IT_p')
    print(f"\n  [{tissue}] ({len(t_df)} gene-lineage combinations)")
    for _, row in t_df.head(20).iterrows():
        print(f"    {row['lineage']:<10} {row['gene']:<12} {row['NL_IT_pct']:>8.1f}% "
              f"p={row['NL_IT_p']:.4f}{row['NL_IT_sig']}")


# ============================================================
# CELL 8: VISUALIZATION — Dot/Box Plots for Top Findings
# ============================================================
print("\n" + "=" * 70)
print("  CELL 8: GENERATING DOT/BOX PLOTS")
print("=" * 70)

def plot_gene_dotbox(liver_df, blood_df, gene, lineage, save_dir,
                     liver_stats_df=None, blood_stats_df=None):
    """
    Create side-by-side dot/box plot for a gene × lineage.
    Left: Liver (red circles), Right: Blood (blue triangles)
    All 5 groups shown, NL always included.
    P-values displayed.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    stage_order = ['NL', 'IT', 'IA', 'AR', 'CR']
    stage_colors = {'NL': '#2ecc71', 'IT': '#e74c3c', 'IA': '#e67e22',
                    'AR': '#3498db', 'CR': '#9b59b6'}

    for ax_idx, (df, tissue, color, marker, label) in enumerate([
        (liver_df, 'Liver', '#e74c3c', 'o', 'Liver'),
        (blood_df, 'Blood', '#3498db', '^', 'Blood')
    ]):
        ax = axes[ax_idx]
        lin_df = df[df['lineage'] == lineage]

        if gene not in lin_df.columns or len(lin_df) == 0:
            ax.set_title(f'{tissue}: No data', fontsize=12)
            continue

        positions = []
        box_data = []

        for i, stage in enumerate(stage_order):
            vals = lin_df[lin_df['Stage'] == stage][gene].dropna().values
            if len(vals) > 0:
                # Jittered dots
                jitter = np.random.uniform(-0.15, 0.15, len(vals))
                ax.scatter([i] * len(vals) + jitter, vals,
                          c=stage_colors[stage], marker=marker, s=80,
                          edgecolors='black', linewidths=0.5, alpha=0.8,
                          zorder=3)
                positions.append(i)
                box_data.append(vals)

                # Mean line
                ax.hlines(np.mean(vals), i - 0.3, i + 0.3,
                         colors='black', linewidths=2, zorder=4)

        # Box plot (transparent)
        if box_data:
            bp = ax.boxplot(box_data, positions=positions, widths=0.5,
                           patch_artist=True, showfliers=False, zorder=2)
            for patch in bp['boxes']:
                patch.set_facecolor('white')
                patch.set_alpha(0.3)

        # P-values for NL comparisons
        stats_df = liver_stats_df if tissue == 'Liver' else blood_stats_df
        if stats_df is not None:
            gene_lin_stats = stats_df[
                (stats_df['gene'] == gene) & (stats_df['lineage'] == lineage)
            ]
            for _, row in gene_lin_stats.iterrows():
                if row['grp1'] == 'NL' and row['p_value'] < 0.10:
                    grp2_idx = stage_order.index(row['grp2'])
                    p_text = f"p={row['p_value']:.3f}{row['sig']}"
                    y_max = ax.get_ylim()[1]
                    ax.annotate(p_text, xy=(grp2_idx, y_max * 0.95),
                              fontsize=8, ha='center', color='red' if row['p_value'] < 0.05 else 'gray')

        ax.set_xticks(range(len(stage_order)))
        ax.set_xticklabels(stage_order, fontsize=11)
        ax.set_title(f'{tissue}', fontsize=14, fontweight='bold', color=color)
        ax.set_ylabel('Expression' if ax_idx == 0 else '', fontsize=11)
        ax.grid(axis='y', alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle(f'{gene} — {lineage}', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()

    fname = f'{save_dir}/fig_C3_{gene}_{lineage}.png'
    fig.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return fname


# Create figures directory
fig_dir = f'{RESULTS_DIR}/figures'
os.makedirs(fig_dir, exist_ok=True)

# Generate plots for all significant findings
sig_genes = set()
for _, row in all_stats[all_stats['p_value'] < 0.05].iterrows():
    sig_genes.add((row['gene'], row['lineage']))

# Also add key v17 genes of interest
key_genes = [
    ('AIM2', 'Myeloid'), ('NLRC4', 'Myeloid'), ('MEFV', 'Myeloid'),
    ('CASP1', 'Myeloid'), ('LGALS9', 'Myeloid'), ('TGFB1', 'Myeloid'),
    ('AICDA', 'B'), ('JCHAIN', 'B'), ('JCHAIN', 'PlasmaB'),
    ('MTOR', 'Myeloid'), ('MTOR', 'CD4_T'), ('MTOR', 'NK'),
    ('JAK1', 'Myeloid'), ('JAK1', 'CD4_T'),
    ('RICTOR', 'Myeloid'), ('RICTOR', 'CD4_T'),
    ('IL1RN', 'CD4_T'), ('IL1RN', 'CD8_T'),
    ('AIM2', 'B'), ('TOX2', 'Myeloid'), ('CD27', 'CD4_T'),
    ('PYCARD', 'Myeloid'), ('PYCARD', 'PlasmaB'),
]
sig_genes.update(key_genes)

print(f"Generating {len(sig_genes)} dot/box plots...")
t0 = time.time()

for i, (gene, lineage) in enumerate(sorted(sig_genes)):
    if gene in valid_genes:
        plot_gene_dotbox(liver_df, blood_df, gene, lineage, fig_dir,
                        liver_stats, blood_stats)
    if (i + 1) % 20 == 0:
        print(f"  [{i+1}/{len(sig_genes)}] plots done")

print(f"✅ All plots generated in {time.time()-t0:.1f}s")


# ============================================================
# CELL 9: V17 COMPARISON — Armed Tolerance Genes (Liver-only)
# ============================================================
print("\n" + "=" * 70)
print("  CELL 9: v17 KEY FINDINGS REPLICATION CHECK (Liver-only)")
print("=" * 70)

armed_tolerance_genes = {
    'Sensors': ['AIM2', 'NLRC4', 'MEFV', 'CASP1'],
    'Suppressors': ['LGALS9', 'TGFB1'],
    'Pan-immune': ['MTOR', 'JAK1', 'RPTOR', 'PYCARD'],
    'B cell': ['AICDA', 'JCHAIN', 'PRDM1'],
    'CR scar': ['RICTOR'],
    'C3b new': ['IL1RN', 'CD27', 'TOX2', 'AIM2'],
}

print("\n--- v17 Key Findings: Liver-only Replication ---")
print(f"{'Category':<15} {'Gene':<10} {'Lineage':<10} {'v18 Liver%':>12} {'p-value':>10} {'Sig':>4} {'Consistency':>12}")
print("-" * 80)

liver_nl_it = liver_stats[liver_stats['comparison'] == 'NL→IT']

for category, genes in armed_tolerance_genes.items():
    for gene in genes:
        gene_data = liver_nl_it[liver_nl_it['gene'] == gene]
        for _, row in gene_data.iterrows():
            print(f"{category:<15} {gene:<10} {row['lineage']:<10} "
                  f"{row['pct_change']:>10.1f}% {row['p_value']:>10.4f} {row['sig']:>4} "
                  f"{row['consistency']:>12}")


# ============================================================
# CELL 10: SUMMARY & STATUS
# ============================================================
print("\n" + "=" * 70)
print("  V18 C3: COMPLETE SUMMARY")
print("=" * 70)

print(f"""
V18-C3: INDIVIDUAL GENE EXPRESSION (Tissue-Separated, GPU-accelerated)
================================================================
Genes analyzed: {len(valid_genes)}
Lineages: {', '.join(LINEAGES)}
Tissues: Liver, Blood (NO combined)
Comparisons: NL→IT, NL→IA, NL→AR, NL→CR, IT→IA, IA→AR, CR→AR

SIGNIFICANT (p<0.05):
  Liver: {n_sig_liver}
  Blood: {n_sig_blood}

PATTERN CLASSIFICATION:
""")

for tissue in TISSUES:
    t_df = pattern_df[pattern_df['tissue'] == tissue]
    counts = t_df['pattern'].value_counts()
    print(f"  [{tissue}]")
    for pat, cnt in counts.items():
        if pat != 'NS':
            print(f"    {pat}: {cnt}")

print(f"""
FILES: {RESULTS_DIR}/
  - C3_gene_list_{len(valid_genes)}genes.csv
  - C3_liver_donor_gene_expression.csv.gz
  - C3_blood_donor_gene_expression.csv.gz
  - C3_all_statistics.csv.gz
  - C3_NL_vs_IT_all.csv
  - C3_NL_vs_IT_significant.csv
  - C3_liver_blood_discrepancy_NL_IT.csv
  - C3_pattern_classification.csv
  - figures/ (dot/box plots for significant genes)

STATUS: ✅ C3 COMPLETE — Ready for C4
""")

print("=" * 70)
print("  v18-C3 COMPLETE")
print("=" * 70)