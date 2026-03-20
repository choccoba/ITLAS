#!/usr/bin/env python3
"""
V18 Supplementary Table S1/S2 Verification
Google Colab notebook — copy each cell block into separate Colab cells.
Reads actual C2/C3/C4/C5/C9 CSV outputs and cross-checks against S1/S2.
"""

# =================================================================
# CELL 1: SETUP & MOUNT
# =================================================================
from google.colab import drive
drive.mount('/content/drive')

import pandas as pd
import numpy as np
import json, os
from collections import Counter

BASE = '/content/drive/MyDrive/ITLAS/results/version18-analysis'
C3_DIR = f'{BASE}/C3_gene_expression'
C4_DIR = f'{BASE}/C4_pathway'
C5_DIR = f'{BASE}/C5_genes'
C9_DIR = f'{BASE}/C9_method_fixes'

print("=" * 70)
print("  S1/S2 VERIFICATION — Reading actual Colab outputs")
print("=" * 70)

# List available files
for folder, label in [(C3_DIR, 'C3'), (C4_DIR, 'C4'), (C5_DIR, 'C5'), (C9_DIR, 'C9')]:
    if os.path.exists(folder):
        files = os.listdir(folder)
        print(f"\n  {label}: {len(files)} files")
        for f in sorted(files)[:10]:
            print(f"    {f}")
    else:
        print(f"\n  ⚠️ {label}: folder not found at {folder}")

print("\n✅ Cell 1 complete")


# =================================================================
# CELL 2: LOAD C3 GENE LIST (196 genes — GROUND TRUTH)
# =================================================================
print("=" * 70)
print("  LOADING C3 GENE LIST (196 genes)")
print("=" * 70)

c3_gene_file = f'{C3_DIR}/C3_gene_list_196genes.csv'
if os.path.exists(c3_gene_file):
    df_c3_genes = pd.read_csv(c3_gene_file)
    print(f"  C3 gene list: {len(df_c3_genes)} rows")
    print(f"  Columns: {list(df_c3_genes.columns)}")
    print(f"\n  First 10 genes:")
    print(df_c3_genes.head(10).to_string(index=False))
    C3_GENES = sorted(df_c3_genes.iloc[:, 0].tolist())  # first column = gene names
    print(f"\n  Total C3 genes: {len(C3_GENES)}")
else:
    print(f"  ⚠️ File not found: {c3_gene_file}")
    print("  Trying alternative: extract from C3_all_statistics.csv.gz")
    c3_stats = f'{C3_DIR}/C3_all_statistics.csv.gz'
    if os.path.exists(c3_stats):
        df_c3 = pd.read_csv(c3_stats)
        C3_GENES = sorted(df_c3['gene'].unique().tolist())
        print(f"  Extracted {len(C3_GENES)} unique genes from C3_all_statistics")
    else:
        C3_GENES = []
        print("  🔴 Cannot find C3 gene list!")

print(f"\n✅ C3 genes loaded: {len(C3_GENES)}")


# =================================================================
# CELL 3: LOAD C5 GENE LIST (148 genes — GROUND TRUTH)
# =================================================================
print("=" * 70)
print("  LOADING C5 GENE LIST (148 genes)")
print("=" * 70)

# Try C4_selected_genes_for_C5.csv first (this defines which genes went into C5)
c4_sel_file = f'{C4_DIR}/C4_selected_genes_for_C5.csv'
if os.path.exists(c4_sel_file):
    df_c4_sel = pd.read_csv(c4_sel_file)
    print(f"  C4 selected genes: {len(df_c4_sel)} rows")
    print(f"  Columns: {list(df_c4_sel.columns)}")
    C5_GENES_from_C4 = sorted(df_c4_sel.iloc[:, 0].tolist())
    print(f"  Total from C4 selection: {len(C5_GENES_from_C4)}")
else:
    C5_GENES_from_C4 = []
    print(f"  ⚠️ {c4_sel_file} not found")

# Also extract from C5 results directly
c5_liver = f'{C5_DIR}/C5_genes_liver.csv'
c5_blood = f'{C5_DIR}/C5_genes_blood.csv'
C5_GENES_from_results = set()
for f in [c5_liver, c5_blood]:
    if os.path.exists(f):
        df = pd.read_csv(f)
        print(f"\n  {os.path.basename(f)}: {len(df)} rows, columns: {list(df.columns)[:8]}")
        # Find gene column
        gene_col = None
        for candidate in ['gene', 'Gene', 'gene_name', 'symbol']:
            if candidate in df.columns:
                gene_col = candidate
                break
        if gene_col is None and len(df.columns) > 0:
            # Try first column
            gene_col = df.columns[0]
        if gene_col:
            genes = df[gene_col].unique().tolist()
            C5_GENES_from_results.update(genes)
            print(f"    Unique genes: {len(genes)}")

C5_GENES_from_results = sorted(C5_GENES_from_results)
print(f"\n  C5 genes from C4 selection: {len(C5_GENES_from_C4)}")
print(f"  C5 genes from C5 results:   {len(C5_GENES_from_results)}")

# Use whichever is available; prefer C4 selection as it's the definitive list
C5_GENES = C5_GENES_from_C4 if C5_GENES_from_C4 else C5_GENES_from_results
print(f"\n✅ C5 genes loaded: {len(C5_GENES)}")


# =================================================================
# CELL 4: LOAD C2/C4 PATHWAY DEFINITIONS (26 original)
# =================================================================
print("=" * 70)
print("  LOADING PATHWAY DEFINITIONS")
print("=" * 70)

# The 26 original pathways are defined in C2 code.
# We reconstruct them here and verify against C4 pathway results.
GENE_SETS_26 = {
    'inflammasome':     ['NLRP3', 'CASP1', 'IL1B', 'IL18', 'PYCARD', 'GSDMD'],
    'cytotoxicity':     ['GZMB', 'GZMA', 'GZMK', 'PRF1', 'GNLY', 'NKG7', 'FASLG'],
    'checkpoint':       ['PDCD1', 'CTLA4', 'HAVCR2', 'LAG3', 'TIGIT', 'TOX'],
    'exhaustion':       ['TOX', 'PDCD1', 'HAVCR2', 'LAG3', 'TIGIT', 'ENTPD1'],
    'nk_function':      ['NCAM1', 'KLRD1', 'KLRC2', 'KLRK1', 'NCR1', 'NCR3'],
    'nk_il15_dual':     ['IL2RB', 'IL2RG', 'GZMB', 'PRF1', 'GNLY', 'FCGR3A', 'KLRD1', 'KLRC2'],
    'il15_mtor':        ['IL2RB', 'IL2RG', 'JAK1', 'JAK3', 'STAT5A', 'STAT5B', 'MTOR', 'RPTOR', 'RPS6KB1', 'EIF4EBP1'],
    'immune_evasion':   ['CD274', 'PDCD1LG2', 'LGALS9', 'IDO1', 'TGFB1', 'IL10', 'HAVCR2', 'VTCN1'],
    'treg':             ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'TNFRSF18'],
    'naive_t':          ['LEF1', 'TCF7', 'CCR7', 'SELL', 'IL7R'],
    'memory_t':         ['IL7R', 'CD44', 'EOMES', 'BCL6', 'ID3', 'PRDM1', 'TCF7', 'GZMK'],
    'tf_programs':      ['TBX21', 'EOMES', 'GATA3', 'RORC', 'BCL6', 'PRDM1', 'IRF4', 'BATF', 'MAF'],
    'tissue_resident':  ['ITGAE', 'CXCR6', 'ZNF683', 'PRDM1', 'RUNX3'],
    'stemness':         ['TCF7', 'LEF1', 'MYB', 'KLF2', 'SELL', 'IL7R', 'BCL2'],
    'glycolysis':       ['HK1', 'HK2', 'PFKM', 'PKM', 'LDHA', 'SLC2A1', 'ENO1', 'GAPDH'],
    'oxphos':           ['ATP5F1A', 'ATP5F1B', 'NDUFA1', 'NDUFB1', 'UQCRC1', 'COX4I1', 'SDHB'],
    'mito_dysfunction': ['MT-ND1', 'MT-ND2', 'MT-CO1', 'MT-CO2', 'MT-ATP6', 'MT-CYB', 'PINK1', 'PRKN'],
    'metabolic_recovery': ['PPARGC1A', 'TFAM', 'NRF1', 'ESRRA', 'SIRT1', 'SIRT3'],
    'cancer_associated': ['ATM', 'BCL2', 'ZEB1', 'SNAI1', 'TWIST1', 'CDH1', 'VIM', 'MYC', 'CCND1'],
    'fibrosis':         ['COL1A1', 'COL3A1', 'ACTA2', 'FAP', 'TGFB1', 'TGFB2', 'TGFBR1', 'TGFBR2', 'CTGF', 'LOX'],
    'senescence':       ['CDKN1A', 'CDKN2A', 'GLB1', 'TP53', 'RB1', 'SERPINE1', 'IGFBP7'],
    'epigenetics':      ['TET2', 'TOX', 'EZH2', 'DNMT1', 'DNMT3A', 'DNMT3B', 'HDAC1', 'KDM6A'],
    'angiogenesis':     ['HIF1A', 'VEGFA', 'VEGFB', 'KDR', 'FLT1', 'ANGPT1', 'ANGPT2', 'TEK', 'PECAM1'],
    'cell_cycle':       ['MKI67', 'TOP2A', 'PCNA', 'CDK1', 'CDK2', 'CDK4', 'CCNB1', 'CCND1', 'CCNE1', 'RB1'],
    'proliferation':    ['MKI67', 'TOP2A', 'PCNA', 'CDK1', 'CCNB1'],
    'apoptosis':        ['BCL2', 'BAX', 'BAK1', 'CASP3', 'CASP8', 'FAS'],
}

# 3 new pathways from C9
GENE_SETS_3NEW = {
    'antigen_presentation': ['HLA-DRA', 'HLA-DRB1', 'HLA-DPB1', 'HLA-DPA1', 'HLA-DQB1', 'CD74', 'B2M', 'TAP1', 'TAP2', 'CIITA'],
    'type1_ifn':            ['MX1', 'ISG15', 'STAT1', 'STAT2', 'IRF3', 'IRF7', 'IFNAR1', 'OAS1', 'IFIT1', 'DDX58'],
    'tgfb_signaling':       ['TGFB1', 'TGFB2', 'TGFBR1', 'TGFBR2', 'SMAD2', 'SMAD3', 'SMAD4', 'SMAD7', 'ACVR1', 'LTBP1'],
}

ALL_29 = {**GENE_SETS_26, **GENE_SETS_3NEW}

# Verify against C4 pathway results
c4_pw_file = f'{C4_DIR}/C4_pathway_liver.csv'
if os.path.exists(c4_pw_file):
    df_c4 = pd.read_csv(c4_pw_file)
    c4_pathways = sorted(df_c4['pathway'].unique().tolist()) if 'pathway' in df_c4.columns else []
    print(f"  C4 liver pathways in CSV: {len(c4_pathways)}")
    print(f"  Our 26 defined: {len(GENE_SETS_26)}")
    
    # Check if C4 CSV pathways match our definitions
    our_26 = set(GENE_SETS_26.keys())
    csv_set = set(c4_pathways)
    if our_26 == csv_set:
        print("  ✅ C4 pathway names MATCH our 26 definitions exactly")
    else:
        missing = our_26 - csv_set
        extra = csv_set - our_26
        if missing: print(f"  ⚠️ In our list but NOT in C4 CSV: {missing}")
        if extra: print(f"  ⚠️ In C4 CSV but NOT in our list: {extra}")

# Unique gene counts
all_pw_genes_26 = set()
for genes in GENE_SETS_26.values():
    all_pw_genes_26.update(genes)

all_pw_genes_29 = set()
for genes in ALL_29.values():
    all_pw_genes_29.update(genes)

# Multi-pathway genes
gene_count = Counter()
for pw, genes in ALL_29.items():
    for g in genes:
        gene_count[g] += 1
multi = {g for g, c in gene_count.items() if c > 1}

print(f"\n  26 original pathways → {len(all_pw_genes_26)} unique genes")
print(f"  29 pathways (with 3 new) → {len(all_pw_genes_29)} unique genes")
print(f"  Genes in ≥2 pathways: {len(multi)}")
print(f"  Multi-pathway genes: {sorted(multi)}")

print(f"\n✅ Cell 4 complete")


# =================================================================
# CELL 5: CROSS-VERIFICATION — C3 vs C5 overlap
# =================================================================
print("=" * 70)
print("  CROSS-VERIFICATION: C3 ∩ C5")
print("=" * 70)

if C3_GENES and C5_GENES:
    c3_set = set(C3_GENES)
    c5_set = set(C5_GENES)
    shared = c3_set & c5_set
    c3_only = c3_set - c5_set
    c5_only = c5_set - c3_set
    
    print(f"  C3: {len(c3_set)} genes")
    print(f"  C5: {len(c5_set)} genes")
    print(f"  Shared (C3∩C5): {len(shared)}")
    print(f"  C3-only: {len(c3_only)}")
    print(f"  C5-only: {len(c5_only)}")
    
    # Expected from manuscript: C3=196, C5=148, Shared=102
    print(f"\n  === MANUSCRIPT CLAIM CHECK ===")
    print(f"  C3 = {len(c3_set)} (expected 196) → {'✅' if len(c3_set)==196 else '❌ MISMATCH'}")
    print(f"  C5 = {len(c5_set)} (expected 148) → {'✅' if len(c5_set)==148 else '❌ MISMATCH'}")
    print(f"  Shared = {len(shared)} (expected 102) → {'✅' if len(shared)==102 else '⚠️ CHECK'}")
    
    # C9B genes check
    C9B = ['MX1','ISG15','STAT2','IRF3','IRF7','SOCS1','SOCS3','AICDA','JCHAIN',
           'HLA-DRA','HLA-DRB1','HLA-DPB1','HLA-DPA1','CD74','B2M','TAP1',
           'IL1RN','STAT1','IL1B']
    c9b_in_c3 = [g for g in C9B if g in c3_set]
    c9b_not_in_c3 = [g for g in C9B if g not in c3_set]
    print(f"\n  C9B genes in C3: {len(c9b_in_c3)}/19")
    if c9b_not_in_c3:
        print(f"  ⚠️ C9B genes NOT in C3: {c9b_not_in_c3}")
    
    # Key C3-only genes for manuscript
    key_c3_only = ['MX1','ISG15','STAT2','IRF3','IRF7','SOCS1','SOCS3',
                   'AICDA','JCHAIN','HLA-DRA','HLA-DRB1','HLA-DPB1','HLA-DPA1',
                   'CD74','TAP1','IL1RN','STAT1','IL1B','LDHA']
    print(f"\n  Key manuscript genes — C3/C5 status:")
    for g in key_c3_only:
        in3 = '✅' if g in c3_set else '❌'
        in5 = '✅' if g in c5_set else '❌'
        print(f"    {g:<12s}  C3:{in3}  C5:{in5}")
else:
    print("  🔴 Cannot verify — gene lists not loaded")

print(f"\n✅ Cell 5 complete")


# =================================================================
# CELL 6: PATHWAY GENE UNIVERSE — Verify 182 unique genes
# =================================================================
print("=" * 70)
print("  PATHWAY GENE UNIVERSE VERIFICATION")
print("=" * 70)

# Check against h5ad to find which genes are actually present
DATA_PATH = '/content/drive/MyDrive/ITLAS/data/processed/GSE182159_gut2021_annotated.h5ad'
try:
    import scanpy as sc
    adata = sc.read_h5ad(DATA_PATH, backed='r')
    h5ad_genes = set(adata.var_names)
    print(f"  h5ad gene universe: {len(h5ad_genes)} genes")
    
    # Check each pathway
    print(f"\n  {'Pathway':<25s} {'Defined':>7s} {'Found':>5s} {'Missing':>7s}")
    print("  " + "-" * 50)
    
    total_defined = 0
    total_found = 0
    all_found_genes = set()
    all_missing_genes = set()
    
    for pw in sorted(ALL_29.keys()):
        genes = ALL_29[pw]
        found = [g for g in genes if g in h5ad_genes]
        missing = [g for g in genes if g not in h5ad_genes]
        total_defined += len(genes)
        total_found += len(found)
        all_found_genes.update(found)
        all_missing_genes.update(missing)
        
        status = '✅' if not missing else '⚠️'
        print(f"  {status} {pw:<23s} {len(genes):>7d} {len(found):>5d} {len(missing):>7d}", end='')
        if missing:
            print(f"  → {missing}", end='')
        print()
    
    print(f"\n  SUMMARY:")
    print(f"  Total gene mentions (with overlap): {total_defined}")
    print(f"  Unique genes across 29 pathways: {len(all_pw_genes_29)}")
    print(f"  Found in h5ad: {len(all_found_genes)}")
    print(f"  Missing from h5ad: {len(all_missing_genes)} → {sorted(all_missing_genes)}")
    print(f"\n  Manuscript claims 182 unique genes → actual: {len(all_found_genes)}")
    print(f"  {'✅ MATCH' if len(all_found_genes)==182 else '⚠️ DISCREPANCY: ' + str(len(all_found_genes))}")
    
    adata.file.close()
    
except Exception as e:
    print(f"  ⚠️ Cannot load h5ad: {e}")
    print(f"  Using definition-only count: {len(all_pw_genes_29)} unique genes")
    print(f"  Manuscript claims 182 → {'✅' if len(all_pw_genes_29)==182 else '⚠️ ' + str(len(all_pw_genes_29))}")

print(f"\n✅ Cell 6 complete")


# =================================================================
# CELL 7: GENERATE CORRECTED S1/S2 (if discrepancies found)
# =================================================================
print("=" * 70)
print("  FINAL REPORT — S1/S2 VERIFICATION SUMMARY")
print("=" * 70)

print(f"""
  S1 (29 Pathway Definitions):
    Pathways: {len(ALL_29)} (26 original + 3 new)
    Unique genes (definition): {len(all_pw_genes_29)}
    Multi-pathway genes: {len(multi)}
    Expected in manuscript: 29 pathways, 182 unique genes
    
  S2 (Gene Candidates):
    C3: {len(C3_GENES)} genes (expected 196)
    C5: {len(C5_GENES)} genes (expected 148)
    Shared: {len(set(C3_GENES) & set(C5_GENES)) if C3_GENES and C5_GENES else 'N/A'}
    C3-only: {len(set(C3_GENES) - set(C5_GENES)) if C3_GENES and C5_GENES else 'N/A'}
    C5-only: {len(set(C5_GENES) - set(C3_GENES)) if C3_GENES and C5_GENES else 'N/A'}
    C9B: 19 genes
    Expected in manuscript: 196 / 148 / 102 shared
""")

# If any mismatch, output corrected gene lists for S2 regeneration
if C3_GENES and len(C3_GENES) != 196:
    print(f"  🔴 C3 MISMATCH: {len(C3_GENES)} vs 196")
    print(f"     Actual C3 genes saved to: C3_verified_genes.txt")
    with open(f'{BASE}/C3_verified_genes.txt', 'w') as f:
        for g in sorted(C3_GENES):
            f.write(g + '\n')

if C5_GENES and len(C5_GENES) != 148:
    print(f"  🔴 C5 MISMATCH: {len(C5_GENES)} vs 148")
    print(f"     Actual C5 genes saved to: C5_verified_genes.txt")
    with open(f'{BASE}/C5_verified_genes.txt', 'w') as f:
        for g in sorted(C5_GENES):
            f.write(g + '\n')

# Save full verification report
report = {
    'C3_count': len(C3_GENES),
    'C5_count': len(C5_GENES),
    'shared_count': len(set(C3_GENES) & set(C5_GENES)) if C3_GENES and C5_GENES else None,
    'pathway_unique_genes': len(all_pw_genes_29),
    'multi_pathway_genes': len(multi),
    'C3_genes': sorted(C3_GENES),
    'C5_genes': sorted(C5_GENES),
}
with open(f'{BASE}/S1_S2_verification_report.json', 'w') as f:
    json.dump(report, f, indent=2)

print(f"\n  Report saved: {BASE}/S1_S2_verification_report.json")
print(f"\n{'=' * 70}")
print(f"  ✅ S1/S2 VERIFICATION COMPLETE")
print(f"  If all ✅ → S1/S2 are correct as generated")
print(f"  If any ❌ → use verified gene lists to regenerate S1/S2")
print(f"{'=' * 70}")
