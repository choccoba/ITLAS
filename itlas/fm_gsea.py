"""FM-GSEA: Foundation Model Gene Set Enrichment Analysis"""

import numpy as np
import pandas as pd
from typing import Dict, Optional

PATHWAYS = {
    'mTOR_signaling': {
        'genes': ['MTOR', 'RPTOR', 'RICTOR', 'AKT1', 'PIK3CA', 'TSC1', 'TSC2',
                  'RHEB', 'EIF4EBP1', 'RPS6KB1', 'ULK1', 'PTEN'],
        'relevance': 'Yu et al. 2025 - HBsAg inhibits mTOR in NK cells'
    },
    'glycolysis': {
        'genes': ['HK1', 'HK2', 'GPI', 'PFKM', 'ALDOA', 'GAPDH', 'PGK1',
                  'ENO1', 'PKM', 'LDHA', 'LDHB', 'SLC2A1', 'SLC2A3'],
        'relevance': 'NK cell metabolic dysfunction in IT phase'
    },
    'oxidative_phosphorylation': {
        'genes': ['MT-ND1', 'MT-ND2', 'MT-ND3', 'MT-ND4', 'MT-ND4L', 'MT-ND5',
                  'MT-CYB', 'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-ATP6', 'MT-ATP8'],
        'relevance': 'Mito-high clusters (21, 23) in IT phase'
    },
    'nk_cell_cytotoxicity': {
        'genes': ['PRF1', 'GZMA', 'GZMB', 'GZMK', 'GNLY', 'FASLG', 'IFNG',
                  'TNF', 'KLRK1', 'KLRD1', 'NCR1', 'NCR3', 'CD226'],
        'relevance': 'NK cell collapse signature (Cluster 15)'
    },
    'il15_signaling': {
        'genes': ['IL15', 'IL15RA', 'IL2RB', 'IL2RG', 'JAK1', 'JAK3',
                  'STAT5A', 'STAT5B', 'SYK', 'LCK'],
        'relevance': 'Yu et al. 2025 - HBsAg binds IL-15Rβ'
    },
    'b_cell_differentiation': {
        'genes': ['PAX5', 'BCL6', 'PRDM1', 'IRF4', 'XBP1', 'CD79A', 'CD79B',
                  'MS4A1', 'CD19', 'CD22', 'JCHAIN', 'MZB1'],
        'relevance': 'B->Plasma differentiation block in IT'
    }
}


def calculate_pathway_scores(adata, pathways=None, method='mean'):
    """Calculate pathway activity scores for each cell."""
    if pathways is None:
        pathways = PATHWAYS
    
    scores = {}
    for name, info in pathways.items():
        genes = info['genes']
        available = [g for g in genes if g in adata.var_names]
        
        if len(available) < 3:
            scores[name] = np.zeros(adata.n_obs)
            continue
        
        expr = adata[:, available].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray()
        scores[name] = np.mean(expr, axis=1)
    
    return pd.DataFrame(scores, index=adata.obs_names)


def compare_pathway_by_stage(adata, stage_col='Stage'):
    """Compare pathway activities across disease stages."""
    from scipy.stats import mannwhitneyu
    
    pathway_scores = calculate_pathway_scores(adata)
    pathway_scores['Stage'] = adata.obs[stage_col].values
    
    results = []
    for pathway in PATHWAYS.keys():
        row = {'Pathway': pathway}
        for stage in pathway_scores['Stage'].unique():
            mask = pathway_scores['Stage'] == stage
            row[f'{stage}_mean'] = pathway_scores.loc[mask, pathway].mean()
        
        if 'IT' in pathway_scores['Stage'].values and 'NL' in pathway_scores['Stage'].values:
            it_scores = pathway_scores.loc[pathway_scores['Stage'] == 'IT', pathway]
            nl_scores = pathway_scores.loc[pathway_scores['Stage'] == 'NL', pathway]
            _, pval = mannwhitneyu(it_scores, nl_scores, alternative='two-sided')
            row['IT_vs_NL_pval'] = pval
        
        results.append(row)
    
    return pd.DataFrame(results)
