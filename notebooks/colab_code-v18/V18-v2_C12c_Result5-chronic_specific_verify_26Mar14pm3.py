# C12c: IT-Oriented Re-verification — Chronic-Specific vs Normal Hepatitis Response
# Key question: Which findings are chronic HBV-specific (IT/IA ≠ AR) vs general hepatitis (IA ≈ AR)?
# Run after C12/C12b notebooks (uses obs, safe_clonality)

from scipy.stats import mannwhitneyu
import numpy as np, pandas as pd
import warnings; warnings.filterwarnings('ignore')

# Ensure obs and functions are available
try:
    _ = obs['donor']
except:
    import scanpy as sc
    from google.colab import drive
    drive.mount('/content/drive')
    DATA_PATH = '/content/drive/MyDrive/ITLAS/data/processed/GSE182159_gut2021_annotated.h5ad'
    adata = sc.read_h5ad(DATA_PATH, backed='r')
    obs = adata.obs.copy()
    obs['donor'] = obs['sample'].astype(str).str.split('_').str[1]

def safe_clonality(cc):
    nu=len(cc); nt=cc.sum()
    if nu<=0 or nt<=0: return 0.0
    if nu==1: return 1.0 if nt>1 else 0.0
    fr=cc.values/nt; fr=fr[fr>0]
    ent=-np.sum(fr*np.log2(fr))
    return 1-(ent/np.log2(nu)) if np.log2(nu)>0 else 0.0

def mw(a, b, label=''):
    a=a.dropna(); b=b.dropna()
    if len(a)<2 or len(b)<2: return None
    stat,p = mannwhitneyu(a, b, alternative='two-sided')
    am,bm = a.mean(),b.mean()
    d = '↑' if bm>am else '↓'
    pct = ((bm-am)/am*100) if am!=0 else float('inf')
    pt=0; pc=0
    for av in a:
        for bv in b:
            pt+=1
            if (bm>am and bv>av) or (bm<=am and bv<=av): pc+=1
    sig = '★' if p<0.05 else '†' if p<0.10 else ' '
    return {'label':label,'A_m':am,'B_m':bm,'d':d,'pct':pct,'p':p,'cons':f'{pc}/{pt}','sig':sig,'nA':len(a),'nB':len(b)}

# ============================================================
# PART 1: B/PlasmaB SUBCLUSTER PROPORTIONS — Full Pairwise
# ============================================================
print('='*70)
print('PART 1: Liver B/PlasmaB Subcluster — Chronic-Specific Analysis')
print('  Key: IA≈AR → normal hepatitis response')
print('       IT/IA ≠ AR → chronic HBV-specific')
print('='*70)

bp = obs[obs['major_lineage'].isin(['B','PlasmaB'])].copy()
sc_list = sorted(bp['gut2021_subcluster_v2'].unique())

sc_rows = []
for (stage, tissue, donor), grp in bp.groupby(['Stage','tissue','donor'], observed=True):
    total = len(grp)
    if total < 5: continue
    counts = grp['gut2021_subcluster_v2'].value_counts()
    for sc in sc_list:
        sc_rows.append({'Stage':stage,'tissue':tissue,'donor':donor,
                        'subcluster':sc,'proportion':counts.get(sc,0)/total*100})
sc_df = pd.DataFrame(sc_rows)

key_subclusters = ['B_c04-COCH', 'B_c07-FCRL5', 'plasmaB_c01-SDC1', 'plasmaB_c03-MKI67']
comparisons = [('NL','IT'), ('NL','IA'), ('NL','AR'), ('IT','IA'), ('IT','AR'), ('IA','AR')]

for sc in key_subclusters:
    print(f'\n{"─"*70}')
    print(f'  {sc} — LIVER')
    print(f'{"─"*70}')
    sub = sc_df[(sc_df['tissue']=='Liver') & (sc_df['subcluster']==sc)]
    
    # Show individual donor values
    for stage in ['NL','IT','IA','AR','CR']:
        s = sub[sub['Stage']==stage]
        if len(s)==0: continue
        vals = sorted(s['proportion'].values, reverse=True)
        print(f'  {stage} (n={len(s)}): mean={s.proportion.mean():.1f}%, '
              f'donors=[{", ".join(f"{v:.1f}" for v in vals)}]')
    
    # All pairwise
    print(f'\n  Pairwise comparisons:')
    for s1, s2 in comparisons:
        a = sub[sub['Stage']==s1]['proportion']
        b = sub[sub['Stage']==s2]['proportion']
        r = mw(a, b, f'{s1}→{s2}')
        if r:
            print(f'    {r["sig"]} {s1}→{s2}: {r["A_m"]:.1f}%→{r["B_m"]:.1f}% '
                  f'({r["d"]}{abs(r["pct"]):.0f}%) p={r["p"]:.4f} [{r["cons"]}]')
        else:
            print(f'      {s1}→{s2}: insufficient data')
    
    # Chronic-specific test: (IT+IA) vs AR
    chronic = sub[sub['Stage'].isin(['IT','IA'])]['proportion']
    ar = sub[sub['Stage']=='AR']['proportion']
    r_chr = mw(chronic, ar, '(IT+IA) vs AR')
    if r_chr:
        print(f'\n    → (IT+IA) vs AR: {r_chr["A_m"]:.1f}%→{r_chr["B_m"]:.1f}% '
              f'({r_chr["d"]}{abs(r_chr["pct"]):.0f}%) p={r_chr["p"]:.4f} [{r_chr["cons"]}]')
    
    # Classification
    nl_it = mw(sub[sub['Stage']=='NL']['proportion'], sub[sub['Stage']=='IT']['proportion'])
    ia_ar = mw(sub[sub['Stage']=='IA']['proportion'], sub[sub['Stage']=='AR']['proportion'])
    it_ar = mw(sub[sub['Stage']=='IT']['proportion'], sub[sub['Stage']=='AR']['proportion'])
    
    nl_it_sig = nl_it and nl_it['p'] < 0.10
    ia_ar_sig = ia_ar and ia_ar['p'] < 0.10
    it_ar_sig = it_ar and it_ar['p'] < 0.10
    
    if nl_it_sig and not ia_ar_sig:
        interpret = '⚠️ NL→IT sig BUT IA≈AR → may be general hepatitis response'
    elif nl_it_sig and ia_ar_sig:
        interpret = '🔥 NL→IT sig AND IA≠AR → CHRONIC HBV-SPECIFIC'
    elif nl_it_sig and it_ar_sig:
        interpret = '🔥 NL→IT sig AND IT≠AR → CHRONIC HBV-SPECIFIC (IT level)'
    else:
        interpret = '  NL→IT NS'
    print(f'\n    INTERPRETATION: {interpret}')


# ============================================================
# PART 2: BCR METRICS — Chronic-Specific Analysis
# ============================================================
print(f'\n\n{"="*70}')
print('PART 2: BCR Metrics — Chronic-Specific Analysis')
print('='*70)

bcr_rows = []
for (stage, tissue, donor), grp in obs.groupby(['Stage','tissue','donor'], observed=True):
    n = len(grp)
    bcr = grp[grp['BCR_clone.id'].notna()]
    nb = len(bcr)
    if nb == 0:
        bcr_rows.append({'Stage':stage,'tissue':tissue,'donor':donor,
                         'n_bcr':0,'clonality':np.nan,'pct_IgD':np.nan,
                         'pct_IgM':np.nan,'pct_switched':np.nan,'top_clone':0,
                         'pct_IGHV3_23':np.nan,'pct_bcr':0})
        continue
    cc=bcr['BCR_clone.id'].value_counts()
    iso=bcr['BCR_CType'].value_counts(); it=iso.sum()
    vg=bcr['BCR_v_gene'].value_counts()
    bcr_rows.append({'Stage':stage,'tissue':tissue,'donor':donor,
                     'n_bcr':nb,'clonality':safe_clonality(cc),
                     'pct_IgD':iso.get('IGHD',0)/it*100,
                     'pct_IgM':iso.get('IGHM',0)/it*100,
                     'pct_switched':(iso.get('IGHG',0)+iso.get('IGHA',0))/it*100,
                     'top_clone':cc.max(),
                     'pct_IGHV3_23':vg.get('IGHV3-23',0)/nb*100,
                     'pct_bcr':nb/n*100})
bcr_df = pd.DataFrame(bcr_rows)

key_bcr = [('Blood', 'pct_IgD', 'IgD% (naive B depletion)'),
           ('Blood', 'top_clone', 'Top clone size'),
           ('Blood', 'clonality', 'BCR clonality'),
           ('Liver', 'pct_IGHV3_23', 'IGHV3-23%')]

for tissue_val, metric, label in key_bcr:
    print(f'\n{"─"*70}')
    print(f'  {tissue_val} {label}')
    print(f'{"─"*70}')
    df = bcr_df[(bcr_df['tissue']==tissue_val) & (bcr_df['n_bcr']>0)]
    
    for stage in ['NL','IT','IA','AR']:
        s = df[df['Stage']==stage]
        if len(s)==0: continue
        vals = s[metric].values
        print(f'  {stage} (n={len(s)}): mean={s[metric].mean():.2f}, '
              f'donors=[{", ".join(f"{v:.2f}" for v in vals)}]')
    
    print(f'\n  Key comparisons:')
    for s1, s2 in [('NL','IT'),('NL','IA'),('IT','IA'),('IA','AR'),('IT','AR')]:
        a = df[df['Stage']==s1][metric]
        b = df[df['Stage']==s2][metric]
        r = mw(a, b)
        if r:
            print(f'    {r["sig"]} {s1}→{s2}: {r["A_m"]:.2f}→{r["B_m"]:.2f} '
                  f'({r["d"]}{abs(r["pct"]):.0f}%) p={r["p"]:.4f} [{r["cons"]}]')


# ============================================================
# PART 3: TCR METRICS — Chronic-Specific Analysis
# ============================================================
print(f'\n\n{"="*70}')
print('PART 3: TCR Metrics — Chronic-Specific Analysis')
print('='*70)

tcr_rows = []
for (stage, tissue, donor), grp in obs.groupby(['Stage','tissue','donor'], observed=True):
    n = len(grp)
    tcr = grp[grp['TCR_clone.id'].notna()]
    nt = len(tcr)
    if nt == 0:
        tcr_rows.append({'Stage':stage,'tissue':tissue,'donor':donor,
                         'n_tcr':0,'clonality':np.nan,'pct_tcr':0})
        continue
    cc=tcr['TCR_clone.id'].value_counts()
    tcr_rows.append({'Stage':stage,'tissue':tissue,'donor':donor,
                     'n_tcr':nt,'clonality':safe_clonality(cc),'pct_tcr':nt/n*100})
tcr_df = pd.DataFrame(tcr_rows)

for tissue_val in ['Blood', 'Liver']:
    print(f'\n{"─"*70}')
    print(f'  {tissue_val} TCR Clonality')
    print(f'{"─"*70}')
    df = tcr_df[(tcr_df['tissue']==tissue_val) & (tcr_df['n_tcr']>0)]
    
    for stage in ['NL','IT','IA','AR']:
        s = df[df['Stage']==stage]
        if len(s)==0: continue
        print(f'  {stage} (n={len(s)}): mean={s.clonality.mean():.4f}, '
              f'donors=[{", ".join(f"{v:.4f}" for v in s.clonality.values)}]')
    
    print(f'\n  Key comparisons:')
    for s1, s2 in [('NL','IT'),('NL','IA'),('IT','IA'),('IA','AR'),('IT','AR')]:
        a = df[df['Stage']==s1]['clonality']
        b = df[df['Stage']==s2]['clonality']
        r = mw(a, b)
        if r:
            print(f'    {r["sig"]} {s1}→{s2}: {r["A_m"]:.4f}→{r["B_m"]:.4f} '
                  f'({r["d"]}{abs(r["pct"]):.0f}%) p={r["p"]:.4f} [{r["cons"]}]')


# ============================================================
# PART 4: SYNTHESIS — Classification Table
# ============================================================
print(f'\n\n{"="*70}')
print('SYNTHESIS: Chronic HBV-Specific vs General Hepatitis Response')
print('='*70)
print('''
Classification logic:
  - NL→IT sig + IA≈AR → General hepatitis response (not chronic-specific)
  - NL→IT sig + IA≠AR or IT≠AR → Chronic HBV-specific
  - NL→IT sig + IA≈AR but IT≠AR → IT-specific (reversed at IA)

Key question for anti-HBs failure:
  Features present in IT/IA/CR but NOT in AR (which makes anti-HBs)
  = candidates for anti-HBs production blockade mechanism
''')

import os
SAVE_DIR = '/content/drive/MyDrive/ITLAS/results/version18-analysis-v2/BCR_TCR'
os.makedirs(SAVE_DIR, exist_ok=True)

# Save all pairwise results
bcr_df.to_csv(f'{SAVE_DIR}/C12c_BCR_donor_level_full.csv', index=False)
tcr_df.to_csv(f'{SAVE_DIR}/C12c_TCR_donor_level_full.csv', index=False)
sc_df.to_csv(f'{SAVE_DIR}/C12c_subcluster_proportions_full.csv', index=False)
print(f'\nSaved to: {SAVE_DIR}')
