# ITLAS: Integrated Transcriptomic Landscape Analysis System

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

ITLAS is a single-cell RNA sequencing reanalysis pipeline for comprehensive characterization of immune landscapes across chronic HBV infection disease stages. Using AUCell pathway scoring combined with individual gene expression analysis across 26 immunological pathways, ITLAS reveals multi-lineage immune tolerance mechanisms in the Immune Tolerant (IT) phase.

**Associated Publication:**
Park YM. *Single-Cell Transcriptomics Reveals Myeloid-Driven Multi-Lineage Immune Tolerance in the IT Phase of Chronic HBV Infection.* Clinical and Molecular Hepatology (submitted).

## Dataset

Reanalysis of **GSE182159** (Zhang et al., *Gut* 2023):
- **243,000 single cells** from peripheral blood
- **59 subclusters** across 7 major lineages
- **23 donors** across 5 disease groups:
  - NL: Normal (healthy controls)
  - IT: Immune Tolerant (vertically-infected chronic HBV)
  - IA: Immune Active (vertically-infected chronic HBV)
  - AR: Acute Resolved (adult-acquired HBV, independent cohort)
  - CR: Chronic Resolved (vertically-infected chronic HBV)

> **Note:** IT > IA > CR represent sequential phases of vertically-infected chronic HBV.
> AR constitutes an independent adult acute HBV recovery cohort -- never part of the chronic continuum.

## Key Findings

### Inflammatory Tolerance (IT Phase)
The IT phase is not immunologically "silent" -- it exhibits active inflammation held in check by brake mechanisms:

1. **Myeloid "Fire + Firewall"**: 9/10 myeloid subclusters simultaneously elevate inflammasome (NLRP3, CASP1, IL1B) and immune evasion genes (LGALS9, IDO1, TNFAIP3, IL1RN)
2. **IT to IA Brake Release**: Transition involves immune evasion gene downregulation while inflammasome activity persists
3. **NK "Quantity not Quality" Decline**: Proportional reduction without per-cell functional impairment
4. **B Cell "Switched but Stuck"**: Class-switch recombination without differentiation to antibody-secreting cells

## Analysis Pipeline

```
Notebook 01: Data QC & annotation (h5ad to 59 subclusters x 5 stages)
Notebook 02: Cell composition (enrichment/depletion analysis)
Notebook 03: AUCell pathway scoring (26 gene sets, pyscenic)
Notebook 04: Global pathway analysis (26 pathways x 5 stages)
Notebook 05: Lineage-stratified analysis (myeloid/NK/CD8T/B deep dive)
Notebook 06: 59-subcluster 3D matrix (subcluster x pathway x stage)
Notebook 07: Individual gene analysis (6 genes/pathway x 26 pathways)
Notebook 08: Donor-level statistics (Mann-Whitney U, per-donor validation)
Notebook 09: Publication figures (Fig 1-5, S1-S10)
Notebook 10: Supplementary tables & data compilation
```

## Installation

```bash
git clone https://github.com/choccoba/ITLAS.git
cd ITLAS
pip install -r requirements.txt
```

## Requirements

- Python >= 3.8
- scanpy >= 1.9
- pyscenic (for AUCell)
- scipy, pandas, numpy, matplotlib, seaborn

## Citation

```bibtex
@article{park2026itlas,
  author = {Park, Young Min},
  title = {Single-Cell Transcriptomics Reveals Myeloid-Driven Multi-Lineage
           Immune Tolerance in the IT Phase of Chronic HBV Infection},
  journal = {Clinical and Molecular Hepatology},
  year = {2026}
}
```

## References

1. Zhang C, et al. (2023). Dissecting the single-cell transcriptome network underlying chronic HBV infection. *Gut*, 72(10), 1903-1914.
2. Yu J, et al. (2025). HBsAg inhibits NK cell function via IL-15Rb-mTOR axis. *Nature*.

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contact

- **Author**: Young Min Park, MD
- **GitHub**: [@choccoba](https://github.com/choccoba)
