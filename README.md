# ITLAS
ITLAS: IT-phase Liver Immune Atlas Scoring - Single-cell analysis toolkit for identifying immune tolerant phase in chronic HBV infection

# README 파일 생성
readme_content = """# ITLAS: IT-phase Liver Immune Atlas Scoring

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

ITLAS is a single-cell RNA sequencing analysis toolkit designed to identify and characterize the **Immune Tolerant (IT) phase** of chronic Hepatitis B Virus (HBV) infection. The toolkit provides scoring systems, pathway analysis, and machine learning classifiers based on immunometabolic signatures.

## Key Features

- **IT Signature Scoring**: Mito-high cluster-based scoring system for IT phase identification
- **FM-GSEA**: Foundation Model Gene Set Enrichment Analysis for pathway activity quantification
- **Phase Classifier**: XGBoost-based classifier for HBV infection stage prediction (AUC=0.719)
- **Tahoe-x1 Integration**: Validated against Tahoe-x1 foundation model clustering

## Background

Chronic HBV infection progresses through distinct immunological phases:
- **IT (Immune Tolerant)**: High viral load, minimal liver inflammation
- **IA (Immune Active)**: Active immune response, liver damage
- **AR (Acute Resolved)**: Viral clearance
- **AC (Acute Chronic)**: Chronic progression

Understanding IT phase immunopathogenesis is crucial for treatment decisions, as patients in this phase typically do not require immediate antiviral therapy.

## Key Findings

### IT Phase Signatures
| Feature | IT vs NL | Significance |
|---------|----------|--------------|
| OXPHOS | ↑ Elevated | p < 1e-300 |
| NK Cytotoxicity | ↓ Reduced | p < 1e-300 |
| mTOR Signaling | ↑ Elevated | p = 9.68e-108 |
| Glycolysis | ↓ Reduced | p = 2.37e-13 |

### IT-Exclusive Clusters (Tahoe-x1)
| Cluster | Log2 OR in IT | Cell Type |
|---------|---------------|-----------|
| 21 | +11.35 | Mito-high/Quiescence |
| 23 | +10.03 | Mito-high |
| 25 | +7.67 | Naive B cells |
| 15 | -5.15 | NK cells (collapsed) |

## Installation
```bash
git clone https://github.com/choccoba/ITLAS.git
cd ITLAS
pip install -r requirements.txt
```

## Usage

### 1. IT Signature Scoring
```python
from itlas.signature_score import calculate_it_signature_score

scores = calculate_it_signature_score(adata)
adata.obs['IT_score'] = scores['IT_signature']
```

### 2. Pathway Analysis
```python
from itlas.fm_gsea import calculate_pathway_scores

pathway_scores = calculate_pathway_scores(adata)
```

### 3. Phase Classification
```python
from itlas.classifier import predict_phase

predictions = predict_phase(adata)
```

## Project Structure
```
ITLAS/
├── data/
│   ├── raw/                    # Original h5ad files
│   ├── processed/              # Processed data with scores
│   └── embeddings/             # Tahoe-x1 embeddings
├── itlas/
│   ├── __init__.py
│   ├── signature_score.py      # IT signature scoring
│   ├── fm_gsea.py              # Pathway analysis
│   ├── classifier.py           # Phase classifier
│   └── utils.py                # Utilities
├── notebooks/
│   ├── 00_ITLAS_Setup.ipynb
│   ├── 01_IT_Signature_Score.ipynb
│   ├── 02_FM_GSEA.ipynb
│   └── 03_Phase_Classifier.ipynb
├── results/
│   ├── figures/                # Visualizations
│   ├── tables/                 # Statistical results
│   └── models/                 # Trained classifiers
└── README.md
```

## Dataset

Analysis performed on **GSE182159** dataset:
- **243,000 single cells** from chronic HBV patients and healthy controls
- 5 disease stages: NL, IT, IA, AR, AC
- Original publication: Zhang et al., Gut 2023

## Model Performance

### Multi-class Classification (5 stages)
- **Accuracy**: 43.4% (vs 20% random baseline)
- **F1 Score**: 42.9%

### Binary Classification (IT vs non-IT)
- **AUC-ROC**: 0.719

### Top Predictive Features
1. IT_nk_collapse (17.3%)
2. PW_oxidative_phosphorylation (12.8%)
3. PW_nk_cell_cytotoxicity (12.2%)

## Biological Significance

This work supports the **HBsAg-IL15Rβ-mTOR axis** model of immune tolerance:

1. HBsAg binds to IL-15Rβ on NK cells
2. mTOR signaling is disrupted
3. Metabolic shift from glycolysis to OXPHOS
4. NK cell functional exhaustion
5. Immune tolerance established

Reference: Yu et al., Cell Death & Disease, 2025

## Citation

If you use ITLAS in your research, please cite:
```bibtex
@software{itlas2025,
  author = {Young Min Park},
  title = {ITLAS: IT-phase Liver Immune Atlas Scoring},
  year = {2025},
  url = {https://github.com/choccoba/ITLAS}
}
```

## References

1. Zhang et al. (2023). Single-cell transcriptomics of chronic HBV infection. *Gut*.
2. Yu et al. (2025). HBsAg inhibits NK cell function via IL-15Rβ-mTOR axis. *Cell Death & Disease*.
3. Tahoe Bio (2025). Tahoe-x1: Perturbation-trained single-cell foundation model.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- **Author**: YoungMin Park
- **GitHub**: [@chocco_ba](https://github.com/choccoba)

---

*ITLAS - Illuminating Immune Tolerance in Chronic HBV Infection*
"""

# 저장
with open(f"{PROJECT_ROOT}/README.md", 'w') as f:
    f.write(readme_content)
print("✓ Saved: README.md")

# requirements.txt 생성
requirements = """# ITLAS Requirements
scanpy>=1.9.0
anndata>=0.8.0
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
scikit-learn>=1.0.0
xgboost>=1.5.0
matplotlib>=3.5.0
seaborn>=0.11.0
"""

with open(f"{PROJECT_ROOT}/requirements.txt", 'w') as f:
    f.write(requirements)
print("✓ Saved: requirements.txt")

# LICENSE 파일 생성
license_text = """MIT License

Copyright (c) 2025 Young Min Park

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

with open(f"{PROJECT_ROOT}/LICENSE", 'w') as f:
    f.write(license_text)
print("✓ Saved: LICENSE")
