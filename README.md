# Seven-Protein Plasma Proteomic Ensemble for PDAC PFS Prediction

Code accompanying:

**"A seven-protein plasma proteomic ensemble predicts progression-free survival in therapy-naive stage IV PDAC"**

Khoshnevis N, Becher H, Kestler HA, Lahusen T, Daiss N, Melzer MK, Seufferlein T, Eiseler T.
*Molecular Cancer* (2026).

## Overview

This repository contains the analysis code for building and validating a 7-protein ensemble classifier combining:

- **ML classifier**: 3-feature Naive Bayes
- **Cox model**: 4-feature Cox PH
- **Ensemble**: Cox-weighted integration P = (P_NB + 2 x P_Cox) / 3

## Repository Structure

```
00_cox_regression/                # Cox PH regression pipeline (Fig 3A-F, S3)
  cox_regression_pipeline.ipynb       # Univariate/multivariate Cox, QC, risk stratification (designed for Google Colab)

01_nested_cv_robustness/          # Monte Carlo robustness & noise injection (Fig 2D-N, S6)
  monte_carlo_robustness_4feature.py    # 4-feature NB robustness analysis
  monte_carlo_robustness_3feature.py    # 3-feature NB robustness analysis
  model_comparison_3vs4_features.py     # Head-to-head 3 vs 4 feature comparison
  advanced_model_analysis.py            # SHAP, DCA, ROC with confidence intervals
  generate_comparison_figures.py        # Individual metric comparison panels

02_cox_vs_ml_comparison/          # Cox vs ML head-to-head (Fig 4A-G)
  cox_vs_ml_comparison.py

03_ensemble_model/                # Ensemble integration & robustness (Fig 4H-R, S7)
  ensemble_nb3_cox_comparison.py        # 7-strategy ensemble evaluation
  ensemble_monte_carlo_robustness.py    # Ensemble noise sensitivity (train/test)
  ensemble_monte_carlo_robustness_cv.py # Ensemble noise sensitivity (cross-validation)

04_transferability/               # Endpoint transferability (Fig S4E-I)
  create_transferability_figures.py

05_nested_cv_pipeline/            # Nested CV feature selection & model training (Fig 2A-C, S5)
  nested_cv_pipeline.ipynb            # Full pipeline (USE_ELASTIC_NET tunable for mode selection)

06_wgcna/                         # WGCNA co-expression network analysis (Fig 5, S8-S10)
  wgcna_Nika_final.R                  # Full pipeline: network construction, module detection,
                                      #   hub genes, pathway enrichment, clinical associations (R 4.x)
```

## Data Availability

Input data files are not included in this repository and are available upon reasonable request from the corresponding author. Place data files in a `data/` subdirectory relative to each script:

| File | Description |
|------|-------------|
| `FinalwinnerML_Cox_dataset.xlsx` | Combined ML + Cox feature matrix (50 patients x 891 proteins) |
| `champion_model_NB_4features.pkl` | Trained 4-feature NB classifier |
| `champion_model_NB_3features.pkl` | Trained 3-feature NB classifier |
| `PFS_ML_train_cof_wo_IG.xlsx` | Training set (n=34) |
| `PFS_ML_test_cof_wo_IG.xlsx` | Holdout test set (n=16) |
| `COMPREHENSIVE_RESULTS.xlsx` | Transferability analysis results |
| `correctedHLA_Rawdata_PFSfiltered.csv` | Expression + clinical data for WGCNA (891 proteins x 50 samples) |
| `mastertable_Signature.csv` | Gene signature database (150 PDAC-specific signatures) |

## Requirements

```
pip install -r requirements.txt
```

See `requirements.txt` for the full list. Key dependencies: Python >= 3.9, scikit-learn >= 1.3, lifelines >= 0.27, shap >= 0.42.

## Execution Order

1. **`05_nested_cv_pipeline/nested_cv_pipeline.ipynb`** — generates `champion_model_NB_4features.pkl`
2. **`00_cox_regression/cox_regression_pipeline.ipynb`** — Cox PH pipeline, generates risk groups
3. **`01_nested_cv_robustness/model_comparison_3vs4_features.py`** — generates `champion_model_NB_3features.pkl` (must run before other robustness scripts)
4. **`01_nested_cv_robustness/`** remaining scripts (any order):
   - `monte_carlo_robustness_4feature.py` — requires 4-feature .pkl
   - `monte_carlo_robustness_3feature.py` — retrains 3-feature NB internally
   - `advanced_model_analysis.py` — requires both .pkl files
   - `generate_comparison_figures.py` — requires both .pkl files
5. **`02_cox_vs_ml_comparison/cox_vs_ml_comparison.py`** — requires trained model files
6. **`03_ensemble_model/ensemble_nb3_cox_comparison.py`** first, then robustness scripts
7. **`04_transferability/create_transferability_figures.py`** — independent
8. **`06_wgcna/wgcna_Nika_final.R`** — independent (requires R 4.x with WGCNA, tidyverse, GSVA)

## Figure Mapping

| Figure | Script | Analysis |
|--------|--------|----------|
| Fig 2A-C | `nested_cv_pipeline.ipynb` | Nested CV, champion selection, calibration |
| Fig 3A-F | `cox_regression_pipeline.ipynb` | Univariate/multivariate Cox, QC, risk stratification |
| Fig 2D-E | `advanced_model_analysis.py` | SHAP importance, Decision Curve Analysis |
| Fig 2F-G | `monte_carlo_robustness_4feature.py` | Feature noise injection, dropout |
| Fig 2H-L | `model_comparison_3vs4_features.py` | 3 vs 4 feature comparison |
| Fig 2M | `monte_carlo_robustness_4feature.py` + `monte_carlo_robustness_3feature.py` | Permutation test overlay |
| Fig 2N | `monte_carlo_robustness_3feature.py` | 3-feature ablation |
| Fig 4A-G | `cox_vs_ml_comparison.py` | Cox vs ML bootstrap CV comparison |
| Fig 4H-R | `ensemble_nb3_cox_comparison.py` | Ensemble strategy evaluation |
| Fig 4P-Q | `ensemble_monte_carlo_robustness_cv.py` | Ensemble noise + permutation test |
| Fig S4E-I | `create_transferability_figures.py` | Endpoint transferability |
| Fig S6A-E | `monte_carlo_robustness_4feature.py` + `monte_carlo_robustness_3feature.py` | Extended robustness |
| Fig S7G-H | `ensemble_monte_carlo_robustness.py` | Ensemble noise sensitivity |

## License

This code is provided for academic and research purposes accompanying the manuscript.

## Contact

For questions about the code or data access, please contact the corresponding author.
