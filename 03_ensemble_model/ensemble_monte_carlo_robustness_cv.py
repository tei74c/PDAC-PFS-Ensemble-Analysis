"""
Monte Carlo Robustness Analysis: Ensemble Model (NB3 + Cox) using Cross-Validation
====================================================================================
Robustness analysis using cross-validation for stable performance estimates.

Analyses:
1. Targeted Gene Dropout (Leave-One-Out Feature Analysis)
2. Random Gene Dropout Analysis
3. Expression Noise Injection
4. Monte Carlo Label Shuffling (Permutation Test)
5. Feature-Specific Noise Sensitivity

CV Methodology: Repeated Stratified K-Fold (5 folds x 10 repeats = 50 evaluations)

Statistical correction: Phipson-Smyth permutation p-values
  p = (r + 1) / (n + 1)
  Reference: Phipson & Smyth (2010), Stat Appl Genet Mol Biol 9:Article39

Reference:
Khoshnevis et al. Molecular Cancer (2026).
"""

import pandas as pd
import numpy as np
from sklearn.metrics import (accuracy_score, f1_score, roc_auc_score,
                             precision_score, recall_score, matthews_corrcoef,
                             brier_score_loss)
from sklearn.naive_bayes import GaussianNB
from sklearn.calibration import CalibratedClassifierCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from scipy.special import expit
import matplotlib.pyplot as plt
from scipy import stats
import warnings
import os

warnings.filterwarnings('ignore')

# Publication-quality figure settings with TrueType fonts for Illustrator
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

np.random.seed(42)

# ============================================================================
# CONFIGURATION
# ============================================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'data')
OUTPUT_DIR = os.path.join(BASE_DIR, 'results', 'robustness_cv')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Feature definitions
NB_FEATURES = ['DYNC2H1', 'ECM2', 'PPIB']
COX_FEATURES = ['TFRC', 'APOF', 'ANG', 'FABP4']
ALL_FEATURES = NB_FEATURES + COX_FEATURES

# Cox model coefficients (from multivariate Cox regression)
COX_BETAS = {
    'TFRC': 0.390057,
    'APOF': 0.391837,
    'ANG': 0.269207,
    'FABP4': 0.280084
}

# Ensemble weighting strategy: (NB + 2*Cox) / 3
ENSEMBLE_STRATEGY = 'weighted_1_2'

# Cross-validation configuration
N_SPLITS = 5
N_REPEATS = 10
TOTAL_FOLDS = N_SPLITS * N_REPEATS

# ============================================================================
# MODEL CLASSES
# ============================================================================

class NBChampionClassifier:
    """3-feature Naive Bayes classifier with calibration."""

    def __init__(self, features=None, calibrate=True):
        self.features = features or NB_FEATURES
        self.calibrate = calibrate
        self.model_ = None
        self.scaler_ = None

    def fit(self, X, y):
        self.scaler_ = StandardScaler()
        X_scaled = self.scaler_.fit_transform(X)
        base_clf = GaussianNB()
        if self.calibrate and len(np.unique(y)) > 1:
            min_class = min(np.sum(y == 0), np.sum(y == 1))
            cv_folds = min(3, min_class)
            if cv_folds >= 2:
                self.model_ = CalibratedClassifierCV(base_clf, method='sigmoid', cv=cv_folds)
            else:
                self.model_ = base_clf
        else:
            self.model_ = base_clf
        self.model_.fit(X_scaled, y)
        return self

    def predict(self, X):
        X_scaled = self.scaler_.transform(X)
        return self.model_.predict(X_scaled)

    def predict_proba(self, X):
        X_scaled = self.scaler_.transform(X)
        return self.model_.predict_proba(X_scaled)


class CoxRiskClassifier:
    """Cox model with fixed coefficients converted to a binary classifier."""

    def __init__(self, betas=None, features=None):
        self.betas = betas or COX_BETAS
        self.features = features or COX_FEATURES
        self.threshold_ = None

    def _compute_risk(self, X):
        risk = np.zeros(len(X))
        for i, feat in enumerate(self.features):
            if feat in self.betas:
                risk += self.betas[feat] * X[:, i]
        return risk

    def fit(self, X, y):
        risk_scores = self._compute_risk(X)
        self.threshold_ = np.median(risk_scores)
        return self

    def predict(self, X):
        proba = self.predict_proba(X)[:, 1]
        return (proba >= 0.5).astype(int)

    def predict_proba(self, X):
        risk_scores = self._compute_risk(X)
        proba_pos = expit(risk_scores - self.threshold_)
        return np.column_stack([1 - proba_pos, proba_pos])


class EnsembleClassifier:
    """Ensemble combining NB and Cox models with configurable weighting."""

    def __init__(self, strategy='weighted_1_2'):
        self.strategy = strategy
        self.nb_model_ = None
        self.cox_model_ = None

    def fit(self, X_nb, X_cox, y):
        self.nb_model_ = NBChampionClassifier()
        self.nb_model_.fit(X_nb, y)
        self.cox_model_ = CoxRiskClassifier()
        self.cox_model_.fit(X_cox, y)
        return self

    def predict(self, X_nb, X_cox):
        proba = self.predict_proba(X_nb, X_cox)[:, 1]
        return (proba >= 0.5).astype(int)

    def predict_proba(self, X_nb, X_cox):
        proba_nb = self.nb_model_.predict_proba(X_nb)[:, 1]
        proba_cox = self.cox_model_.predict_proba(X_cox)[:, 1]
        if self.strategy == 'weighted_1_2':
            proba_pos = (proba_nb + 2 * proba_cox) / 3
        else:
            proba_pos = 0.5 * proba_nb + 0.5 * proba_cox
        return np.column_stack([1 - proba_pos, proba_pos])


# ============================================================================
# LOAD DATA
# ============================================================================
print("=" * 70)
print("MONTE CARLO ROBUSTNESS ANALYSIS FOR ENSEMBLE MODEL (CV VERSION)")
print("=" * 70)

data_path = os.path.join(DATA_DIR, "FinalwinnerML_Cox_dataset.xlsx")
df = pd.read_excel(data_path)

# Encode labels: L (Long PFS) = 0, S (Short PFS) = 1
label_map = {'L': 0, 'S': 1}
df['y'] = df['PFS_group'].map(label_map)

X_nb = df[NB_FEATURES].values
X_cox = df[COX_FEATURES].values
y = df['y'].values

print(f"\nEnsemble Strategy: {ENSEMBLE_STRATEGY}")
print(f"NB Features: {NB_FEATURES}")
print(f"Cox Features: {COX_FEATURES}")
print(f"\nTotal samples: {len(y)} (Short PFS: {sum(y)}, Long PFS: {len(y)-sum(y)})")
print(f"CV Strategy: {N_SPLITS}-fold x {N_REPEATS} repeats = {TOTAL_FOLDS} evaluations")

# ============================================================================
# BASELINE PERFORMANCE (Cross-Validated)
# ============================================================================
print(f"\n{'='*70}")
print("BASELINE PERFORMANCE (Cross-Validated)")
print("=" * 70)

cv = RepeatedStratifiedKFold(n_splits=N_SPLITS, n_repeats=N_REPEATS, random_state=42)

cv_aucs = []
cv_accs = []
cv_f1s = []
cv_mccs = []

for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X_nb, y)):
    X_nb_train, X_nb_test = X_nb[train_idx], X_nb[test_idx]
    X_cox_train, X_cox_test = X_cox[train_idx], X_cox[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    model = EnsembleClassifier(strategy=ENSEMBLE_STRATEGY)
    model.fit(X_nb_train, X_cox_train, y_train)

    pred = model.predict(X_nb_test, X_cox_test)
    proba = model.predict_proba(X_nb_test, X_cox_test)[:, 1]

    cv_aucs.append(roc_auc_score(y_test, proba))
    cv_accs.append(accuracy_score(y_test, pred))
    cv_f1s.append(f1_score(y_test, pred))
    cv_mccs.append(matthews_corrcoef(y_test, pred))

baseline_metrics = {
    'ROC-AUC': np.mean(cv_aucs),
    'ROC-AUC_std': np.std(cv_aucs),
    'ROC-AUC_CI_low': np.percentile(cv_aucs, 2.5),
    'ROC-AUC_CI_high': np.percentile(cv_aucs, 97.5),
    'Accuracy': np.mean(cv_accs),
    'Accuracy_std': np.std(cv_accs),
    'F1-Score': np.mean(cv_f1s),
    'F1-Score_std': np.std(cv_f1s),
    'MCC': np.mean(cv_mccs),
    'MCC_std': np.std(cv_mccs)
}

print(f"\n  ROC-AUC: {baseline_metrics['ROC-AUC']:.4f} +/- {baseline_metrics['ROC-AUC_std']:.4f}")
print(f"  95% CI: [{baseline_metrics['ROC-AUC_CI_low']:.4f}, {baseline_metrics['ROC-AUC_CI_high']:.4f}]")
print(f"  Accuracy: {baseline_metrics['Accuracy']:.4f} +/- {baseline_metrics['Accuracy_std']:.4f}")
print(f"  F1-Score: {baseline_metrics['F1-Score']:.4f} +/- {baseline_metrics['F1-Score_std']:.4f}")
print(f"  MCC: {baseline_metrics['MCC']:.4f} +/- {baseline_metrics['MCC_std']:.4f}")

# ============================================================================
# 1. TARGETED GENE DROPOUT (CV)
# ============================================================================
print(f"\n{'='*70}")
print("1. TARGETED GENE DROPOUT ANALYSIS (CV)")
print("=" * 70)

dropout_results = []

# NB Feature Dropout
print("\n  NB Feature Dropout:")
for feature in NB_FEATURES:
    remaining = [f for f in NB_FEATURES if f != feature]
    remaining_idx = [NB_FEATURES.index(f) for f in remaining]

    cv_aucs_drop = []
    for train_idx, test_idx in cv.split(X_nb, y):
        X_nb_train = X_nb[train_idx][:, remaining_idx]
        X_nb_test = X_nb[test_idx][:, remaining_idx]
        X_cox_train, X_cox_test = X_cox[train_idx], X_cox[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        nb_drop = NBChampionClassifier(features=remaining)
        nb_drop.fit(X_nb_train, y_train)
        cox_model = CoxRiskClassifier()
        cox_model.fit(X_cox_train, y_train)

        proba_nb = nb_drop.predict_proba(X_nb_test)[:, 1]
        proba_cox = cox_model.predict_proba(X_cox_test)[:, 1]
        proba_ens = (proba_nb + 2 * proba_cox) / 3

        cv_aucs_drop.append(roc_auc_score(y_test, proba_ens))

    mean_auc = np.mean(cv_aucs_drop)
    auc_drop = baseline_metrics['ROC-AUC'] - mean_auc

    dropout_results.append({
        'Dropped_Feature': feature, 'Model': 'NB',
        'AUC_Without': mean_auc, 'AUC_Without_std': np.std(cv_aucs_drop),
        'AUC_Drop': auc_drop
    })
    print(f"    Without {feature}: AUC = {mean_auc:.4f} +/- {np.std(cv_aucs_drop):.4f} (drop: {auc_drop:+.4f})")

# Cox Feature Dropout
print("\n  Cox Feature Dropout:")
for feature in COX_FEATURES:
    remaining = [f for f in COX_FEATURES if f != feature]
    remaining_betas = {k: v for k, v in COX_BETAS.items() if k != feature}
    remaining_idx = [COX_FEATURES.index(f) for f in remaining]

    cv_aucs_drop = []
    for train_idx, test_idx in cv.split(X_nb, y):
        X_nb_train, X_nb_test = X_nb[train_idx], X_nb[test_idx]
        X_cox_train = X_cox[train_idx][:, remaining_idx]
        X_cox_test = X_cox[test_idx][:, remaining_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        nb_model = NBChampionClassifier()
        nb_model.fit(X_nb_train, y_train)
        cox_drop = CoxRiskClassifier(betas=remaining_betas, features=remaining)
        cox_drop.fit(X_cox_train, y_train)

        proba_nb = nb_model.predict_proba(X_nb_test)[:, 1]
        proba_cox = cox_drop.predict_proba(X_cox_test)[:, 1]
        proba_ens = (proba_nb + 2 * proba_cox) / 3

        cv_aucs_drop.append(roc_auc_score(y_test, proba_ens))

    mean_auc = np.mean(cv_aucs_drop)
    auc_drop = baseline_metrics['ROC-AUC'] - mean_auc

    dropout_results.append({
        'Dropped_Feature': feature, 'Model': 'Cox',
        'AUC_Without': mean_auc, 'AUC_Without_std': np.std(cv_aucs_drop),
        'AUC_Drop': auc_drop
    })
    print(f"    Without {feature}: AUC = {mean_auc:.4f} +/- {np.std(cv_aucs_drop):.4f} (drop: {auc_drop:+.4f})")

targeted_df = pd.DataFrame(dropout_results)
ranked = targeted_df.sort_values('AUC_Drop', ascending=False)

print(f"\n  Feature Importance Ranking:")
for _, row in ranked.iterrows():
    print(f"    {row['Dropped_Feature']} ({row['Model']}): {row['AUC_Drop']:+.4f}")

# ============================================================================
# 2. RANDOM GENE DROPOUT (CV)
# ============================================================================
print(f"\n{'='*70}")
print("2. RANDOM GENE DROPOUT ANALYSIS (CV)")
print("=" * 70)

n_iterations = 100
random_dropout_results = {1: [], 2: [], 3: []}

for n_drop in [1, 2, 3]:
    print(f"\n  Dropping {n_drop} feature(s) randomly ({n_iterations} iterations)...")

    for _ in range(n_iterations):
        drop_features = np.random.choice(ALL_FEATURES, n_drop, replace=False)
        nb_remain = [f for f in NB_FEATURES if f not in drop_features]
        cox_remain = [f for f in COX_FEATURES if f not in drop_features]

        if len(nb_remain) == 0 or len(cox_remain) == 0:
            continue

        nb_idx = [NB_FEATURES.index(f) for f in nb_remain]
        cox_idx = [COX_FEATURES.index(f) for f in cox_remain]
        cox_betas_remain = {k: v for k, v in COX_BETAS.items() if k in cox_remain}

        fold_aucs = []
        for train_idx, test_idx in StratifiedKFold(n_splits=5, shuffle=True, random_state=None).split(X_nb, y):
            try:
                X_nb_train = X_nb[train_idx][:, nb_idx]
                X_nb_test = X_nb[test_idx][:, nb_idx]
                X_cox_train = X_cox[train_idx][:, cox_idx]
                X_cox_test = X_cox[test_idx][:, cox_idx]
                y_train, y_test = y[train_idx], y[test_idx]

                nb_drop = NBChampionClassifier(features=nb_remain)
                nb_drop.fit(X_nb_train, y_train)
                cox_drop = CoxRiskClassifier(betas=cox_betas_remain, features=cox_remain)
                cox_drop.fit(X_cox_train, y_train)

                proba_nb = nb_drop.predict_proba(X_nb_test)[:, 1]
                proba_cox = cox_drop.predict_proba(X_cox_test)[:, 1]
                proba_ens = (proba_nb + 2 * proba_cox) / 3

                fold_aucs.append(roc_auc_score(y_test, proba_ens))
            except:
                continue

        if fold_aucs:
            random_dropout_results[n_drop].append({'auc': np.mean(fold_aucs)})

    if random_dropout_results[n_drop]:
        aucs = [r['auc'] for r in random_dropout_results[n_drop]]
        print(f"    AUC: {np.mean(aucs):.4f} +/- {np.std(aucs):.4f}")

# ============================================================================
# 3. NOISE INJECTION (CV)
# ============================================================================
print(f"\n{'='*70}")
print("3. EXPRESSION NOISE INJECTION ANALYSIS (CV)")
print("=" * 70)

noise_levels = [0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.50]
n_noise_iter = 50
noise_results = {level: {'auc': []} for level in noise_levels}

# Feature standard deviations from full dataset
nb_stds = np.std(X_nb, axis=0)
cox_stds = np.std(X_cox, axis=0)

for noise_level in noise_levels:
    print(f"\n  Noise level: {int(noise_level*100)}%...")

    for _ in range(n_noise_iter):
        fold_aucs = []
        for train_idx, test_idx in StratifiedKFold(n_splits=5, shuffle=True).split(X_nb, y):
            X_nb_train, X_nb_test = X_nb[train_idx], X_nb[test_idx].copy()
            X_cox_train, X_cox_test = X_cox[train_idx], X_cox[test_idx].copy()
            y_train, y_test = y[train_idx], y[test_idx]

            # Add noise to test set only
            X_nb_test += np.random.normal(0, nb_stds * noise_level, X_nb_test.shape)
            X_cox_test += np.random.normal(0, cox_stds * noise_level, X_cox_test.shape)

            model = EnsembleClassifier(strategy=ENSEMBLE_STRATEGY)
            model.fit(X_nb_train, X_cox_train, y_train)
            proba = model.predict_proba(X_nb_test, X_cox_test)[:, 1]

            try:
                fold_aucs.append(roc_auc_score(y_test, proba))
            except:
                continue

        if fold_aucs:
            noise_results[noise_level]['auc'].append(np.mean(fold_aucs))

    mean_auc = np.mean(noise_results[noise_level]['auc'])
    std_auc = np.std(noise_results[noise_level]['auc'])
    degradation = (baseline_metrics['ROC-AUC'] - mean_auc) / baseline_metrics['ROC-AUC'] * 100
    print(f"    AUC: {mean_auc:.4f} +/- {std_auc:.4f} ({degradation:.1f}% degradation)")

# ============================================================================
# 4. MONTE CARLO LABEL SHUFFLING (Permutation Test with CV)
# ============================================================================
print(f"\n{'='*70}")
print("4. MONTE CARLO LABEL SHUFFLING (Permutation Test - CV)")
print("=" * 70)
print("Shuffling TRAINING labels only within each CV fold")
print("Evaluating against TRUE test fold labels")

n_permutations = 500
permutation_aucs = []
permutation_accs = []
permutation_f1s = []

print(f"\n  Running {n_permutations} permutations...")
for perm_idx in range(n_permutations):
    if (perm_idx + 1) % 100 == 0:
        print(f"    Completed {perm_idx+1}/{n_permutations} permutations...")

    fold_aucs = []
    fold_accs = []
    fold_f1s = []

    for train_idx, test_idx in StratifiedKFold(n_splits=5, shuffle=True).split(X_nb, y):
        X_nb_train, X_nb_test = X_nb[train_idx], X_nb[test_idx]
        X_cox_train, X_cox_test = X_cox[train_idx], X_cox[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Shuffle ONLY training labels (correct methodology)
        y_train_shuffled = np.random.permutation(y_train)

        try:
            model = EnsembleClassifier(strategy=ENSEMBLE_STRATEGY)
            model.fit(X_nb_train, X_cox_train, y_train_shuffled)

            pred = model.predict(X_nb_test, X_cox_test)
            proba = model.predict_proba(X_nb_test, X_cox_test)[:, 1]

            # Evaluate against TRUE test labels (NOT shuffled)
            fold_aucs.append(roc_auc_score(y_test, proba))
            fold_accs.append(accuracy_score(y_test, pred))
            fold_f1s.append(f1_score(y_test, pred, zero_division=0))
        except:
            continue

    if fold_aucs:
        permutation_aucs.append(np.mean(fold_aucs))
        permutation_accs.append(np.mean(fold_accs))
        permutation_f1s.append(np.mean(fold_f1s))

# Phipson-Smyth corrected empirical p-values: p = (r + 1) / (n + 1)
# Reference: Phipson & Smyth (2010), Stat Appl Genet Mol Biol 9:Article39
# This correction avoids zero p-values and provides valid confidence bounds
n_perm = len(permutation_aucs)
empirical_p_auc = (np.sum(np.array(permutation_aucs) >= baseline_metrics['ROC-AUC']) + 1) / (n_perm + 1)
empirical_p_acc = (np.sum(np.array(permutation_accs) >= baseline_metrics['Accuracy']) + 1) / (n_perm + 1)
empirical_p_f1 = (np.sum(np.array(permutation_f1s) >= baseline_metrics['F1-Score']) + 1) / (n_perm + 1)

cohens_d_auc = (baseline_metrics['ROC-AUC'] - np.mean(permutation_aucs)) / np.std(permutation_aucs)

print(f"\n  Null Distribution Statistics:")
print(f"    AUC:  mean={np.mean(permutation_aucs):.4f} +/- {np.std(permutation_aucs):.4f}")
print(f"    Accuracy: mean={np.mean(permutation_accs):.4f} +/- {np.std(permutation_accs):.4f}")
print(f"    F1-Score: mean={np.mean(permutation_f1s):.4f} +/- {np.std(permutation_f1s):.4f}")

print(f"\n  Observed vs Null:")
print(f"    AUC:  Observed={baseline_metrics['ROC-AUC']:.4f}, Null={np.mean(permutation_aucs):.4f}, p={empirical_p_auc:.4f}")
print(f"    Effect Size (Cohen's d): {cohens_d_auc:.2f}")
print(f"    Interpretation: {'SIGNIFICANT (p < 0.05)' if empirical_p_auc < 0.05 else 'Not significant'}")

# ============================================================================
# 5. FEATURE-SPECIFIC NOISE SENSITIVITY (CV)
# ============================================================================
print(f"\n{'='*70}")
print("5. FEATURE-SPECIFIC NOISE SENSITIVITY (CV)")
print("=" * 70)

feature_noise_sensitivity = {}

for feat_type, features, X_data, stds in [('NB', NB_FEATURES, X_nb, nb_stds), ('Cox', COX_FEATURES, X_cox, cox_stds)]:
    print(f"\n  {feat_type} Features:")

    for feat_idx, feature in enumerate(features):
        sensitivities = []

        for noise_level in [0.1, 0.2, 0.3]:
            fold_aucs = []

            for _ in range(20):
                for train_idx, test_idx in StratifiedKFold(n_splits=5, shuffle=True).split(X_nb, y):
                    X_nb_train, X_nb_test = X_nb[train_idx], X_nb[test_idx].copy()
                    X_cox_train, X_cox_test = X_cox[train_idx], X_cox[test_idx].copy()
                    y_train, y_test = y[train_idx], y[test_idx]

                    # Add noise to single feature only
                    if feat_type == 'NB':
                        X_nb_test[:, feat_idx] += np.random.normal(0, stds[feat_idx] * noise_level, len(test_idx))
                    else:
                        X_cox_test[:, feat_idx] += np.random.normal(0, stds[feat_idx] * noise_level, len(test_idx))

                    model = EnsembleClassifier(strategy=ENSEMBLE_STRATEGY)
                    model.fit(X_nb_train, X_cox_train, y_train)
                    proba = model.predict_proba(X_nb_test, X_cox_test)[:, 1]

                    try:
                        fold_aucs.append(roc_auc_score(y_test, proba))
                    except:
                        continue

            sensitivities.append(baseline_metrics['ROC-AUC'] - np.mean(fold_aucs))

        feature_noise_sensitivity[feature] = {
            'model': feat_type,
            'mean_degradation': np.mean(sensitivities),
            'sensitivities': sensitivities
        }
        print(f"    {feature}: Mean AUC degradation = {np.mean(sensitivities):.4f}")

# ============================================================================
# 6. COMPREHENSIVE VISUALIZATION
# ============================================================================
print(f"\n{'='*70}")
print("6. GENERATING COMPREHENSIVE VISUALIZATIONS")
print("=" * 70)

fig = plt.figure(figsize=(16, 12))

# Panel A: Feature Importance
ax1 = fig.add_subplot(2, 3, 1)
features_list = targeted_df['Dropped_Feature'].tolist()
auc_drops = targeted_df['AUC_Drop'].tolist()
models = targeted_df['Model'].tolist()
colors = ['#d62728' if d > 0.02 else '#2ca02c' if d < 0.01 else '#ff7f0e' for d in auc_drops]

bars = ax1.bar(range(len(features_list)), auc_drops, color=colors, edgecolor='black', linewidth=0.5)
ax1.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
ax1.set_xticks(range(len(features_list)))
ax1.set_xticklabels([f"{f}\n({m})" for f, m in zip(features_list, models)], rotation=45, ha='right', fontsize=8)
ax1.set_ylabel('AUC Drop')
ax1.set_title('A. Feature Importance (CV)')

# Panel B: Random Dropout
ax2 = fig.add_subplot(2, 3, 2)
dropout_data = [[r['auc'] for r in random_dropout_results[n]] for n in [1, 2, 3]]
bp = ax2.boxplot(dropout_data, patch_artist=True, labels=['1', '2', '3'])
for patch, color in zip(bp['boxes'], ['#1f77b4', '#ff7f0e', '#d62728']):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax2.axhline(y=baseline_metrics['ROC-AUC'], color='green', linestyle='--', linewidth=2,
            label=f'Baseline ({baseline_metrics["ROC-AUC"]:.3f})')
ax2.axhline(y=0.5, color='red', linestyle=':', linewidth=1.5)
ax2.set_ylabel('ROC-AUC')
ax2.set_xlabel('Features Dropped')
ax2.set_title('B. Random Dropout (CV)')
ax2.legend(loc='lower left', fontsize=8)

# Panel C: Noise Injection
ax3 = fig.add_subplot(2, 3, 3)
noise_pct = [int(n * 100) for n in noise_levels]
noise_means = [np.mean(noise_results[n]['auc']) for n in noise_levels]
noise_stds = [np.std(noise_results[n]['auc']) for n in noise_levels]

ax3.errorbar(noise_pct, noise_means, yerr=noise_stds, fmt='o-', color='#1f77b4',
             capsize=4, markersize=8, linewidth=2)
ax3.fill_between(noise_pct, [m-s for m, s in zip(noise_means, noise_stds)],
                 [m+s for m, s in zip(noise_means, noise_stds)], alpha=0.2)
ax3.axhline(y=baseline_metrics['ROC-AUC'], color='green', linestyle='--', linewidth=2)
ax3.axhline(y=0.5, color='red', linestyle=':', linewidth=1.5)
ax3.set_xlabel('Noise Level (%)')
ax3.set_ylabel('ROC-AUC')
ax3.set_title('C. Noise Robustness (CV)')

# Panel D: Permutation Test
ax4 = fig.add_subplot(2, 3, 4)
ax4.hist(permutation_aucs, bins=30, color='lightblue', edgecolor='black', alpha=0.7)
ax4.axvline(x=baseline_metrics['ROC-AUC'], color='red', linestyle='-', linewidth=2.5,
            label=f'Observed ({baseline_metrics["ROC-AUC"]:.3f})')
ax4.axvline(x=np.mean(permutation_aucs), color='blue', linestyle='--', linewidth=1.5,
            label=f'Null ({np.mean(permutation_aucs):.3f})')

# Add Phipson-Smyth corrected p-value annotation
p_text = f'p = {empirical_p_auc:.3f}' if empirical_p_auc >= 0.001 else 'p < 0.001'
ax4.text(0.95, 0.95, f'{p_text}\nd = {cohens_d_auc:.2f}', transform=ax4.transAxes,
         fontsize=10, va='top', ha='right', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
ax4.set_xlabel('ROC-AUC')
ax4.set_ylabel('Frequency')
ax4.set_title('D. Permutation Test (CV)')
ax4.legend(loc='upper left', fontsize=8)

# Panel E: Summary Table
ax5 = fig.add_subplot(2, 3, 5)
ax5.axis('off')

summary_data = [
    ['Metric', 'Observed', 'Null Mean', 'p-value', "Cohen's d"],
    ['ROC-AUC', f'{baseline_metrics["ROC-AUC"]:.4f}', f'{np.mean(permutation_aucs):.4f}',
     f'{empirical_p_auc:.4f}', f'{cohens_d_auc:.2f}'],
    ['Accuracy', f'{baseline_metrics["Accuracy"]:.4f}', f'{np.mean(permutation_accs):.4f}',
     f'{empirical_p_acc:.4f}', f'{(baseline_metrics["Accuracy"]-np.mean(permutation_accs))/np.std(permutation_accs):.2f}'],
    ['F1-Score', f'{baseline_metrics["F1-Score"]:.4f}', f'{np.mean(permutation_f1s):.4f}',
     f'{empirical_p_f1:.4f}', f'{(baseline_metrics["F1-Score"]-np.mean(permutation_f1s))/np.std(permutation_f1s):.2f}']
]

table = ax5.table(cellText=summary_data, loc='center', cellLoc='center',
                  colWidths=[0.2, 0.18, 0.18, 0.15, 0.15])
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.8)
for i in range(5):
    table[(0, i)].set_facecolor('#4472C4')
    table[(0, i)].set_text_props(color='white', weight='bold')
ax5.set_title('E. Statistical Summary (CV)', pad=20)

# Panel F: Radar Chart
ax6 = fig.add_subplot(2, 3, 6, polar=True)
categories = ['Feature\nImportance', 'Noise\nTolerance', 'Statistical\nSignificance', 'Model\nStability']
N = len(categories)

# Calculate noise tolerance (highest noise level where AUC >= 0.75)
noise_tol = 0
for level in noise_levels:
    if np.mean(noise_results[level]['auc']) >= 0.75:
        noise_tol = level
noise_tol = noise_tol / max(noise_levels)

# Stability from random dropout
stability = np.mean([np.mean([r['auc'] for r in random_dropout_results[1]]),
                     np.mean([r['auc'] for r in random_dropout_results[2]])]) / baseline_metrics['ROC-AUC']

values = [np.mean(auc_drops) / baseline_metrics['ROC-AUC'],
          noise_tol, 1 - empirical_p_auc, stability]
values += values[:1]

angles = [n / float(N) * 2 * np.pi for n in range(N)]
angles += angles[:1]

ax6.plot(angles, values, 'o-', linewidth=2, color='#2ca02c')
ax6.fill(angles, values, alpha=0.25, color='#2ca02c')
ax6.set_xticks(angles[:-1])
ax6.set_xticklabels(categories, size=9)
ax6.set_ylim(0, 1)
ax6.set_title('F. Robustness Profile', pad=15)

plt.suptitle(f'Monte Carlo Robustness Analysis (CV)\nEnsemble Model ({ENSEMBLE_STRATEGY})',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout(rect=[0, 0, 1, 0.96])

for fmt in ['png', 'pdf', 'svg']:
    plt.savefig(os.path.join(OUTPUT_DIR, f'Ensemble_Monte_Carlo_Robustness_CV.{fmt}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("\n  Saved: Ensemble_Monte_Carlo_Robustness_CV.png/pdf/svg")

# ============================================================================
# 7. SAVE RESULTS TO EXCEL
# ============================================================================
print(f"\n{'='*70}")
print("7. SAVING RESULTS TO EXCEL")
print("=" * 70)

with pd.ExcelWriter(os.path.join(OUTPUT_DIR, 'Ensemble_Robustness_Results_CV.xlsx'), engine='openpyxl') as writer:
    pd.DataFrame([baseline_metrics]).to_excel(writer, sheet_name='Baseline_CV', index=False)
    targeted_df.to_excel(writer, sheet_name='Feature_Dropout', index=False)

    random_summary = [{'N_Dropped': n, 'Mean_AUC': np.mean([r['auc'] for r in random_dropout_results[n]]),
                       'Std_AUC': np.std([r['auc'] for r in random_dropout_results[n]])}
                      for n in [1, 2, 3]]
    pd.DataFrame(random_summary).to_excel(writer, sheet_name='Random_Dropout', index=False)

    noise_summary = [{'Noise_Pct': int(n*100), 'Mean_AUC': np.mean(noise_results[n]['auc']),
                      'Std_AUC': np.std(noise_results[n]['auc'])} for n in noise_levels]
    pd.DataFrame(noise_summary).to_excel(writer, sheet_name='Noise_Injection', index=False)

    # Permutation summary with Phipson-Smyth corrected p-values
    perm_summary = pd.DataFrame([{
        'Observed_AUC': baseline_metrics['ROC-AUC'],
        'Null_Mean': np.mean(permutation_aucs), 'Null_Std': np.std(permutation_aucs),
        'p_value': empirical_p_auc, 'Cohens_d': cohens_d_auc,
        'N_Permutations': len(permutation_aucs), 'CV_Folds': 5
    }])
    perm_summary.to_excel(writer, sheet_name='Permutation_Test', index=False)

print("  Saved: Ensemble_Robustness_Results_CV.xlsx")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print(f"\n{'='*70}")
print("FINAL SUMMARY: ENSEMBLE MODEL ROBUSTNESS (CV)")
print("=" * 70)

print(f"""
MODEL: Ensemble ({ENSEMBLE_STRATEGY}) - CROSS-VALIDATED
CV Strategy: {N_SPLITS}-fold x {N_REPEATS} repeats = {TOTAL_FOLDS} evaluations

BASELINE PERFORMANCE (CV Mean +/- Std):
   - ROC-AUC: {baseline_metrics['ROC-AUC']:.4f} +/- {baseline_metrics['ROC-AUC_std']:.4f}
   - 95% CI: [{baseline_metrics['ROC-AUC_CI_low']:.4f}, {baseline_metrics['ROC-AUC_CI_high']:.4f}]
   - Accuracy: {baseline_metrics['Accuracy']:.4f} +/- {baseline_metrics['Accuracy_std']:.4f}

STATISTICAL SIGNIFICANCE (Phipson-Smyth Corrected Permutation Test):
   - p-value: {empirical_p_auc:.4f} ({'SIGNIFICANT' if empirical_p_auc < 0.05 else 'Not significant'})
   - Cohen's d: {cohens_d_auc:.2f} ({'Large' if cohens_d_auc > 0.8 else 'Medium' if cohens_d_auc > 0.5 else 'Small'} effect)

FEATURE IMPORTANCE (CV):
   - Most important: {ranked.iloc[0]['Dropped_Feature']} ({ranked.iloc[0]['Model']}) = {ranked.iloc[0]['AUC_Drop']:+.4f}
   - Least important: {ranked.iloc[-1]['Dropped_Feature']} ({ranked.iloc[-1]['Model']}) = {ranked.iloc[-1]['AUC_Drop']:+.4f}

NOISE ROBUSTNESS (CV):
   - 10% noise: AUC = {np.mean(noise_results[0.10]['auc']):.4f}
   - 30% noise: AUC = {np.mean(noise_results[0.30]['auc']):.4f}

CONCLUSIONS:
   1. CV provides more stable and reliable estimates
   2. Model {'SIGNIFICANTLY' if empirical_p_auc < 0.05 else 'does not significantly'} outperform chance
   3. Effect size is {'LARGE' if cohens_d_auc > 0.8 else 'MODERATE' if cohens_d_auc > 0.5 else 'SMALL'}
""")

print("=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
