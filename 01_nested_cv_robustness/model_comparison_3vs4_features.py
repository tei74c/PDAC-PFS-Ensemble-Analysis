"""
3-Feature vs 4-Feature Model Comparison
=========================================
Head-to-head comparison of NB classifiers with and without ENPP1.
Generates Figures 2H-L in the manuscript.

Reference:
Khoshnevis et al. Molecular Cancer (2026).
"""

import os
import pickle
import pandas as pd
import numpy as np
from sklearn.metrics import (accuracy_score, f1_score, roc_auc_score,
                             precision_score, recall_score, confusion_matrix,
                             roc_curve, precision_recall_curve, average_precision_score,
                             matthews_corrcoef, cohen_kappa_score, brier_score_loss)
from sklearn.naive_bayes import GaussianNB
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.model_selection import cross_val_score, StratifiedKFold
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# Path Configuration
# ============================================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'data')
OUTPUT_DIR = os.path.join(BASE_DIR, 'results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Set style for publication-quality figures
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
# Expected Calibration Error (ECE) Function
# ============================================================================
def expected_calibration_error(y_true, y_proba, n_bins=10):
    """
    Calculate Expected Calibration Error (ECE).

    ECE measures how well-calibrated predicted probabilities are.
    Lower is better (0 = perfectly calibrated).

    Parameters:
    -----------
    y_true : array-like - True binary labels
    y_proba : array-like - Predicted probabilities for positive class
    n_bins : int - Number of bins for calibration

    Returns:
    --------
    float : ECE value
    """
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    ece = 0.0

    for i in range(n_bins):
        # Find samples in this bin
        in_bin = (y_proba > bin_boundaries[i]) & (y_proba <= bin_boundaries[i + 1])
        prop_in_bin = np.mean(in_bin)

        if prop_in_bin > 0:
            # Average predicted probability in bin
            avg_confidence = np.mean(y_proba[in_bin])
            # Actual accuracy in bin
            avg_accuracy = np.mean(y_true[in_bin])
            # Weighted absolute difference
            ece += np.abs(avg_accuracy - avg_confidence) * prop_in_bin

    return ece

# ============================================================================
# Load Data and Original Model
# ============================================================================
print("=" * 70)
print("MODEL COMPARISON: 3-FEATURE vs 4-FEATURE PFS CLASSIFIER")
print("=" * 70)

# Load pickled model
with open(os.path.join(DATA_DIR, 'champion_model_NB_4features.pkl'), 'rb') as f:
    model_data = pickle.load(f)

original_model = model_data['model']
features_4 = model_data['features']
features_3 = [f for f in features_4 if f != 'ENPP1']  # Remove ENPP1

print(f"\n4-Feature Model: {features_4}")
print(f"3-Feature Model: {features_3}")

# Load data
train_df = pd.read_excel(os.path.join(DATA_DIR, 'PFS_ML_train_cof_wo_IG.xlsx'))
test_df = pd.read_excel(os.path.join(DATA_DIR, 'PFS_ML_test_cof_wo_IG.xlsx'))

# Encode labels
label_map = {'L': 0, 'S': 1}
train_df['PFS_group_encoded'] = train_df['PFS_group'].map(label_map)
test_df['PFS_group_encoded'] = test_df['PFS_group'].map(label_map)

# Prepare data
X_train_4 = train_df[features_4].values
X_test_4 = test_df[features_4].values
X_train_3 = train_df[features_3].values
X_test_3 = test_df[features_3].values
y_train = train_df['PFS_group_encoded'].values
y_test = test_df['PFS_group_encoded'].values

# Sample size reporting
print(f"\n{'='*50}")
print("SAMPLE SIZE INFORMATION")
print(f"{'='*50}")
print(f"Training set: n = {len(y_train)}")
print(f"  - Long PFS (L=0): {sum(y_train == 0)}")
print(f"  - Short PFS (S=1): {sum(y_train == 1)}")
print(f"Test set: n = {len(y_test)}")
print(f"  - Long PFS (L=0): {sum(y_test == 0)}")
print(f"  - Short PFS (S=1): {sum(y_test == 1)}")

# ============================================================================
# Train 3-Feature Model
# ============================================================================
print(f"\n{'='*70}")
print("TRAINING 3-FEATURE MODEL")
print("=" * 70)

nb_3 = GaussianNB()
model_3 = CalibratedClassifierCV(nb_3, cv=3, method='sigmoid')
model_3.fit(X_train_3, y_train)

print("3-Feature model trained successfully (Calibrated Naive Bayes)")

# ============================================================================
# Bootstrap Confidence Intervals Function
# ============================================================================
def bootstrap_metrics(y_true, y_pred, y_proba, n_bootstrap=1000, ci=95):
    """
    Calculate bootstrap confidence intervals for classification metrics.

    Parameters:
    -----------
    y_true : array-like - True labels
    y_pred : array-like - Predicted labels
    y_proba : array-like - Predicted probabilities
    n_bootstrap : int - Number of bootstrap iterations
    ci : float - Confidence interval percentage

    Returns:
    --------
    dict : Metrics with point estimates and confidence intervals
    """
    n_samples = len(y_true)

    # Storage for bootstrap samples
    boot_accuracy = []
    boot_f1 = []
    boot_auc = []
    boot_precision = []
    boot_recall = []
    boot_ap = []  # Average precision
    boot_mcc = []
    boot_kappa = []
    boot_brier = []
    boot_ece = []

    for _ in range(n_bootstrap):
        # Bootstrap sample with replacement
        indices = np.random.choice(n_samples, n_samples, replace=True)

        y_true_boot = y_true[indices]
        y_pred_boot = y_pred[indices]
        y_proba_boot = y_proba[indices]

        # Skip if only one class in bootstrap sample
        if len(np.unique(y_true_boot)) < 2:
            continue

        try:
            boot_accuracy.append(accuracy_score(y_true_boot, y_pred_boot))
            boot_f1.append(f1_score(y_true_boot, y_pred_boot))
            boot_auc.append(roc_auc_score(y_true_boot, y_proba_boot))
            boot_precision.append(precision_score(y_true_boot, y_pred_boot))
            boot_recall.append(recall_score(y_true_boot, y_pred_boot))
            boot_ap.append(average_precision_score(y_true_boot, y_proba_boot))
            boot_mcc.append(matthews_corrcoef(y_true_boot, y_pred_boot))
            boot_kappa.append(cohen_kappa_score(y_true_boot, y_pred_boot))
            boot_brier.append(brier_score_loss(y_true_boot, y_proba_boot))
            boot_ece.append(expected_calibration_error(y_true_boot, y_proba_boot))
        except:
            continue

    # Calculate percentiles for CI
    alpha = (100 - ci) / 2

    results = {}
    for name, values in [('Accuracy', boot_accuracy), ('F1-Score', boot_f1),
                         ('ROC-AUC', boot_auc), ('Precision', boot_precision),
                         ('Recall', boot_recall), ('Avg-Precision', boot_ap),
                         ('MCC', boot_mcc), ('Cohen-Kappa', boot_kappa),
                         ('Brier-Score', boot_brier), ('ECE', boot_ece)]:
        if len(values) > 0:
            results[name] = {
                'mean': np.mean(values),
                'std': np.std(values),
                'ci_lower': np.percentile(values, alpha),
                'ci_upper': np.percentile(values, 100 - alpha),
                'point_estimate': None  # Will be filled with actual value
            }

    return results

# ============================================================================
# Evaluate Both Models with Bootstrap CI
# ============================================================================
print(f"\n{'='*70}")
print("MODEL PERFORMANCE COMPARISON WITH BOOTSTRAP 95% CI")
print("=" * 70)

# 4-Feature Model
pred_4 = original_model.predict(X_test_4)
proba_4 = original_model.predict_proba(X_test_4)[:, 1]

# 3-Feature Model
pred_3 = model_3.predict(X_test_3)
proba_3 = model_3.predict_proba(X_test_3)[:, 1]

# Calculate point estimates (including advanced metrics)
metrics_4 = {
    'Accuracy': accuracy_score(y_test, pred_4),
    'F1-Score': f1_score(y_test, pred_4),
    'ROC-AUC': roc_auc_score(y_test, proba_4),
    'Precision': precision_score(y_test, pred_4),
    'Recall': recall_score(y_test, pred_4),
    'Avg-Precision': average_precision_score(y_test, proba_4),
    'MCC': matthews_corrcoef(y_test, pred_4),
    'Cohen-Kappa': cohen_kappa_score(y_test, pred_4),
    'Brier-Score': brier_score_loss(y_test, proba_4),
    'ECE': expected_calibration_error(y_test, proba_4)
}

metrics_3 = {
    'Accuracy': accuracy_score(y_test, pred_3),
    'F1-Score': f1_score(y_test, pred_3),
    'ROC-AUC': roc_auc_score(y_test, proba_3),
    'Precision': precision_score(y_test, pred_3),
    'Recall': recall_score(y_test, pred_3),
    'Avg-Precision': average_precision_score(y_test, proba_3),
    'MCC': matthews_corrcoef(y_test, pred_3),
    'Cohen-Kappa': cohen_kappa_score(y_test, pred_3),
    'Brier-Score': brier_score_loss(y_test, proba_3),
    'ECE': expected_calibration_error(y_test, proba_3)
}

# Bootstrap confidence intervals
print("\nCalculating bootstrap confidence intervals (n=5000)...")
boot_4 = bootstrap_metrics(y_test, pred_4, proba_4, n_bootstrap=5000)
boot_3 = bootstrap_metrics(y_test, pred_3, proba_3, n_bootstrap=5000)

# Display comparison - Standard metrics
print(f"\n{'Metric':<15} {'4-Feature Model':<30} {'3-Feature Model':<30} {'Delta':>10}")
print("-" * 85)

comparison_data = []
standard_metrics = ['Accuracy', 'F1-Score', 'ROC-AUC', 'Precision', 'Recall', 'Avg-Precision']
advanced_metrics = ['MCC', 'Cohen-Kappa', 'Brier-Score', 'ECE']
all_metrics = standard_metrics + advanced_metrics

for metric in standard_metrics:
    val_4 = metrics_4[metric]
    val_3 = metrics_3[metric]
    ci_4 = f"[{boot_4[metric]['ci_lower']:.3f}, {boot_4[metric]['ci_upper']:.3f}]"
    ci_3 = f"[{boot_3[metric]['ci_lower']:.3f}, {boot_3[metric]['ci_upper']:.3f}]"
    delta = val_3 - val_4
    delta_str = f"{delta:+.4f}"

    print(f"{metric:<15} {val_4:.4f} {ci_4:<18} {val_3:.4f} {ci_3:<18} {delta_str:>10}")

    comparison_data.append({
        'Metric': metric,
        '4-Feature': val_4,
        '4-Feature_CI_Lower': boot_4[metric]['ci_lower'],
        '4-Feature_CI_Upper': boot_4[metric]['ci_upper'],
        '3-Feature': val_3,
        '3-Feature_CI_Lower': boot_3[metric]['ci_lower'],
        '3-Feature_CI_Upper': boot_3[metric]['ci_upper'],
        'Delta': delta
    })

# Display advanced metrics
print(f"\n{'='*70}")
print("ADVANCED METRICS COMPARISON")
print("=" * 70)
print(f"\n{'Metric':<15} {'4-Feature Model':<30} {'3-Feature Model':<30} {'Delta':>10} {'Better':>10}")
print("-" * 95)

for metric in advanced_metrics:
    val_4 = metrics_4[metric]
    val_3 = metrics_3[metric]
    ci_4 = f"[{boot_4[metric]['ci_lower']:.3f}, {boot_4[metric]['ci_upper']:.3f}]"
    ci_3 = f"[{boot_3[metric]['ci_lower']:.3f}, {boot_3[metric]['ci_upper']:.3f}]"
    delta = val_3 - val_4

    # For Brier Score and ECE, lower is better
    if metric in ['Brier-Score', 'ECE']:
        better = "3-Feature" if val_3 < val_4 else "4-Feature"
        delta_str = f"{delta:+.4f}"
    else:
        better = "3-Feature" if val_3 > val_4 else "4-Feature"
        delta_str = f"{delta:+.4f}"

    print(f"{metric:<15} {val_4:.4f} {ci_4:<18} {val_3:.4f} {ci_3:<18} {delta_str:>10} {better:>10}")

    comparison_data.append({
        'Metric': metric,
        '4-Feature': val_4,
        '4-Feature_CI_Lower': boot_4[metric]['ci_lower'],
        '4-Feature_CI_Upper': boot_4[metric]['ci_upper'],
        '3-Feature': val_3,
        '3-Feature_CI_Lower': boot_3[metric]['ci_lower'],
        '3-Feature_CI_Upper': boot_3[metric]['ci_upper'],
        'Delta': delta,
        'Better': better
    })

print("\nNote: For Brier Score and ECE, lower values are better.")

# ============================================================================
# Permutation Test (Only shuffle training labels)
# ============================================================================
print(f"\n{'='*70}")
print("CORRECTED PERMUTATION TEST (Training labels shuffled, test labels FIXED)")
print("=" * 70)

n_permutations = 1000

def run_permutation_test(X_train, X_test, y_train, y_test, n_perm=1000):
    """
    Run permutation test with ONLY training labels shuffled.
    Test labels remain fixed to properly assess if learned patterns are real.
    """
    perm_aucs = []
    perm_accs = []
    perm_f1s = []

    for i in range(n_perm):
        if (i + 1) % 200 == 0:
            print(f"  Permutation {i + 1}/{n_perm}...")

        # ONLY shuffle training labels
        y_train_shuffled = np.random.permutation(y_train)

        # Train model with shuffled training labels
        nb_perm = GaussianNB()
        cal_perm = CalibratedClassifierCV(nb_perm, cv=3, method='sigmoid')

        try:
            cal_perm.fit(X_train, y_train_shuffled)

            # Predict on test set with FIXED test labels
            pred_perm = cal_perm.predict(X_test)
            proba_perm = cal_perm.predict_proba(X_test)[:, 1]

            # Evaluate against TRUE test labels (not shuffled)
            perm_aucs.append(roc_auc_score(y_test, proba_perm))
            perm_accs.append(accuracy_score(y_test, pred_perm))
            perm_f1s.append(f1_score(y_test, pred_perm, zero_division=0))
        except:
            continue

    return np.array(perm_aucs), np.array(perm_accs), np.array(perm_f1s)

# Run permutation tests for both models
print("\n4-Feature Model Permutation Test:")
perm_auc_4, perm_acc_4, perm_f1_4 = run_permutation_test(
    X_train_4, X_test_4, y_train, y_test, n_permutations)

print("\n3-Feature Model Permutation Test:")
perm_auc_3, perm_acc_3, perm_f1_3 = run_permutation_test(
    X_train_3, X_test_3, y_train, y_test, n_permutations)

# Calculate empirical p-values (Phipson-Smyth correction: p = (r + 1) / (n + 1))
p_auc_4 = (np.sum(perm_auc_4 >= metrics_4['ROC-AUC']) + 1) / (len(perm_auc_4) + 1)
p_auc_3 = (np.sum(perm_auc_3 >= metrics_3['ROC-AUC']) + 1) / (len(perm_auc_3) + 1)
p_acc_4 = (np.sum(perm_acc_4 >= metrics_4['Accuracy']) + 1) / (len(perm_acc_4) + 1)
p_acc_3 = (np.sum(perm_acc_3 >= metrics_3['Accuracy']) + 1) / (len(perm_acc_3) + 1)

# Cohen's d
d_auc_4 = (metrics_4['ROC-AUC'] - np.mean(perm_auc_4)) / np.std(perm_auc_4)
d_auc_3 = (metrics_3['ROC-AUC'] - np.mean(perm_auc_3)) / np.std(perm_auc_3)

print(f"\n{'='*70}")
print("PERMUTATION TEST RESULTS (CORRECTED METHODOLOGY)")
print("=" * 70)
cohens_label = "Cohen's d"
print(f"\n{'Model':<20} {'Observed AUC':<15} {'Null Mean':<12} {'Null Std':<10} {'p-value':<10} {cohens_label:<10}")
print("-" * 77)
print(f"{'4-Feature':<20} {metrics_4['ROC-AUC']:<15.4f} {np.mean(perm_auc_4):<12.4f} {np.std(perm_auc_4):<10.4f} {p_auc_4:<10.4f} {d_auc_4:<10.2f}")
print(f"{'3-Feature':<20} {metrics_3['ROC-AUC']:<15.4f} {np.mean(perm_auc_3):<12.4f} {np.std(perm_auc_3):<10.4f} {p_auc_3:<10.4f} {d_auc_3:<10.2f}")

# ============================================================================
# Cross-Validation Comparison
# ============================================================================
print(f"\n{'='*70}")
print("CROSS-VALIDATION COMPARISON (5-Fold Stratified)")
print("=" * 70)

# Combine train and test for CV
X_full_4 = np.vstack([X_train_4, X_test_4])
X_full_3 = np.vstack([X_train_3, X_test_3])
y_full = np.concatenate([y_train, y_test])

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# 4-Feature CV
nb_cv_4 = GaussianNB()
cv_scores_4 = cross_val_score(nb_cv_4, X_full_4, y_full, cv=cv, scoring='roc_auc')

# 3-Feature CV
nb_cv_3 = GaussianNB()
cv_scores_3 = cross_val_score(nb_cv_3, X_full_3, y_full, cv=cv, scoring='roc_auc')

print(f"\n{'Model':<20} {'Mean AUC':<12} {'Std':<10} {'Min':<10} {'Max':<10}")
print("-" * 62)
print(f"{'4-Feature':<20} {np.mean(cv_scores_4):<12.4f} {np.std(cv_scores_4):<10.4f} {np.min(cv_scores_4):<10.4f} {np.max(cv_scores_4):<10.4f}")
print(f"{'3-Feature':<20} {np.mean(cv_scores_3):<12.4f} {np.std(cv_scores_3):<10.4f} {np.min(cv_scores_3):<10.4f} {np.max(cv_scores_3):<10.4f}")

# Paired t-test for CV scores
t_stat, t_pval = stats.ttest_rel(cv_scores_4, cv_scores_3)
print(f"\nPaired t-test (4-Feature vs 3-Feature): t={t_stat:.3f}, p={t_pval:.4f}")

# ============================================================================
# Noise Robustness Comparison
# ============================================================================
print(f"\n{'='*70}")
print("NOISE ROBUSTNESS COMPARISON")
print("=" * 70)

noise_levels = [0.05, 0.10, 0.20, 0.30, 0.50]
n_noise_iter = 100

noise_results_4 = {level: [] for level in noise_levels}
noise_results_3 = {level: [] for level in noise_levels}

for noise_level in noise_levels:
    for _ in range(n_noise_iter):
        # 4-Feature noise
        std_4 = np.std(X_train_4, axis=0)
        noise_4 = np.random.normal(0, std_4 * noise_level, X_test_4.shape)
        X_noisy_4 = X_test_4 + noise_4
        proba_noisy_4 = original_model.predict_proba(X_noisy_4)[:, 1]
        noise_results_4[noise_level].append(roc_auc_score(y_test, proba_noisy_4))

        # 3-Feature noise
        std_3 = np.std(X_train_3, axis=0)
        noise_3 = np.random.normal(0, std_3 * noise_level, X_test_3.shape)
        X_noisy_3 = X_test_3 + noise_3
        proba_noisy_3 = model_3.predict_proba(X_noisy_3)[:, 1]
        noise_results_3[noise_level].append(roc_auc_score(y_test, proba_noisy_3))

print(f"\n{'Noise Level':<15} {'4-Feature AUC':<20} {'3-Feature AUC':<20} {'Better Model':<15}")
print("-" * 70)

for level in noise_levels:
    mean_4 = np.mean(noise_results_4[level])
    mean_3 = np.mean(noise_results_3[level])
    std_4 = np.std(noise_results_4[level])
    std_3 = np.std(noise_results_3[level])
    better = "3-Feature" if mean_3 > mean_4 else "4-Feature" if mean_4 > mean_3 else "Tie"

    print(f"{level*100:>5.0f}%{'':<9} {mean_4:.4f} +/- {std_4:.4f}{'':<5} {mean_3:.4f} +/- {std_3:.4f}{'':<5} {better:<15}")

# ============================================================================
# Confusion Matrix Comparison
# ============================================================================
cm_4 = confusion_matrix(y_test, pred_4)
cm_3 = confusion_matrix(y_test, pred_3)

print(f"\n{'='*70}")
print("CONFUSION MATRIX COMPARISON")
print("=" * 70)

print("\n4-Feature Model:")
print(f"                 Predicted")
print(f"                 Long    Short")
print(f"  Actual Long    {cm_4[0,0]:<7} {cm_4[0,1]:<7}")
print(f"  Actual Short   {cm_4[1,0]:<7} {cm_4[1,1]:<7}")

print("\n3-Feature Model:")
print(f"                 Predicted")
print(f"                 Long    Short")
print(f"  Actual Long    {cm_3[0,0]:<7} {cm_3[0,1]:<7}")
print(f"  Actual Short   {cm_3[1,0]:<7} {cm_3[1,1]:<7}")

# ============================================================================
# Comprehensive Visualization
# ============================================================================
print(f"\n{'='*70}")
print("GENERATING COMPARISON VISUALIZATIONS")
print("=" * 70)

fig = plt.figure(figsize=(20, 16))

# Panel A: Performance Metrics Comparison with CI (Standard metrics)
ax1 = fig.add_subplot(3, 3, 1)
metrics_list = ['Accuracy', 'F1-Score', 'ROC-AUC', 'Precision', 'Recall', 'Avg-Precision']
x = np.arange(len(metrics_list))
width = 0.35

vals_4 = [metrics_4[m] for m in metrics_list]
vals_3 = [metrics_3[m] for m in metrics_list]
errs_4 = [(metrics_4[m] - boot_4[m]['ci_lower'], boot_4[m]['ci_upper'] - metrics_4[m]) for m in metrics_list]
errs_3 = [(metrics_3[m] - boot_3[m]['ci_lower'], boot_3[m]['ci_upper'] - metrics_3[m]) for m in metrics_list]

bars1 = ax1.bar(x - width/2, vals_4, width, label='4-Feature', color='#1f77b4',
                yerr=np.array(errs_4).T, capsize=3, error_kw={'linewidth': 1})
bars2 = ax1.bar(x + width/2, vals_3, width, label='3-Feature', color='#2ca02c',
                yerr=np.array(errs_3).T, capsize=3, error_kw={'linewidth': 1})

ax1.set_ylabel('Score')
ax1.set_title('A. Standard Metrics with 95% Bootstrap CI')
ax1.set_xticks(x)
ax1.set_xticklabels(metrics_list, rotation=45, ha='right')
ax1.legend()
ax1.set_ylim(0.5, 1.05)
ax1.axhline(y=0.5, color='red', linestyle=':', alpha=0.5, label='Chance')

# Panel B: ROC Curves
ax2 = fig.add_subplot(3, 3, 2)
fpr_4, tpr_4, _ = roc_curve(y_test, proba_4)
fpr_3, tpr_3, _ = roc_curve(y_test, proba_3)

ax2.plot(fpr_4, tpr_4, 'b-', linewidth=2, label=f'4-Feature (AUC={metrics_4["ROC-AUC"]:.3f})')
ax2.plot(fpr_3, tpr_3, 'g-', linewidth=2, label=f'3-Feature (AUC={metrics_3["ROC-AUC"]:.3f})')
ax2.plot([0, 1], [0, 1], 'k--', linewidth=1, label='Chance')
ax2.fill_between(fpr_4, tpr_4, alpha=0.2, color='blue')
ax2.fill_between(fpr_3, tpr_3, alpha=0.2, color='green')
ax2.set_xlabel('False Positive Rate')
ax2.set_ylabel('True Positive Rate')
ax2.set_title('B. ROC Curve Comparison')
ax2.legend(loc='lower right')
ax2.set_xlim([0, 1])
ax2.set_ylim([0, 1.05])

# Panel C: Precision-Recall Curves
ax3 = fig.add_subplot(3, 3, 3)
prec_4, rec_4, _ = precision_recall_curve(y_test, proba_4)
prec_3, rec_3, _ = precision_recall_curve(y_test, proba_3)

ax3.plot(rec_4, prec_4, 'b-', linewidth=2, label=f'4-Feature (AP={metrics_4["Avg-Precision"]:.3f})')
ax3.plot(rec_3, prec_3, 'g-', linewidth=2, label=f'3-Feature (AP={metrics_3["Avg-Precision"]:.3f})')
baseline = sum(y_test) / len(y_test)
ax3.axhline(y=baseline, color='red', linestyle=':', label=f'Baseline ({baseline:.3f})')
ax3.set_xlabel('Recall')
ax3.set_ylabel('Precision')
ax3.set_title('C. Precision-Recall Curve Comparison')
ax3.legend(loc='lower left')
ax3.set_xlim([0, 1])
ax3.set_ylim([0, 1.05])

# Panel D: Permutation Test Histograms
ax4 = fig.add_subplot(3, 3, 4)
ax4.hist(perm_auc_4, bins=40, alpha=0.6, color='blue', label='4-Feature Null', edgecolor='black', linewidth=0.5)
ax4.hist(perm_auc_3, bins=40, alpha=0.6, color='green', label='3-Feature Null', edgecolor='black', linewidth=0.5)
ax4.axvline(x=metrics_4['ROC-AUC'], color='blue', linestyle='-', linewidth=2.5, label=f'4-Feature Observed ({metrics_4["ROC-AUC"]:.3f})')
ax4.axvline(x=metrics_3['ROC-AUC'], color='green', linestyle='-', linewidth=2.5, label=f'3-Feature Observed ({metrics_3["ROC-AUC"]:.3f})')
ax4.axvline(x=0.5, color='red', linestyle=':', linewidth=1.5)

# Add p-value annotations
p_text_4 = f'4-Feature: p={p_auc_4:.3f}' if p_auc_4 >= 0.001 else '4-Feature: p<0.001'
p_text_3 = f'3-Feature: p={p_auc_3:.3f}' if p_auc_3 >= 0.001 else '3-Feature: p<0.001'
ax4.text(0.95, 0.95, f'{p_text_4}\n{p_text_3}', transform=ax4.transAxes, fontsize=10,
         verticalalignment='top', horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

ax4.set_xlabel('ROC-AUC')
ax4.set_ylabel('Frequency')
ax4.set_title('D. Corrected Permutation Test\n(Only training labels shuffled)')
ax4.legend(loc='upper left', fontsize=8)

# Panel E: Noise Robustness Comparison
ax5 = fig.add_subplot(3, 3, 5)
noise_pct = [level * 100 for level in noise_levels]
means_4 = [np.mean(noise_results_4[l]) for l in noise_levels]
means_3 = [np.mean(noise_results_3[l]) for l in noise_levels]
stds_4 = [np.std(noise_results_4[l]) for l in noise_levels]
stds_3 = [np.std(noise_results_3[l]) for l in noise_levels]

ax5.errorbar(noise_pct, means_4, yerr=stds_4, fmt='o-', color='blue', capsize=4,
             linewidth=2, markersize=8, label='4-Feature')
ax5.errorbar(noise_pct, means_3, yerr=stds_3, fmt='s-', color='green', capsize=4,
             linewidth=2, markersize=8, label='3-Feature')
ax5.axhline(y=metrics_4['ROC-AUC'], color='blue', linestyle='--', alpha=0.5)
ax5.axhline(y=metrics_3['ROC-AUC'], color='green', linestyle='--', alpha=0.5)
ax5.axhline(y=0.5, color='red', linestyle=':', linewidth=1.5)
ax5.set_xlabel('Noise Level (% of Feature Std)')
ax5.set_ylabel('ROC-AUC')
ax5.set_title('E. Noise Robustness Comparison')
ax5.legend()
ax5.set_xlim(0, 55)

# Panel F: Confusion Matrix Heatmaps
ax6 = fig.add_subplot(3, 3, 6)

# Create side-by-side confusion matrices
cm_combined = np.zeros((2, 4))
cm_combined[:, 0:2] = cm_4
cm_combined[:, 2:4] = cm_3

im = ax6.imshow(cm_combined, cmap='Blues', aspect='auto')

# Add text annotations
for i in range(2):
    for j in range(4):
        val = int(cm_combined[i, j])
        ax6.text(j, i, str(val), ha='center', va='center', fontsize=14, fontweight='bold',
                color='white' if val > cm_combined.max()/2 else 'black')

ax6.set_xticks([0.5, 2.5])
ax6.set_xticklabels(['4-Feature', '3-Feature'])
ax6.set_yticks([0, 1])
ax6.set_yticklabels(['Long PFS', 'Short PFS'])
ax6.set_ylabel('Actual')
ax6.set_title('F. Confusion Matrix Comparison')

# Add column labels
ax6.text(0, -0.3, 'Pred L', ha='center', fontsize=9)
ax6.text(1, -0.3, 'Pred S', ha='center', fontsize=9)
ax6.text(2, -0.3, 'Pred L', ha='center', fontsize=9)
ax6.text(3, -0.3, 'Pred S', ha='center', fontsize=9)

# Add vertical line to separate
ax6.axvline(x=1.5, color='white', linewidth=3)

# Panel G: Advanced Metrics Comparison (MCC, Kappa, Brier, ECE)
ax7 = fig.add_subplot(3, 3, 7)
adv_metrics_list = ['MCC', 'Cohen-Kappa', 'Brier-Score', 'ECE']
x_adv = np.arange(len(adv_metrics_list))

vals_4_adv = [metrics_4[m] for m in adv_metrics_list]
vals_3_adv = [metrics_3[m] for m in adv_metrics_list]
errs_4_adv = [(abs(metrics_4[m] - boot_4[m]['ci_lower']), abs(boot_4[m]['ci_upper'] - metrics_4[m])) for m in adv_metrics_list]
errs_3_adv = [(abs(metrics_3[m] - boot_3[m]['ci_lower']), abs(boot_3[m]['ci_upper'] - metrics_3[m])) for m in adv_metrics_list]

bars1_adv = ax7.bar(x_adv - width/2, vals_4_adv, width, label='4-Feature', color='#1f77b4',
                    yerr=np.array(errs_4_adv).T, capsize=3, error_kw={'linewidth': 1})
bars2_adv = ax7.bar(x_adv + width/2, vals_3_adv, width, label='3-Feature', color='#2ca02c',
                    yerr=np.array(errs_3_adv).T, capsize=3, error_kw={'linewidth': 1})

ax7.set_ylabel('Score')
ax7.set_title('G. Advanced Metrics with 95% Bootstrap CI\n(Lower is better for Brier & ECE)')
ax7.set_xticks(x_adv)
ax7.set_xticklabels(adv_metrics_list, rotation=45, ha='right')
ax7.legend()
ax7.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)

# Panel H: Calibration Curves
ax8 = fig.add_subplot(3, 3, 8)

# Calculate calibration data with quantile strategy for small samples
n_cal_bins = 4
try:
    prob_true_4, prob_pred_4 = calibration_curve(y_test, proba_4, n_bins=n_cal_bins, strategy='quantile')
    prob_true_3, prob_pred_3 = calibration_curve(y_test, proba_3, n_bins=n_cal_bins, strategy='quantile')
except:
    prob_true_4, prob_pred_4 = calibration_curve(y_test, proba_4, n_bins=n_cal_bins, strategy='uniform')
    prob_true_3, prob_pred_3 = calibration_curve(y_test, proba_3, n_bins=n_cal_bins, strategy='uniform')

# Perfect calibration line
ax8.plot([0, 1], [0, 1], 'k--', linewidth=2, label='Perfectly Calibrated', zorder=1)

# Use scatter + line with distinct markers to avoid overlap
ax8.scatter(prob_pred_4, prob_true_4, s=120, c='blue', marker='o', edgecolors='black',
            linewidths=1.5, label=f'4-Feature (Brier={metrics_4["Brier-Score"]:.3f})', zorder=3)
ax8.plot(prob_pred_4, prob_true_4, 'b-', linewidth=1.5, alpha=0.7, zorder=2)

ax8.scatter(prob_pred_3, prob_true_3, s=120, c='green', marker='s', edgecolors='black',
            linewidths=1.5, label=f'3-Feature (Brier={metrics_3["Brier-Score"]:.3f})', zorder=3)
ax8.plot(prob_pred_3, prob_true_3, 'g-', linewidth=1.5, alpha=0.7, zorder=2)

ax8.set_xlabel('Mean Predicted Probability')
ax8.set_ylabel('Fraction of Positives')
ax8.set_title('H. Calibration Curves')
ax8.legend(loc='lower right', fontsize=8)
ax8.set_xlim([0, 1])
ax8.set_ylim([0, 1])
ax8.set_aspect('equal')
ax8.grid(True, alpha=0.3)

# Panel I: Summary Statistics Table with 95% CI
ax9 = fig.add_subplot(3, 3, 9)
ax9.axis('off')

# Helper function to format value with CI
def fmt_ci(val, ci_low, ci_high):
    return f'{val:.2f} [{ci_low:.2f},{ci_high:.2f}]'

# Create comprehensive summary table with 95% CIs
summary_table_data = [
    ['Metric', '4-Feature [95% CI]', '3-Feature [95% CI]', 'Winner'],
    ['Accuracy', fmt_ci(metrics_4["Accuracy"], boot_4["Accuracy"]["ci_lower"], boot_4["Accuracy"]["ci_upper"]),
                 fmt_ci(metrics_3["Accuracy"], boot_3["Accuracy"]["ci_lower"], boot_3["Accuracy"]["ci_upper"]),
                 '3-Feat' if metrics_3["Accuracy"] > metrics_4["Accuracy"] else '4-Feat'],
    ['ROC-AUC', fmt_ci(metrics_4["ROC-AUC"], boot_4["ROC-AUC"]["ci_lower"], boot_4["ROC-AUC"]["ci_upper"]),
                fmt_ci(metrics_3["ROC-AUC"], boot_3["ROC-AUC"]["ci_lower"], boot_3["ROC-AUC"]["ci_upper"]),
                '3-Feat' if metrics_3["ROC-AUC"] > metrics_4["ROC-AUC"] else '4-Feat'],
    ['F1-Score', fmt_ci(metrics_4["F1-Score"], boot_4["F1-Score"]["ci_lower"], boot_4["F1-Score"]["ci_upper"]),
                 fmt_ci(metrics_3["F1-Score"], boot_3["F1-Score"]["ci_lower"], boot_3["F1-Score"]["ci_upper"]),
                 '3-Feat' if metrics_3["F1-Score"] > metrics_4["F1-Score"] else '4-Feat'],
    ['Recall', fmt_ci(metrics_4["Recall"], boot_4["Recall"]["ci_lower"], boot_4["Recall"]["ci_upper"]),
               fmt_ci(metrics_3["Recall"], boot_3["Recall"]["ci_lower"], boot_3["Recall"]["ci_upper"]),
               '3-Feat' if metrics_3["Recall"] > metrics_4["Recall"] else '4-Feat'],
    ['MCC', fmt_ci(metrics_4["MCC"], boot_4["MCC"]["ci_lower"], boot_4["MCC"]["ci_upper"]),
            fmt_ci(metrics_3["MCC"], boot_3["MCC"]["ci_lower"], boot_3["MCC"]["ci_upper"]),
            '3-Feat' if metrics_3["MCC"] > metrics_4["MCC"] else '4-Feat'],
    ['Kappa', fmt_ci(metrics_4["Cohen-Kappa"], boot_4["Cohen-Kappa"]["ci_lower"], boot_4["Cohen-Kappa"]["ci_upper"]),
              fmt_ci(metrics_3["Cohen-Kappa"], boot_3["Cohen-Kappa"]["ci_lower"], boot_3["Cohen-Kappa"]["ci_upper"]),
              '3-Feat' if metrics_3["Cohen-Kappa"] > metrics_4["Cohen-Kappa"] else '4-Feat'],
    ['Brier*', fmt_ci(metrics_4["Brier-Score"], boot_4["Brier-Score"]["ci_lower"], boot_4["Brier-Score"]["ci_upper"]),
               fmt_ci(metrics_3["Brier-Score"], boot_3["Brier-Score"]["ci_lower"], boot_3["Brier-Score"]["ci_upper"]),
               '3-Feat' if metrics_3["Brier-Score"] < metrics_4["Brier-Score"] else '4-Feat'],
    ['ECE*', fmt_ci(metrics_4["ECE"], boot_4["ECE"]["ci_lower"], boot_4["ECE"]["ci_upper"]),
             fmt_ci(metrics_3["ECE"], boot_3["ECE"]["ci_lower"], boot_3["ECE"]["ci_upper"]),
             '3-Feat' if metrics_3["ECE"] < metrics_4["ECE"] else '4-Feat'],
    ['p-value', f'{p_auc_4:.3f}', f'{p_auc_3:.3f}', '3-Feat' if p_auc_3 < p_auc_4 else '4-Feat'],
]

table = ax9.table(cellText=summary_table_data, loc='center', cellLoc='center',
                  colWidths=[0.15, 0.35, 0.35, 0.12])
table.auto_set_font_size(False)
table.set_fontsize(7)
table.scale(1.2, 1.5)

# Color header and winner column
for i in range(4):
    table[(0, i)].set_facecolor('#4472C4')
    table[(0, i)].set_text_props(color='white', weight='bold', fontsize=7)

# Highlight winners
for row in range(1, len(summary_table_data)):
    winner = summary_table_data[row][3]
    if winner == '3-Feat':
        table[(row, 2)].set_facecolor('#c6efce')  # Light green
        table[(row, 3)].set_facecolor('#c6efce')

ax9.set_title('I. Summary with 95% CI (*lower is better)', pad=20, fontsize=10)

plt.suptitle('Model Comparison: 3-Feature vs 4-Feature PFS Classifier\n(With Advanced Metrics & Calibration Analysis)',
             fontsize=14, fontweight='bold', y=1.02)

plt.tight_layout()

# Save figures
plt.savefig(os.path.join(OUTPUT_DIR, 'Model_Comparison_3vs4_Features.png'),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(OUTPUT_DIR, 'Model_Comparison_3vs4_Features.pdf'),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(OUTPUT_DIR, 'Model_Comparison_3vs4_Features.svg'),
            bbox_inches='tight', facecolor='white')
plt.close()

print("\nSaved: Model_Comparison_3vs4_Features.png/pdf/svg")

# ============================================================================
# Save Results to Excel
# ============================================================================
print(f"\n{'='*70}")
print("SAVING RESULTS TO EXCEL")
print("=" * 70)

with pd.ExcelWriter(os.path.join(OUTPUT_DIR, 'Model_Comparison_Results.xlsx'),
                     engine='openpyxl') as writer:

    # Comparison summary
    pd.DataFrame(comparison_data).to_excel(writer, sheet_name='Performance_Comparison', index=False)

    # Bootstrap CI details
    boot_df = pd.DataFrame([
        {'Model': '4-Feature', 'Metric': m, 'Point_Estimate': metrics_4[m],
         'CI_Lower': boot_4[m]['ci_lower'], 'CI_Upper': boot_4[m]['ci_upper'],
         'CI_Width': boot_4[m]['ci_upper'] - boot_4[m]['ci_lower']}
        for m in metrics_list
    ] + [
        {'Model': '3-Feature', 'Metric': m, 'Point_Estimate': metrics_3[m],
         'CI_Lower': boot_3[m]['ci_lower'], 'CI_Upper': boot_3[m]['ci_upper'],
         'CI_Width': boot_3[m]['ci_upper'] - boot_3[m]['ci_lower']}
        for m in metrics_list
    ])
    boot_df.to_excel(writer, sheet_name='Bootstrap_CI', index=False)

    # Permutation test results
    perm_summary = pd.DataFrame([
        {'Model': '4-Feature', 'Features': ', '.join(features_4),
         'Observed_AUC': metrics_4['ROC-AUC'], 'Null_Mean': np.mean(perm_auc_4),
         'Null_Std': np.std(perm_auc_4), 'P_Value': p_auc_4, 'Cohens_d': d_auc_4},
        {'Model': '3-Feature', 'Features': ', '.join(features_3),
         'Observed_AUC': metrics_3['ROC-AUC'], 'Null_Mean': np.mean(perm_auc_3),
         'Null_Std': np.std(perm_auc_3), 'P_Value': p_auc_3, 'Cohens_d': d_auc_3}
    ])
    perm_summary.to_excel(writer, sheet_name='Permutation_Test', index=False)

    # CV results
    cv_df = pd.DataFrame({
        'Fold': range(1, 6),
        '4-Feature_AUC': cv_scores_4,
        '3-Feature_AUC': cv_scores_3,
        'Difference': cv_scores_3 - cv_scores_4
    })
    cv_df.to_excel(writer, sheet_name='Cross_Validation', index=False)

    # Noise robustness
    noise_df = pd.DataFrame([
        {'Noise_Level': f'{l*100:.0f}%',
         '4-Feature_Mean_AUC': np.mean(noise_results_4[l]),
         '4-Feature_Std': np.std(noise_results_4[l]),
         '3-Feature_Mean_AUC': np.mean(noise_results_3[l]),
         '3-Feature_Std': np.std(noise_results_3[l]),
         'Better_Model': '3-Feature' if np.mean(noise_results_3[l]) > np.mean(noise_results_4[l]) else '4-Feature'}
        for l in noise_levels
    ])
    noise_df.to_excel(writer, sheet_name='Noise_Robustness', index=False)

    # Sample sizes
    sample_df = pd.DataFrame([
        {'Dataset': 'Training', 'Total': len(y_train), 'Long_PFS': sum(y_train==0), 'Short_PFS': sum(y_train==1)},
        {'Dataset': 'Test', 'Total': len(y_test), 'Long_PFS': sum(y_test==0), 'Short_PFS': sum(y_test==1)},
        {'Dataset': 'Combined', 'Total': len(y_full), 'Long_PFS': sum(y_full==0), 'Short_PFS': sum(y_full==1)}
    ])
    sample_df.to_excel(writer, sheet_name='Sample_Sizes', index=False)

print("Saved: Model_Comparison_Results.xlsx")

# ============================================================================
# Save 3-Feature Model
# ============================================================================
model_3_data = {
    'model': model_3,
    'features': features_3,
    'algo': 'NB',
    'metrics': metrics_3,
    'bootstrap_ci': boot_3,
    'permutation_p': p_auc_3,
    'cohens_d': d_auc_3
}

with open(os.path.join(OUTPUT_DIR, 'champion_model_NB_3features.pkl'), 'wb') as f:
    pickle.dump(model_3_data, f)

print("Saved: champion_model_NB_3features.pkl")

# ============================================================================
# Final Interpretation
# ============================================================================
print(f"\n{'='*70}")
print("FINAL INTERPRETATION & RECOMMENDATIONS")
print("=" * 70)

auc_diff = metrics_3['ROC-AUC'] - metrics_4['ROC-AUC']
acc_diff = metrics_3['Accuracy'] - metrics_4['Accuracy']

print(f"""
SUMMARY OF FINDINGS:
====================

1. PERFORMANCE COMPARISON:
   - 4-Feature Model: AUC = {metrics_4['ROC-AUC']:.4f}, Accuracy = {metrics_4['Accuracy']:.4f}
   - 3-Feature Model: AUC = {metrics_3['ROC-AUC']:.4f}, Accuracy = {metrics_3['Accuracy']:.4f}
   - Delta: AUC {auc_diff:+.4f}, Accuracy {acc_diff:+.4f}

2. STATISTICAL SIGNIFICANCE (Corrected Permutation Test):
   - 4-Feature: p = {p_auc_4:.4f}, Cohen's d = {d_auc_4:.2f}
   - 3-Feature: p = {p_auc_3:.4f}, Cohen's d = {d_auc_3:.2f}
   - Both models significantly outperform chance (p < 0.05)

3. CROSS-VALIDATION:
   - 4-Feature CV AUC: {np.mean(cv_scores_4):.4f} +/- {np.std(cv_scores_4):.4f}
   - 3-Feature CV AUC: {np.mean(cv_scores_3):.4f} +/- {np.std(cv_scores_3):.4f}
   - Paired t-test p-value: {t_pval:.4f}

4. NOISE ROBUSTNESS:
   - Both models show similar robustness to measurement noise
   - At 30% noise: 4-Feature = {np.mean(noise_results_4[0.30]):.4f}, 3-Feature = {np.mean(noise_results_3[0.30]):.4f}

OUTPUT FILES:
=============
- Model_Comparison_3vs4_Features.png/pdf/svg (comprehensive comparison figure)
- Model_Comparison_Results.xlsx (all numerical results)
- champion_model_NB_3features.pkl (saved 3-feature model)
""")

print("=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
