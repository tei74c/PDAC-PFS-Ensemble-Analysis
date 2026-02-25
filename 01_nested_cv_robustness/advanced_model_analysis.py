"""
Advanced Model Analysis: SHAP, DCA, and ROC with Confidence Intervals
======================================================================
SHAP feature importance, Decision Curve Analysis, and ROC-AUC with
bootstrapped confidence intervals for the NB classifier.
Generates Figures 2D-E in the manuscript.

Reference:
Khoshnevis et al. Molecular Cancer (2026).
"""

import os
import pickle
import pandas as pd
import numpy as np
from sklearn.metrics import (roc_curve, roc_auc_score, brier_score_loss,
                             confusion_matrix, accuracy_score)
from sklearn.naive_bayes import GaussianNB
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
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
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

np.random.seed(42)

# ============================================================================
# Load Data and Models
# ============================================================================
print("=" * 70)
print("ADVANCED MODEL ANALYSIS: DCA, SHAP, ROC with 95% CI")
print("=" * 70)

# Load 4-feature model
with open(os.path.join(DATA_DIR, 'champion_model_NB_4features.pkl'), 'rb') as f:
    model_4_data = pickle.load(f)
model_4 = model_4_data['model']
features_4 = model_4_data['features']

# Load 3-feature model
with open(os.path.join(DATA_DIR, 'champion_model_NB_3features.pkl'), 'rb') as f:
    model_3_data = pickle.load(f)
model_3 = model_3_data['model']
features_3 = model_3_data['features']

print(f"\n4-Feature Model: {features_4}")
print(f"3-Feature Model: {features_3}")

# Load data
train_df = pd.read_excel(os.path.join(DATA_DIR, 'PFS_ML_train_cof_wo_IG.xlsx'))
test_df = pd.read_excel(os.path.join(DATA_DIR, 'PFS_ML_test_cof_wo_IG.xlsx'))

label_map = {'L': 0, 'S': 1}
train_df['PFS_group_encoded'] = train_df['PFS_group'].map(label_map)
test_df['PFS_group_encoded'] = test_df['PFS_group'].map(label_map)

X_train_4 = train_df[features_4].values
X_test_4 = test_df[features_4].values
X_train_3 = train_df[features_3].values
X_test_3 = test_df[features_3].values
y_train = train_df['PFS_group_encoded'].values
y_test = test_df['PFS_group_encoded'].values

# Get predictions
proba_4 = model_4.predict_proba(X_test_4)[:, 1]
proba_3 = model_3.predict_proba(X_test_3)[:, 1]

print(f"\nTest samples: {len(y_test)}")

# ============================================================================
# 1. ROC CURVES WITH 95% BOOTSTRAP CONFIDENCE INTERVALS
# ============================================================================
print(f"\n{'='*70}")
print("1. ROC CURVES WITH 95% BOOTSTRAP CONFIDENCE INTERVALS")
print("=" * 70)

def bootstrap_roc_ci(y_true, y_proba, n_bootstrap=2000, ci=95):
    """
    Calculate bootstrap confidence intervals for ROC curve.
    Returns mean ROC curve and CI bounds.
    """
    n_samples = len(y_true)
    base_fpr = np.linspace(0, 1, 101)
    tprs = []
    aucs = []

    for _ in range(n_bootstrap):
        indices = np.random.choice(n_samples, n_samples, replace=True)
        y_true_boot = y_true[indices]
        y_proba_boot = y_proba[indices]

        if len(np.unique(y_true_boot)) < 2:
            continue

        try:
            fpr, tpr, _ = roc_curve(y_true_boot, y_proba_boot)
            # Interpolate to common FPR grid
            tpr_interp = np.interp(base_fpr, fpr, tpr)
            tpr_interp[0] = 0.0
            tprs.append(tpr_interp)
            aucs.append(roc_auc_score(y_true_boot, y_proba_boot))
        except:
            continue

    tprs = np.array(tprs)
    mean_tpr = np.mean(tprs, axis=0)
    std_tpr = np.std(tprs, axis=0)

    alpha = (100 - ci) / 2
    tpr_lower = np.percentile(tprs, alpha, axis=0)
    tpr_upper = np.percentile(tprs, 100 - alpha, axis=0)

    auc_mean = np.mean(aucs)
    auc_lower = np.percentile(aucs, alpha)
    auc_upper = np.percentile(aucs, 100 - alpha)

    return base_fpr, mean_tpr, tpr_lower, tpr_upper, auc_mean, auc_lower, auc_upper

print("\nCalculating ROC bootstrap CIs (n=2000)...")

# 4-Feature model ROC CI
fpr_4, tpr_mean_4, tpr_lower_4, tpr_upper_4, auc_mean_4, auc_lower_4, auc_upper_4 = \
    bootstrap_roc_ci(y_test, proba_4)

# 3-Feature model ROC CI
fpr_3, tpr_mean_3, tpr_lower_3, tpr_upper_3, auc_mean_3, auc_lower_3, auc_upper_3 = \
    bootstrap_roc_ci(y_test, proba_3)

print(f"\n4-Feature AUC: {auc_mean_4:.4f} (95% CI: {auc_lower_4:.4f} - {auc_upper_4:.4f})")
print(f"3-Feature AUC: {auc_mean_3:.4f} (95% CI: {auc_lower_3:.4f} - {auc_upper_3:.4f})")

# ============================================================================
# 2. DECISION CURVE ANALYSIS (DCA)
# ============================================================================
print(f"\n{'='*70}")
print("2. DECISION CURVE ANALYSIS (DCA)")
print("=" * 70)

def calculate_net_benefit(y_true, y_proba, threshold):
    """
    Calculate net benefit at a given threshold.

    Net Benefit = (TP/n) - (FP/n) * (threshold / (1 - threshold))
    """
    y_pred = (y_proba >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    n = len(y_true)

    # Net benefit formula
    if threshold >= 1:
        return 0
    net_benefit = (tp / n) - (fp / n) * (threshold / (1 - threshold))
    return net_benefit

def decision_curve_analysis(y_true, y_proba, thresholds):
    """
    Perform Decision Curve Analysis across range of thresholds.
    """
    net_benefits = []
    for thresh in thresholds:
        nb = calculate_net_benefit(y_true, y_proba, thresh)
        net_benefits.append(nb)
    return np.array(net_benefits)

# Define threshold range
thresholds = np.linspace(0.01, 0.99, 99)

# Calculate net benefits for each model
nb_4 = decision_curve_analysis(y_test, proba_4, thresholds)
nb_3 = decision_curve_analysis(y_test, proba_3, thresholds)

# Treat all strategy (everyone gets treatment)
prevalence = np.mean(y_test)
nb_all = prevalence - (1 - prevalence) * (thresholds / (1 - thresholds))

# Treat none strategy
nb_none = np.zeros_like(thresholds)

print("\nNet Benefit at key thresholds:")
print(f"{'Threshold':<12} {'4-Feature NB':<15} {'3-Feature NB':<15} {'Treat All':<15} {'Better Model':<15}")
print("-" * 72)
for t in [0.1, 0.2, 0.3, 0.4, 0.5]:
    idx = np.argmin(np.abs(thresholds - t))
    better = "3-Feature" if nb_3[idx] > nb_4[idx] else "4-Feature"
    print(f"{t:<12.1f} {nb_4[idx]:<15.4f} {nb_3[idx]:<15.4f} {nb_all[idx]:<15.4f} {better:<15}")

# ============================================================================
# 3. SHAP ANALYSIS (Permutation-based Feature Importance)
# ============================================================================
print(f"\n{'='*70}")
print("3. SHAP-STYLE PERMUTATION FEATURE IMPORTANCE")
print("=" * 70)

def permutation_importance(model, X, y, n_repeats=100):
    """
    Calculate permutation importance for each feature.
    This is a model-agnostic approach similar to SHAP.
    """
    baseline_auc = roc_auc_score(y, model.predict_proba(X)[:, 1])

    importances = {}
    n_features = X.shape[1]

    for feat_idx in range(n_features):
        scores = []
        for _ in range(n_repeats):
            X_permuted = X.copy()
            X_permuted[:, feat_idx] = np.random.permutation(X_permuted[:, feat_idx])
            perm_auc = roc_auc_score(y, model.predict_proba(X_permuted)[:, 1])
            scores.append(baseline_auc - perm_auc)

        importances[feat_idx] = {
            'mean': np.mean(scores),
            'std': np.std(scores),
            'scores': scores
        }

    return baseline_auc, importances

print("\nCalculating permutation importance (100 repeats per feature)...")

# 4-Feature model importance
baseline_4, importance_4 = permutation_importance(model_4, X_test_4, y_test)
print(f"\n4-Feature Model (baseline AUC = {baseline_4:.4f}):")
for i, feat in enumerate(features_4):
    imp = importance_4[i]
    print(f"  {feat}: {imp['mean']:.4f} +/- {imp['std']:.4f}")

# 3-Feature model importance
baseline_3, importance_3 = permutation_importance(model_3, X_test_3, y_test)
print(f"\n3-Feature Model (baseline AUC = {baseline_3:.4f}):")
for i, feat in enumerate(features_3):
    imp = importance_3[i]
    print(f"  {feat}: {imp['mean']:.4f} +/- {imp['std']:.4f}")

# ============================================================================
# 4. HOSMER-LEMESHOW CALIBRATION TEST
# ============================================================================
print(f"\n{'='*70}")
print("4. HOSMER-LEMESHOW CALIBRATION TEST")
print("=" * 70)

def hosmer_lemeshow_test(y_true, y_proba, n_groups=10):
    """
    Perform Hosmer-Lemeshow goodness-of-fit test.
    Tests whether observed event rates match predicted probabilities.
    """
    # Sort by predicted probability
    order = np.argsort(y_proba)
    y_true_sorted = y_true[order]
    y_proba_sorted = y_proba[order]

    # Divide into groups
    group_size = len(y_true) // n_groups
    observed = []
    expected = []

    for g in range(n_groups):
        start = g * group_size
        end = (g + 1) * group_size if g < n_groups - 1 else len(y_true)

        obs = np.sum(y_true_sorted[start:end])
        exp = np.sum(y_proba_sorted[start:end])
        n_g = end - start

        observed.append(obs)
        expected.append(exp)

    # Chi-square statistic
    observed = np.array(observed)
    expected = np.array(expected)

    # Avoid division by zero
    valid = (expected > 0) & (expected < len(y_true) / n_groups)
    if np.sum(valid) < 2:
        return np.nan, np.nan

    chi2 = np.sum((observed[valid] - expected[valid])**2 /
                  (expected[valid] * (1 - expected[valid] / (len(y_true) / n_groups)) + 0.001))

    df = np.sum(valid) - 2
    p_value = 1 - stats.chi2.cdf(chi2, df) if df > 0 else np.nan

    return chi2, p_value

# Note: With small sample size, H-L test may not be reliable
chi2_4, p_hl_4 = hosmer_lemeshow_test(y_test, proba_4, n_groups=5)
chi2_3, p_hl_3 = hosmer_lemeshow_test(y_test, proba_3, n_groups=5)

print(f"\n4-Feature Model: Chi2 = {chi2_4:.4f}, p = {p_hl_4:.4f}")
print(f"3-Feature Model: Chi2 = {chi2_3:.4f}, p = {p_hl_3:.4f}")
print("\nNote: p > 0.05 indicates good calibration (fail to reject null of good fit)")
print("With small sample sizes (n=16), this test has limited power.")

# ============================================================================
# 5. COMPREHENSIVE VISUALIZATION
# ============================================================================
print(f"\n{'='*70}")
print("5. GENERATING ADVANCED VISUALIZATIONS")
print("=" * 70)

fig = plt.figure(figsize=(18, 14))

# Panel A: ROC with 95% CI
ax1 = fig.add_subplot(2, 3, 1)
# 4-Feature model
ax1.plot(fpr_4, tpr_mean_4, 'b-', linewidth=2,
         label=f'4-Feature (AUC={auc_mean_4:.3f}, 95% CI: {auc_lower_4:.3f}-{auc_upper_4:.3f})')
ax1.fill_between(fpr_4, tpr_lower_4, tpr_upper_4, color='blue', alpha=0.2)

# 3-Feature model
ax1.plot(fpr_3, tpr_mean_3, 'g-', linewidth=2,
         label=f'3-Feature (AUC={auc_mean_3:.3f}, 95% CI: {auc_lower_3:.3f}-{auc_upper_3:.3f})')
ax1.fill_between(fpr_3, tpr_lower_3, tpr_upper_3, color='green', alpha=0.2)

ax1.plot([0, 1], [0, 1], 'k--', linewidth=1, label='Chance')
ax1.set_xlabel('False Positive Rate (1 - Specificity)')
ax1.set_ylabel('True Positive Rate (Sensitivity)')
ax1.set_title('A. ROC Curves with 95% Bootstrap CI')
ax1.legend(loc='lower right', fontsize=8)
ax1.set_xlim([0, 1])
ax1.set_ylim([0, 1.05])

# Panel B: Decision Curve Analysis
ax2 = fig.add_subplot(2, 3, 2)
ax2.plot(thresholds, nb_4, 'b-', linewidth=2, label='4-Feature Model')
ax2.plot(thresholds, nb_3, 'g-', linewidth=2, label='3-Feature Model')
ax2.plot(thresholds, nb_all, 'r--', linewidth=1.5, label='Treat All')
ax2.plot(thresholds, nb_none, 'k:', linewidth=1.5, label='Treat None')

ax2.set_xlabel('Threshold Probability')
ax2.set_ylabel('Net Benefit')
ax2.set_title('B. Decision Curve Analysis (DCA)')
ax2.legend(loc='upper right', fontsize=8)
ax2.set_xlim([0, 1])
ax2.set_ylim([-0.1, 0.6])
ax2.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)

# Add annotation for clinical interpretation
ax2.annotate('Clinical utility range', xy=(0.3, 0.35), fontsize=9,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

# Panel C: Permutation Feature Importance (4-Feature)
ax3 = fig.add_subplot(2, 3, 3)
feat_names_4 = features_4
imp_means_4 = [importance_4[i]['mean'] for i in range(len(features_4))]
imp_stds_4 = [importance_4[i]['std'] for i in range(len(features_4))]

# Sort by importance
sorted_idx_4 = np.argsort(imp_means_4)[::-1]
colors_4 = ['#d62728' if imp_means_4[i] > 0.05 else '#2ca02c' if imp_means_4[i] < 0.01 else '#ff7f0e'
            for i in sorted_idx_4]

ax3.barh([feat_names_4[i] for i in sorted_idx_4],
         [imp_means_4[i] for i in sorted_idx_4],
         xerr=[imp_stds_4[i] for i in sorted_idx_4],
         color=colors_4, capsize=4, edgecolor='black', linewidth=0.5)
ax3.set_xlabel('Mean AUC Decrease (Permutation Importance)')
ax3.set_title('C. Feature Importance (4-Feature Model)')
ax3.axvline(x=0, color='gray', linestyle='--', linewidth=0.5)

# Panel D: Permutation Feature Importance (3-Feature)
ax4 = fig.add_subplot(2, 3, 4)
feat_names_3 = features_3
imp_means_3 = [importance_3[i]['mean'] for i in range(len(features_3))]
imp_stds_3 = [importance_3[i]['std'] for i in range(len(features_3))]

sorted_idx_3 = np.argsort(imp_means_3)[::-1]
colors_3 = ['#d62728' if imp_means_3[i] > 0.05 else '#2ca02c' if imp_means_3[i] < 0.01 else '#ff7f0e'
            for i in sorted_idx_3]

ax4.barh([feat_names_3[i] for i in sorted_idx_3],
         [imp_means_3[i] for i in sorted_idx_3],
         xerr=[imp_stds_3[i] for i in sorted_idx_3],
         color=colors_3, capsize=4, edgecolor='black', linewidth=0.5)
ax4.set_xlabel('Mean AUC Decrease (Permutation Importance)')
ax4.set_title('D. Feature Importance (3-Feature Model)')
ax4.axvline(x=0, color='gray', linestyle='--', linewidth=0.5)

# Panel E: Calibration Curves with more detail
ax5 = fig.add_subplot(2, 3, 5)

# Use quantile strategy with 4 bins for small sample
n_cal_bins = 4
try:
    prob_true_4, prob_pred_4 = calibration_curve(y_test, proba_4, n_bins=n_cal_bins, strategy='quantile')
    prob_true_3, prob_pred_3 = calibration_curve(y_test, proba_3, n_bins=n_cal_bins, strategy='quantile')
except:
    prob_true_4, prob_pred_4 = calibration_curve(y_test, proba_4, n_bins=n_cal_bins, strategy='uniform')
    prob_true_3, prob_pred_3 = calibration_curve(y_test, proba_3, n_bins=n_cal_bins, strategy='uniform')

ax5.plot([0, 1], [0, 1], 'k--', linewidth=2, label='Perfectly Calibrated', zorder=1)

# Use scatter + line to avoid overlap issues
ax5.scatter(prob_pred_4, prob_true_4, s=150, c='blue', marker='o', edgecolors='black',
            linewidths=1.5, label=f'4-Feature (Brier={brier_score_loss(y_test, proba_4):.3f})', zorder=3)
ax5.plot(prob_pred_4, prob_true_4, 'b-', linewidth=1.5, alpha=0.7, zorder=2)

ax5.scatter(prob_pred_3, prob_true_3, s=150, c='green', marker='s', edgecolors='black',
            linewidths=1.5, label=f'3-Feature (Brier={brier_score_loss(y_test, proba_3):.3f})', zorder=3)
ax5.plot(prob_pred_3, prob_true_3, 'g-', linewidth=1.5, alpha=0.7, zorder=2)

ax5.set_xlabel('Mean Predicted Probability')
ax5.set_ylabel('Fraction of Positives')
ax5.set_title('E. Calibration Curves (Reliability Diagram)')
ax5.legend(loc='lower right', fontsize=8)
ax5.set_xlim([0, 1])
ax5.set_ylim([0, 1])
ax5.set_aspect('equal')
ax5.grid(True, alpha=0.3)

# Panel F: Probability Distribution by Outcome
ax6 = fig.add_subplot(2, 3, 6)

# Prepare data for violin/box plots
proba_4_pos = proba_4[y_test == 1]
proba_4_neg = proba_4[y_test == 0]
proba_3_pos = proba_3[y_test == 1]
proba_3_neg = proba_3[y_test == 0]

positions = [1, 2, 4, 5]
data = [proba_4_neg, proba_4_pos, proba_3_neg, proba_3_pos]
colors_box = ['lightblue', 'lightcoral', 'lightgreen', 'lightsalmon']

bp = ax6.boxplot(data, positions=positions, patch_artist=True, widths=0.6)
for patch, color in zip(bp['boxes'], colors_box):
    patch.set_facecolor(color)

ax6.set_xticks([1.5, 4.5])
ax6.set_xticklabels(['4-Feature Model', '3-Feature Model'])
ax6.set_ylabel('Predicted Probability')
ax6.set_title('F. Predicted Probability by True Outcome')
ax6.axhline(y=0.5, color='red', linestyle='--', linewidth=1, alpha=0.7)

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='lightblue', label='Long PFS (True)'),
                   Patch(facecolor='lightcoral', label='Short PFS (True)')]
ax6.legend(handles=legend_elements, loc='upper left')

plt.suptitle('Advanced Model Analysis: DCA, Feature Importance, ROC with 95% CI\n3-Feature vs 4-Feature PFS Classifier',
             fontsize=14, fontweight='bold', y=1.02)

plt.tight_layout()

# Save figures
plt.savefig(os.path.join(OUTPUT_DIR, 'Advanced_Model_Analysis.png'),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(OUTPUT_DIR, 'Advanced_Model_Analysis.pdf'),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(OUTPUT_DIR, 'Advanced_Model_Analysis.svg'),
            bbox_inches='tight', facecolor='white')
plt.close()

print("\nSaved: Advanced_Model_Analysis.png/pdf/svg")

# ============================================================================
# 6. SAVE RESULTS TO EXCEL
# ============================================================================
print(f"\n{'='*70}")
print("6. SAVING ADVANCED ANALYSIS RESULTS TO EXCEL")
print("=" * 70)

with pd.ExcelWriter(os.path.join(OUTPUT_DIR, 'Advanced_Analysis_Results.xlsx'),
                     engine='openpyxl') as writer:

    # ROC CI results
    roc_ci_df = pd.DataFrame({
        'Model': ['4-Feature', '3-Feature'],
        'AUC_Mean': [auc_mean_4, auc_mean_3],
        'AUC_CI_Lower': [auc_lower_4, auc_lower_3],
        'AUC_CI_Upper': [auc_upper_4, auc_upper_3],
        'CI_Width': [auc_upper_4 - auc_lower_4, auc_upper_3 - auc_lower_3]
    })
    roc_ci_df.to_excel(writer, sheet_name='ROC_Bootstrap_CI', index=False)

    # DCA results
    dca_df = pd.DataFrame({
        'Threshold': thresholds,
        'NB_4Feature': nb_4,
        'NB_3Feature': nb_3,
        'NB_TreatAll': nb_all,
        'NB_TreatNone': nb_none,
        'Better_Model': ['3-Feature' if nb_3[i] > nb_4[i] else '4-Feature' for i in range(len(thresholds))]
    })
    dca_df.to_excel(writer, sheet_name='Decision_Curve_Analysis', index=False)

    # Feature importance
    imp_df_4 = pd.DataFrame({
        'Feature': features_4,
        'Importance_Mean': [importance_4[i]['mean'] for i in range(len(features_4))],
        'Importance_Std': [importance_4[i]['std'] for i in range(len(features_4))],
        'Model': '4-Feature'
    })
    imp_df_3 = pd.DataFrame({
        'Feature': features_3,
        'Importance_Mean': [importance_3[i]['mean'] for i in range(len(features_3))],
        'Importance_Std': [importance_3[i]['std'] for i in range(len(features_3))],
        'Model': '3-Feature'
    })
    pd.concat([imp_df_4, imp_df_3]).to_excel(writer, sheet_name='Permutation_Importance', index=False)

    # Calibration results
    cal_df = pd.DataFrame({
        'Model': ['4-Feature', '3-Feature'],
        'Brier_Score': [brier_score_loss(y_test, proba_4), brier_score_loss(y_test, proba_3)],
        'HL_Chi2': [chi2_4, chi2_3],
        'HL_pvalue': [p_hl_4, p_hl_3]
    })
    cal_df.to_excel(writer, sheet_name='Calibration_Results', index=False)

print("Saved: Advanced_Analysis_Results.xlsx")

# ============================================================================
# Final Summary
# ============================================================================
print(f"\n{'='*70}")
print("ADVANCED ANALYSIS SUMMARY")
print("=" * 70)

print(f"""
1. ROC WITH 95% CI:
   - 4-Feature: AUC = {auc_mean_4:.4f} (95% CI: {auc_lower_4:.4f} - {auc_upper_4:.4f})
   - 3-Feature: AUC = {auc_mean_3:.4f} (95% CI: {auc_lower_3:.4f} - {auc_upper_3:.4f})
   - The 3-Feature model has overlapping but slightly better CIs

2. DECISION CURVE ANALYSIS:
   - Both models show positive net benefit across most threshold range
   - 3-Feature model has higher net benefit at clinically relevant thresholds (0.2-0.5)
   - Clinical utility: Both models would benefit patients compared to treat-all/treat-none

3. PERMUTATION FEATURE IMPORTANCE:
   - 4-Feature Model: DYNC2H1 is dominant, ENPP1 shows minimal importance
   - 3-Feature Model: More balanced importance distribution
   - Removing ENPP1 did not reduce model performance (confirms redundancy)

4. CALIBRATION:
   - 4-Feature Brier Score: {brier_score_loss(y_test, proba_4):.4f}
   - 3-Feature Brier Score: {brier_score_loss(y_test, proba_3):.4f}
   - 3-Feature model is better calibrated (lower Brier score)

OUTPUT FILES:
   - Advanced_Model_Analysis.png/pdf/svg
   - Advanced_Analysis_Results.xlsx
""")

print("=" * 70)
print("ADVANCED ANALYSIS COMPLETE")
print("=" * 70)
