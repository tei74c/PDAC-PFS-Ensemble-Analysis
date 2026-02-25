"""
Individual Metric Comparison Panels
=====================================
Generates individual comparison panels for 3-feature vs 4-feature models.

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
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# Path Configuration
# ============================================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'data')
OUTPUT_DIR = os.path.join(BASE_DIR, 'results', 'Individual_Comparison_Figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Set style for publication-quality figures
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

np.random.seed(42)

# ============================================================================
# Load Data and Models
# ============================================================================
print("=" * 70)
print("GENERATING INDIVIDUAL HEAD-TO-HEAD COMPARISON FIGURES")
print("=" * 70)

# Load models
with open(os.path.join(DATA_DIR, 'champion_model_NB_4features.pkl'), 'rb') as f:
    model_4_data = pickle.load(f)
model_4 = model_4_data['model']
features_4 = model_4_data['features']

with open(os.path.join(DATA_DIR, 'champion_model_NB_3features.pkl'), 'rb') as f:
    model_3_data = pickle.load(f)
model_3 = model_3_data['model']
features_3 = model_3_data['features']

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
pred_4 = model_4.predict(X_test_4)
proba_4 = model_4.predict_proba(X_test_4)[:, 1]
pred_3 = model_3.predict(X_test_3)
proba_3 = model_3.predict_proba(X_test_3)[:, 1]

print(f"\n4-Feature Model: {features_4}")
print(f"3-Feature Model: {features_3}")
print(f"Test samples: {len(y_test)}")

# ============================================================================
# ECE Function
# ============================================================================
def expected_calibration_error(y_true, y_proba, n_bins=10):
    """Calculate Expected Calibration Error (ECE)."""
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    for i in range(n_bins):
        in_bin = (y_proba > bin_boundaries[i]) & (y_proba <= bin_boundaries[i + 1])
        prop_in_bin = np.mean(in_bin)
        if prop_in_bin > 0:
            avg_confidence = np.mean(y_proba[in_bin])
            avg_accuracy = np.mean(y_true[in_bin])
            ece += np.abs(avg_accuracy - avg_confidence) * prop_in_bin
    return ece

# ============================================================================
# Bootstrap CI Function
# ============================================================================
def bootstrap_metric(y_true, y_pred, y_proba, metric_func, n_bootstrap=5000, ci=95):
    """Bootstrap confidence interval for a single metric."""
    n_samples = len(y_true)
    boot_values = []

    for _ in range(n_bootstrap):
        indices = np.random.choice(n_samples, n_samples, replace=True)
        y_true_boot = y_true[indices]
        y_pred_boot = y_pred[indices]
        y_proba_boot = y_proba[indices]

        if len(np.unique(y_true_boot)) < 2:
            continue

        try:
            if metric_func.__name__ in ['roc_auc_score', 'average_precision_score', 'brier_score_loss']:
                val = metric_func(y_true_boot, y_proba_boot)
            elif metric_func.__name__ == 'expected_calibration_error':
                val = metric_func(y_true_boot, y_proba_boot)
            else:
                val = metric_func(y_true_boot, y_pred_boot)
            boot_values.append(val)
        except:
            continue

    alpha = (100 - ci) / 2
    return {
        'mean': np.mean(boot_values),
        'std': np.std(boot_values),
        'ci_lower': np.percentile(boot_values, alpha),
        'ci_upper': np.percentile(boot_values, 100 - alpha),
        'values': boot_values
    }

# ============================================================================
# Calculate All Metrics with Bootstrap CI
# ============================================================================
print("\nCalculating bootstrap confidence intervals (n=5000)...")

metrics_config = {
    'Accuracy': {'func': accuracy_score, 'higher_better': True, 'unit': '%', 'multiply': 100},
    'F1-Score': {'func': f1_score, 'higher_better': True, 'unit': '', 'multiply': 1},
    'ROC-AUC': {'func': roc_auc_score, 'higher_better': True, 'unit': '', 'multiply': 1},
    'Precision': {'func': precision_score, 'higher_better': True, 'unit': '%', 'multiply': 100},
    'Recall': {'func': recall_score, 'higher_better': True, 'unit': '%', 'multiply': 100},
    'Avg-Precision': {'func': average_precision_score, 'higher_better': True, 'unit': '', 'multiply': 1},
    'MCC': {'func': matthews_corrcoef, 'higher_better': True, 'unit': '', 'multiply': 1},
    'Cohen-Kappa': {'func': cohen_kappa_score, 'higher_better': True, 'unit': '', 'multiply': 1},
    'Brier-Score': {'func': brier_score_loss, 'higher_better': False, 'unit': '', 'multiply': 1},
    'ECE': {'func': expected_calibration_error, 'higher_better': False, 'unit': '', 'multiply': 1},
}

results = {}
for metric_name, config in metrics_config.items():
    print(f"  Bootstrapping {metric_name}...")
    results[metric_name] = {
        '4-Feature': bootstrap_metric(y_test, pred_4, proba_4, config['func']),
        '3-Feature': bootstrap_metric(y_test, pred_3, proba_3, config['func']),
        'config': config
    }

# ============================================================================
# Generate Individual Metric Figures
# ============================================================================
print(f"\n{'='*70}")
print("GENERATING INDIVIDUAL METRIC FIGURES")
print("=" * 70)

def create_metric_figure(metric_name, data, config, filename):
    """Create a single metric comparison figure with 95% CI."""
    fig, ax = plt.subplots(figsize=(8, 6))

    val_4 = data['4-Feature']['mean'] * config['multiply']
    val_3 = data['3-Feature']['mean'] * config['multiply']
    ci_4_lower = data['4-Feature']['ci_lower'] * config['multiply']
    ci_4_upper = data['4-Feature']['ci_upper'] * config['multiply']
    ci_3_lower = data['3-Feature']['ci_lower'] * config['multiply']
    ci_3_upper = data['3-Feature']['ci_upper'] * config['multiply']

    # Determine winner
    if config['higher_better']:
        winner = '3-Feature' if val_3 > val_4 else '4-Feature'
    else:
        winner = '3-Feature' if val_3 < val_4 else '4-Feature'

    # Colors
    color_4 = '#1f77b4'  # Blue
    color_3 = '#2ca02c'  # Green
    winner_color = '#ffd700'  # Gold highlight

    # Bar positions
    x = np.array([0, 1])
    width = 0.5

    # Error bars
    err_4 = [[val_4 - ci_4_lower], [ci_4_upper - val_4]]
    err_3 = [[val_3 - ci_3_lower], [ci_3_upper - val_3]]

    # Plot bars
    bar1 = ax.bar(0, val_4, width, color=color_4, edgecolor='black', linewidth=1.5,
                  yerr=err_4, capsize=8, error_kw={'linewidth': 2, 'capthick': 2})
    bar2 = ax.bar(1, val_3, width, color=color_3, edgecolor='black', linewidth=1.5,
                  yerr=err_3, capsize=8, error_kw={'linewidth': 2, 'capthick': 2})

    # Highlight winner
    if winner == '4-Feature':
        bar1[0].set_edgecolor(winner_color)
        bar1[0].set_linewidth(3)
    else:
        bar2[0].set_edgecolor(winner_color)
        bar2[0].set_linewidth(3)

    # Labels
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['4-Feature Model\n(Original)', '3-Feature Model\n(Optimized)'], fontsize=12)

    unit_str = config['unit']
    ax.set_ylabel(f'{metric_name} {unit_str}', fontsize=14)

    # Add value labels on bars
    for i, (val, ci_l, ci_u) in enumerate([(val_4, ci_4_lower, ci_4_upper),
                                            (val_3, ci_3_lower, ci_3_upper)]):
        if config['multiply'] == 100:
            label = f'{val:.1f}%\n[{ci_l:.1f}-{ci_u:.1f}]'
        else:
            label = f'{val:.3f}\n[{ci_l:.3f}-{ci_u:.3f}]'
        ax.text(i, val + (ci_u - val) + 0.02 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
                label, ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Title with winner indication
    better_text = "Higher is better" if config['higher_better'] else "Lower is better"
    ax.set_title(f'{metric_name} Comparison\n({better_text})', fontsize=14, fontweight='bold')

    # Add winner annotation
    delta = val_3 - val_4
    if config['multiply'] == 100:
        delta_text = f'Delta: {delta:+.1f}%'
    else:
        delta_text = f'Delta: {delta:+.4f}'

    ax.text(0.5, 0.02, f'Winner: {winner} | {delta_text}',
            transform=ax.transAxes, ha='center', fontsize=11,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # Set y-axis limits
    all_vals = [val_4, val_3, ci_4_lower, ci_4_upper, ci_3_lower, ci_3_upper]
    y_min = min(all_vals) * 0.9 if min(all_vals) > 0 else min(all_vals) * 1.1
    y_max = max(all_vals) * 1.15
    ax.set_ylim(y_min, y_max)

    # Add reference line for metrics where 0.5 is chance
    if metric_name in ['Accuracy', 'ROC-AUC', 'Avg-Precision'] and config['multiply'] == 1:
        ax.axhline(y=0.5, color='red', linestyle='--', linewidth=1, alpha=0.7, label='Chance')

    plt.tight_layout()

    # Save
    for ext in ['png', 'pdf', 'svg']:
        plt.savefig(os.path.join(OUTPUT_DIR, f'{filename}.{ext}'),
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    return winner

# Generate figures for each metric
for metric_name, data in results.items():
    filename = f'Comparison_{metric_name.replace("-", "_").replace(" ", "_")}'
    winner = create_metric_figure(metric_name, data, data['config'], filename)
    print(f"  {metric_name}: Winner = {winner} -> {filename}.png")

# ============================================================================
# ROC Curves with 95% CI
# ============================================================================
print(f"\n{'='*70}")
print("GENERATING ROC CURVES WITH 95% CI")
print("=" * 70)

def bootstrap_roc_ci(y_true, y_proba, n_bootstrap=2000):
    """Calculate bootstrap CI for ROC curve."""
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
            tpr_interp = np.interp(base_fpr, fpr, tpr)
            tpr_interp[0] = 0.0
            tprs.append(tpr_interp)
            aucs.append(roc_auc_score(y_true_boot, y_proba_boot))
        except:
            continue

    tprs = np.array(tprs)
    return {
        'fpr': base_fpr,
        'tpr_mean': np.mean(tprs, axis=0),
        'tpr_lower': np.percentile(tprs, 2.5, axis=0),
        'tpr_upper': np.percentile(tprs, 97.5, axis=0),
        'auc_mean': np.mean(aucs),
        'auc_lower': np.percentile(aucs, 2.5),
        'auc_upper': np.percentile(aucs, 97.5)
    }

print("  Calculating ROC bootstrap CIs...")
roc_4 = bootstrap_roc_ci(y_test, proba_4)
roc_3 = bootstrap_roc_ci(y_test, proba_3)

# Combined ROC with CI
fig, ax = plt.subplots(figsize=(10, 8))

# 4-Feature model
ax.plot(roc_4['fpr'], roc_4['tpr_mean'], 'b-', linewidth=2.5,
        label=f'4-Feature (AUC={roc_4["auc_mean"]:.3f}, 95% CI: {roc_4["auc_lower"]:.3f}-{roc_4["auc_upper"]:.3f})')
ax.fill_between(roc_4['fpr'], roc_4['tpr_lower'], roc_4['tpr_upper'],
                color='blue', alpha=0.2, label='4-Feature 95% CI')

# 3-Feature model
ax.plot(roc_3['fpr'], roc_3['tpr_mean'], 'g-', linewidth=2.5,
        label=f'3-Feature (AUC={roc_3["auc_mean"]:.3f}, 95% CI: {roc_3["auc_lower"]:.3f}-{roc_3["auc_upper"]:.3f})')
ax.fill_between(roc_3['fpr'], roc_3['tpr_lower'], roc_3['tpr_upper'],
                color='green', alpha=0.2, label='3-Feature 95% CI')

# Diagonal
ax.plot([0, 1], [0, 1], 'k--', linewidth=1.5, label='Chance (AUC=0.5)')

ax.set_xlabel('False Positive Rate (1 - Specificity)', fontsize=14)
ax.set_ylabel('True Positive Rate (Sensitivity)', fontsize=14)
ax.set_title('ROC Curves with 95% Bootstrap Confidence Intervals\n4-Feature vs 3-Feature Model Comparison',
             fontsize=14, fontweight='bold')
ax.legend(loc='lower right', fontsize=10)
ax.set_xlim([0, 1])
ax.set_ylim([0, 1.02])
ax.grid(True, alpha=0.3)

plt.tight_layout()
for ext in ['png', 'pdf', 'svg']:
    plt.savefig(os.path.join(OUTPUT_DIR, f'ROC_Comparison_with_95CI.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: ROC_Comparison_with_95CI.png/pdf/svg")

# Individual ROC curves (side by side)
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, roc_data, name, color in [(axes[0], roc_4, '4-Feature Model (Original)', 'blue'),
                                   (axes[1], roc_3, '3-Feature Model (Optimized)', 'green')]:
    ax.plot(roc_data['fpr'], roc_data['tpr_mean'], color=color, linewidth=2.5,
            label=f'AUC = {roc_data["auc_mean"]:.3f}')
    ax.fill_between(roc_data['fpr'], roc_data['tpr_lower'], roc_data['tpr_upper'],
                    color=color, alpha=0.3,
                    label=f'95% CI: [{roc_data["auc_lower"]:.3f}, {roc_data["auc_upper"]:.3f}]')
    ax.plot([0, 1], [0, 1], 'k--', linewidth=1.5, label='Chance')

    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title(name, fontsize=13, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1.02])
    ax.grid(True, alpha=0.3)

plt.suptitle('Individual ROC Curves with 95% Bootstrap CI',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout(rect=[0, 0, 1, 0.96])
for ext in ['png', 'pdf', 'svg']:
    plt.savefig(os.path.join(OUTPUT_DIR, f'ROC_Individual_with_95CI.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: ROC_Individual_with_95CI.png/pdf/svg")

# ============================================================================
# Confusion Matrices
# ============================================================================
print(f"\n{'='*70}")
print("GENERATING CONFUSION MATRICES")
print("=" * 70)

cm_4 = confusion_matrix(y_test, pred_4)
cm_3 = confusion_matrix(y_test, pred_3)

# Calculate metrics from CM
def cm_metrics(cm):
    tn, fp, fn, tp = cm.ravel()
    sensitivity = tp / (tp + fn) * 100
    specificity = tn / (tn + fp) * 100
    ppv = tp / (tp + fp) * 100 if (tp + fp) > 0 else 0
    npv = tn / (tn + fn) * 100 if (tn + fn) > 0 else 0
    return {'TN': tn, 'FP': fp, 'FN': fn, 'TP': tp,
            'Sensitivity': sensitivity, 'Specificity': specificity,
            'PPV': ppv, 'NPV': npv}

cm_metrics_4 = cm_metrics(cm_4)
cm_metrics_3 = cm_metrics(cm_3)

# Side-by-side confusion matrices
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, cm, metrics, name, cmap in [(axes[0], cm_4, cm_metrics_4, '4-Feature Model', 'Blues'),
                                     (axes[1], cm_3, cm_metrics_3, '3-Feature Model', 'Greens')]:

    # Create annotated heatmap
    sns.heatmap(cm, annot=False, cmap=cmap, ax=ax, cbar=False,
                linewidths=2, linecolor='white', square=True)

    # Add custom annotations
    labels = [['TN', 'FP'], ['FN', 'TP']]
    for i in range(2):
        for j in range(2):
            val = cm[i, j]
            label = labels[i][j]
            color = 'white' if val > cm.max() / 2 else 'black'
            ax.text(j + 0.5, i + 0.5, f'{label}\n{val}',
                    ha='center', va='center', fontsize=16, fontweight='bold', color=color)

    ax.set_xlabel('Predicted Label', fontsize=12)
    ax.set_ylabel('True Label', fontsize=12)
    ax.set_xticklabels(['Long PFS', 'Short PFS'], fontsize=11)
    ax.set_yticklabels(['Long PFS', 'Short PFS'], fontsize=11, rotation=0)

    # Add metrics below
    metrics_text = (f"Sens: {metrics['Sensitivity']:.1f}% | "
                    f"Spec: {metrics['Specificity']:.1f}%\n"
                    f"PPV: {metrics['PPV']:.1f}% | NPV: {metrics['NPV']:.1f}%")
    ax.text(0.5, -0.15, metrics_text, transform=ax.transAxes, ha='center', fontsize=11,
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    ax.set_title(name, fontsize=13, fontweight='bold')

plt.suptitle('Confusion Matrix Comparison', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout(rect=[0, 0, 1, 0.96])
for ext in ['png', 'pdf', 'svg']:
    plt.savefig(os.path.join(OUTPUT_DIR, f'Confusion_Matrix_Comparison.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: Confusion_Matrix_Comparison.png/pdf/svg")

# Individual CMs
for cm, metrics, name, cmap, fname in [(cm_4, cm_metrics_4, '4-Feature Model (Original)', 'Blues', 'CM_4Feature'),
                                        (cm_3, cm_metrics_3, '3-Feature Model (Optimized)', 'Greens', 'CM_3Feature')]:
    fig, ax = plt.subplots(figsize=(8, 7))

    sns.heatmap(cm, annot=False, cmap=cmap, ax=ax, cbar=True,
                linewidths=2, linecolor='white', square=True)

    labels = [['TN', 'FP'], ['FN', 'TP']]
    for i in range(2):
        for j in range(2):
            val = cm[i, j]
            label = labels[i][j]
            color = 'white' if val > cm.max() / 2 else 'black'
            ax.text(j + 0.5, i + 0.5, f'{label}\n{val}',
                    ha='center', va='center', fontsize=20, fontweight='bold', color=color)

    ax.set_xlabel('Predicted Label', fontsize=14)
    ax.set_ylabel('True Label', fontsize=14)
    ax.set_xticklabels(['Long PFS', 'Short PFS'], fontsize=12)
    ax.set_yticklabels(['Long PFS', 'Short PFS'], fontsize=12, rotation=0)
    ax.set_title(f'Confusion Matrix\n{name}', fontsize=14, fontweight='bold')

    # Metrics table
    metrics_text = (f"Sensitivity: {metrics['Sensitivity']:.1f}%\n"
                    f"Specificity: {metrics['Specificity']:.1f}%\n"
                    f"PPV: {metrics['PPV']:.1f}%\n"
                    f"NPV: {metrics['NPV']:.1f}%")
    ax.text(1.15, 0.5, metrics_text, transform=ax.transAxes, fontsize=12,
            verticalalignment='center',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

    plt.tight_layout()
    for ext in ['png', 'pdf', 'svg']:
        plt.savefig(os.path.join(OUTPUT_DIR, f'{fname}.{ext}'),
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {fname}.png/pdf/svg")

# ============================================================================
# Calibration Curves
# ============================================================================
print(f"\n{'='*70}")
print("GENERATING CALIBRATION CURVES")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Use 4 bins for small sample size (16 samples / 4 = 4 per bin average)
n_cal_bins = 4

# Individual calibration curves for each model
for ax, proba, name, color, marker in [
    (axes[0], proba_4, '4-Feature Model', 'blue', 'o'),
    (axes[1], proba_3, '3-Feature Model', 'green', 's')
]:
    try:
        prob_true, prob_pred = calibration_curve(y_test, proba, n_bins=n_cal_bins, strategy='quantile')
    except:
        prob_true, prob_pred = calibration_curve(y_test, proba, n_bins=n_cal_bins, strategy='uniform')

    brier = brier_score_loss(y_test, proba)
    ece = expected_calibration_error(y_test, proba)

    # Perfect calibration line
    ax.plot([0, 1], [0, 1], 'k--', linewidth=2, label='Perfectly Calibrated', zorder=1)

    # Calibration curve
    ax.plot(prob_pred, prob_true, f'{color[0]}{marker}-', markersize=14, linewidth=2.5,
            label=f'Model (Brier={brier:.3f})', markeredgecolor='black', markeredgewidth=1.5, zorder=3)

    # Add confidence region (simplified)
    ax.fill_between([0, 1], [0, 1], [0.1, 1.1], alpha=0.1, color='gray', label='Ideal zone')

    # Add histogram of predictions at bottom
    ax_hist = ax.inset_axes([0.55, 0.05, 0.4, 0.15])
    ax_hist.hist(proba, bins=10, color=color, alpha=0.7, edgecolor='black')
    ax_hist.set_xlim([0, 1])
    ax_hist.set_xlabel('Pred. Prob.', fontsize=8)
    ax_hist.set_ylabel('Count', fontsize=8)
    ax_hist.tick_params(labelsize=7)

    ax.set_xlabel('Mean Predicted Probability', fontsize=12)
    ax.set_ylabel('Fraction of Positives', fontsize=12)
    ax.set_title(f'{name}\nBrier Score = {brier:.3f}, ECE = {ece:.3f}', fontsize=12, fontweight='bold')
    ax.legend(loc='upper left', fontsize=9)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

plt.suptitle('Calibration Curves (Reliability Diagrams)',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout(rect=[0, 0, 1, 0.96])
for ext in ['png', 'pdf', 'svg']:
    plt.savefig(os.path.join(OUTPUT_DIR, f'Calibration_Curves_Comparison.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: Calibration_Curves_Comparison.png/pdf/svg")

# Combined calibration plot with offset for visibility
fig, ax = plt.subplots(figsize=(10, 8))

# Perfect calibration line
ax.plot([0, 1], [0, 1], 'k--', linewidth=2.5, label='Perfectly Calibrated', zorder=1)

# 4-Feature with slight offset for visibility
prob_true_4, prob_pred_4 = calibration_curve(y_test, proba_4, n_bins=n_cal_bins, strategy='quantile')
ax.scatter(prob_pred_4, prob_true_4, s=200, c='blue', marker='o', edgecolors='black',
           linewidths=2, label=f'4-Feature (Brier={brier_score_loss(y_test, proba_4):.3f})', zorder=3)
ax.plot(prob_pred_4, prob_true_4, 'b-', linewidth=2, alpha=0.7, zorder=2)

# 3-Feature
prob_true_3, prob_pred_3 = calibration_curve(y_test, proba_3, n_bins=n_cal_bins, strategy='quantile')
ax.scatter(prob_pred_3, prob_true_3, s=200, c='green', marker='s', edgecolors='black',
           linewidths=2, label=f'3-Feature (Brier={brier_score_loss(y_test, proba_3):.3f})', zorder=3)
ax.plot(prob_pred_3, prob_true_3, 'g-', linewidth=2, alpha=0.7, zorder=2)

ax.set_xlabel('Mean Predicted Probability', fontsize=14)
ax.set_ylabel('Fraction of Positives (Observed)', fontsize=14)
ax.set_title('Calibration Curves Comparison\n4-Feature vs 3-Feature Model',
             fontsize=14, fontweight='bold')
ax.legend(loc='lower right', fontsize=11)
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
ax.grid(True, alpha=0.3)
ax.set_aspect('equal')

plt.tight_layout()
for ext in ['png', 'pdf', 'svg']:
    plt.savefig(os.path.join(OUTPUT_DIR, f'Calibration_Combined.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: Calibration_Combined.png/pdf/svg")

# ============================================================================
# Comprehensive Summary Figure
# ============================================================================
print(f"\n{'='*70}")
print("GENERATING COMPREHENSIVE SUMMARY FIGURE")
print("=" * 70)

fig, axes = plt.subplots(2, 5, figsize=(24, 10))
axes = axes.flatten()

# Plot all metrics
metric_names = list(results.keys())
for idx, metric_name in enumerate(metric_names):
    ax = axes[idx]
    data = results[metric_name]
    config = data['config']

    val_4 = data['4-Feature']['mean']
    val_3 = data['3-Feature']['mean']
    err_4_low = val_4 - data['4-Feature']['ci_lower']
    err_4_high = data['4-Feature']['ci_upper'] - val_4
    err_3_low = val_3 - data['3-Feature']['ci_lower']
    err_3_high = data['3-Feature']['ci_upper'] - val_3

    x = [0, 1]
    bars = ax.bar(x, [val_4, val_3], color=['#1f77b4', '#2ca02c'],
                  edgecolor='black', linewidth=1, width=0.6,
                  yerr=[[err_4_low, err_3_low], [err_4_high, err_3_high]],
                  capsize=5, error_kw={'capthick': 1.5})

    # Highlight winner
    if config['higher_better']:
        winner_idx = 1 if val_3 > val_4 else 0
    else:
        winner_idx = 1 if val_3 < val_4 else 0
    bars[winner_idx].set_edgecolor('gold')
    bars[winner_idx].set_linewidth(3)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['4-Feat', '3-Feat'], fontsize=9)
    ax.set_title(metric_name, fontsize=11, fontweight='bold')

    # Add values
    for i, val in enumerate([val_4, val_3]):
        ax.text(i, val + 0.02, f'{val:.3f}', ha='center', fontsize=9)

plt.suptitle('Complete Metrics Comparison: 4-Feature vs 3-Feature Model\n'
             '(Gold border = Winner, Error bars = 95% CI)',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout(rect=[0, 0, 1, 0.96])
for ext in ['png', 'pdf', 'svg']:
    plt.savefig(os.path.join(OUTPUT_DIR, f'All_Metrics_Summary.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved: All_Metrics_Summary.png/pdf/svg")

# ============================================================================
# Summary
# ============================================================================
print(f"\n{'='*70}")
print("FIGURE GENERATION COMPLETE")
print("=" * 70)

print(f"\nAll figures saved to: {OUTPUT_DIR}/")
print("\nFiles generated:")
for metric_name in results.keys():
    fname = f'Comparison_{metric_name.replace("-", "_").replace(" ", "_")}'
    print(f"  - {fname}.png/pdf/svg")

print(f"  - ROC_Comparison_with_95CI.png/pdf/svg")
print(f"  - ROC_Individual_with_95CI.png/pdf/svg")
print(f"  - Confusion_Matrix_Comparison.png/pdf/svg")
print(f"  - CM_4Feature.png/pdf/svg")
print(f"  - CM_3Feature.png/pdf/svg")
print(f"  - Calibration_Curves_Comparison.png/pdf/svg")
print(f"  - All_Metrics_Summary.png/pdf/svg")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
