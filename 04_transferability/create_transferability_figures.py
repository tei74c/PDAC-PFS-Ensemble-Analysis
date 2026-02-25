"""
Transferability Analysis Figures
================================
Publication-quality figures for signature transferability across endpoints:
PFS, Treatment Response, and Treatment Assignment.

Reference:
Khoshnevis et al. Molecular Cancer (2026).
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
out_dir = os.path.join(BASE_DIR, "results")
base = os.path.join(BASE_DIR, "data")
os.makedirs(out_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Publication-quality figure settings
# ---------------------------------------------------------------------------
plt.rcParams['pdf.fonttype'] = 42   # TrueType (editable in Illustrator)
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 11
plt.style.use('seaborn-v0_8-whitegrid')

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
df_summary = pd.read_excel(
    os.path.join(base, 'COMPREHENSIVE_RESULTS.xlsx'), sheet_name='Summary')

feat_data = {}
for endpoint in ['PFS', 'Treatment', 'Response']:
    feat_data[endpoint] = pd.read_excel(
        os.path.join(base, 'COMPREHENSIVE_RESULTS.xlsx'),
        sheet_name=f'{endpoint}_Features')

# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------
colors = {'PFS': '#2E86AB', 'Treatment': '#E94F37', 'Response': '#2A9D8F'}
endpoints = ['PFS', 'Treatment', 'Response']


def parse_ci(ci_str):
    """Parse confidence interval string '[lower-upper]' into a tuple."""
    ci_str = ci_str.replace('[', '').replace(']', '')
    parts = ci_str.split('-')
    return float(parts[0]), float(parts[1])


# =============================================================================
# FIGURE 1: Performance Summary Comparison
# =============================================================================
fig = plt.figure(figsize=(16, 14))
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.30)

# Panel A: ROC-AUC Comparison (CV vs Test) --------------------------------
ax = fig.add_subplot(gs[0, 0])
x = np.arange(3)
width = 0.35

cv_auc = [df_summary[df_summary['Dataset'] == e]['CV_ROC_AUC_Mean'].values[0]
          for e in endpoints]
test_auc = [df_summary[df_summary['Dataset'] == e]['Test_ROC_AUC'].values[0]
            for e in endpoints]

cv_ci = [parse_ci(df_summary[df_summary['Dataset'] == e]['CV_ROC_AUC_CI95'].values[0])
         for e in endpoints]
test_ci = [parse_ci(df_summary[df_summary['Dataset'] == e]['Test_ROC_AUC_CI95'].values[0])
           for e in endpoints]

cv_err = [[cv_auc[i] - cv_ci[i][0] for i in range(3)],
          [cv_ci[i][1] - cv_auc[i] for i in range(3)]]
test_err = [[test_auc[i] - test_ci[i][0] for i in range(3)],
            [test_ci[i][1] - test_auc[i] for i in range(3)]]

bars1 = ax.bar(x - width / 2, cv_auc, width, label='Cross-Validation',
               color=[colors[e] for e in endpoints],
               yerr=cv_err, capsize=5, alpha=0.8, edgecolor='black', linewidth=1)
bars2 = ax.bar(x + width / 2, test_auc, width, label='Holdout Test',
               color=[colors[e] for e in endpoints],
               yerr=test_err, capsize=5, alpha=0.5, edgecolor='black',
               linewidth=1, hatch='//')

ax.axhline(y=0.5, color='red', linestyle='--', linewidth=1.5, alpha=0.7,
           label='Chance (0.5)')
ax.set_ylabel('ROC-AUC', fontsize=12, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(endpoints, fontsize=11, fontweight='bold')
ax.set_ylim([0.4, 1.08])
ax.set_title('(A) Discrimination: ROC-AUC', fontsize=14, fontweight='bold',
             pad=15)
ax.legend(loc='lower right', fontsize=9)

# Significance markers
for i, e in enumerate(endpoints):
    if test_ci[i][0] > 0.5:
        ax.annotate('*', xy=(i + width / 2, test_auc[i] + test_err[1][i] + 0.02),
                    ha='center', fontsize=16, fontweight='bold', color='green')
    else:
        ax.annotate('n.s.', xy=(i + width / 2, test_auc[i] + test_err[1][i] + 0.02),
                    ha='center', fontsize=9, color='red')

# Panel B: Feature Significance Heatmap ------------------------------------
ax = fig.add_subplot(gs[0, 1])
genes = ['DYNC2H1', 'ECM2', 'ENPP1', 'PPIB']
f_scores = np.zeros((4, 3))
p_values = np.zeros((4, 3))
significance = np.zeros((4, 3))

for j, endpoint in enumerate(endpoints):
    df_feat = feat_data[endpoint]
    for i, gene in enumerate(genes):
        row = df_feat[df_feat['Feature'] == gene]
        if len(row) > 0:
            f_scores[i, j] = row['F_Score'].values[0]
            p_values[i, j] = row['P_Value_FDR'].values[0]
            significance[i, j] = 1 if row['Significant_FDR'].values[0] else 0

im = ax.imshow(f_scores, cmap='YlOrRd', aspect='auto', vmin=0, vmax=12)

for i in range(4):
    for j in range(3):
        p = p_values[i, j]
        stars = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        text_color = 'white' if f_scores[i, j] > 6 else 'black'
        ax.text(j, i, f'{f_scores[i, j]:.1f}{stars}', ha='center', va='center',
                fontsize=11, fontweight='bold', color=text_color)

ax.set_xticks(range(3))
ax.set_xticklabels(endpoints, fontsize=11, fontweight='bold')
ax.set_yticks(range(4))
ax.set_yticklabels(genes, fontsize=11, fontweight='bold')
ax.set_title('(B) Feature Relevance (ANOVA F-score)', fontsize=14,
             fontweight='bold', pad=15)
cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('F-score', fontsize=10)

# Panel C: Calibration Comparison (Brier & ECE) ----------------------------
ax = fig.add_subplot(gs[1, 0])
metrics = ['Brier Score', 'ECE']
x = np.arange(len(metrics))
width = 0.25

brier_test = [df_summary[df_summary['Dataset'] == e]['Test_Brier'].values[0]
              for e in endpoints]
ece_test = [df_summary[df_summary['Dataset'] == e]['Test_ECE'].values[0]
            for e in endpoints]

for i, endpoint in enumerate(endpoints):
    vals = [brier_test[i], ece_test[i]]
    ax.bar(x + i * width, vals, width, label=endpoint, color=colors[endpoint],
           edgecolor='black', linewidth=1)

ax.set_ylabel('Score (lower is better)', fontsize=12, fontweight='bold')
ax.set_xticks(x + width)
ax.set_xticklabels(metrics, fontsize=11, fontweight='bold')
ax.set_title('(C) Calibration Quality (Holdout Test)', fontsize=14,
             fontweight='bold', pad=15)
ax.legend(loc='upper left', fontsize=9)
ax.set_ylim([0, 0.35])

ax.axhline(y=0.25, color='red', linestyle='--', alpha=0.5, linewidth=1)
ax.text(1.7, 0.26, 'Poor', fontsize=8, color='red')

# Panel D: Summary Assessment -----------------------------------------------
ax = fig.add_subplot(gs[1, 1])
ax.axis('off')

table_data = [
    ['Endpoint', 'Test AUC', 'Features Sig.', 'Transfer?'],
    ['PFS', '0.922', '3/4', 'Reference'],
    ['Treatment', '0.859*', '0/4', 'NO'],
    ['Response', '0.889', '3/4', 'YES']
]

cell_colors = [['lightgray'] * 4,
               ['lightblue'] * 4,
               ['lightcoral'] * 4,
               ['lightgreen'] * 4]

table = ax.table(cellText=table_data, cellColours=cell_colors,
                 loc='center', cellLoc='center',
                 colWidths=[0.25, 0.25, 0.25, 0.25])
table.auto_set_font_size(False)
table.set_fontsize(12)
table.scale(1.3, 2.2)

ax.set_title('(D) Transferability Summary', fontsize=14, fontweight='bold',
             y=0.88, pad=15)
ax.text(0.5, 0.12, '*95% CI includes 0.5 (chance level)', ha='center',
        fontsize=10, style='italic', transform=ax.transAxes)

plt.savefig(os.path.join(out_dir, 'PAPER_Figure1_Transferability_Summary.png'),
            dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(os.path.join(out_dir, 'PAPER_Figure1_Transferability_Summary.pdf'),
            dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()

print('Created Figure 1: Transferability Summary')

# =============================================================================
# FIGURE 2: Detailed Feature Analysis
# =============================================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
plt.subplots_adjust(wspace=0.35, top=0.85, bottom=0.15)

for idx, endpoint in enumerate(endpoints):
    ax = axes[idx]
    df_feat = feat_data[endpoint]

    df_feat = df_feat.sort_values('F_Score', ascending=True)

    colors_bar = ['#2ecc71' if sig else '#e74c3c'
                  for sig in df_feat['Significant_FDR']]

    bars = ax.barh(df_feat['Feature'], df_feat['F_Score'], color=colors_bar,
                   edgecolor='black', linewidth=1, height=0.6)

    for i, (_, row) in enumerate(df_feat.iterrows()):
        p = row['P_Value_FDR']
        p_text = f'p={p:.3f}' if p >= 0.001 else 'p<0.001'
        ax.text(row['F_Score'] + 0.4, i, p_text, va='center', fontsize=9)

    ax.set_xlabel('ANOVA F-score', fontsize=12, fontweight='bold')
    ax.set_title(f'{endpoint}', fontsize=14, fontweight='bold',
                 color=colors[endpoint], pad=12)
    ax.axvline(x=4.0, color='gray', linestyle='--', alpha=0.5,
               label='F=4 threshold')
    ax.set_xlim([0, 15])

legend_elements = [
    Patch(facecolor='#2ecc71', edgecolor='black', label='Significant (FDR<0.05)'),
    Patch(facecolor='#e74c3c', edgecolor='black', label='Not Significant')]
fig.legend(handles=legend_elements, loc='upper center', ncol=2, fontsize=11,
           bbox_to_anchor=(0.5, 0.98))

plt.savefig(os.path.join(out_dir, 'PAPER_Figure2_Feature_Significance.png'),
            dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(os.path.join(out_dir, 'PAPER_Figure2_Feature_Significance.pdf'),
            dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()

print('Created Figure 2: Feature Significance')

# =============================================================================
# FIGURE 3: Treatment Model Comparison (PFS Signature vs Purpose-Built)
# =============================================================================
fig = plt.figure(figsize=(18, 12))
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.40, wspace=0.35,
                       height_ratios=[1, 1])

# Data for comparison
pfs_sig_auc = 0.859
pfs_sig_auc_ci = (0.500, 1.000)
treatment_model_auc = 1.000
treatment_model_auc_ci = (1.000, 1.000)
pfs_for_pfs_auc = 0.922

pfs_sig_f1 = 0.875
treatment_model_f1 = 0.941

pfs_sig_brier = 0.146
treatment_model_brier = 0.044
pfs_for_pfs_brier = 0.122

# SHAP values for Treatment model
shap_features = ['Age', 'MMRN2', 'DYNC2H1']
shap_values = [0.035, 0.195, 0.279]

# Panel A: Discrimination Performance ---------------------------------------
ax = fig.add_subplot(gs[0, 0])
models = ['PFS Signature\n(for PFS)', 'PFS Signature\n(for Treatment)',
          'Treatment-Specific\nModel']
aucs = [pfs_for_pfs_auc, pfs_sig_auc, treatment_model_auc]
model_colors = ['#2E86AB', '#E94F37', '#27ae60']

auc_lower = [0.656, 0.500, 1.000]
auc_upper = [1.000, 1.000, 1.000]
errors = [[aucs[i] - auc_lower[i] for i in range(3)],
          [auc_upper[i] - aucs[i] for i in range(3)]]

bars = ax.bar(range(3), aucs, color=model_colors, edgecolor='black',
              linewidth=1.5, yerr=errors, capsize=8, width=0.6)
ax.axhline(y=0.5, color='red', linestyle='--', linewidth=2, alpha=0.7)
ax.set_ylabel('Test ROC-AUC', fontsize=13, fontweight='bold')
ax.set_xticks(range(3))
ax.set_xticklabels(models, fontsize=10, fontweight='bold')
ax.set_ylim([0.4, 1.12])
ax.set_title('(A) Discrimination Performance', fontsize=14, fontweight='bold',
             pad=15)

for i, (val, bar) in enumerate(zip(aucs, bars)):
    ax.text(bar.get_x() + bar.get_width() / 2, val + errors[1][i] + 0.025,
            f'{val:.3f}', ha='center', va='bottom', fontsize=11,
            fontweight='bold')

ax.annotate('*CI incl. 0.5', xy=(1, 0.53), fontsize=9, color='red',
            ha='center')

# Panel B: Feature Comparison Table -----------------------------------------
ax = fig.add_subplot(gs[0, 1])
ax.axis('off')

table_data = [
    ['', 'PFS Signature', 'Treatment Model'],
    ['Features', 'DYNC2H1, ECM2,\nENPP1, PPIB', 'Age, DYNC2H1,\nMMRN2'],
    ['N Features', '4', '3'],
    ['Shared', '', 'DYNC2H1'],
    ['Unique', 'ECM2, ENPP1, PPIB', 'Age, MMRN2'],
    ['For Treatment', '0/4 significant', '3/3 relevant']
]

cell_colors = [['white', '#fadbd8', '#d5f5e3'],
               ['#f5f5f5', '#fadbd8', '#d5f5e3'],
               ['#f5f5f5', '#fadbd8', '#d5f5e3'],
               ['#f5f5f5', '#fdebd0', '#fdebd0'],
               ['#f5f5f5', '#fadbd8', '#d5f5e3'],
               ['#f5f5f5', '#f1948a', '#82e0aa']]

table = ax.table(cellText=table_data, cellColours=cell_colors,
                 loc='center', cellLoc='center',
                 colWidths=[0.28, 0.36, 0.36])
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1.15, 2.0)

for j in range(3):
    table[(0, j)].set_text_props(fontweight='bold')

ax.set_title('(B) Feature Comparison', fontsize=14, fontweight='bold',
             y=0.92, pad=15)

# Panel C: Treatment Model SHAP --------------------------------------------
ax = fig.add_subplot(gs[0, 2])
bar_colors = ['#9b59b6', '#9b59b6', '#27ae60']
bars = ax.barh(shap_features, shap_values, color=bar_colors,
               edgecolor='black', linewidth=1.5, height=0.5)

for i, (val, bar) in enumerate(zip(shap_values, bars)):
    ax.text(val + 0.008, bar.get_y() + bar.get_height() / 2,
            f'{val:.3f}', va='center', fontsize=11, fontweight='bold')

ax.set_xlabel('Mean |SHAP|', fontsize=13, fontweight='bold')
ax.set_xlim([0, 0.35])
ax.set_title('(C) Treatment Model SHAP', fontsize=14, fontweight='bold',
             pad=15)

legend_elements = [
    Patch(facecolor='#27ae60', edgecolor='black', label='Shared with PFS'),
    Patch(facecolor='#9b59b6', edgecolor='black', label='Unique to Treatment')]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9)

# Panel D: Calibration Quality ----------------------------------------------
ax = fig.add_subplot(gs[1, 0])
models_short = ['PFS\u2192PFS', 'PFS\u2192Treatment', 'Treatment Model']
brier_scores = [pfs_for_pfs_brier, pfs_sig_brier, treatment_model_brier]

bars = ax.bar(range(3), brier_scores, color=model_colors, edgecolor='black',
              linewidth=1.5, width=0.6)
ax.set_ylabel('Brier Score (lower=better)', fontsize=13, fontweight='bold')
ax.set_xticks(range(3))
ax.set_xticklabels(models_short, fontsize=11, fontweight='bold')
ax.set_ylim([0, 0.20])
ax.set_title('(D) Calibration Quality', fontsize=14, fontweight='bold',
             pad=15)

for i, (val, bar) in enumerate(zip(brier_scores, bars)):
    ax.text(bar.get_x() + bar.get_width() / 2, val + 0.005,
            f'{val:.3f}', ha='center', va='bottom', fontsize=11,
            fontweight='bold')

# Panel E: Head-to-Head Comparison Table ------------------------------------
ax = fig.add_subplot(gs[1, 1:])
ax.axis('off')

comparison_data = [
    ['Metric', 'PFS Signature\n(for Treatment)', 'Treatment-Specific\nModel',
     'Improvement'],
    ['Test ROC-AUC', '0.859 [0.500-1.000]*', '1.000 [1.000-1.000]', '+0.141'],
    ['Test F1', '0.875 [0.615-1.000]', '0.941 [0.842-1.000]', '+0.066'],
    ['Brier Score', '0.146', '0.044', '-0.102 (better)'],
    ['Features Sig.', '0/4', '3/3', 'All relevant'],
    ['Conclusion', 'FAILS to transfer', 'PURPOSE-BUILT', '']
]

comp_colors = [['#d4e6f1', '#d4e6f1', '#d4e6f1', '#d4e6f1'],
               ['white', '#fadbd8', '#d5f5e3', '#d5f5e3'],
               ['white', '#fadbd8', '#d5f5e3', '#d5f5e3'],
               ['white', '#fadbd8', '#d5f5e3', '#d5f5e3'],
               ['white', '#f1948a', '#82e0aa', '#d5f5e3'],
               ['white', '#f1948a', '#82e0aa', 'white']]

table = ax.table(cellText=comparison_data, cellColours=comp_colors,
                 loc='center', cellLoc='center',
                 colWidths=[0.22, 0.28, 0.28, 0.22])
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1.2, 2.0)

for j in range(4):
    table[(0, j)].set_text_props(fontweight='bold')

ax.set_title('(E) Head-to-Head Comparison: Transferred vs Purpose-Built Model',
             fontsize=14, fontweight='bold', y=0.92, pad=15)
ax.text(0.5, 0.05, '*95% CI includes 0.5 (chance level)', ha='center',
        fontsize=10, style='italic', transform=ax.transAxes)

fig.suptitle('Treatment Prediction: Transferred PFS Signature vs '
             'Treatment-Specific Model',
             fontsize=16, fontweight='bold', y=0.98)

plt.savefig(os.path.join(out_dir, 'PAPER_Figure3_Treatment_Model_Comparison.png'),
            dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(os.path.join(out_dir, 'PAPER_Figure3_Treatment_Model_Comparison.pdf'),
            dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()

print('Created Figure 3: Treatment Model Comparison')

# =============================================================================
# TABLE 1: Comprehensive Metrics
# =============================================================================
fig, ax = plt.subplots(figsize=(18, 7))
ax.axis('off')

table_data = [
    ['Metric', 'PFS (Reference)', 'Treatment', 'Response'],
    ['', 'CV / Test', 'CV / Test', 'CV / Test'],
    ['ROC-AUC', '0.917 / 0.922', '0.983 / 0.859*', '0.911 / 0.889'],
    ['F1 Score', '0.871 / 0.875', '0.939 / 0.875', '0.826 / 0.800'],
    ['Balanced Acc.', '0.873 / 0.875', '0.942 / 0.875', '0.853 / 0.817'],
    ['Brier Score', '0.103 / 0.122', '0.045 / 0.146', '0.109 / 0.140'],
    ['ECE', '0.072 / 0.185', '0.041 / 0.250', '0.068 / 0.197'],
    ['Permutation p', '0.002', '0.002', '0.005'],
    ['Features Sig.', '3/4', '0/4', '3/4'],
    ['Transfers?', 'Reference', 'NO (CI incl. 0.5)', 'YES']
]

cell_colors = [['#d4e6f1'] * 4,
               ['#eaecee'] * 4,
               ['white', '#d5f5e3', '#fadbd8', '#d5f5e3'],
               ['white', '#d5f5e3', '#fadbd8', '#d5f5e3'],
               ['white', '#d5f5e3', '#fadbd8', '#d5f5e3'],
               ['white', '#d5f5e3', '#fadbd8', '#d5f5e3'],
               ['white', '#d5f5e3', '#fadbd8', '#d5f5e3'],
               ['white', '#d5f5e3', '#d5f5e3', '#d5f5e3'],
               ['white', '#d5f5e3', '#fadbd8', '#d5f5e3'],
               ['white', '#aed6f1', '#f1948a', '#82e0aa']]

table = ax.table(cellText=table_data, cellColours=cell_colors,
                 loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1.4, 2.0)

for i in range(4):
    table[(0, i)].set_text_props(fontweight='bold')
    table[(1, i)].set_text_props(style='italic', fontsize=10)

ax.set_title('Table 1. Transferability of 4-Gene PFS Signature to '
             'Alternative Clinical Endpoints\n',
             fontsize=15, fontweight='bold')
ax.text(0.5, -0.02, '*95% CI includes 0.5 (chance level); '
        'Green=Good, Red=Poor',
        ha='center', fontsize=10, style='italic', transform=ax.transAxes)

plt.savefig(os.path.join(out_dir, 'PAPER_Table1_Comprehensive_Metrics.png'),
            dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(os.path.join(out_dir, 'PAPER_Table1_Comprehensive_Metrics.pdf'),
            dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()

print('Created Table 1: Comprehensive Metrics')
print(f'\nAll figures saved to: {out_dir}')
