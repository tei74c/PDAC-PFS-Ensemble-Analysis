# -*- coding: utf-8 -*-
"""
Ensemble NB3 + Cox Comparison Analysis
=======================================
Systematic evaluation of ensemble strategies combining 3-feature Naive Bayes
and 4-feature Cox models for PFS group prediction.

Ensemble Strategies:
1. Simple Average: (NB + Cox) / 2
2. Weighted 2:1: (2*NB + Cox) / 3  (NB-weighted)
3. Weighted 1:2: (NB + 2*Cox) / 3  (Cox-weighted) -- CHAMPION
4. Stacking (Logistic): Meta-learner on NB/Cox probabilities

ML Features: DYNC2H1, ECM2, PPIB
Cox Features: TFRC, APOF, ANG, FABP4

Methodology:
- 200x Bootstrap resampling with stratified 20x5-fold CV
- Seven ensemble integration strategies evaluated
- Performance heatmap across nine classification metrics

Reference:
Khoshnevis et al. A seven-protein plasma proteomic ensemble predicts
progression-free survival in therapy-naive stage IV PDAC. Molecular Cancer (2026).
"""

# =============================================================================
# CONFIGURATION
# =============================================================================

import os

# -----------------------------------------------------------------------------
# FILE PATHS
# -----------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FULL_DATASET = os.path.join(BASE_DIR, "data", "FinalwinnerML_Cox_dataset.xlsx")
NB_PICKLE = os.path.join(BASE_DIR, "data", "champion_model_NB_3features.pkl")
OUTPUT_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# -----------------------------------------------------------------------------
# COX MODEL PARAMETERS (from Cox model_winner.xlsx)
# -----------------------------------------------------------------------------
COX_FEATURES = ['TFRC', 'APOF', 'ANG', 'FABP4']
COX_BETAS = {
    'TFRC': 0.390057,
    'APOF': 0.391837,
    'ANG': 0.269207,
    'FABP4': 0.280084
}

# -----------------------------------------------------------------------------
# ML MODEL PARAMETERS - 3-Feature Model
# -----------------------------------------------------------------------------
NB_FEATURES = ['DYNC2H1', 'ECM2', 'PPIB']

# -----------------------------------------------------------------------------
# CLASS CONFIGURATION
# -----------------------------------------------------------------------------
LABEL_COLUMN = 'PFS_group'
CLASS_MAP = {'L': 0, 'S': 1}
CLASS_NAMES = ['Long PFS (L)', 'Short PFS (S)']
POSITIVE_CLASS = 'S'
POSITIVE_LABEL = 1

# Survival columns
TIME_COLUMN = 'PFS'
EVENT_COLUMN = 'Event'

# -----------------------------------------------------------------------------
# BOOTSTRAP CV PARAMETERS
# -----------------------------------------------------------------------------
N_BOOTSTRAP = 200
CV_N_REPEATS = 20
CV_N_FOLDS = 5
TOTAL_EVALUATIONS = N_BOOTSTRAP * CV_N_REPEATS * CV_N_FOLDS

# -----------------------------------------------------------------------------
# STATISTICAL PARAMETERS
# -----------------------------------------------------------------------------
CI_LEVEL = 0.95
FDR_ALPHA = 0.05
N_PERMUTATIONS = 500
DELONG_ALPHA = 0.05

# -----------------------------------------------------------------------------
# ENSEMBLE PARAMETERS
# -----------------------------------------------------------------------------
ENSEMBLE_STRATEGIES = [
    'simple_avg',      # (NB + Cox) / 2
    'weighted_2_1',    # (2*NB + Cox) / 3
    'weighted_1_2',    # (NB + 2*Cox) / 3
    'stacking_lr',     # Logistic regression meta-learner
    'stacking_ridge',  # Ridge regression meta-learner
    'max_confidence',  # Use more confident model
    'optimized',       # Grid search optimal weights
]

# Grid for optimized weights
WEIGHT_GRID = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

# -----------------------------------------------------------------------------
# VISUALIZATION PARAMETERS
# -----------------------------------------------------------------------------
FIGURE_DPI = 300
SAVE_FORMATS = ['png', 'pdf', 'svg']
COLORMAP = {
    'NB': '#2E86AB',           # Steel blue (matches COX_vs_ML)
    'Cox': '#A23B72',          # Magenta (matches COX_vs_ML)
    'Ensemble': '#A23B72',     # Magenta (same as Cox for consistency)
    'simple_avg': '#28A745',   # Green
    'weighted_2_1': '#17A2B8', # Cyan
    'weighted_1_2': '#A23B72', # Magenta (best ensemble - same as Cox for comparison)
    'stacking_lr': '#DC3545',  # Red
    'stacking_ridge': '#6F42C1', # Purple
    'max_confidence': '#FD7E14', # Orange
    'optimized': '#20C997',    # Teal
    'ZeroR': '#6C757D',        # Gray
}

# Font settings
FONT_SIZE = 11
TITLE_SIZE = 14
LABEL_SIZE = 12
TICK_SIZE = 10

# TrueType fonts for Adobe Illustrator
TRUETYPE_SETTINGS = {
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
    'axes.unicode_minus': False,
}

# -----------------------------------------------------------------------------
# COMPUTATIONAL PARAMETERS
# -----------------------------------------------------------------------------
RANDOM_STATE = 42
N_JOBS = -1
VERBOSE = True

# -----------------------------------------------------------------------------
# METRICS
# -----------------------------------------------------------------------------
CORE_METRICS = [
    'Accuracy', 'Balanced_Accuracy', 'F1', 'Precision', 'Recall',
    'Specificity', 'NPV', 'ROC_AUC', 'PRC_AUC', 'MCC'
]

PLOT_METRICS = [
    'Accuracy', 'Balanced_Accuracy', 'F1', 'Precision', 'Recall',
    'Specificity', 'NPV', 'ROC_AUC', 'PRC_AUC'
]

# =============================================================================
# IMPORTS
# =============================================================================

import sys
import json
import pickle
import warnings
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any, Union
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
from matplotlib.gridspec import GridSpec
import matplotlib.transforms as transforms
import seaborn as sns
from scipy import stats
from scipy.stats import ttest_rel, wilcoxon
from scipy.special import ndtri, expit
from scipy.optimize import minimize_scalar

from sklearn.model_selection import (
    RepeatedStratifiedKFold, StratifiedKFold,
    cross_val_predict, permutation_test_score
)
from sklearn.preprocessing import StandardScaler
from sklearn.naive_bayes import GaussianNB
from sklearn.dummy import DummyClassifier
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score,
    precision_score, recall_score, roc_auc_score,
    average_precision_score, confusion_matrix,
    roc_curve, precision_recall_curve, matthews_corrcoef,
    cohen_kappa_score, brier_score_loss, log_loss
)
from sklearn.base import clone, BaseEstimator, ClassifierMixin
from sklearn.utils import resample
from sklearn.isotonic import IsotonicRegression

warnings.filterwarnings('ignore')

# Optional imports
try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    HAS_LIFELINES = True
except ImportError:
    HAS_LIFELINES = False
    print("[WARNING] lifelines not installed.")

try:
    from tqdm.auto import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    def tqdm(x, **kwargs): return x

# Set plot defaults with TrueType fonts
plt.rcParams.update({
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
    'axes.unicode_minus': False,
    'font.size': FONT_SIZE,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'axes.labelsize': LABEL_SIZE,
    'axes.titlesize': TITLE_SIZE,
    'axes.linewidth': 1.0,
    'xtick.labelsize': TICK_SIZE,
    'ytick.labelsize': TICK_SIZE,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'legend.fontsize': TICK_SIZE,
    'figure.titlesize': TITLE_SIZE,
    'figure.dpi': 100,
    'savefig.dpi': FIGURE_DPI,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white'
})

# Define numpy-dependent variables
ROC_MEAN_FPR = np.linspace(0, 1, 100)
PRC_MEAN_RECALL = np.linspace(0, 1, 100)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

LOG_PATH = os.path.join(OUTPUT_DIR, "analysis_log.txt")

# Global storage
GRAPH_DATA_STORAGE = {}


def log(msg: str, level: str = "INFO"):
    """Log message to console and file."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] [{level}] {msg}"
    if VERBOSE:
        try:
            print(line)
        except UnicodeEncodeError:
            print(line.encode('ascii', 'replace').decode('ascii'))
    try:
        with open(LOG_PATH, "a", encoding='utf-8') as f:
            f.write(line + "\n")
    except:
        pass


def savefig(name: str, fig=None, subdir: str = None):
    """Save figure in multiple formats."""
    if fig is None:
        fig = plt.gcf()
    save_dir = OUTPUT_DIR if subdir is None else os.path.join(OUTPUT_DIR, subdir)
    os.makedirs(save_dir, exist_ok=True)
    for fmt in SAVE_FORMATS:
        path = os.path.join(save_dir, f"{name}.{fmt}")
        fig.savefig(path, bbox_inches='tight', dpi=FIGURE_DPI, facecolor='white')
    plt.close(fig)
    log(f"Saved figure: {name}")


def save_excel(path: str, sheets: Dict[str, pd.DataFrame]):
    """Save multiple DataFrames to Excel."""
    try:
        with pd.ExcelWriter(path, engine='openpyxl') as writer:
            for name, df in sheets.items():
                if isinstance(df, dict):
                    df = pd.DataFrame(df)
                df.to_excel(writer, sheet_name=str(name)[:31], index=False)
        log(f"Saved Excel: {os.path.basename(path)}")
        return path
    except Exception as e:
        log(f"Excel save failed: {e}", "ERROR")
        return None


def store_graph_data(graph_name: str, data: Dict):
    """Store underlying data for a graph."""
    GRAPH_DATA_STORAGE[graph_name] = data


# =============================================================================
# MODEL CLASSES
# =============================================================================

class NBChampionClassifier(BaseEstimator, ClassifierMixin):
    """Wrapper for the 3-feature Naive Bayes champion model.

    Uses StandardScaler for feature scaling and sigmoid calibration
    with adaptive CV folds based on minimum class count.
    """

    def __init__(self, features: List[str] = None, calibrate: bool = True,
                 calibration_cv: int = 3):
        self.features = features or NB_FEATURES
        self.calibrate = calibrate
        self.calibration_cv = calibration_cv
        self.model_ = None
        self.scaler_ = None
        self.classes_ = np.array([0, 1])

    def fit(self, X, y):
        # Scale features (matching COX_vs_ML)
        self.scaler_ = StandardScaler()
        X_scaled = self.scaler_.fit_transform(X)

        base_clf = GaussianNB()

        if self.calibrate and len(np.unique(y)) > 1:
            min_class_count = min(np.sum(y == 0), np.sum(y == 1))
            cv_folds = min(self.calibration_cv, min_class_count)

            if cv_folds >= 2:
                self.model_ = CalibratedClassifierCV(
                    base_clf, method='sigmoid', cv=cv_folds
                )
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


class CoxRiskClassifier(BaseEstimator, ClassifierMixin):
    """Cox proportional hazards risk classifier with fixed coefficients."""

    def __init__(self, betas: Dict[str, float], features: List[str],
                 threshold_method: str = 'median'):
        self.betas = betas
        self.features = features
        self.threshold_method = threshold_method
        self.threshold_ = None
        self.classes_ = np.array([0, 1])
        self.scaler_ = None
        self.mean_risk_ = None
        self.std_risk_ = None

    def _compute_risk_scores(self, X):
        """Compute linear risk scores."""
        risk = np.zeros(len(X))
        for i, feat in enumerate(self.features):
            if feat in self.betas:
                risk += self.betas[feat] * X[:, i]
        return risk

    def fit(self, X, y):
        risk_scores = self._compute_risk_scores(X)
        self.mean_risk_ = np.mean(risk_scores)
        self.std_risk_ = np.std(risk_scores) + 1e-8

        if self.threshold_method == 'median':
            self.threshold_ = np.median(risk_scores)
        elif self.threshold_method == 'youden':
            proba = expit(risk_scores - np.median(risk_scores))
            fpr, tpr, thresholds = roc_curve(y, proba)
            j_scores = tpr - fpr
            optimal_idx = np.argmax(j_scores)
            self.threshold_ = thresholds[optimal_idx]
        else:
            self.threshold_ = 0.5
        return self

    def predict(self, X):
        proba = self.predict_proba(X)[:, 1]
        return (proba >= 0.5).astype(int)

    def predict_proba(self, X):
        risk_scores = self._compute_risk_scores(X)
        proba_pos = expit(risk_scores - self.threshold_)
        proba_neg = 1 - proba_pos
        return np.column_stack([proba_neg, proba_pos])


class EnsembleClassifier(BaseEstimator, ClassifierMixin):
    """Ensemble classifier combining NB and Cox models."""

    def __init__(self, strategy: str = 'simple_avg',
                 nb_weight: float = 0.5,
                 meta_learner: str = 'logistic'):
        self.strategy = strategy
        self.nb_weight = nb_weight  # Weight for NB (Cox weight = 1 - nb_weight)
        self.meta_learner = meta_learner
        self.nb_model_ = None
        self.cox_model_ = None
        self.meta_model_ = None
        self.optimal_weight_ = None
        self.classes_ = np.array([0, 1])

    def fit(self, X_nb, X_cox, y):
        """Fit ensemble with both feature sets."""
        # Fit base models
        self.nb_model_ = NBChampionClassifier(features=NB_FEATURES, calibrate=True)
        self.nb_model_.fit(X_nb, y)

        self.cox_model_ = CoxRiskClassifier(betas=COX_BETAS, features=COX_FEATURES)
        self.cox_model_.fit(X_cox, y)

        # Strategy-specific fitting
        if self.strategy in ['stacking_lr', 'stacking_ridge']:
            # Get base model probabilities
            proba_nb = self.nb_model_.predict_proba(X_nb)[:, 1]
            proba_cox = self.cox_model_.predict_proba(X_cox)[:, 1]
            X_meta = np.column_stack([proba_nb, proba_cox])

            if self.strategy == 'stacking_lr':
                self.meta_model_ = LogisticRegression(random_state=RANDOM_STATE)
            else:
                self.meta_model_ = LogisticRegression(
                    penalty='l2', C=1.0, random_state=RANDOM_STATE
                )
            self.meta_model_.fit(X_meta, y)

        elif self.strategy == 'optimized':
            # Grid search for optimal weight
            best_auc = 0
            best_weight = 0.5

            proba_nb = self.nb_model_.predict_proba(X_nb)[:, 1]
            proba_cox = self.cox_model_.predict_proba(X_cox)[:, 1]

            for w in WEIGHT_GRID:
                proba_ens = w * proba_nb + (1 - w) * proba_cox
                try:
                    auc = roc_auc_score(y, proba_ens)
                    if auc > best_auc:
                        best_auc = auc
                        best_weight = w
                except:
                    pass

            self.optimal_weight_ = best_weight
            self.nb_weight = best_weight

        return self

    def predict(self, X_nb, X_cox):
        proba = self.predict_proba(X_nb, X_cox)[:, 1]
        return (proba >= 0.5).astype(int)

    def predict_proba(self, X_nb, X_cox):
        proba_nb = self.nb_model_.predict_proba(X_nb)[:, 1]
        proba_cox = self.cox_model_.predict_proba(X_cox)[:, 1]

        if self.strategy == 'simple_avg':
            proba_pos = 0.5 * proba_nb + 0.5 * proba_cox

        elif self.strategy == 'weighted_2_1':
            proba_pos = (2 * proba_nb + proba_cox) / 3

        elif self.strategy == 'weighted_1_2':
            proba_pos = (proba_nb + 2 * proba_cox) / 3

        elif self.strategy in ['stacking_lr', 'stacking_ridge']:
            X_meta = np.column_stack([proba_nb, proba_cox])
            proba_pos = self.meta_model_.predict_proba(X_meta)[:, 1]

        elif self.strategy == 'max_confidence':
            # Use prediction from more confident model
            conf_nb = np.abs(proba_nb - 0.5)
            conf_cox = np.abs(proba_cox - 0.5)
            proba_pos = np.where(conf_nb >= conf_cox, proba_nb, proba_cox)

        elif self.strategy == 'optimized':
            w = self.optimal_weight_ if self.optimal_weight_ is not None else self.nb_weight
            proba_pos = w * proba_nb + (1 - w) * proba_cox

        else:
            # Default to simple average
            proba_pos = 0.5 * proba_nb + 0.5 * proba_cox

        proba_neg = 1 - proba_pos
        return np.column_stack([proba_neg, proba_pos])


# =============================================================================
# METRICS COMPUTATION
# =============================================================================

def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray,
                   y_proba: np.ndarray) -> Dict[str, float]:
    """Compute all classification metrics."""
    metrics = {}

    try:
        metrics['Accuracy'] = accuracy_score(y_true, y_pred)
        metrics['Balanced_Accuracy'] = balanced_accuracy_score(y_true, y_pred)
        metrics['F1'] = f1_score(y_true, y_pred, zero_division=0)
        metrics['Precision'] = precision_score(y_true, y_pred, zero_division=0)
        metrics['Recall'] = recall_score(y_true, y_pred, zero_division=0)

        tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
        metrics['Specificity'] = tn / (tn + fp) if (tn + fp) > 0 else 0
        metrics['NPV'] = tn / (tn + fn) if (tn + fn) > 0 else 0

        if len(np.unique(y_true)) > 1:
            metrics['ROC_AUC'] = roc_auc_score(y_true, y_proba)
            metrics['PRC_AUC'] = average_precision_score(y_true, y_proba)
        else:
            metrics['ROC_AUC'] = 0.5
            metrics['PRC_AUC'] = 0.5

        metrics['MCC'] = matthews_corrcoef(y_true, y_pred)
        metrics['Brier'] = brier_score_loss(y_true, y_proba)

    except Exception as e:
        log(f"Metrics error: {e}", "WARNING")
        for m in CORE_METRICS:
            if m not in metrics:
                metrics[m] = np.nan

    return metrics


def compute_bca_ci(data: np.ndarray, alpha: float = 0.05) -> Tuple[float, float]:
    """Compute BCa confidence interval."""
    n = len(data)
    if n == 0:
        return np.nan, np.nan

    theta_hat = np.mean(data)

    # Bias correction
    z0 = ndtri(np.mean(data < theta_hat))
    if np.isinf(z0):
        z0 = 0

    # Acceleration (jackknife)
    jackknife_means = np.array([
        np.mean(np.delete(data, i)) for i in range(n)
    ])
    jack_mean = np.mean(jackknife_means)
    num = np.sum((jack_mean - jackknife_means) ** 3)
    den = 6 * (np.sum((jack_mean - jackknife_means) ** 2) ** 1.5)
    a = num / den if den != 0 else 0

    # Adjusted percentiles
    z_alpha = ndtri(alpha / 2)
    z_1_alpha = ndtri(1 - alpha / 2)

    alpha1 = stats.norm.cdf(z0 + (z0 + z_alpha) / (1 - a * (z0 + z_alpha)))
    alpha2 = stats.norm.cdf(z0 + (z0 + z_1_alpha) / (1 - a * (z0 + z_1_alpha)))

    alpha1 = np.clip(alpha1, 0.001, 0.999)
    alpha2 = np.clip(alpha2, 0.001, 0.999)

    ci_lower = np.percentile(data, alpha1 * 100)
    ci_upper = np.percentile(data, alpha2 * 100)

    return ci_lower, ci_upper


def delong_test(y_true: np.ndarray, y_proba1: np.ndarray,
                y_proba2: np.ndarray) -> Dict:
    """DeLong test for comparing two AUCs."""
    n = len(y_true)
    pos_idx = y_true == 1
    neg_idx = y_true == 0
    n_pos = np.sum(pos_idx)
    n_neg = np.sum(neg_idx)

    if n_pos == 0 or n_neg == 0:
        return {'auc1': np.nan, 'auc2': np.nan, 'difference': np.nan,
                'z_score': np.nan, 'p_value': 1.0, 'significant': False}

    auc1 = roc_auc_score(y_true, y_proba1)
    auc2 = roc_auc_score(y_true, y_proba2)

    # Placement values
    V10_1 = np.array([np.mean(y_proba1[pos_idx] > p) + 0.5 * np.mean(y_proba1[pos_idx] == p)
                      for p in y_proba1[neg_idx]])
    V01_1 = np.array([np.mean(y_proba1[neg_idx] < p) + 0.5 * np.mean(y_proba1[neg_idx] == p)
                      for p in y_proba1[pos_idx]])

    V10_2 = np.array([np.mean(y_proba2[pos_idx] > p) + 0.5 * np.mean(y_proba2[pos_idx] == p)
                      for p in y_proba2[neg_idx]])
    V01_2 = np.array([np.mean(y_proba2[neg_idx] < p) + 0.5 * np.mean(y_proba2[neg_idx] == p)
                      for p in y_proba2[pos_idx]])

    # Variances
    S10_1 = np.var(V10_1, ddof=1) if len(V10_1) > 1 else 0
    S01_1 = np.var(V01_1, ddof=1) if len(V01_1) > 1 else 0
    S10_2 = np.var(V10_2, ddof=1) if len(V10_2) > 1 else 0
    S01_2 = np.var(V01_2, ddof=1) if len(V01_2) > 1 else 0

    # Covariances
    S10_12 = np.cov(V10_1, V10_2)[0, 1] if len(V10_1) > 1 else 0
    S01_12 = np.cov(V01_1, V01_2)[0, 1] if len(V01_1) > 1 else 0

    var1 = S10_1 / n_neg + S01_1 / n_pos
    var2 = S10_2 / n_neg + S01_2 / n_pos
    cov12 = S10_12 / n_neg + S01_12 / n_pos

    var_diff = var1 + var2 - 2 * cov12
    se_diff = np.sqrt(max(var_diff, 1e-10))

    z_score = (auc1 - auc2) / se_diff
    p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))

    return {
        'auc1': auc1,
        'auc2': auc2,
        'difference': auc1 - auc2,
        'z_score': z_score,
        'p_value': p_value,
        'significant': p_value < DELONG_ALPHA
    }


def compute_nri_idi(y_true: np.ndarray, y_proba_old: np.ndarray,
                    y_proba_new: np.ndarray, threshold: float = 0.5) -> Dict:
    """Compute Net Reclassification Improvement and Integrated Discrimination Improvement."""
    n = len(y_true)
    events = y_true == 1
    non_events = y_true == 0

    # Categorical NRI
    pred_old = (y_proba_old >= threshold).astype(int)
    pred_new = (y_proba_new >= threshold).astype(int)

    # Events
    events_up = np.sum((pred_new > pred_old) & events)
    events_down = np.sum((pred_new < pred_old) & events)
    n_events = np.sum(events)

    # Non-events
    non_events_up = np.sum((pred_new > pred_old) & non_events)
    non_events_down = np.sum((pred_new < pred_old) & non_events)
    n_non_events = np.sum(non_events)

    nri_events = (events_up - events_down) / n_events if n_events > 0 else 0
    nri_non_events = (non_events_down - non_events_up) / n_non_events if n_non_events > 0 else 0
    nri = nri_events + nri_non_events

    # NRI SE and p-value
    nri_se = np.sqrt(
        (events_up + events_down) / (n_events ** 2) +
        (non_events_up + non_events_down) / (n_non_events ** 2)
    ) if n_events > 0 and n_non_events > 0 else 0
    nri_z = nri / nri_se if nri_se > 0 else 0
    nri_p = 2 * (1 - stats.norm.cdf(abs(nri_z)))

    # IDI
    p_new_events = np.mean(y_proba_new[events]) if np.sum(events) > 0 else 0.5
    p_new_non_events = np.mean(y_proba_new[non_events]) if np.sum(non_events) > 0 else 0.5
    p_old_events = np.mean(y_proba_old[events]) if np.sum(events) > 0 else 0.5
    p_old_non_events = np.mean(y_proba_old[non_events]) if np.sum(non_events) > 0 else 0.5

    idi = (p_new_events - p_new_non_events) - (p_old_events - p_old_non_events)

    # IDI SE
    var_new = (np.var(y_proba_new[events]) / n_events if n_events > 0 else 0) + \
              (np.var(y_proba_new[non_events]) / n_non_events if n_non_events > 0 else 0)
    var_old = (np.var(y_proba_old[events]) / n_events if n_events > 0 else 0) + \
              (np.var(y_proba_old[non_events]) / n_non_events if n_non_events > 0 else 0)
    idi_se = np.sqrt(var_new + var_old)
    idi_z = idi / idi_se if idi_se > 0 else 0
    idi_p = 2 * (1 - stats.norm.cdf(abs(idi_z)))

    return {
        'NRI': nri,
        'NRI_events': nri_events,
        'NRI_non_events': nri_non_events,
        'NRI_SE': nri_se,
        'NRI_z': nri_z,
        'NRI_p': nri_p,
        'IDI': idi,
        'IDI_SE': idi_se,
        'IDI_z': idi_z,
        'IDI_p': idi_p,
        'events_up': events_up,
        'events_down': events_down,
        'non_events_up': non_events_up,
        'non_events_down': non_events_down
    }


# =============================================================================
# COMPREHENSIVE STATISTICAL COMPARISON
# =============================================================================

def compute_cohens_d(group1: np.ndarray, group2: np.ndarray) -> float:
    """Compute Cohen's d effect size for paired samples."""
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    if pooled_std == 0:
        return 0.0
    return (np.mean(group1) - np.mean(group2)) / pooled_std


def compute_pairwise_statistics(dfs: Dict, model1: str, model2: str,
                                 metrics: List[str] = None) -> Dict:
    """
    Compute comprehensive pairwise statistics between two models.

    Returns paired t-test p-values with Bonferroni correction and Cohen's d
    for all specified metrics.
    """
    if metrics is None:
        metrics = PLOT_METRICS

    results = {}
    n_comparisons = len(metrics)
    bonferroni_alpha = FDR_ALPHA / n_comparisons

    for metric in metrics:
        scores1 = dfs[model1][metric].values
        scores2 = dfs[model2][metric].values

        # Paired t-test
        try:
            stat, p_val = ttest_rel(scores1, scores2)
        except:
            stat, p_val = np.nan, 1.0

        # Cohen's d
        cohens_d = compute_cohens_d(scores1, scores2)

        # Determine significance with Bonferroni correction
        if p_val < 0.001 / n_comparisons:
            sig_level = '***'
        elif p_val < 0.01 / n_comparisons:
            sig_level = '**'
        elif p_val < 0.05 / n_comparisons:
            sig_level = '*'
        else:
            sig_level = 'ns'

        results[metric] = {
            'mean_1': np.mean(scores1),
            'mean_2': np.mean(scores2),
            'difference': np.mean(scores1) - np.mean(scores2),
            't_stat': stat,
            'p_value': p_val,
            'p_bonferroni': min(p_val * n_comparisons, 1.0),
            'significant': p_val < bonferroni_alpha,
            'sig_level': sig_level,
            'cohens_d': cohens_d,
            'effect_interpretation': (
                'large' if abs(cohens_d) >= 0.8 else
                'medium' if abs(cohens_d) >= 0.5 else
                'small' if abs(cohens_d) >= 0.2 else
                'negligible'
            )
        }

    return results


def plot_head_to_head_comparison(dfs: Dict, summaries: Dict,
                                  model1: str, model2: str,
                                  pairwise_stats: Dict):
    """
    Head-to-head comparison plot with significance brackets and Cohen's d.
    Matches the style of COX_vs_ML_Comparison script.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    metrics = list(pairwise_stats.keys())
    n_metrics = len(metrics)

    # Panel A: Bar chart with significance brackets
    ax = axes[0]
    x = np.arange(n_metrics)
    width = 0.35

    means1 = [summaries[model1][m]['mean'] for m in metrics]
    means2 = [summaries[model2][m]['mean'] for m in metrics]
    ci1 = [(summaries[model1][m]['ci_upper'] - summaries[model1][m]['ci_lower'])/2 for m in metrics]
    ci2 = [(summaries[model2][m]['ci_upper'] - summaries[model2][m]['ci_lower'])/2 for m in metrics]

    bars1 = ax.bar(x - width/2, means1, width, label=model1,
                   color=COLORMAP.get(model1, '#2E86AB'), edgecolor='white', linewidth=0.5)
    bars2 = ax.bar(x + width/2, means2, width, label=model2,
                   color=COLORMAP.get(model2, '#A23B72'), edgecolor='white', linewidth=0.5)

    ax.errorbar(x - width/2, means1, yerr=ci1, fmt='none', color='black', capsize=3, linewidth=1)
    ax.errorbar(x + width/2, means2, yerr=ci2, fmt='none', color='black', capsize=3, linewidth=1)

    # Add significance brackets
    for i, m in enumerate(metrics):
        sig = pairwise_stats[m]['sig_level']
        if sig != 'ns':
            max_y = max(means1[i] + ci1[i], means2[i] + ci2[i])
            bracket_y = max_y + 0.03
            ax.plot([i - width/2, i - width/2, i + width/2, i + width/2],
                   [max_y + 0.01, bracket_y, bracket_y, max_y + 0.01],
                   color='black', linewidth=1)
            ax.text(i, bracket_y + 0.01, sig, ha='center', fontsize=10, fontweight='bold')

    ax.set_ylabel('Score')
    ax.set_title(f'Performance Comparison: {model1} vs {model2}\n'
                 f'Bonferroni-corrected: *** p<0.001/{n_metrics}, ** p<0.01/{n_metrics}, * p<0.05/{n_metrics}')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, rotation=45, ha='right')
    ax.legend(loc='lower right')
    ax.set_ylim(0, 1.15)
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)

    # Panel B: Effect size (Cohen's d) horizontal bar chart
    ax = axes[1]
    cohens_d_vals = [pairwise_stats[m]['cohens_d'] for m in metrics]
    colors = ['#2E86AB' if d > 0 else '#A23B72' for d in cohens_d_vals]

    y_pos = np.arange(n_metrics)
    ax.barh(y_pos, cohens_d_vals, color=colors, edgecolor='white', linewidth=0.5, height=0.6)

    # Add value labels
    for i, (d, m) in enumerate(zip(cohens_d_vals, metrics)):
        label_x = d + 0.02 if d >= 0 else d - 0.02
        ha = 'left' if d >= 0 else 'right'
        ax.text(label_x, i, f'd={d:.2f}', va='center', ha=ha, fontsize=9)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(metrics)
    ax.set_xlabel(f'Difference ({model1} - {model2})\n(Positive = {model1} better)')
    ax.set_title(f'Performance Difference\n(Cohen\'s d effect size)')
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)

    # Add effect size reference lines
    for threshold, label in [(0.2, 'small'), (0.5, 'medium'), (0.8, 'large')]:
        ax.axvline(x=threshold, color='gray', linestyle=':', alpha=0.5)
        ax.axvline(x=-threshold, color='gray', linestyle=':', alpha=0.5)

    plt.suptitle(f'Head-to-Head: {model1} vs {model2}', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    savefig('head_to_head_comparison', fig)

    # Store underlying data
    store_graph_data('Head_to_Head_Comparison', {
        'Metric': metrics,
        f'{model1}_Mean': means1,
        f'{model2}_Mean': means2,
        'Difference': [pairwise_stats[m]['difference'] for m in metrics],
        'Cohen_d': cohens_d_vals,
        'P_Value': [pairwise_stats[m]['p_value'] for m in metrics],
        'P_Bonferroni': [pairwise_stats[m]['p_bonferroni'] for m in metrics],
        'Significant': [pairwise_stats[m]['sig_level'] for m in metrics]
    })

    return fig


# =============================================================================
# BOOTSTRAP CV EVALUATION
# =============================================================================

def run_bootstrap_cv_all_ensembles(X_nb: np.ndarray, X_cox: np.ndarray, y: np.ndarray,
                                   patient_ids: np.ndarray = None,
                                   n_bootstrap: int = 200,
                                   n_repeats: int = 20,
                                   n_folds: int = 5,
                                   seed: int = 42) -> Dict:
    """Run bootstrap CV for NB, all ensemble strategies, and ZeroR."""

    log(f"\nRunning Bootstrap CV: {n_bootstrap} bootstrap x {n_repeats}x{n_folds} CV")
    log(f"Total evaluations: {n_bootstrap * n_repeats * n_folds}")
    log(f"Testing {len(ENSEMBLE_STRATEGIES)} ensemble strategies")

    rng = np.random.RandomState(seed)
    n_samples = len(y)

    if patient_ids is None:
        patient_ids = np.arange(n_samples)

    # Storage for all models
    all_metrics = {
        'NB': defaultdict(list),
        'ZeroR': defaultdict(list)
    }
    for strategy in ENSEMBLE_STRATEGIES:
        all_metrics[strategy] = defaultdict(list)

    # Prediction storage
    all_y_true = []
    all_y_proba = {'NB': []}
    for strategy in ENSEMBLE_STRATEGIES:
        all_y_proba[strategy] = []

    # Optimal weights tracker for 'optimized' strategy
    optimal_weights = []

    # ROC/PRC curve storage for CI visualization
    bootstrap_roc = {'NB': []}  # List of (fpr, tpr) tuples per bootstrap
    bootstrap_prc = {'NB': []}  # List of (recall, precision) tuples per bootstrap
    for strategy in ENSEMBLE_STRATEGIES:
        bootstrap_roc[strategy] = []
        bootstrap_prc[strategy] = []

    for boot_idx in tqdm(range(n_bootstrap), desc="Bootstrap iterations"):
        boot_indices = rng.choice(n_samples, size=n_samples, replace=True)

        X_nb_boot = X_nb[boot_indices]
        X_cox_boot = X_cox[boot_indices]
        y_boot = y[boot_indices]

        cv = RepeatedStratifiedKFold(n_splits=n_folds, n_repeats=n_repeats,
                                     random_state=seed + boot_idx)

        fold_metrics = {name: defaultdict(list) for name in all_metrics.keys()}

        boot_y_true = []
        boot_y_proba = {name: [] for name in all_y_proba.keys()}

        for train_idx, val_idx in cv.split(X_nb_boot, y_boot):
            X_nb_tr, X_nb_val = X_nb_boot[train_idx], X_nb_boot[val_idx]
            X_cox_tr, X_cox_val = X_cox_boot[train_idx], X_cox_boot[val_idx]
            y_tr, y_val = y_boot[train_idx], y_boot[val_idx]

            if len(np.unique(y_val)) < 2:
                continue

            # NB Model
            try:
                nb_clf = NBChampionClassifier(features=NB_FEATURES, calibrate=True)
                nb_clf.fit(X_nb_tr, y_tr)
                y_pred_nb = nb_clf.predict(X_nb_val)
                y_proba_nb = nb_clf.predict_proba(X_nb_val)[:, 1]
            except:
                y_pred_nb = np.zeros(len(y_val))
                y_proba_nb = np.full(len(y_val), 0.5)

            metrics_nb = compute_metrics(y_val, y_pred_nb, y_proba_nb)
            for k, v in metrics_nb.items():
                fold_metrics['NB'][k].append(v)

            boot_y_true.extend(y_val)
            boot_y_proba['NB'].extend(y_proba_nb)

            # All ensemble strategies
            for strategy in ENSEMBLE_STRATEGIES:
                try:
                    ens_clf = EnsembleClassifier(strategy=strategy)
                    ens_clf.fit(X_nb_tr, X_cox_tr, y_tr)

                    # Get uncalibrated predictions on training set for calibrator
                    y_proba_train_uncal = ens_clf.predict_proba(X_nb_tr, X_cox_tr)[:, 1]

                    # Fit isotonic calibrator on training predictions
                    calibrator = IsotonicRegression(out_of_bounds='clip')
                    calibrator.fit(y_proba_train_uncal, y_tr)

                    # Get uncalibrated predictions on validation set
                    y_proba_val_uncal = ens_clf.predict_proba(X_nb_val, X_cox_val)[:, 1]

                    # Apply isotonic calibration to validation predictions
                    y_proba_ens = calibrator.predict(y_proba_val_uncal)

                    # Clip to valid probability range
                    y_proba_ens = np.clip(y_proba_ens, 0.001, 0.999)

                    y_pred_ens = (y_proba_ens >= 0.5).astype(int)

                    if strategy == 'optimized' and ens_clf.optimal_weight_ is not None:
                        optimal_weights.append(ens_clf.optimal_weight_)
                except Exception as e:
                    y_pred_ens = np.zeros(len(y_val))
                    y_proba_ens = np.full(len(y_val), 0.5)

                metrics_ens = compute_metrics(y_val, y_pred_ens, y_proba_ens)
                for k, v in metrics_ens.items():
                    fold_metrics[strategy][k].append(v)

                boot_y_proba[strategy].extend(y_proba_ens)

            # ZeroR baseline
            zeror = DummyClassifier(strategy='most_frequent')
            zeror.fit(X_nb_tr, y_tr)
            y_pred_zeror = zeror.predict(X_nb_val)
            y_proba_zeror = np.full(len(y_val), np.mean(y_tr))

            metrics_zeror = compute_metrics(y_val, y_pred_zeror, y_proba_zeror)
            for k, v in metrics_zeror.items():
                fold_metrics['ZeroR'][k].append(v)

        # Aggregate fold metrics for this bootstrap
        for model_name in all_metrics.keys():
            for metric_name, values in fold_metrics[model_name].items():
                if values:
                    all_metrics[model_name][metric_name].append(np.mean(values))

        # Store predictions
        all_y_true.extend(boot_y_true)
        for model_name in all_y_proba.keys():
            all_y_proba[model_name].extend(boot_y_proba[model_name])

        # Compute and store ROC/PRC curves for this bootstrap (for CI visualization)
        boot_y_true_arr = np.array(boot_y_true)
        if len(np.unique(boot_y_true_arr)) >= 2:
            for model_name in bootstrap_roc.keys():
                boot_proba = np.array(boot_y_proba[model_name])
                try:
                    fpr, tpr, _ = roc_curve(boot_y_true_arr, boot_proba)
                    bootstrap_roc[model_name].append((fpr, tpr))

                    prec, rec, _ = precision_recall_curve(boot_y_true_arr, boot_proba)
                    bootstrap_prc[model_name].append((rec, prec))
                except:
                    pass

    # Convert to arrays
    all_y_true = np.array(all_y_true)
    for model_name in all_y_proba.keys():
        all_y_proba[model_name] = np.array(all_y_proba[model_name])

    # Compute summaries with BCa CI
    summaries = {}
    for model_name, metrics_dict in all_metrics.items():
        summary = {}
        for metric_name, values in metrics_dict.items():
            values = np.array(values)
            ci_low, ci_high = compute_bca_ci(values)
            summary[metric_name] = {
                'mean': np.mean(values),
                'std': np.std(values),
                'ci_lower': ci_low,
                'ci_upper': ci_high
            }
        summaries[model_name] = summary

    # Convert metrics to DataFrames
    dfs = {}
    for model_name, metrics_dict in all_metrics.items():
        dfs[model_name] = pd.DataFrame(dict(metrics_dict))

    return {
        'summaries': summaries,
        'dataframes': dfs,
        'all_y_true': all_y_true,
        'all_y_proba': all_y_proba,
        'optimal_weights': np.array(optimal_weights) if optimal_weights else None,
        'bootstrap_roc': bootstrap_roc,
        'bootstrap_prc': bootstrap_prc
    }


# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def plot_ensemble_comparison_bars(summaries: Dict, best_ensemble: str):
    """Bar chart comparing NB vs best ensemble across metrics."""
    fig, ax = plt.subplots(figsize=(14, 8))

    metrics = PLOT_METRICS
    x = np.arange(len(metrics))
    width = 0.35

    nb_means = [summaries['NB'][m]['mean'] for m in metrics]
    nb_cis = [(summaries['NB'][m]['mean'] - summaries['NB'][m]['ci_lower'],
               summaries['NB'][m]['ci_upper'] - summaries['NB'][m]['mean']) for m in metrics]

    ens_means = [summaries[best_ensemble][m]['mean'] for m in metrics]
    ens_cis = [(summaries[best_ensemble][m]['mean'] - summaries[best_ensemble][m]['ci_lower'],
                summaries[best_ensemble][m]['ci_upper'] - summaries[best_ensemble][m]['mean']) for m in metrics]

    bars1 = ax.bar(x - width/2, nb_means, width, label='NB (3-feature)',
                   color=COLORMAP['NB'], edgecolor='white', linewidth=0.5)
    ax.errorbar(x - width/2, nb_means, yerr=np.array(nb_cis).T, fmt='none',
                color='black', capsize=3, linewidth=1)

    bars2 = ax.bar(x + width/2, ens_means, width, label=f'Best Ensemble ({best_ensemble})',
                   color=COLORMAP.get(best_ensemble, COLORMAP['Ensemble']),
                   edgecolor='white', linewidth=0.5)
    ax.errorbar(x + width/2, ens_means, yerr=np.array(ens_cis).T, fmt='none',
                color='black', capsize=3, linewidth=1)

    ax.set_ylabel('Score')
    ax.set_title('NB (3-feature) vs Best Ensemble Comparison\n(95% BCa CI)')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, rotation=45, ha='right')
    ax.legend(loc='lower right')
    ax.set_ylim(0, 1.05)
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='Chance')

    # Add significance markers
    for i, m in enumerate(metrics):
        nb_val = summaries['NB'][m]['mean']
        ens_val = summaries[best_ensemble][m]['mean']
        nb_ci = (summaries['NB'][m]['ci_lower'], summaries['NB'][m]['ci_upper'])
        ens_ci = (summaries[best_ensemble][m]['ci_lower'], summaries[best_ensemble][m]['ci_upper'])

        # Check if CIs overlap
        if nb_ci[1] < ens_ci[0] or ens_ci[1] < nb_ci[0]:
            max_y = max(nb_val + nb_cis[i][1], ens_val + ens_cis[i][1])
            ax.text(i, max_y + 0.02, '*', ha='center', fontsize=14, color='red')

    plt.tight_layout()
    savefig('ensemble_comparison_bars', fig)

    store_graph_data('Ensemble_Comparison_Bars', {
        'Metric': metrics,
        'NB_Mean': nb_means,
        'NB_CI_Low': [summaries['NB'][m]['ci_lower'] for m in metrics],
        'NB_CI_High': [summaries['NB'][m]['ci_upper'] for m in metrics],
        'Ensemble_Mean': ens_means,
        'Ensemble_CI_Low': [summaries[best_ensemble][m]['ci_lower'] for m in metrics],
        'Ensemble_CI_High': [summaries[best_ensemble][m]['ci_upper'] for m in metrics],
    })


def plot_comprehensive_metrics_comparison(summaries: Dict, dfs: Dict, best_ensemble: str):
    """Multi-panel bar chart comparing metrics with significance brackets (matches COX_vs_ML style)."""
    metrics = PLOT_METRICS
    n_metrics = len(metrics)

    # Calculate grid dimensions
    n_cols = 3
    n_rows = (n_metrics + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4 * n_rows))
    axes = axes.flatten()

    # Compute pairwise statistics
    pairwise_stats = compute_pairwise_statistics(dfs, best_ensemble, 'NB', metrics)

    for idx, metric in enumerate(metrics):
        ax = axes[idx]

        # Get means and CIs
        nb_mean = summaries['NB'][metric]['mean']
        nb_ci_low = summaries['NB'][metric]['ci_lower']
        nb_ci_high = summaries['NB'][metric]['ci_upper']

        ens_mean = summaries[best_ensemble][metric]['mean']
        ens_ci_low = summaries[best_ensemble][metric]['ci_lower']
        ens_ci_high = summaries[best_ensemble][metric]['ci_upper']

        # Bar positions
        x_pos = [0, 1]
        means = [nb_mean, ens_mean]
        colors = [COLORMAP['NB'], COLORMAP.get(best_ensemble, COLORMAP['Ensemble'])]
        ci_errors = [[nb_mean - nb_ci_low, ens_mean - ens_ci_low],
                     [nb_ci_high - nb_mean, ens_ci_high - ens_mean]]

        # Plot bars
        bars = ax.bar(x_pos, means, color=colors, edgecolor='white', linewidth=0.5, width=0.6)
        ax.errorbar(x_pos, means, yerr=ci_errors, fmt='none', color='black', capsize=5, linewidth=1.5)

        # Add mean value labels
        for i, (bar, mean) in enumerate(zip(bars, means)):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                    f'{mean:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

        # Add 95% CI text
        ax.text(0, nb_ci_low - 0.05, f'95% CI: {nb_ci_low:.3f}-{nb_ci_high:.3f}',
                ha='center', va='top', fontsize=8, color='black', fontweight='bold')
        ax.text(1, ens_ci_low - 0.05, f'95% CI: {ens_ci_low:.3f}-{ens_ci_high:.3f}',
                ha='center', va='top', fontsize=8, color='black', fontweight='bold')

        # Add significance bracket
        stats_result = pairwise_stats.get(metric, {})
        if stats_result:
            sig_level = stats_result.get('sig_level', '')
            if sig_level:
                max_y = max(nb_ci_high, ens_ci_high) + 0.08
                # Draw bracket
                ax.plot([0, 0, 1, 1], [max_y - 0.02, max_y, max_y, max_y - 0.02],
                        color='black', linewidth=1)
                ax.text(0.5, max_y + 0.01, sig_level, ha='center', va='bottom',
                        fontsize=12, fontweight='bold')

        ax.set_xticks(x_pos)
        ax.set_xticklabels(['NB\n(3-feat)', f'{best_ensemble}'], fontsize=9)
        ax.set_ylabel('Score')
        ax.set_title(metric.replace('_', ' '), fontweight='bold')
        ax.set_ylim(0, 1.15)

    # Remove empty subplots
    for idx in range(len(metrics), len(axes)):
        fig.delaxes(axes[idx])

    n_comparisons = len(metrics)
    fig.suptitle(f'Bootstrap CV: NB vs {best_ensemble}\n'
                 f'Bonferroni-corrected: *** p<0.001/{n_comparisons}, '
                 f'** p<0.01/{n_comparisons}, * p<0.05/{n_comparisons}, ns: not significant',
                 fontsize=11, y=1.02)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    savefig('comprehensive_metrics_comparison', fig)


def plot_all_ensemble_strategies(summaries: Dict):
    """Heatmap comparing all ensemble strategies across metrics with mean and std."""
    fig, ax = plt.subplots(figsize=(16, 12))

    models = ['NB'] + ENSEMBLE_STRATEGIES
    metrics = PLOT_METRICS

    data = np.zeros((len(models), len(metrics)))
    stds = np.zeros((len(models), len(metrics)))
    for i, model in enumerate(models):
        for j, metric in enumerate(metrics):
            data[i, j] = summaries[model][metric]['mean']
            stds[i, j] = summaries[model][metric]['std']

    im = ax.imshow(data, cmap='RdYlGn', aspect='auto', vmin=0.5, vmax=1.0)

    ax.set_xticks(np.arange(len(metrics)))
    ax.set_yticks(np.arange(len(models)))
    ax.set_xticklabels(metrics, rotation=45, ha='right', fontsize=10)

    model_labels = ['NB (3-feat)'] + [s.replace('_', ' ').title() for s in ENSEMBLE_STRATEGIES]
    ax.set_yticklabels(model_labels, fontsize=10)

    # Add text annotations with mean +/- std
    for i in range(len(models)):
        for j in range(len(metrics)):
            val = data[i, j]
            std = stds[i, j]
            color = 'white' if val < 0.7 else 'black'
            # Show mean on first line, +/- std on second line
            ax.text(j, i - 0.15, f'{val:.3f}', ha='center', va='center',
                   color=color, fontsize=9, fontweight='bold')
            ax.text(j, i + 0.15, f'\u00b1{std:.3f}', ha='center', va='center',
                   color=color, fontsize=7)

    plt.colorbar(im, ax=ax, label='Score')
    ax.set_title('All Ensemble Strategies Comparison\n(Mean \u00b1 Std)', fontsize=14, fontweight='bold')

    plt.tight_layout()
    savefig('all_ensemble_strategies_heatmap', fig)

    # Store data with CI
    df_data = {'Model': model_labels}
    for j, metric in enumerate(metrics):
        df_data[f'{metric}_mean'] = data[:, j].tolist()
        df_data[f'{metric}_std'] = stds[:, j].tolist()
    store_graph_data('All_Ensemble_Strategies', df_data)


def plot_ensemble_ranking(summaries: Dict):
    """Rank all models by key metrics."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    key_metrics = ['ROC_AUC', 'PRC_AUC', 'F1', 'Balanced_Accuracy']
    models = ['NB'] + ENSEMBLE_STRATEGIES

    for ax, metric in zip(axes.flat, key_metrics):
        values = [(m, summaries[m][metric]['mean'],
                   summaries[m][metric]['ci_lower'],
                   summaries[m][metric]['ci_upper']) for m in models]
        values.sort(key=lambda x: x[1], reverse=True)

        model_names = [v[0] for v in values]
        means = [v[1] for v in values]
        ci_lows = [v[2] for v in values]
        ci_highs = [v[3] for v in values]

        colors = [COLORMAP.get(m, '#999999') for m in model_names]

        y_pos = np.arange(len(model_names))
        bars = ax.barh(y_pos, means, color=colors, edgecolor='white', height=0.6)
        ax.errorbar(means, y_pos, xerr=[np.array(means) - np.array(ci_lows),
                                         np.array(ci_highs) - np.array(means)],
                    fmt='none', color='black', capsize=3)

        ax.set_yticks(y_pos)
        labels = ['NB (3-feat)' if m == 'NB' else m.replace('_', ' ').title()
                  for m in model_names]
        ax.set_yticklabels(labels)
        ax.set_xlabel('Score')
        ax.set_title(metric.replace('_', ' '))
        ax.set_xlim(0.5, 1.0)
        ax.axvline(x=0.5, color='gray', linestyle='--', alpha=0.5)

        # Add value labels
        for i, (m, ci_h) in enumerate(zip(means, ci_highs)):
            ax.text(ci_h + 0.01, i, f'{m:.3f}', va='center', fontsize=9)

    plt.suptitle('Model Ranking by Key Metrics\n(95% BCa CI)', fontsize=14, y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    savefig('ensemble_ranking', fig)


def plot_bootstrap_distributions_ensemble(dfs: Dict, best_ensemble: str):
    """Violin plots of bootstrap distributions for NB vs best ensemble."""
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))

    for ax, metric in zip(axes.flat, PLOT_METRICS):
        data = pd.DataFrame({
            'NB': dfs['NB'][metric].values,
            best_ensemble: dfs[best_ensemble][metric].values
        })
        data_melted = data.melt(var_name='Model', value_name='Score')

        sns.violinplot(x='Model', y='Score', data=data_melted, ax=ax,
                       palette=[COLORMAP['NB'], COLORMAP.get(best_ensemble, COLORMAP['Ensemble'])],
                       inner='box')

        ax.set_title(metric.replace('_', ' '))
        ax.set_xlabel('')
        ax.set_ylim(0, 1)

    plt.suptitle(f'Bootstrap Distributions: NB vs {best_ensemble}', fontsize=14, y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    savefig('bootstrap_distributions_ensemble', fig)


def plot_roc_comparison_ensemble(all_y_true: np.ndarray, all_y_proba: Dict,
                                  best_ensemble: str, bootstrap_roc: Dict = None):
    """ROC curves comparing NB vs best ensemble with DeLong test and 95% CI bands.

    Uses pre-computed bootstrap curves with mean +/- 1.96*std CI bands.
    """
    fig, ax = plt.subplots(figsize=(10, 10))

    # Common FPR grid for interpolation
    mean_fpr = np.linspace(0, 1, 100)

    # Store results for legend
    legend_handles = []
    auc_results = {}

    for model_name, proba_key, color in [('NB (3-feat)', 'NB', COLORMAP['NB']),
                                          (best_ensemble, best_ensemble,
                                           COLORMAP.get(best_ensemble, COLORMAP['Ensemble']))]:
        proba = all_y_proba[proba_key]

        # Main ROC curve and AUC from all predictions
        auc_val = roc_auc_score(all_y_true, proba)

        # Process stored bootstrap ROC curves
        tprs = []
        if bootstrap_roc and proba_key in bootstrap_roc:
            for fpr_boot, tpr_boot in bootstrap_roc[proba_key]:
                tpr_interp = np.interp(mean_fpr, fpr_boot, tpr_boot)
                tpr_interp[0] = 0.0
                tprs.append(tpr_interp)

        if tprs:
            tprs = np.array(tprs)
            mean_tpr = np.mean(tprs, axis=0)
            mean_tpr[-1] = 1.0
            std_tpr = np.std(tprs, axis=0)
            tpr_upper = np.clip(mean_tpr + 1.96 * std_tpr, 0, 1)
            tpr_lower = np.clip(mean_tpr - 1.96 * std_tpr, 0, 1)

            # Plot CI band first
            ci_band = ax.fill_between(mean_fpr, tpr_lower, tpr_upper,
                                       color=color, alpha=0.2, label=f'{model_name} 95% CI')

            # Plot mean curve
            line, = ax.plot(mean_fpr, mean_tpr, color=color, lw=2.5,
                           label=f'{model_name} AUC={auc_val:.3f}')

            auc_results[proba_key] = {'auc': auc_val, 'mean_tpr': mean_tpr, 'std_tpr': std_tpr}
        else:
            # Fallback: plot raw ROC curve without CI
            fpr, tpr, _ = roc_curve(all_y_true, proba)
            line, = ax.plot(fpr, tpr, color=color, lw=2.5,
                           label=f'{model_name} AUC={auc_val:.3f}')
            auc_results[proba_key] = {'auc': auc_val}

    # Diagonal
    ax.plot([0, 1], [0, 1], 'k--', lw=1, label='Random (AUC=0.5)')

    # DeLong test
    delong = delong_test(all_y_true, all_y_proba[best_ensemble], all_y_proba['NB'])

    # Add DeLong test box
    delong_text = (f"DeLong Test\n"
                   f"dAUC={delong['difference']:.3f}\n"
                   f"z={delong['z_score']:.2f}\n"
                   f"p={delong['p_value']:.4f}")
    ax.text(0.55, 0.15, delong_text, transform=ax.transAxes, fontsize=11,
            verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))

    ax.set_xlabel('False Positive Rate (1 - Specificity)', fontsize=12)
    ax.set_ylabel('True Positive Rate (Sensitivity)', fontsize=12)
    ax.set_title('ROC Curves with 95% Bootstrap Confidence Bands', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_aspect('equal')

    plt.tight_layout()
    savefig('roc_comparison_ensemble', fig)

    return delong


def plot_prc_comparison_ensemble(all_y_true: np.ndarray, all_y_proba: Dict,
                                  best_ensemble: str, bootstrap_prc: Dict = None):
    """PRC curves comparing NB vs best ensemble with 95% CI bands.

    Uses pre-computed bootstrap curves with mean +/- 1.96*std CI bands.
    """
    fig, ax = plt.subplots(figsize=(10, 10))

    # Common recall grid for interpolation
    mean_recall = np.linspace(0, 1, 100)

    for model_name, proba_key, color in [('NB (3-feat)', 'NB', COLORMAP['NB']),
                                          (best_ensemble, best_ensemble,
                                           COLORMAP.get(best_ensemble, COLORMAP['Ensemble']))]:
        proba = all_y_proba[proba_key]

        # Main AP from all predictions
        ap_val = average_precision_score(all_y_true, proba)

        # Process stored bootstrap PRC curves
        precs = []
        if bootstrap_prc and proba_key in bootstrap_prc:
            for rec_boot, prec_boot in bootstrap_prc[proba_key]:
                # Reverse for interpolation (recall is decreasing in original)
                prec_interp = np.interp(mean_recall, rec_boot[::-1], prec_boot[::-1])
                precs.append(prec_interp)

        if precs:
            precs = np.array(precs)
            mean_prec = np.mean(precs, axis=0)
            std_prec = np.std(precs, axis=0)
            prec_upper = np.clip(mean_prec + 1.96 * std_prec, 0, 1)
            prec_lower = np.clip(mean_prec - 1.96 * std_prec, 0, 1)

            # Plot CI band first
            ax.fill_between(mean_recall, prec_lower, prec_upper,
                           color=color, alpha=0.2, label=f'{model_name} 95% CI')

            # Plot mean curve
            ax.plot(mean_recall, mean_prec, color=color, lw=2.5,
                   label=f'{model_name} AP={ap_val:.3f}')
        else:
            # Fallback: plot raw PRC curve without CI
            prec, rec, _ = precision_recall_curve(all_y_true, proba)
            ax.plot(rec, prec, color=color, lw=2.5,
                   label=f'{model_name} AP={ap_val:.3f}')

    # Baseline (no-skill classifier)
    baseline = np.mean(all_y_true)
    ax.axhline(y=baseline, color='gray', linestyle='--', linewidth=1,
              label=f'No Skill (Prevalence={baseline:.2f})')

    ax.set_xlabel('Recall (Sensitivity)', fontsize=12)
    ax.set_ylabel('Precision (PPV)', fontsize=12)
    ax.set_title('Precision-Recall Curves with 95% Bootstrap Confidence Bands',
                fontsize=14, fontweight='bold')
    ax.legend(loc='lower left', fontsize=10)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

    plt.tight_layout()
    savefig('prc_comparison_ensemble', fig)


def plot_calibration_comparison(all_y_true: np.ndarray, all_y_proba: Dict,
                                 best_ensemble: str):
    """Calibration curves for NB vs best ensemble."""
    fig, ax = plt.subplots(figsize=(10, 9))

    # NB calibration
    prob_true_nb, prob_pred_nb = calibration_curve(all_y_true, all_y_proba['NB'], n_bins=10)
    ax.plot(prob_pred_nb, prob_true_nb, 's-', color=COLORMAP['NB'], lw=2,
            label='NB (3-feat)', markersize=8)

    # Best ensemble calibration
    prob_true_ens, prob_pred_ens = calibration_curve(all_y_true, all_y_proba[best_ensemble], n_bins=10)
    ax.plot(prob_pred_ens, prob_true_ens, 'o-',
            color=COLORMAP.get(best_ensemble, COLORMAP['Ensemble']), lw=2,
            label=f'{best_ensemble}', markersize=8)

    # Perfect calibration
    ax.plot([0, 1], [0, 1], 'k--', lw=1, label='Perfect calibration')

    # Brier scores
    brier_nb = brier_score_loss(all_y_true, all_y_proba['NB'])
    brier_ens = brier_score_loss(all_y_true, all_y_proba[best_ensemble])

    ax.set_xlabel('Mean Predicted Probability')
    ax.set_ylabel('Fraction of Positives')
    ax.set_title(f'Calibration Curves\nBrier: NB = {brier_nb:.4f}, {best_ensemble} = {brier_ens:.4f}')
    ax.legend(loc='lower right')
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)

    plt.tight_layout()
    savefig('calibration_comparison_ensemble', fig)


def plot_decision_curve_analysis(all_y_true: np.ndarray, all_y_proba: Dict,
                                  best_ensemble: str):
    """Decision curve analysis for NB vs best ensemble."""
    fig, ax = plt.subplots(figsize=(10, 8))

    thresholds = np.linspace(0.01, 0.99, 100)
    n = len(all_y_true)
    prevalence = np.mean(all_y_true)

    net_benefit_nb = []
    net_benefit_ens = []
    net_benefit_all = []
    net_benefit_none = []

    for thresh in thresholds:
        # Treat all
        nb_all = prevalence - (1 - prevalence) * thresh / (1 - thresh)
        net_benefit_all.append(nb_all)
        net_benefit_none.append(0)

        # NB
        pred_nb = (all_y_proba['NB'] >= thresh).astype(int)
        tp_nb = np.sum((pred_nb == 1) & (all_y_true == 1))
        fp_nb = np.sum((pred_nb == 1) & (all_y_true == 0))
        nb_nb = tp_nb / n - fp_nb / n * thresh / (1 - thresh)
        net_benefit_nb.append(nb_nb)

        # Ensemble
        pred_ens = (all_y_proba[best_ensemble] >= thresh).astype(int)
        tp_ens = np.sum((pred_ens == 1) & (all_y_true == 1))
        fp_ens = np.sum((pred_ens == 1) & (all_y_true == 0))
        nb_ens = tp_ens / n - fp_ens / n * thresh / (1 - thresh)
        net_benefit_ens.append(nb_ens)

    ax.plot(thresholds, net_benefit_nb, color=COLORMAP['NB'], lw=2, label='NB (3-feat)')
    ax.plot(thresholds, net_benefit_ens,
            color=COLORMAP.get(best_ensemble, COLORMAP['Ensemble']), lw=2,
            label=f'{best_ensemble}')
    ax.plot(thresholds, net_benefit_all, 'k--', lw=1, label='Treat All')
    ax.plot(thresholds, net_benefit_none, 'k:', lw=1, label='Treat None')

    ax.set_xlabel('Threshold Probability')
    ax.set_ylabel('Net Benefit')
    ax.set_title('Decision Curve Analysis')
    ax.legend(loc='upper right')
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.1, max(prevalence, max(net_benefit_nb), max(net_benefit_ens)) + 0.05)

    plt.tight_layout()
    savefig('decision_curve_ensemble', fig)

    store_graph_data('Decision_Curve_Ensemble', {
        'Threshold': thresholds.tolist(),
        'NB_Net_Benefit': net_benefit_nb,
        'Ensemble_Net_Benefit': net_benefit_ens,
        'Treat_All': net_benefit_all
    })


def plot_nri_idi_analysis(nri_idi: Dict, best_ensemble: str):
    """Visualize NRI and IDI results."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # NRI
    ax = axes[0]
    categories = ['Events\nReclassified', 'Non-events\nReclassified', 'Total NRI']
    values = [nri_idi['NRI_events'], nri_idi['NRI_non_events'], nri_idi['NRI']]
    colors = ['#28A745' if v > 0 else '#DC3545' for v in values]

    bars = ax.bar(categories, values, color=colors, edgecolor='white', width=0.6)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax.set_ylabel('NRI')
    ax.set_title(f'Net Reclassification Improvement\n{best_ensemble} vs NB\np = {nri_idi["NRI_p"]:.4f}')

    for bar, val in zip(bars, values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.01 * np.sign(height),
                f'{val:.3f}', ha='center', va='bottom' if height > 0 else 'top')

    # IDI
    ax = axes[1]
    ax.bar(['IDI'], [nri_idi['IDI']],
           color='#28A745' if nri_idi['IDI'] > 0 else '#DC3545',
           edgecolor='white', width=0.4)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax.set_ylabel('IDI')
    ax.set_title(f'Integrated Discrimination Improvement\n{best_ensemble} vs NB\np = {nri_idi["IDI_p"]:.4f}')
    ax.text(0, nri_idi['IDI'] + 0.005 * np.sign(nri_idi['IDI']),
            f'{nri_idi["IDI"]:.4f}', ha='center')

    plt.tight_layout()
    savefig('nri_idi_ensemble', fig)


def plot_optimal_weights_distribution(optimal_weights: np.ndarray):
    """Distribution of optimal weights from grid search."""
    if optimal_weights is None or len(optimal_weights) == 0:
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(optimal_weights, bins=len(WEIGHT_GRID), color=COLORMAP['optimized'],
            edgecolor='white', alpha=0.8)
    ax.axvline(x=np.mean(optimal_weights), color='red', linestyle='--', lw=2,
               label=f'Mean: {np.mean(optimal_weights):.3f}')
    ax.axvline(x=np.median(optimal_weights), color='blue', linestyle=':', lw=2,
               label=f'Median: {np.median(optimal_weights):.3f}')

    ax.set_xlabel('NB Weight (Cox weight = 1 - NB weight)')
    ax.set_ylabel('Frequency')
    ax.set_title('Distribution of Optimal NB:Cox Weights\n(Grid Search across Bootstrap CV)')
    ax.legend()

    plt.tight_layout()
    savefig('optimal_weights_distribution', fig)

    store_graph_data('Optimal_Weights', {
        'weights': optimal_weights.tolist(),
        'mean': np.mean(optimal_weights),
        'median': np.median(optimal_weights),
        'std': np.std(optimal_weights)
    })


def plot_confusion_matrices(all_y_true: np.ndarray, all_y_proba: Dict, best_ensemble: str):
    """Side-by-side confusion matrices."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, (name, proba) in zip(axes, [('NB (3-feat)', all_y_proba['NB']),
                                          (best_ensemble, all_y_proba[best_ensemble])]):
        y_pred = (proba >= 0.5).astype(int)
        cm = confusion_matrix(all_y_true, y_pred)
        cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)

        sns.heatmap(cm_norm, annot=True, fmt='.2%', cmap='Blues', ax=ax,
                    xticklabels=CLASS_NAMES, yticklabels=CLASS_NAMES,
                    cbar=True, vmin=0, vmax=1)

        # Add raw counts
        for i in range(2):
            for j in range(2):
                ax.text(j + 0.5, i + 0.7, f'(n={cm[i, j]})', ha='center', va='center',
                        fontsize=9, color='gray')

        ax.set_xlabel('Predicted')
        ax.set_ylabel('Actual')
        ax.set_title(name)

    plt.suptitle('Confusion Matrices (Row-Normalized)', fontsize=14, y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    savefig('confusion_matrices_ensemble', fig)


def plot_km_risk_stratification(data: pd.DataFrame, X_nb: np.ndarray, X_cox: np.ndarray,
                                 y: np.ndarray, best_ensemble: str, optimal_weights: np.ndarray = None):
    """Kaplan-Meier curves for risk stratification by NB vs best ensemble.

    Trains final models on full data to get proper per-patient predictions
    aligned with survival data, instead of using misaligned CV predictions.
    """
    if not HAS_LIFELINES:
        log("Lifelines not available for KM analysis", "WARNING")
        return None

    # Train final models on full data to get per-patient predictions
    log("  Training final models for KM risk stratification...")

    # NB model (already calibrated internally)
    nb_clf = NBChampionClassifier(features=NB_FEATURES, calibrate=True)
    nb_clf.fit(X_nb, y)
    proba_nb = nb_clf.predict_proba(X_nb)[:, 1]

    # Best ensemble model with isotonic calibration
    ens_clf = EnsembleClassifier(strategy=best_ensemble)
    ens_clf.fit(X_nb, X_cox, y)

    # For optimized strategy, use mean optimal weights from CV
    if best_ensemble == 'optimized' and optimal_weights is not None:
        mean_weights = np.mean(optimal_weights, axis=0)
        ens_clf.optimal_weight_ = mean_weights[0]  # NB weight

    # Get uncalibrated ensemble predictions
    proba_ens_uncal = ens_clf.predict_proba(X_nb, X_cox)[:, 1]

    # Apply isotonic calibration to ensemble predictions
    calibrator = IsotonicRegression(out_of_bounds='clip')
    calibrator.fit(proba_ens_uncal, y)
    proba_ens = calibrator.predict(proba_ens_uncal)
    proba_ens = np.clip(proba_ens, 0.001, 0.999)

    log(f"  NB predictions: mean={np.mean(proba_nb):.3f}, std={np.std(proba_nb):.3f}")
    log(f"  {best_ensemble} predictions (calibrated): mean={np.mean(proba_ens):.3f}, std={np.std(proba_ens):.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    results = {}
    predictions = {'NB': proba_nb, best_ensemble: proba_ens}

    for ax, (name, proba_key) in zip(axes, [('NB (3-feat)', 'NB'),
                                              (best_ensemble, best_ensemble)]):
        proba = predictions[proba_key]

        # Stratify by 0.5 threshold (consistent with COX_vs_ML comparison)
        threshold = 0.5
        high_risk = proba >= threshold
        low_risk = ~high_risk

        # Convert boolean index to integer index for iloc
        high_idx = np.where(high_risk)[0]
        low_idx = np.where(low_risk)[0]

        kmf_high = KaplanMeierFitter()
        kmf_low = KaplanMeierFitter()

        # High risk
        time_high = data.iloc[high_idx][TIME_COLUMN].values
        event_high = data.iloc[high_idx][EVENT_COLUMN].values
        kmf_high.fit(time_high, event_high, label='High Risk')
        kmf_high.plot_survival_function(ax=ax, ci_show=True, color='#DC3545')

        # Low risk
        time_low = data.iloc[low_idx][TIME_COLUMN].values
        event_low = data.iloc[low_idx][EVENT_COLUMN].values
        kmf_low.fit(time_low, event_low, label='Low Risk')
        kmf_low.plot_survival_function(ax=ax, ci_show=True, color='#28A745')

        # Log-rank test
        lr_result = logrank_test(time_high, time_low, event_high, event_low)

        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Survival Probability')
        ax.set_title(f'{name}\nLog-rank p = {lr_result.p_value:.4f}')
        ax.legend(loc='lower left')

        # Add hazard ratio if possible
        try:
            cox_df = pd.DataFrame({
                'time': np.concatenate([time_high, time_low]),
                'event': np.concatenate([event_high, event_low]),
                'high_risk': np.concatenate([np.ones(len(time_high)), np.zeros(len(time_low))])
            })
            cph = CoxPHFitter()
            cph.fit(cox_df, duration_col='time', event_col='event')
            hr = np.exp(cph.params_['high_risk'])
            hr_ci = np.exp(cph.confidence_intervals_.loc['high_risk'])
            ax.text(0.95, 0.95, f'HR = {hr:.2f}\n(95% CI: {hr_ci.iloc[0]:.2f}-{hr_ci.iloc[1]:.2f})',
                   transform=ax.transAxes, fontsize=10, va='top', ha='right',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        except Exception as e:
            log(f"  Warning: Could not compute HR for {name}: {e}")

        results[proba_key] = {
            'n_high': np.sum(high_risk),
            'n_low': np.sum(low_risk),
            'median_high': kmf_high.median_survival_time_,
            'median_low': kmf_low.median_survival_time_,
            'logrank_p': lr_result.p_value
        }

        log(f"  {name}: n_high={results[proba_key]['n_high']}, n_low={results[proba_key]['n_low']}, "
            f"logrank p={lr_result.p_value:.4f}")

    plt.suptitle('Risk Stratification: Kaplan-Meier Survival Curves', fontsize=14, y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    savefig('km_risk_stratification_ensemble', fig)

    return results


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Main execution function."""

    # Clear log file
    with open(LOG_PATH, 'w') as f:
        f.write("")

    log("=" * 70)
    log("ENSEMBLE NB3 + COX COMPARISON ANALYSIS")
    log("=" * 70)
    log(f"Output directory: {OUTPUT_DIR}")
    log(f"Bootstrap iterations: {N_BOOTSTRAP}")
    log(f"CV configuration: {CV_N_REPEATS}x{CV_N_FOLDS}-fold")
    log(f"Total evaluations: {TOTAL_EVALUATIONS}")
    log(f"Ensemble strategies: {ENSEMBLE_STRATEGIES}")

    # -------------------------------------------------------------------------
    # 1. LOAD DATA
    # -------------------------------------------------------------------------
    log("\n[1/8] Loading data...")

    data = pd.read_excel(FULL_DATASET)
    log(f"  Loaded dataset: {data.shape}")

    # Prepare features
    X_nb = data[NB_FEATURES].values
    X_cox = data[COX_FEATURES].values
    y = data[LABEL_COLUMN].map(CLASS_MAP).values
    patient_ids = data.index.values

    log(f"  NB features: {NB_FEATURES}")
    log(f"  Cox features: {COX_FEATURES}")
    log(f"  Class distribution: {np.sum(y==0)} Long, {np.sum(y==1)} Short")

    # -------------------------------------------------------------------------
    # 2. BOOTSTRAP CV EVALUATION
    # -------------------------------------------------------------------------
    log("\n[2/8] Running Bootstrap CV evaluation for all models...")

    cv_results = run_bootstrap_cv_all_ensembles(
        X_nb, X_cox, y, patient_ids,
        n_bootstrap=N_BOOTSTRAP,
        n_repeats=CV_N_REPEATS,
        n_folds=CV_N_FOLDS,
        seed=RANDOM_STATE
    )

    summaries = cv_results['summaries']
    dfs = cv_results['dataframes']
    all_y_true = cv_results['all_y_true']
    all_y_proba = cv_results['all_y_proba']
    optimal_weights = cv_results['optimal_weights']
    bootstrap_roc = cv_results['bootstrap_roc']
    bootstrap_prc = cv_results['bootstrap_prc']

    # -------------------------------------------------------------------------
    # 3. DETERMINE BEST ENSEMBLE
    # -------------------------------------------------------------------------
    log("\n[3/8] Determining best ensemble strategy...")

    # Rank by ROC-AUC
    ensemble_aucs = {s: summaries[s]['ROC_AUC']['mean'] for s in ENSEMBLE_STRATEGIES}
    best_ensemble = max(ensemble_aucs, key=ensemble_aucs.get)

    log(f"  Ensemble ROC-AUC rankings:")
    for s, auc in sorted(ensemble_aucs.items(), key=lambda x: x[1], reverse=True):
        log(f"    {s}: {auc:.4f}")
    log(f"  Best ensemble: {best_ensemble} (AUC = {ensemble_aucs[best_ensemble]:.4f})")
    log(f"  NB (3-feat) AUC: {summaries['NB']['ROC_AUC']['mean']:.4f}")

    # -------------------------------------------------------------------------
    # 4. STATISTICAL COMPARISONS
    # -------------------------------------------------------------------------
    log("\n[4/8] Running statistical comparisons...")

    # DeLong test: Best ensemble vs NB
    delong_result = delong_test(all_y_true, all_y_proba[best_ensemble], all_y_proba['NB'])
    log(f"  DeLong test ({best_ensemble} vs NB):")
    log(f"    ΔAUC = {delong_result['difference']:.4f}")
    log(f"    p-value = {delong_result['p_value']:.6f}")
    log(f"    Significant: {delong_result['significant']}")

    # NRI/IDI: Best ensemble vs NB (ensemble as "new", NB as "old")
    nri_idi_result = compute_nri_idi(all_y_true, all_y_proba['NB'],
                                      all_y_proba[best_ensemble])
    log(f"  NRI/IDI ({best_ensemble} vs NB):")
    log(f"    NRI = {nri_idi_result['NRI']:.4f} (p = {nri_idi_result['NRI_p']:.4f})")
    log(f"    IDI = {nri_idi_result['IDI']:.4f} (p = {nri_idi_result['IDI_p']:.4f})")

    # Comprehensive pairwise statistics
    log(f"\n  Paired t-test with Bonferroni correction ({best_ensemble} vs NB):")
    pairwise_stats = compute_pairwise_statistics(dfs, best_ensemble, 'NB', PLOT_METRICS)
    for metric, stats_dict in pairwise_stats.items():
        sig_marker = stats_dict['sig_level']
        d = stats_dict['cohens_d']
        p = stats_dict['p_value']
        log(f"    {metric}: d={stats_dict['difference']:.4f}, d={d:.3f} ({stats_dict['effect_interpretation']}), "
            f"p={p:.4f} {sig_marker}")

    # -------------------------------------------------------------------------
    # 5. GENERATE VISUALIZATIONS
    # -------------------------------------------------------------------------
    log("\n[5/8] Generating visualizations...")

    # Main comparison plots
    plot_ensemble_comparison_bars(summaries, best_ensemble)
    plot_comprehensive_metrics_comparison(summaries, dfs, best_ensemble)
    plot_all_ensemble_strategies(summaries)
    plot_ensemble_ranking(summaries)
    plot_bootstrap_distributions_ensemble(dfs, best_ensemble)

    # Head-to-head comparison with significance brackets
    plot_head_to_head_comparison(dfs, summaries, best_ensemble, 'NB', pairwise_stats)

    # ROC/PRC curves with bootstrap CI bands
    plot_roc_comparison_ensemble(all_y_true, all_y_proba, best_ensemble, bootstrap_roc)
    plot_prc_comparison_ensemble(all_y_true, all_y_proba, best_ensemble, bootstrap_prc)

    # Calibration and decision curves
    plot_calibration_comparison(all_y_true, all_y_proba, best_ensemble)
    plot_decision_curve_analysis(all_y_true, all_y_proba, best_ensemble)

    # NRI/IDI visualization
    plot_nri_idi_analysis(nri_idi_result, best_ensemble)

    # Confusion matrices
    plot_confusion_matrices(all_y_true, all_y_proba, best_ensemble)

    # Optimal weights distribution
    if optimal_weights is not None:
        plot_optimal_weights_distribution(optimal_weights)

    # KM risk stratification using per-patient final model predictions
    km_results = plot_km_risk_stratification(data, X_nb, X_cox, y, best_ensemble, optimal_weights)

    # -------------------------------------------------------------------------
    # 6. CREATE SUMMARY TABLES
    # -------------------------------------------------------------------------
    log("\n[6/8] Creating summary tables...")

    # Summary metrics table
    summary_rows = []
    for model in ['NB'] + ENSEMBLE_STRATEGIES:
        row = {'Model': model}
        for metric in PLOT_METRICS:
            row[f'{metric}_Mean'] = summaries[model][metric]['mean']
            row[f'{metric}_CI_Low'] = summaries[model][metric]['ci_lower']
            row[f'{metric}_CI_High'] = summaries[model][metric]['ci_upper']
        summary_rows.append(row)
    summary_df = pd.DataFrame(summary_rows)

    # DeLong results for all ensembles vs NB
    delong_rows = []
    for strategy in ENSEMBLE_STRATEGIES:
        dl = delong_test(all_y_true, all_y_proba[strategy], all_y_proba['NB'])
        delong_rows.append({
            'Ensemble': strategy,
            'AUC_Ensemble': dl['auc1'],
            'AUC_NB': dl['auc2'],
            'Delta_AUC': dl['difference'],
            'Z_Score': dl['z_score'],
            'P_Value': dl['p_value'],
            'Significant': dl['significant']
        })
    delong_df = pd.DataFrame(delong_rows)

    # -------------------------------------------------------------------------
    # 7. EXPORT RESULTS
    # -------------------------------------------------------------------------
    log("\n[7/8] Exporting results...")

    sheets = {
        '1_Configuration': pd.DataFrame({
            'Parameter': ['N_BOOTSTRAP', 'CV_N_REPEATS', 'CV_N_FOLDS',
                          'TOTAL_EVALUATIONS', 'NB_FEATURES', 'COX_FEATURES',
                          'BEST_ENSEMBLE'],
            'Value': [N_BOOTSTRAP, CV_N_REPEATS, CV_N_FOLDS, TOTAL_EVALUATIONS,
                      str(NB_FEATURES), str(COX_FEATURES), best_ensemble]
        }),
        '2_Summary_Metrics': summary_df,
        '3_DeLong_Results': delong_df,
        '4_NRI_IDI': pd.DataFrame([nri_idi_result]),
        '5_Optimal_Weights': pd.DataFrame({
            'statistic': ['mean', 'median', 'std', 'min', 'max'],
            'value': [np.mean(optimal_weights), np.median(optimal_weights),
                      np.std(optimal_weights), np.min(optimal_weights),
                      np.max(optimal_weights)] if optimal_weights is not None else [np.nan]*5
        }),
        '6_Pairwise_Statistics': pd.DataFrame([
            {
                'Metric': metric,
                f'{best_ensemble}_Mean': stats_val['mean_1'],
                'NB_Mean': stats_val['mean_2'],
                'Difference': stats_val['difference'],
                'Cohens_d': stats_val['cohens_d'],
                'Effect_Size': stats_val['effect_interpretation'],
                'T_Statistic': stats_val['t_stat'],
                'P_Value': stats_val['p_value'],
                'P_Bonferroni': stats_val['p_bonferroni'],
                'Significant': stats_val['sig_level']
            }
            for metric, stats_val in pairwise_stats.items()
        ])
    }

    # Add bootstrap DataFrames
    for model in ['NB'] + ENSEMBLE_STRATEGIES:
        sheets[f'Bootstrap_{model[:20]}'] = dfs[model]

    # Add graph data
    for name, data_dict in GRAPH_DATA_STORAGE.items():
        try:
            if isinstance(data_dict, dict):
                sheets[f'G_{name[:25]}'] = pd.DataFrame(data_dict)
            elif isinstance(data_dict, pd.DataFrame):
                sheets[f'G_{name[:25]}'] = data_dict
        except:
            pass

    excel_path = os.path.join(OUTPUT_DIR, 'ENSEMBLE_COMPARISON_RESULTS.xlsx')
    save_excel(excel_path, sheets)

    # -------------------------------------------------------------------------
    # 8. FINAL SUMMARY
    # -------------------------------------------------------------------------
    log("\n[8/8] Final Summary")
    log("=" * 70)
    log(f"\nBest Ensemble Strategy: {best_ensemble}")
    log(f"\nPerformance Comparison (ROC-AUC):")
    log(f"  NB (3-feat):    {summaries['NB']['ROC_AUC']['mean']:.4f} "
        f"({summaries['NB']['ROC_AUC']['ci_lower']:.4f}-{summaries['NB']['ROC_AUC']['ci_upper']:.4f})")
    log(f"  {best_ensemble}: {summaries[best_ensemble]['ROC_AUC']['mean']:.4f} "
        f"({summaries[best_ensemble]['ROC_AUC']['ci_lower']:.4f}-{summaries[best_ensemble]['ROC_AUC']['ci_upper']:.4f})")
    log(f"\nDeLong Test: dAUC = {delong_result['difference']:.4f}, p = {delong_result['p_value']:.4f}")
    log(f"NRI = {nri_idi_result['NRI']:.4f}, IDI = {nri_idi_result['IDI']:.4f}")

    if km_results:
        log(f"\nKaplan-Meier Risk Stratification:")
        for model, res in km_results.items():
            log(f"  {model}: High={res['n_high']}({res['median_high']:.0f}d), "
                f"Low={res['n_low']}({res['median_low']:.0f}d), p={res['logrank_p']:.4f}")

    log("\n" + "=" * 70)
    log("Analysis complete!")
    log(f"Total visualizations generated: 12+")
    log(f"Excel file: ENSEMBLE_COMPARISON_RESULTS.xlsx")
    log("=" * 70)


if __name__ == "__main__":
    main()
