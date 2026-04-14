#!/usr/bin/env python3
"""
03_survival.py

Kaplan-Meier survival analysis for multi-gene expression groups.
Splits patients into High/Low for each gene, creates 2^N combo groups,
and plots KM curves with risk tables. Supports 2 or more genes.
"""

import os
import sys
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D
import argparse
import yaml
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests

matplotlib.use("Agg")


# ============================================================
# Configuration
# ============================================================

EXPRESSION_FILE = os.path.join("results", "expression_clean.tsv")
SURVIVAL_FILE = os.path.join("data", "TCGA_master_clinical_survival.csv")
RESULTS_DIR = "results"

# ============================================================
# Helper functions
# ============================================================

def get_group_colors(groups):
    """
    Generate a color dict for N groups using matplotlib's Dark2 palette.
    Supports up to 8 groups (2^3).
    """
    cmap = plt.get_cmap("Dark2")
    return {group: cmap(i % 8) for i, group in enumerate(groups)}


def generate_all_groups(n_genes):
    """
    Return all 2^n possible Low/High combinations as strings.
    E.g., n=2 → ['Low/Low', 'Low/High', 'High/Low', 'High/High']
    """
    from itertools import product
    return ["/".join(combo) for combo in product(["Low", "High"], repeat=n_genes)]

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", required=True,
                        help="Comma-separated gene symbols")
    return parser.parse_args()

with open("config.yaml") as f:
    config = yaml.safe_load(f)

AML_PREFIX = config["aml_prefix"]

def load_and_merge(expression_file, survival_file, genes, sample_filter=None):
    """
    Load expression and survival data, merge them by Sample_ID.
    Returns a DataFrame with columns: Sample_ID, OS.time, OS, expr1, expr2, ..., exprN
    """
    tpm = pd.read_csv(expression_file, sep="\t", index_col=0)
    expr = tpm.loc[genes].T
    expr.columns = [f"expr{i+1}" for i in range(len(genes))]
    expr.index.name = "Sample_ID"
    expr = expr.reset_index()

    surv = pd.read_csv(survival_file)
    if sample_filter:
        surv = surv[surv["Sample_ID"].str.startswith(sample_filter)]
    surv = surv[["Sample_ID", "OS.time", "OS"]].dropna()
    surv = surv[surv["OS.time"] > 0]
    surv["OS"] = surv["OS"].astype(int)

    dat = pd.merge(surv, expr, on="Sample_ID", how="inner")
    return dat

def find_optimal_cutpoint(dat, expr_col, time_col="OS.time", event_col="OS", minprop=0.1):
    """
    Find the expression cutpoint that gives the most significant
    log-rank test split (maximally selected rank statistics).
    """
    values = dat[expr_col].values
    times = dat[time_col].values
    events = dat[event_col].values
    n_total = len(values)

    # Use percentiles of ALL values, not unique values
    low_bound = np.percentile(values, minprop * 100)
    high_bound = np.percentile(values, (1 - minprop) * 100)

    # Get candidate cutpoints within the allowed range
    sorted_vals = np.sort(np.unique(values))
    candidates = sorted_vals[(sorted_vals >= low_bound) & (sorted_vals <= high_bound)]

    best_stat = -1
    best_cut = np.median(values)

    for cut in candidates:
        high = values > cut
        low = ~high

        # Enforce minimum proportion in each group
        # R equivalent: minprop = 0.1 means at least 10% in each group
        if high.sum() < n_total * minprop or low.sum() < n_total * minprop:
            continue

        result = logrank_test(
            times[high], times[low],
            event_observed_A=events[high],
            event_observed_B=events[low],
        )

        if result.test_statistic > best_stat:
            best_stat = result.test_statistic
            best_cut = cut

    return best_cut


def assign_groups(dat, genes, split=None):
    """
    Split patients into High/Low for each gene and create combo groups.
    Works for any number of genes.
    """
    if split is None:
        split = config["split_method"]

    cuts = []
    for i, gene in enumerate(genes):
        expr_col = f"expr{i+1}"
        if split == "optimal":
            cut = find_optimal_cutpoint(dat, expr_col)
        else:
            cut = dat[expr_col].median()
        cuts.append(cut)
        print(f"  {'Optimal' if split == 'optimal' else 'Median'} cutpoint: {gene} = {cut:.3f}")

    # Assign High/Low per gene
    for i in range(len(genes)):
        dat[f"G{i+1}"] = np.where(dat[f"expr{i+1}"] > cuts[i], "High", "Low")

    # Create combo group: "High/Low/High" etc.
    g_cols = [f"G{i+1}" for i in range(len(genes))]
    dat["Combo"] = dat[g_cols].agg("/".join, axis=1)

    print(f"  Group sizes: {dat['Combo'].value_counts().to_dict()}")

    # Warn about underpowered groups (Peduzzi et al. 1995, J Clin Epidemiol)
    MIN_EVENTS = 10
    group_events = dat.groupby("Combo")["OS"].sum()
    underpowered = group_events[group_events < MIN_EVENTS]
    if len(underpowered) > 0:
        print(f"  ⚠ WARNING: Groups with fewer than {MIN_EVENTS} events detected:")
        for group, n_events in underpowered.items():
            print(f"    {group}: {int(n_events)} events")
        print(f"  Results for these groups may be unreliable")
        
    return dat

def compute_logrank_pvalue(dat):
    """
    Compute overall log-rank p-value across all 4 groups.
    """
    from lifelines.statistics import multivariate_logrank_test
    result = multivariate_logrank_test(
        dat["OS.time"], dat["Combo"], dat["OS"]
    )
    return result.p_value


def format_pvalue(pval):
    """Format p-value for display on plot."""
    if pval < 1e-16:
        return "p < 1e-16"
    elif pval < 0.001:
        return f"p = {pval:.2e}"
    else:
        return f"p = {pval:.4f}"


def make_km_plot(dat, genes, title_prefix, output_path,
                 xlim_days=365, break_time=73):
    """
    Create Kaplan-Meier plot with risk table
    """
    groups = generate_all_groups(len(genes))
    group_colors = get_group_colors(groups)

    # --- Compute p-value ---
    pval = compute_logrank_pvalue(dat)
    pval_str = format_pvalue(pval)
    print(f"  Log-rank {pval_str}")

    # --- Set up figure: KM plot on top, risk table on bottom ---
    fig, (ax_km, ax_risk) = plt.subplots(
        2, 1,
        figsize=(10, 11),
        gridspec_kw={"height_ratios": [3.5, 1]},
    )

    # --- Fit and plot KM curves per group ---
    for group in groups:
        mask = dat["Combo"] == group
        group_data = dat[mask]

        if len(group_data) == 0:
            continue

        kmf = KaplanMeierFitter()
        kmf.fit(
            group_data["OS.time"],
            event_observed=group_data["OS"],
            label=group,
        )

        # Plot the KM step curve manually for full control
        times_plot = kmf.survival_function_.index.values
        surv_plot = kmf.survival_function_.iloc[:, 0].values

        # Clip to xlim
        mask_time = times_plot <= xlim_days
        times_plot = times_plot[mask_time]
        surv_plot = surv_plot[mask_time]

        ax_km.step(
            times_plot, surv_plot,
            where="post",
            color=group_colors[group],
            linewidth=2.5,
            label=group,
            zorder=2,
        )

        # Add censor tick marks
        censor_times = group_data.loc[group_data["OS"] == 0, "OS.time"].values
        censor_times = censor_times[censor_times <= xlim_days]

        for ct in censor_times:
            idx = np.searchsorted(times_plot, ct, side="right") - 1
            if idx >= 0 and idx < len(surv_plot):
                surv_at_censor = surv_plot[idx]
                ax_km.plot(
                    ct, surv_at_censor, "|",
                    color=group_colors[group],
                    markersize=8,
                    markeredgewidth=1.5,
                    zorder=3,
                )

    # --- Format KM plot ---
    ax_km.set_ylim(0, 1.02)
    ax_km.set_xlim(0, xlim_days)
    ax_km.set_ylabel("Percent survival", fontsize=25, fontweight="bold")
    ax_km.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax_km.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize = 25)
    ax_km.set_xticks(np.arange(0, xlim_days + 1, break_time))
    ax_km.set_xticklabels(ax_km.get_xticks().astype(int), fontsize=25)
    ax_km.set_xlabel("")

    # Title
    fig.text(
        0.05, 0.94, title_prefix,
        fontsize=32, fontweight="bold",
        verticalalignment="top",
    )
    
    # P-value annotation
    ax_km.text(
        xlim_days * 0.25, 0.92,
        f"Log-rank P {'< 0.0001' if pval < 0.0001 else '= ' + f'{pval:.4f}'}",
        fontsize=18, fontweight="bold",
    )

    # Legend
    legend_elements = [
        Line2D([0], [0], color=group_colors[g], linewidth=2.5, label=g)
        for g in groups
        if (dat["Combo"] == g).any()
    ]
    leg = ax_km.legend(
        handles=legend_elements,
        title=f"$\\it{{{'/'.join(genes)}}}$",
        loc="upper right",
        fontsize=18,
        title_fontsize=19,
        framealpha=0.9,
        edgecolor="none",
    )
    leg.get_title().set_fontstyle("italic")

    # Border
    for spine in ax_km.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(1)

    # --- Risk table ---
    time_points = np.arange(0, xlim_days + 1, break_time)

    # Title for risk table
    ax_risk.text(
        -0.02, 1.25, "Number at risk",
        transform=ax_risk.transAxes,
        fontsize=26, fontweight="bold",
        verticalalignment="top",
    )
                     
    # Gene-order header above group labels
    ax_risk.text(
        -0.06, len(groups) - 0.3,
        "/".join(genes),
        ha="right", va="center",
        fontsize=16, fontweight="bold", fontstyle="italic",
        color="black",
        transform=ax_risk.get_yaxis_transform(),
    )                 

    # Risk table: group labels and counts
    for i, group in enumerate(groups):
        mask = dat["Combo"] == group
        group_data = dat[mask]

        if len(group_data) == 0:
            continue

        # Group label (colored, no dash)
        ax_risk.text(
            -0.06, i, group,
            ha="right", va="center",
            fontsize=17, fontweight="bold",
            color=group_colors[group],
            transform=ax_risk.get_yaxis_transform(),
        )

        # Count patients at risk at each time point
        at_risk = []
        for t in time_points:
            n = (group_data["OS.time"] >= t).sum()
            at_risk.append(n)

        # Place counts
        for j, (t, n) in enumerate(zip(time_points, at_risk)):
            ax_risk.text(
                t, i, str(n),
                ha="center", va="center",
                fontsize=18, fontweight="bold",
                color="black",
            )

    # Format risk table
    ax_risk.set_ylim(-0.5, len(groups) - 0.5)
    ax_risk.set_xlim(0, xlim_days)
    ax_risk.set_yticks(range(len(groups)))
    ax_risk.set_yticklabels([""] * len(groups))
    ax_risk.set_xlabel("Time (days)", fontsize=24, fontweight="bold", labelpad=23)
    ax_risk.set_xticks(time_points)
    ax_risk.tick_params(labelsize=25, colors="black", width=1.5, length=6)

    # Clean up risk table borders
    ax_risk.spines["top"].set_visible(False)
    ax_risk.spines["right"].set_visible(False)
    ax_risk.spines["left"].set_visible(False)
    ax_risk.tick_params(left=False, bottom=False)

    plt.subplots_adjust(hspace=0.30)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    github_path = output_path.replace(".png", "_github.png")
    plt.savefig(github_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved plot: {output_path}")

def run_cox_fdr(dat, genes):
    """
    Univariate Cox regression per gene with BH-FDR correction.
    Returns a DataFrame sorted by q-value.
    """
    # Check event count (Peduzzi 1995 — need ≥10 events for Cox to be reliable)
    n_events = int(dat["OS"].sum())
    if n_events < 10:
        print(f"  ⚠ WARNING: Only {n_events} events in this cohort (<10).")
        print(f"  Cox regression results may be unreliable.")

    results = []
    for i, gene in enumerate(genes):
        expr_col = f"expr{i+1}"
        cox_data = dat[["OS.time", "OS", expr_col]].rename(columns={expr_col: "expression"})

        cph = CoxPHFitter()
        try:
            cph.fit(cox_data, duration_col="OS.time", event_col="OS")
            summary = cph.summary.loc["expression"]
            results.append({
                "gene": gene,
                "HR": summary["exp(coef)"],
                "HR_lower_95": summary["exp(coef) lower 95%"],
                "HR_upper_95": summary["exp(coef) upper 95%"],
                "p_value": summary["p"],
                "n_events": n_events,
            })
        except Exception as e:
            print(f"  ⚠ Cox fit failed for {gene}: {e}")
            results.append({
                "gene": gene,
                "HR": np.nan,
                "HR_lower_95": np.nan,
                "HR_upper_95": np.nan,
                "p_value": np.nan,
                "n_events": n_events,
            })

    df = pd.DataFrame(results)

    # BH-FDR correction on non-NaN p-values
    valid = df["p_value"].notna()
    if valid.sum() > 0:
        _, qvals, _, _ = multipletests(df.loc[valid, "p_value"], method="fdr_bh")
        df.loc[valid, "q_value"] = qvals
    else:
        df["q_value"] = np.nan

    # Sort by q-value ascending (most significant first)
    df = df.sort_values("q_value", na_position="last").reset_index(drop=True)
    return df

def make_forest_plot(df, title_prefix, output_path):
    """
    Horizontal forest plot of HRs with 95% CI error bars.
    Genes sorted by q-value (most significant at top).
    HR axis on log scale.
    """
    df_plot = df.dropna(subset=["HR"]).copy()
    if len(df_plot) == 0:
        print(f"  No valid Cox results to plot for {title_prefix}.")
        return

    # Sort for display: most significant (smallest q) at top
    df_plot = df_plot.sort_values("q_value", ascending=True).reset_index(drop=True)
    df_plot = df_plot.iloc[::-1].reset_index(drop=True)  # reverse so smallest q is at top of plot

    n = len(df_plot)
    fig_height = max(4, 0.6 * n + 2)
    fig, ax = plt.subplots(figsize=(9, fig_height))

    y_positions = np.arange(n)

    # Error bars (95% CI)
    for y, (_, row) in zip(y_positions, df_plot.iterrows()):
        ax.plot([row["HR_lower_95"], row["HR_upper_95"]], [y, y],
                color="black", linewidth=1.5, zorder=2)
        ax.plot(row["HR"], y, "s", color="#4477AA",
                markersize=10, markeredgecolor="black", zorder=3)

    # Reference line at HR=1
    ax.axvline(1, color="gray", linestyle="--", linewidth=1, zorder=1)

    # Gene labels on y-axis (italic)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(df_plot["gene"], fontsize=14, fontstyle="italic")

    # Log-scale x-axis with explicit tick marks
    from matplotlib.ticker import FixedLocator, FixedFormatter
    ax.set_xscale("log")
    tick_values = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0]
    ax.xaxis.set_major_locator(FixedLocator(tick_values))
    ax.xaxis.set_major_formatter(FixedFormatter([str(t) for t in tick_values]))
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel("Hazard Ratio (95% CI)", fontsize=14, fontweight="bold")
    ax.tick_params(labelsize=12)

    # Annotations placed to the right of the plot box in axes coordinates
    # (x > 1.0 means outside the right edge of the data area)
    header_y = n - 0.4

    # Column x-positions in axes fraction (0 = left edge, 1 = right edge of plot box)
    col_hr_x = 1.04
    col_ci_x = 1.14
    col_p_x = 1.34
    col_q_x = 1.52

    # Column headers above the top row
    ax.text(col_hr_x, header_y, "HR",
            transform=ax.get_yaxis_transform(),
            va="bottom", ha="left",
            fontsize=12, fontweight="bold", family="monospace")
    ax.text(col_ci_x, header_y, "95% CI",
            transform=ax.get_yaxis_transform(),
            va="bottom", ha="left",
            fontsize=12, fontweight="bold", family="monospace")
    ax.text(col_p_x, header_y, "P-value (Wald)",
            transform=ax.get_yaxis_transform(),
            va="bottom", ha="left",
            fontsize=12, fontweight="bold", family="monospace")
    ax.text(col_q_x, header_y, "P_adj (B-H)",
            transform=ax.get_yaxis_transform(),
            va="bottom", ha="left",
            fontsize=12, fontweight="bold", family="monospace")

    # Data rows
    for y, (_, row) in zip(y_positions, df_plot.iterrows()):
        hr_text = f"{row['HR']:.2f}"
        ci_text = f"({row['HR_lower_95']:.2f}–{row['HR_upper_95']:.2f})"
        p_text = f"{row['p_value']:.3g}" if pd.notna(row["p_value"]) else "NA"
        q_text = f"{row['q_value']:.3g}" if pd.notna(row["q_value"]) else "NA"
        ax.text(col_hr_x, y, hr_text,
                transform=ax.get_yaxis_transform(),
                va="center", ha="left",
                fontsize=11, family="monospace")
        ax.text(col_ci_x, y, ci_text,
                transform=ax.get_yaxis_transform(),
                va="center", ha="left",
                fontsize=11, family="monospace")
        ax.text(col_p_x, y, p_text,
                transform=ax.get_yaxis_transform(),
                va="center", ha="left",
                fontsize=11, family="monospace")
        ax.text(col_q_x, y, q_text,
                transform=ax.get_yaxis_transform(),
                va="center", ha="left",
                fontsize=11, family="monospace")

    # Keep x-axis limits focused on the HR data only (not the annotations)
    ax.set_xlim(left=min(df_plot["HR_lower_95"].min() * 0.8, 0.5),
                right=4.5)

    # Title
    fig.text(0.05, 0.95, title_prefix,
             fontsize=18, fontweight="bold",
             verticalalignment="top")

    # Clean spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    github_path = output_path.replace(".png", "_github.png")
    plt.savefig(github_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved forest plot: {output_path}")

# ============================================================
# Main analysis
# ============================================================

def run_2gene_survival(genes):
    print("=== Pan-Cancer Survival ===")
    dat_pan = load_and_merge(EXPRESSION_FILE, SURVIVAL_FILE, genes)
    print(f"  Merged samples: {len(dat_pan)}")
    dat_pan = assign_groups(dat_pan, genes)
    make_km_plot(
        dat_pan, genes,
        title_prefix="Pan-Cancer",
        output_path=os.path.join(RESULTS_DIR, "survival_pancancer.png"),
        xlim_days=config["pancancer_xlim_days"],
        break_time=config["pancancer_break_time"],
    )

    print("\n=== AML Survival ===")
    dat_aml = load_and_merge(
        EXPRESSION_FILE, SURVIVAL_FILE, genes,
        sample_filter=AML_PREFIX,
    )
    print(f"  Merged AML samples: {len(dat_aml)}")
    dat_aml = assign_groups(dat_aml, genes)
    make_km_plot(
        dat_aml, genes,
        title_prefix="Acute Myeloid Leukemia",
        output_path=os.path.join(RESULTS_DIR, "survival_aml.png"),
        xlim_days=config["aml_xlim_days"],
        break_time=config["aml_break_time"],
    )

def run_multigene_cox(genes):
    print(f"=== Pan-Cancer Cox Analysis ({len(genes)} genes) ===")
    dat_pan = load_and_merge(EXPRESSION_FILE, SURVIVAL_FILE, genes)
    print(f"  Merged samples: {len(dat_pan)}")
    df_pan = run_cox_fdr(dat_pan, genes)
    df_pan.to_csv(os.path.join(RESULTS_DIR, "cox_pancancer.csv"), index=False)
    print(f"  Saved table: results/cox_pancancer.csv")
    make_forest_plot(df_pan, "Pan-Cancer",
                     os.path.join(RESULTS_DIR, "cox_forest_pancancer.png"))

    print(f"\n=== AML Cox Analysis ({len(genes)} genes) ===")
    dat_aml = load_and_merge(EXPRESSION_FILE, SURVIVAL_FILE, genes,
                             sample_filter=AML_PREFIX)
    print(f"  Merged AML samples: {len(dat_aml)}")
    df_aml = run_cox_fdr(dat_aml, genes)
    df_aml.to_csv(os.path.join(RESULTS_DIR, "cox_aml.csv"), index=False)
    print(f"  Saved table: results/cox_aml.csv")
    make_forest_plot(df_aml, "Acute Myeloid Leukemia",
                     os.path.join(RESULTS_DIR, "cox_forest_aml.png"))


def main():
    args = parse_args()
    genes = [g.strip() for g in args.genes.split(",")]
    if len(genes) < 2:
        print("Error: provide at least 2 genes.")
        sys.exit(1)

    if not os.path.exists(EXPRESSION_FILE):
        print(f"Error: {EXPRESSION_FILE} not found. Run 01_data_prep.py first.")
        sys.exit(1)
    if not os.path.exists(SURVIVAL_FILE):
        print(f"Error: {SURVIVAL_FILE} not found. Run 01_data_prep.py first.")
        sys.exit(1)

    os.makedirs(RESULTS_DIR, exist_ok=True)

    if len(genes) == 2:
        run_2gene_survival(genes)
    else:
        run_multigene_cox(genes)

    print("\nDone.")


if __name__ == "__main__":
    main()
