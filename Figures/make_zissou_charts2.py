#!/usr/bin/env python3
"""
make_zissou_charts2.py

A second set of figures from cpg_brain_peripheral_correlation_database.csv and
buccal_brain_high_corr.csv, in the Wes Anderson "Zissou1" palette
(https://github.com/karthik/wesanderson), highlighting findings not covered
by the first Zissou set (brain-region focus) or the earlier default-palette
figures (tissue-pair / chromosome / gene-count / genomic-context overviews).

Zissou1 = ["#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"]

Figures:
  1. zissou_fig5_mqtl_corsiv.png     - genetic (mQTL) vs environmental (CoRSIV)
                                       drivers behind buccal-brain correlated CpGs
  2. zissou_fig6_context_by_source.png - genomic context differs between the
                                       blood-brain and buccal-brain datasets
  3. zissou_fig7_island_relation.png  - CpG-island relation of qualifying sites
  4. zissou_fig8_corr_distribution.png - overall correlation-value distribution,
                                       Pearson (blood) vs Spearman (buccal)
"""

import csv
import textwrap
from collections import Counter

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

plt.rcdefaults()
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'figure.dpi': 150,
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.titleweight': 'bold',
    'font.family': 'DejaVu Sans',
})

ZISSOU1 = ["#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"]

MASTER_CSV = "/home/user/workspace/cpg_brain_peripheral_correlation_database.csv"
BUCCAL_CSV = "/home/user/workspace/buccal_brain_high_corr.csv"

SOURCE_NOTE_MASTER = ("Data: Hannon et al. 2015, Epigenetics 10(11):1024-32 "
                       "(https://pmc.ncbi.nlm.nih.gov/articles/PMC4844197/) and Sommerer et al. 2022, "
                       "Clinical Epigenetics 14:118 (https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-022-01357-w). "
                       "Combined database, |r|/|rho| > 0.6 threshold, n=9,127 CpG x tissue-pair rows.")
SOURCE_NOTE_BUCCAL = ("Data: Sommerer et al. 2022, Clinical Epigenetics 14:118 "
                       "(https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-022-01357-w), "
                       "Supplementary Table S6, buccal vs prefrontal cortex, |Spearman rho| > 0.6, n=1,875 CpGs.")

with open(MASTER_CSV, newline="", encoding="utf-8") as f:
    master_rows = list(csv.DictReader(f))
    for r in master_rows:
        r["corr_value"] = float(r["corr_value"])

with open(BUCCAL_CSV, newline="", encoding="utf-8") as f:
    buccal_rows = list(csv.DictReader(f))

print(f"Loaded {len(master_rows):,} master rows, {len(buccal_rows):,} buccal rows")


# ===========================================================================
# FIGURE 5 — mQTL / CoRSIV / IMAGE-CpG status of buccal-brain correlated CpGs
# A key mechanistic question for cross-tissue epigenetic markers: is the
# correlation genetically anchored (mQTL, i.e. driven by underlying DNA
# sequence variation) or does it reflect stable environmentally-influenced
# variation (CoRSIV)? This matters directly for causal interpretation.
# ===========================================================================
n = len(buccal_rows)
categories = [
    ("Blood mQTL", sum(1 for r in buccal_rows if r["blood_mqtl"] == "yes")),
    ("Buccal mQTL", sum(1 for r in buccal_rows if r["buccal_mqtl"] == "yes")),
    ("Brain mQTL", sum(1 for r in buccal_rows if r["brain_mqtl"] == "yes")),
    ("CoRSIV region", sum(1 for r in buccal_rows if r["corsiv"] == "yes")),
    ("Prior IMAGE-CpG\n(Braun et al.)", sum(1 for r in buccal_rows if r["in_braun_IMAGE_CpG"] == "yes")),
]

fig, ax = plt.subplots(figsize=(9.5, 6), layout='constrained')
labels = [c[0] for c in categories]
values = [c[1] for c in categories]
pcts = [v / n * 100 for v in values]
colors = [ZISSOU1[0], ZISSOU1[1], ZISSOU1[2], ZISSOU1[3], ZISSOU1[4]]

bars = ax.barh(labels, pcts, color=colors, edgecolor='white', linewidth=1.2, height=0.62)
ax.invert_yaxis()
for bar, val, pct in zip(bars, values, pcts):
    ax.annotate(f"{val:,} ({pct:.0f}%)", (bar.get_width(), bar.get_y() + bar.get_height() / 2),
                xytext=(8, 0), textcoords='offset points', va='center', fontsize=11, fontweight='bold',
                color='#28251D')

ax.set_xlim(0, 100)
ax.set_xlabel(f"% of {n:,} qualifying buccal-brain CpGs (|Spearman rho| > 0.6)")
ax.set_title("Most cross-tissue correlated CpGs sit at genetically-influenced (mQTL) sites",
             loc='left', fontsize=14, fontweight='bold')
ax.text(0, 1.09, "mQTL and CoRSIV annotation from Sommerer et al. 2022 Supplementary Table S6",
        transform=ax.transAxes, fontsize=10.5, color='#7A7974')
ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f"{int(x)}%"))

note = textwrap.fill(SOURCE_NOTE_BUCCAL, 115)
fig.text(0.01, -0.06, note, ha='left', va='top', fontsize=8, color='#7A7974')

fig.savefig("/home/user/workspace/zissou_fig5_mqtl_corsiv.png", bbox_inches='tight', dpi=150)
plt.close(fig)
print("Saved zissou_fig5_mqtl_corsiv.png")


# ===========================================================================
# FIGURE 6 — Genomic context distribution: does the blood-brain dataset skew
# differently from the buccal-brain dataset in terms of promoter / gene-body
# / intergenic location?
# ===========================================================================
CONTEXT_ORDER = ["Promoter-associated", "Gene body/intragenic", "Unannotated/Intergenic"]
SOURCE_GROUPS = {
    "Blood vs brain\n(Hannon 2015)": [r for r in master_rows if r["corr_type"] == "Pearson r"],
    "Buccal vs brain\n(Sommerer 2022)": [r for r in master_rows if r["corr_type"] == "Spearman rho"],
}

fig, ax = plt.subplots(figsize=(9, 6), layout='constrained')
x = np.arange(len(SOURCE_GROUPS))
width = 0.25
for i, ctx in enumerate(CONTEXT_ORDER):
    pcts = []
    for group_rows in SOURCE_GROUPS.values():
        cnt = sum(1 for r in group_rows if r["context"] == ctx)
        pcts.append(cnt / len(group_rows) * 100)
    offset = (i - 1) * width
    bars = ax.bar(x + offset, pcts, width=width, label=ctx, color=ZISSOU1[i * 2 if i < 2 else 4],
                  edgecolor='white', linewidth=1)
    for bar, pct in zip(bars, pcts):
        ax.annotate(f"{pct:.0f}%", (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                    xytext=(0, 4), textcoords='offset points', ha='center', fontsize=9.5, fontweight='bold',
                    color='#28251D')

ax.set_xticks(x)
ax.set_xticklabels(SOURCE_GROUPS.keys(), fontsize=11)
ax.set_ylabel("% of qualifying CpGs in dataset")
ax.set_ylim(0, 95)
ax.set_title("Buccal-brain correlated CpGs skew more intergenic than blood-brain ones",
             loc='left', fontsize=14, fontweight='bold')
ax.text(0, 1.06, "Genomic context (UCSC RefGene group) of high-correlation CpGs, by source dataset",
        transform=ax.transAxes, fontsize=10.5, color='#7A7974')
legend = ax.legend(loc='upper center', ncol=1, frameon=True, fontsize=9.5)
legend.get_frame().set_facecolor('#F9F8F5')
legend.get_frame().set_edgecolor('#D4D1CA')

note = textwrap.fill(SOURCE_NOTE_MASTER, 130)
fig.text(0.01, -0.06, note, ha='left', va='top', fontsize=8, color='#7A7974')

fig.savefig("/home/user/workspace/zissou_fig6_context_by_source.png", bbox_inches='tight', dpi=150)
plt.close(fig)
print("Saved zissou_fig6_context_by_source.png")


# ===========================================================================
# FIGURE 7 — CpG-island relation of qualifying sites (Island / Shore / Shelf
# / Open sea), across the full combined database.
# ===========================================================================
ISLAND_ORDER = ["Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", ""]
ISLAND_LABELS = {"Island": "CpG Island", "N_Shore": "N Shore", "S_Shore": "S Shore",
                  "N_Shelf": "N Shelf", "S_Shelf": "S Shelf", "": "Open sea /\nunannotated"}

counts = Counter(r["island_relation"] for r in master_rows)
values = [counts.get(k, 0) for k in ISLAND_ORDER]
total = sum(values)
labels = [ISLAND_LABELS[k] for k in ISLAND_ORDER]

fig, ax = plt.subplots(figsize=(9.5, 6), layout='constrained')
colors = [ZISSOU1[4], ZISSOU1[3], ZISSOU1[2], ZISSOU1[1], ZISSOU1[0], "#BAB9B4"]
bars = ax.bar(labels, values, color=colors, edgecolor='white', linewidth=1.2, width=0.62)
for bar, val in zip(bars, values):
    pct = val / total * 100
    ax.annotate(f"{val:,}\n({pct:.0f}%)", (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                xytext=(0, 6), textcoords='offset points', ha='center', va='bottom', fontsize=10, fontweight='bold',
                color='#28251D')

ax.set_ylim(0, max(values) * 1.25)
ax.set_ylabel("Qualifying CpG x tissue-pair rows")
ax.set_title("Most high-correlation CpGs fall outside annotated CpG islands",
             loc='left', fontsize=14, fontweight='bold')
ax.text(0, 1.06, f"CpG-island relation across all {total:,} qualifying rows in the combined database",
        transform=ax.transAxes, fontsize=10.5, color='#7A7974')

note = textwrap.fill(SOURCE_NOTE_MASTER, 130)
fig.text(0.01, -0.06, note, ha='left', va='top', fontsize=8, color='#7A7974')

fig.savefig("/home/user/workspace/zissou_fig7_island_relation.png", bbox_inches='tight', dpi=150)
plt.close(fig)
print("Saved zissou_fig7_island_relation.png")


# ===========================================================================
# FIGURE 8 — Overall correlation-value distribution: Pearson r (blood-brain)
# vs Spearman rho (buccal-brain), highlighting how the two source studies'
# statistics compare above the shared |r|/|rho| > 0.6 threshold.
# ===========================================================================
pearson_vals = [r["corr_value"] for r in master_rows if r["corr_type"] == "Pearson r"]
spearman_vals = [r["corr_value"] for r in master_rows if r["corr_type"] == "Spearman rho"]

fig, ax = plt.subplots(figsize=(9.5, 6), layout='constrained')
bins = np.linspace(0.6, 1.0, 41)
ax.hist(pearson_vals, bins=bins, color=ZISSOU1[0], alpha=0.75, label=f"Pearson r, blood-brain (n={len(pearson_vals):,})",
        edgecolor='white', linewidth=0.4)
ax.hist(spearman_vals, bins=bins, color=ZISSOU1[4], alpha=0.75, label=f"Spearman rho, buccal-brain (n={len(spearman_vals):,})",
        edgecolor='white', linewidth=0.4)

ax.set_ylim(0, 460)
ax.axvline(np.median(pearson_vals), color=ZISSOU1[0], linestyle='--', linewidth=1.5)
ax.axvline(np.median(spearman_vals), color=ZISSOU1[4], linestyle='--', linewidth=1.5)
ax.annotate(f"median {np.median(pearson_vals):.3f}", (np.median(pearson_vals), 250),
            xytext=(-8, 0), textcoords='offset points', ha='right', va='center', fontsize=9.5,
            fontweight='bold', color=ZISSOU1[0])
ax.annotate(f"median {np.median(spearman_vals):.3f}", (np.median(spearman_vals), 250),
            xytext=(8, 0), textcoords='offset points', ha='left', va='center', fontsize=9.5,
            fontweight='bold', color='#A13544')

ax.set_xlabel("Absolute correlation value (|r| or |rho|)")
ax.set_ylabel("Number of CpG x tissue-pair rows")
ax.set_title("Blood-brain correlations skew higher and tighter than buccal-brain ones",
             loc='left', fontsize=14, fontweight='bold')
ax.text(0, 1.06, "Distribution of qualifying correlation values by source dataset and statistic type",
        transform=ax.transAxes, fontsize=10.5, color='#7A7974')
legend = ax.legend(loc='upper left', frameon=True, fontsize=10)
legend.get_frame().set_facecolor('#F9F8F5')
legend.get_frame().set_edgecolor('#D4D1CA')

note = textwrap.fill(SOURCE_NOTE_MASTER, 130)
fig.text(0.01, -0.06, note, ha='left', va='top', fontsize=8, color='#7A7974')

fig.savefig("/home/user/workspace/zissou_fig8_corr_distribution.png", bbox_inches='tight', dpi=150)
plt.close(fig)
print("Saved zissou_fig8_corr_distribution.png")

print("\nAll 4 additional figures generated.")
