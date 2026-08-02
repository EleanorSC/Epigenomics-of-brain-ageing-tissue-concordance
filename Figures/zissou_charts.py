#!/usr/bin/env python3
"""
make_zissou_charts.py

Brain-region-focused figures from cpg_brain_peripheral_correlation_database.csv,
styled with the Wes Anderson "Zissou1" color palette
(https://github.com/karthik/wesanderson).

Zissou1 = ["#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"]

Figures:
  1. zissou_fig1_corr_by_region.png   - distribution of Pearson r by brain region
  2. zissou_fig2_region_counts.png    - qualifying CpG counts per brain region
  3. zissou_fig3_manhattan_regions.png- genome-wide CpG positions colored by region
  4. zissou_fig4_region_specificity.png - how many regions each CpG correlates in
"""

import csv
import textwrap
from collections import Counter, defaultdict

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

DATA_PATH = "/home/user/workspace/cpg_brain_peripheral_correlation_database.csv"

REGION_ORDER = ["Blood vs EC (brain)", "Blood vs STG (brain)",
                "Blood vs PFC (brain)", "Blood vs CER (brain)"]
REGION_LABELS = {"Blood vs EC (brain)": "EC", "Blood vs STG (brain)": "STG",
                  "Blood vs PFC (brain)": "PFC", "Blood vs CER (brain)": "CER"}
REGION_COLOR = dict(zip(REGION_ORDER, ZISSOU1[:4]))  # blue -> yellow -> orange -> red across 4 regions

SOURCE_NOTE = ("Data: Hannon et al. 2015, Epigenetics 10(11):1024-32 "
               "(https://pmc.ncbi.nlm.nih.gov/articles/PMC4844197/) — blood vs brain (EC/STG/PFC/CER), "
               "Pearson r; |r| > 0.6 threshold.")

# ---------------------------------------------------------------------------
# Load data (blood-brain rows only — buccal excluded since it has no region split)
# ---------------------------------------------------------------------------
rows = []
with open(DATA_PATH, newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        if r["tissue_pair"] in REGION_ORDER:
            r["corr_value"] = float(r["corr_value"])
            r["pos"] = int(r["pos"]) if r["pos"] else None
            r["chr"] = r["chr"] if r["chr"] else None
            rows.append(r)

print(f"Loaded {len(rows):,} blood-brain rows across 4 regions")


# ===========================================================================
# FIGURE 1 — Violin/box distribution of correlation strength by brain region
# ===========================================================================
fig, ax = plt.subplots(figsize=(9, 6), layout='constrained')

data_by_region = [[r["corr_value"] for r in rows if r["tissue_pair"] == reg] for reg in REGION_ORDER]
positions = range(len(REGION_ORDER))

parts = ax.violinplot(data_by_region, positions=positions, showmedians=True, widths=0.75)
for i, body in enumerate(parts['bodies']):
    body.set_facecolor(ZISSOU1[i])
    body.set_edgecolor(ZISSOU1[i])
    body.set_alpha(0.75)
for key in ['cbars', 'cmins', 'cmaxes', 'cmedians']:
    parts[key].set_color('#28251D')
    parts[key].set_linewidth(1.2)

# overlay median labels
for i, vals in enumerate(data_by_region):
    med = np.median(vals)
    ax.annotate(f"{med:.3f}", (i, med), xytext=(12, 0), textcoords='offset points',
                fontsize=10, fontweight='bold', va='center', color='#28251D')

ax.set_xticks(list(positions))
ax.set_xticklabels([REGION_LABELS[r] for r in REGION_ORDER], fontsize=12)
ax.set_ylabel("Pearson r (blood vs brain region)")
ax.set_ylim(0.55, 1.02)
ax.set_title("Cortical regions show near-identical blood-methylation correlation strength",
             loc='left', fontsize=14, fontweight='bold')
ax.text(0, 1.06, "Distribution of Pearson r values for 1,813 CpGs per region, Hannon et al. 2015 blood-brain dataset",
        transform=ax.transAxes, fontsize=10.5, color='#7A7974')

note = textwrap.fill(SOURCE_NOTE, 110)
fig.text(0.01, -0.05, note, ha='left', va='top', fontsize=8, color='#7A7974')

fig.savefig("/home/user/workspace/zissou_fig1_corr_by_region.png", bbox_inches='tight', dpi=150)
plt.close(fig)
print("Saved zissou_fig1_corr_by_region.png")


# ===========================================================================
# FIGURE 2 — Bar chart: qualifying CpG counts per brain region (all equal at
# the raw count, so instead show count of CpGs uniquely vs shared across all
# regions per bar, stacked) -- simpler: show mean |r| and count as dual bars
# ===========================================================================
fig, ax = plt.subplots(figsize=(9, 5.5), layout='constrained')

counts = [len(v) for v in data_by_region]
means = [np.mean(v) for v in data_by_region]

bars = ax.bar([REGION_LABELS[r] for r in REGION_ORDER], counts,
              color=[REGION_COLOR[r] for r in REGION_ORDER], width=0.6, edgecolor='white', linewidth=1.5)

for bar, mean_r in zip(bars, means):
    h = bar.get_height()
    ax.annotate(f"{h:,}\nmean r={mean_r:.3f}", (bar.get_x() + bar.get_width() / 2, h),
                xytext=(0, 6), textcoords='offset points', ha='center', va='bottom',
                fontsize=10.5, fontweight='bold', color='#28251D')

ax.set_ylim(0, max(counts) * 1.22)
ax.set_ylabel("Qualifying CpGs (|Pearson r| > 0.6)")
ax.set_title("All four brain regions carry the same 1,813 high-correlation CpGs",
             loc='left', fontsize=14, fontweight='bold')
ax.text(0, 1.06, "Hannon et al.'s 'highly variable' probe set qualifies uniformly across EC, STG, PFC, and CER",
        transform=ax.transAxes, fontsize=10.5, color='#7A7974')
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, p: f"{int(x):,}"))

note = textwrap.fill(SOURCE_NOTE, 110)
fig.text(0.01, -0.05, note, ha='left', va='top', fontsize=8, color='#7A7974')

fig.savefig("/home/user/workspace/zissou_fig2_region_counts.png", bbox_inches='tight', dpi=150)
plt.close(fig)
print("Saved zissou_fig2_region_counts.png")


# ===========================================================================
# FIGURE 3 — Manhattan-style genome-wide plot of CpG sites, colored by brain
# region, y-axis = correlation strength
# ===========================================================================
CHR_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y"]
chr_sizes_hg19 = {  # approximate hg19 chromosome lengths (bp), for cumulative x-axis
    "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260,
    "6": 171115067, "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747,
    "11": 135006516, "12": 133851895, "13": 115169878, "14": 107349540, "15": 102531392,
    "16": 90354753, "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
    "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
}
present_chrs = [c for c in CHR_ORDER if c in chr_sizes_hg19]
cum_offset = {}
offset = 0
for c in present_chrs:
    cum_offset[c] = offset
    offset += chr_sizes_hg19[c]

fig, ax = plt.subplots(figsize=(14, 6), layout='constrained')

for reg in REGION_ORDER:
    xs, ys = [], []
    for r in rows:
        if r["tissue_pair"] != reg or r["chr"] is None or r["pos"] is None:
            continue
        if r["chr"] not in cum_offset:
            continue
        xs.append(cum_offset[r["chr"]] + r["pos"])
        ys.append(r["corr_value"])
    ax.scatter(xs, ys, s=10, alpha=0.55, color=REGION_COLOR[reg], label=REGION_LABELS[reg], linewidths=0)

# chromosome boundary lines + labels
tick_pos, tick_labels = [], []
for c in present_chrs:
    start = cum_offset[c]
    ax.axvline(start, color='#D4D1CA', linewidth=0.5, zorder=0)
    tick_pos.append(start + chr_sizes_hg19[c] / 2)
    tick_labels.append(c)

ax.set_xticks(tick_pos)
ax.set_xticklabels(tick_labels, fontsize=8, rotation=90)
ax.set_xlim(0, offset)
ax.set_ylim(0.55, 1.02)
ax.set_ylabel("Pearson r (blood vs brain)")
ax.set_xlabel("Chromosome (hg19 coordinates)")
ax.set_title("High-correlation CpGs are distributed across the genome, not clustered by region",
             loc='left', fontsize=14, fontweight='bold')
ax.text(0, 1.05, "7,252 blood-brain CpG x region pairs with |Pearson r| > 0.6, colored by brain region compared",
        transform=ax.transAxes, fontsize=10.5, color='#7A7974')
legend = ax.legend(loc='lower right', frameon=True, fontsize=10, title="Brain region")
legend.get_frame().set_facecolor('#F9F8F5')
legend.get_frame().set_edgecolor('#D4D1CA')

note = textwrap.fill(SOURCE_NOTE, 150)
fig.text(0.01, -0.09, note, ha='left', va='top', fontsize=8, color='#7A7974')

fig.savefig("/home/user/workspace/zissou_fig3_manhattan_regions.png", bbox_inches='tight', dpi=150)
plt.close(fig)
print("Saved zissou_fig3_manhattan_regions.png")


# ===========================================================================
# FIGURE 4 — Pairwise consistency of correlation strength across brain
# regions: for each CpG present in two regions, plot r(region A) vs
# r(region B). Near-perfect diagonal clustering shows a CpG's blood-brain
# correlation strength is a property of the SITE, not the specific brain
# region it is compared against.
# ===========================================================================
cpg_region_r = defaultdict(dict)
for r in rows:
    cpg_region_r[r["cpg"]][r["tissue_pair"]] = r["corr_value"]

REGION_PAIRS = [
    ("Blood vs EC (brain)", "Blood vs STG (brain)"),
    ("Blood vs EC (brain)", "Blood vs PFC (brain)"),
    ("Blood vs EC (brain)", "Blood vs CER (brain)"),
    ("Blood vs STG (brain)", "Blood vs PFC (brain)"),
    ("Blood vs STG (brain)", "Blood vs CER (brain)"),
    ("Blood vs PFC (brain)", "Blood vs CER (brain)"),
]

fig, axes = plt.subplots(2, 3, figsize=(13, 8.5), layout='constrained')
for ax, (reg_a, reg_b), color in zip(axes.flat, REGION_PAIRS, ZISSOU1 + [ZISSOU1[0]]):
    xs, ys = [], []
    for cpg, region_map in cpg_region_r.items():
        if reg_a in region_map and reg_b in region_map:
            xs.append(region_map[reg_a])
            ys.append(region_map[reg_b])
    xs, ys = np.array(xs), np.array(ys)
    ax.scatter(xs, ys, s=9, alpha=0.45, color=color, linewidths=0)
    ax.plot([0.55, 1.0], [0.55, 1.0], color='#28251D', linewidth=1, linestyle='--', alpha=0.6)
    meta_r = np.corrcoef(xs, ys)[0, 1] if len(xs) > 1 else float('nan')
    ax.set_xlim(0.55, 1.02)
    ax.set_ylim(0.55, 1.02)
    ax.set_xlabel(f"r, {REGION_LABELS[reg_a]}", fontsize=10)
    ax.set_ylabel(f"r, {REGION_LABELS[reg_b]}", fontsize=10)
    ax.set_title(f"{REGION_LABELS[reg_a]} vs {REGION_LABELS[reg_b]}  (n={len(xs):,}, meta-r={meta_r:.3f})",
                 fontsize=11, fontweight='bold')
    ax.tick_params(labelsize=9)

fig.suptitle("A CpG's blood-brain correlation strength is a property of the site, not the region\nPairwise comparison of Pearson r across the 4 brain regions, for CpGs shared between each region pair",
             fontsize=14.5, fontweight='bold', x=0.01, ha='left')

note = textwrap.fill(SOURCE_NOTE, 150)
fig.text(0.01, -0.03, note, ha='left', va='top', fontsize=8, color='#7A7974')

fig.savefig("/home/user/workspace/zissou_fig4_region_specificity.png", bbox_inches='tight', dpi=150)
plt.close(fig)
print("Saved zissou_fig4_region_specificity.png")

print("\nAll 4 figures generated.")
