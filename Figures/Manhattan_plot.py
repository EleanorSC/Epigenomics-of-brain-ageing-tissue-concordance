import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import textwrap

plt.rcdefaults()
plt.rcParams.update({
    'figure.dpi': 150,
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.titleweight': 'bold',
    'font.family': 'sans-serif',
    'axes.grid': False,
})

TEAL = '#20808D'
DARK_TEAL = '#1B474D'
RUST = '#A84B2F'
TEXT = '#28251D'
MUTED = '#7A7974'
BAND_LIGHT = '#F2F1EC'

np.random.seed(42)

df = pd.read_excel('/home/user/workspace/CpG_Brain_Peripheral_Correlation_Database.xlsx', sheet_name='CpG_Database')

# Collapse to one row per unique CpG site: take the strongest correlation observed across tissue pairs
site = (
    df.sort_values('Correlation Value', ascending=False)
      .groupby('CpG ID', as_index=False)
      .first()[['CpG ID', 'Chromosome', 'Position (hg19)', 'Gene', 'Correlation Value']]
)
site['Gene_label'] = site['Gene'].str.split(';').str[-1].str.strip()

chrom_order = [f'chr{i}' for i in range(1, 23)]
site = site[site['Chromosome'].isin(chrom_order)].copy()
site['chrom_idx'] = site['Chromosome'].map({c: i for i, c in enumerate(chrom_order)})

# Equal-width chromosome slots with tight geom_jitter-style horizontal scatter within each
# (mirrors the dense, clustered "towers" typical of real EWAS Manhattan plots, where hundreds of
# thousands of probes per chromosome pack points into a narrow band rather than spreading them
# thinly across true base-pair coordinates)
SLOT_WIDTH = 1.0
JITTER_WIDTH = 0.34  # tight grouping around each chromosome's center
site['genome_x'] = site['chrom_idx'] * SLOT_WIDTH + SLOT_WIDTH / 2 + np.random.uniform(
    -JITTER_WIDTH, JITTER_WIDTH, size=len(site)
)

cum = len(chrom_order) * SLOT_WIDTH

fig, ax = plt.subplots(figsize=(14, 5.5), layout='constrained')
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

# Alternating light background bands per chromosome (classic Manhattan-plot convention)
for i, c in enumerate(chrom_order):
    if i % 2 == 1:
        ax.axvspan(i * SLOT_WIDTH, (i + 1) * SLOT_WIDTH, color=BAND_LIGHT, zorder=0, linewidth=0)

colors_alt = [TEAL, DARK_TEAL]
tick_positions = []
tick_labels = []
for i, c in enumerate(chrom_order):
    sub = site[site['Chromosome'] == c]
    ax.scatter(sub['genome_x'], sub['Correlation Value'], s=15, color=colors_alt[i % 2],
               alpha=0.8, linewidths=0, zorder=2)
    tick_positions.append(i * SLOT_WIDTH + SLOT_WIDTH / 2)
    tick_labels.append(c.replace('chr', ''))

# Genome-wide correlation threshold line (classic red dashed significance-line convention)
ax.axhline(0.6, color=RUST, linestyle='--', linewidth=1.4, zorder=3)
ax.text(cum * 0.005, 0.615, 'Inclusion threshold: r/\u03c1 = 0.6', ha='left', va='bottom',
        fontsize=9.5, color=RUST, fontweight='bold', zorder=4,
        bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='none', alpha=0.88))

# Highlight and label the top 10 strongest single-site correlations, in red
top_sites = site.sort_values('Correlation Value', ascending=False).head(10).sort_values('genome_x').reset_index(drop=True)
for _, row in top_sites.iterrows():
    ax.scatter(row['genome_x'], row['Correlation Value'], s=58, color=RUST, zorder=5,
               edgecolors='white', linewidths=0.9)

# Stagger label y-offsets in a repeating pattern so 10 closely spaced labels stay legible
offset_pattern = [16, 34, 16, 34, 52, 16, 34, 16, 34, 16]
for k, (_, row) in enumerate(top_sites.iterrows()):
    dy = offset_pattern[k % len(offset_pattern)]
    ax.annotate(
        f"{row['Gene_label']}",
        xy=(row['genome_x'], row['Correlation Value']),
        xytext=(0, dy), textcoords='offset points',
        fontsize=8.7, color=RUST, ha='center', fontweight='bold', zorder=6,
        arrowprops=dict(arrowstyle='-', color=RUST, lw=0.5, alpha=0.6),
    )

ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels, fontsize=9.5)
ax.set_xlim(0, cum)
ax.set_ylim(0.55, 1.14)
ax.set_yticks([0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_xlabel('Chromosome', fontsize=11)
ax.set_ylabel('Strongest cross-tissue correlation (r or \u03c1)', fontsize=11)
ax.set_title(
    "970 individual CpG sites show cross-tissue methylation concordance genome-wide\n"
    "Manhattan plot of strongest blood/buccal-brain correlation per CpG site, by chromosome",
    loc='left'
)
ax.yaxis.grid(True, alpha=0.18, linewidth=0.7)
ax.xaxis.grid(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color(MUTED)
ax.spines['bottom'].set_color(MUTED)
ax.tick_params(colors=TEXT)

note = textwrap.fill(
    "Source: curated cross-tissue CpG correlation database (970 unique CpG sites underlying 2,937 CpG-tissue-pair records). "
    "Each point is one CpG site, jittered within its chromosome band and plotted at its strongest observed correlation "
    "across all tested tissue pairs; horizontal position is not to genomic scale. Top 10 sites by correlation strength "
    "labeled in red with gene symbol.",
    120
)
fig.text(0.01, -0.06, note, ha='left', va='top', fontsize=9, color=TEXT)

plt.savefig('/home/user/workspace/fig6_cpg_sites_manhattan.png', bbox_inches='tight', dpi=150, facecolor='white')
plt.close()
print("done")
