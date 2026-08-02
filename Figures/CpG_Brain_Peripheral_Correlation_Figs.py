import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import textwrap

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'figure.dpi': 150,
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.titleweight': 'bold',
    'font.family': 'sans-serif',
})

TEAL = '#20808D'
RUST = '#A84B2F'
DARK_TEAL = '#1B474D'
LIGHT_CYAN = '#BCE2E7'
MAUVE = '#944454'
GOLD = '#FFC553'
OLIVE = '#848456'
BROWN = '#6E522B'
PALETTE = [TEAL, RUST, DARK_TEAL, MAUVE, GOLD, OLIVE, BROWN, LIGHT_CYAN]

TEXT = '#28251D'
MUTED = '#7A7974'
BORDER = '#D4D1CA'

df = pd.read_excel('/home/user/workspace/CpG_Brain_Peripheral_Correlation_Database.xlsx', sheet_name='CpG_Database')
gs = pd.read_excel('/home/user/workspace/CpG_Brain_Peripheral_Correlation_Database.xlsx', sheet_name='Gene_Summary')

print(df.shape, gs.shape)
print(df.columns.tolist())
print(df['Tissue Pair'].value_counts())
print(df['Genomic Context'].value_counts())
print(df['Correlation Value'].describe())

# ---------- Figure 1: Distribution of correlation values ----------
fig, ax = plt.subplots(figsize=(9, 5.5), layout='constrained')
vals = df['Correlation Value'].dropna()
ax.hist(vals, bins=40, color=TEAL, edgecolor='white', linewidth=0.4)
median_v = vals.median()
ax.axvline(median_v, color=RUST, linestyle='--', linewidth=1.6)
ax.text(median_v + 0.01, ax.get_ylim()[1]*0.95, f'Median r = {median_v:.2f}', color=RUST, fontsize=10, va='top')
ax.set_title(f"Cross-tissue CpGs are concentrated at high correlation strength\nDistribution of {len(vals):,} qualifying blood/buccal-brain correlations (r or ρ > 0.6)", loc='left')
ax.set_xlabel('Correlation value (Pearson r or Spearman ρ)')
ax.set_ylabel('Number of CpG-tissue-pair records')
ax.spines[['top', 'right']].set_visible(False)
note = textwrap.fill("Source: curated CpG-level database compiled from Hannon et al. 2015 (Blood-Brain DNA Methylation Comparison Tool) and Sommerer et al. 2022 (buccal-brain correlation map); n=9,127 pre-filter master rows, 2,937 gene-annotated rows shown here.", 105)
fig.text(0.01, -0.06, note, ha='left', va='top', fontsize=8, color=MUTED)
plt.savefig('/home/user/workspace/fig1_correlation_distribution.png', bbox_inches='tight', dpi=150)
plt.close()

# ---------- Figure 2: Correlation strength by tissue pair (box plot) ----------
fig, ax = plt.subplots(figsize=(10, 6), layout='constrained')
tissue_order = df.groupby('Tissue Pair')['Correlation Value'].median().sort_values(ascending=False).index.tolist()
data_by_tissue = [df.loc[df['Tissue Pair'] == t, 'Correlation Value'].dropna().values for t in tissue_order]
bp = ax.boxplot(data_by_tissue, vert=False, patch_artist=True, showfliers=False, widths=0.6)
for i, box in enumerate(bp['boxes']):
    box.set_facecolor(PALETTE[i % len(PALETTE)])
    box.set_alpha(0.85)
    box.set_edgecolor(DARK_TEAL)
for median in bp['medians']:
    median.set_color(TEXT)
    median.set_linewidth(1.5)
ax.set_yticks(range(1, len(tissue_order) + 1))
ax.set_yticklabels(tissue_order, fontsize=10)
ax.set_xlabel('Correlation value (r or ρ)')
ax.set_title("Blood-brain pairs show the strongest and most consistent concordance\nCorrelation strength distribution by tissue-pair comparison", loc='left')
ax.spines[['top', 'right']].set_visible(False)
n_counts = df['Tissue Pair'].value_counts()
for i, t in enumerate(tissue_order):
    ax.text(1.005, i + 1, f'n={n_counts[t]:,}', transform=ax.get_yaxis_transform(), fontsize=8.5, color=MUTED, va='center')
note = textwrap.fill("Source: curated cross-tissue CpG correlation database (2,937 gene-annotated rows, r/ρ > 0.6 threshold). EC = entorhinal cortex; STG = superior temporal gyrus; PFC = prefrontal cortex; CER = cerebellum.", 110)
fig.text(0.01, -0.06, note, ha='left', va='top', fontsize=8, color=MUTED)
plt.savefig('/home/user/workspace/fig2_correlation_by_tissue_pair.png', bbox_inches='tight', dpi=150)
plt.close()

# ---------- Figure 3: Qualifying CpGs by chromosome ----------
chrom_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
chrom_counts = df['Chromosome'].value_counts().reindex(chrom_order).dropna()
fig, ax = plt.subplots(figsize=(11, 5.5), layout='constrained')
colors = [TEAL if c not in ('chrX', 'chrY') else RUST for c in chrom_counts.index]
ax.bar(chrom_counts.index, chrom_counts.values, color=colors, edgecolor='white', linewidth=0.4)
ax.set_title(f"Qualifying CpGs are distributed across all autosomes, with chr1 and chr6 leading\nCross-tissue correlated CpGs (r/ρ > 0.6) by chromosome, n={int(chrom_counts.sum()):,}", loc='left')
ax.set_ylabel('Number of CpG-tissue-pair records')
ax.set_xlabel('Chromosome')
ax.tick_params(axis='x', rotation=45)
ax.spines[['top', 'right']].set_visible(False)
note = textwrap.fill("Source: curated cross-tissue CpG correlation database. No qualifying CpGs (r/ρ > 0.6) were annotated to sex chromosomes in this dataset.", 110)
fig.text(0.01, -0.05, note, ha='left', va='top', fontsize=8, color=MUTED)
plt.savefig('/home/user/workspace/fig3_cpgs_by_chromosome.png', bbox_inches='tight', dpi=150)
plt.close()

# ---------- Figure 4: Genomic context breakdown ----------
ctx_counts = df['Genomic Context'].value_counts()
ctx_counts = ctx_counts[ctx_counts.index != 'n.a.'].sort_values(ascending=True)
fig, ax = plt.subplots(figsize=(9.5, 5), layout='constrained')
bars = ax.barh(ctx_counts.index, ctx_counts.values, color=TEAL, edgecolor='white', linewidth=0.4)
total = ctx_counts.sum()
for bar, val in zip(bars, ctx_counts.values):
    pct = val / total * 100
    ax.text(val + total*0.01, bar.get_y() + bar.get_height()/2, f'{val:,} ({pct:.0f}%)', va='center', fontsize=9.5, color=TEXT)
ax.set_title("Most cross-tissue correlated CpGs fall in gene bodies rather than promoters\nGenomic context of qualifying CpGs (annotated rows only)", loc='left')
ax.set_xlabel('Number of CpG-tissue-pair records')
ax.spines[['top', 'right']].set_visible(False)
ax.set_xlim(0, max(ctx_counts.values) * 1.22)
note = textwrap.fill("Source: curated cross-tissue CpG correlation database; genomic context derived from Illumina array UCSC RefGene annotation. Rows without a resolvable context (n.a.) excluded.", 110)
fig.text(0.01, -0.06, note, ha='left', va='top', fontsize=8, color=MUTED)
plt.savefig('/home/user/workspace/fig4_genomic_context.png', bbox_inches='tight', dpi=150)
plt.close()

# ---------- Figure 5: Top genes by number of qualifying CpGs ----------
top_genes = gs.sort_values('# Qualifying CpGs (r>0.6)', ascending=False).head(15).iloc[::-1]
fig, ax = plt.subplots(figsize=(9.5, 7), layout='constrained')
bars = ax.barh(top_genes['Gene'], top_genes['# Qualifying CpGs (r>0.6)'], color=TEAL, edgecolor='white', linewidth=0.4)
for bar, gene, maxcorr in zip(bars, top_genes['Gene'], top_genes['Max |Correlation|']):
    ax.text(bar.get_width() + 0.15, bar.get_y() + bar.get_height()/2, f"max |r|={maxcorr:.2f}", va='center', fontsize=8.5, color=MUTED)
ax.set_title("A small set of genes carry disproportionately strong cross-tissue signal\nTop 15 genes by number of qualifying CpGs (r/ρ > 0.6)", loc='left')
ax.set_xlabel('Number of qualifying CpGs')
ax.spines[['top', 'right']].set_visible(False)
ax.set_xlim(0, top_genes['# Qualifying CpGs (r>0.6)'].max() * 1.35)
note = textwrap.fill("Source: curated cross-tissue CpG correlation database, Gene_Summary sheet (815 genes total). Max |r| = strongest single correlation observed for that gene across all tissue-pair comparisons.", 110)
fig.text(0.01, -0.05, note, ha='left', va='top', fontsize=8, color=MUTED)
plt.savefig('/home/user/workspace/fig5_top_genes.png', bbox_inches='tight', dpi=150)
plt.close()

print("done")
