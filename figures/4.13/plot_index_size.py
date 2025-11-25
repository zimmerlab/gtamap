import pandas as pd
from scipy import stats
import seaborn as sns
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import os
import sys

# sns.set_theme(style="ticks", context="paper", palette="Set2")
# sns.set_theme(style="ticks", context="paper", palette="colorblind")

sns.set_theme(style="ticks", context="paper", palette="husl")

tool_order_a = ["GTAMap", "Minimap2", "HISAT2", "STAR", "BWA", "Bowtie2"]
base_colors = sns.color_palette("husl", n_colors=6)


tool_colors = dict(zip(tool_order_a, base_colors))
tool_colors["GTAMap"] = sns.color_palette("Reds", 6)[4]  # bright red tone

sns.set_context(
    "paper",
    rc={
        "font.size": 15,  # base font size
        "axes.titlesize": 17,  # title
        "axes.labelsize": 15,  # axis labels
        "xtick.labelsize": 14,  # x tick labels
        "ytick.labelsize": 14,  # y tick labels
        "legend.fontsize": 12,  # legend
        "lines.linewidth": 1.5,  # line width
        "lines.markersize": 6,  # marker size
    },
)

plt.rcParams["figure.figsize"] = (6.5, 4.5)  # width x height in inches

# --- Colors ---
fill_color = "#CFD5E8"
edge_color = "black"


# --- Input files ---
tsv_file = sys.argv[1]  # Your mapper file size TSV
genes_meta_file = sys.argv[2]  # File with gene lengths
out_dir = sys.argv[3] if len(sys.argv) > 3 else "plots"
os.makedirs(out_dir, exist_ok=True)

# --- Read and process file size data ---
# Assuming your TSV has columns: mapper, gene, size_MB
df_sizes = pd.read_csv(tsv_file, sep="\t")


# Clean the file_size_MB column (remove 'K', 'M' and convert to numeric)
def clean_size(size_str):
    """Convert du -sh output to MB (handles European decimal comma)"""
    if not isinstance(size_str, str):
        return size_str

    size_str = size_str.strip()

    # Replace comma with period for decimal point
    size_str = size_str.replace(",", ".")

    # Extract number and unit
    import re

    match = re.match(r"([\d.]+)([KMGT]?)", size_str)
    if not match:
        return float(size_str)

    value = float(match.group(1))
    unit = match.group(2)

    # Convert to MB
    conversions = {
        "": 1 / (1024 * 1024),  # bytes to MB
        "K": 1 / 1024,  # KB to MB
        "M": 1,  # MB to MB
        "G": 1024,  # GB to MB
        "T": 1024 * 1024,  # TB to MB
    }

    return value * conversions.get(unit, 1)


df_sizes["size_MB"] = df_sizes["size_MB"].apply(clean_size)
df_sizes["size_KB"] = df_sizes["size_MB"] * 1024
df_sizes["size_GB"] = df_sizes["size_MB"] / 1024

df_gb = df_sizes.groupby(["mapper", "gene"], as_index=False)["size_GB"].sum()
df_kb = df_sizes.groupby(["mapper", "gene"], as_index=False)["size_KB"].sum()

df_sizes = df_sizes.groupby(["mapper", "gene"], as_index=False)["size_MB"].sum()


print(df_sizes)

# --- Read gene metadata ---
# Assuming genes_meta has columns: gene_id, gene_length
df_genes = pd.read_csv(
    genes_meta_file,
    sep="\t",
    names=["gene_id", "gene_start", "gene_length", "strand", "gene_name", "chr"],
)

# Merge with file size data
df_merged = pd.merge(df_sizes, df_genes, left_on="gene", right_on="gene_id", how="left")
df_merged = df_merged.dropna(subset=["gene_length"])


# For each mapper, get the index of the row with the minimum size
idx = df_kb.groupby("mapper")["size_KB"].idxmin()

# Use those indices to get the full rows, including gene_name
min_sizes_with_gene = df_sizes.loc[idx, ["mapper", "size_MB", "gene"]].reset_index(
    drop=True
)

print(min_sizes_with_gene)

# For each mapper, get the index of the row with the minimum size
idx = df_kb.groupby("mapper")["size_KB"].idxmax()

# Use those indices to get the full rows, including gene_name
max_sizes_with_gene = df_sizes.loc[idx, ["mapper", "size_MB", "gene"]].reset_index(
    drop=True
)
print(max_sizes_with_gene)


# --- Your existing plotting functions ---
def get_whiskers(series):
    """Compute Tukey-style whiskers (1.5*IQR)."""
    q1 = series.quantile(0.25)
    q3 = series.quantile(0.75)
    iqr = q3 - q1
    low = q1 - 1.5 * iqr
    high = q3 + 1.5 * iqr
    return low, high


def boxplot_metric_zoom(
    data, x, y, title, ylabel, filename, out_dir, figsize=(6.5, 4.5)
):
    plt.figure(figsize=figsize)

    ax = sns.boxplot(
        data=data,
        x=x,
        y=y,
        showfliers=False,  # Hide outliers
        width=0.6,
        boxprops=dict(edgecolor=edge_color, linewidth=1.5),
        whiskerprops=dict(color=edge_color, linewidth=1.5),
        capprops=dict(color=edge_color, linewidth=1.5),
        medianprops=dict(color=edge_color, linewidth=1.5),
        hue="mapper",
        palette=tool_colors,
        hue_order=tool_order_a,
        order=tool_order_a,
    )

    # Compute whiskers for zoom
    whisker_lows, whisker_highs = [], []
    for _, group in data.groupby(x):
        low, high = get_whiskers(group[y].dropna())
        whisker_lows.append(low)
        whisker_highs.append(high)
    global_low = min(whisker_lows)
    global_high = max(whisker_highs)
    padding = 0.002 * (global_high - global_low)
    y_start = min(0, global_low - padding)
    ax.set_ylim(0, global_high + padding)

    ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    ax.xaxis.grid(False)  # optional: turn off vertical grid lines

    # Styling
    ax.set_title(title, weight="bold", pad=15)
    plt.xlabel("Mapper")
    ax.set_ylabel(ylabel)
    sns.despine()
    plt.tight_layout()

    # Save plot
    plt.savefig(os.path.join(out_dir, filename), dpi=300)
    plt.close()


################################### COMMENTED OUT FOR SPEED #####################################
#
# --- Create boxplots ---
# Boxplot with y-axis in MB
boxplot_metric_zoom(
    data=df_sizes,
    x="mapper",
    y="size_MB",
    title="Index File Sizes by Mapper",
    ylabel="File Size (MB)",
    filename="index_sizes_MB_box.png",
    out_dir=out_dir,
)
#
# plt.figure(figsize=(9, 5))
# ax = sns.lineplot(
#     data=df_merged,
#     x="gene_length",
#     y="size_MB",
#     hue="mapper",
#     hue_order=tool_order_a,
#     linewidth=1.5,
#     marker="o",  # add points
#     markersize=5,
# )

#
#
# # --- Create scatter plot: gene length vs file size ---
# plt.figure(figsize=(9, 5))
# ax = sns.scatterplot(
#     data=df_merged,
#     x="gene_length",  # from genes_meta
#     y="size_MB",
#     hue="mapper",
#     s=30,
#     hue_order=tool_order_a,
#     edgecolor="black",  # Add black edges
#     linewidths=0.5,
# )
#
# plt.title("Index File Size vs Gene Length", weight="bold", pad=15)
# plt.xlabel("Gene Length [log(length)]")
# plt.ylabel("Size [log(MB)]")
# # plt.legend(
# #     title="Mapper",
# #     fontsize=12,
# #     title_fontsize=14,
# #     loc="center left",  # anchor legend to the left side of the box
# #     bbox_to_anchor=(1.02, 0.5),
# # )
# plt.xscale("log")
# plt.yscale("log")
# sns.despine()
# plt.tight_layout()
# plt.savefig(os.path.join(out_dir, "gene_length_vs_size_log.png"), dpi=300)
# legend = plt.legend(
#     title="Mapper",
#     fontsize=12,
#     title_fontsize=14,
#     loc="center left",
#     bbox_to_anchor=(1.02, 0.5),
# )
#
# fig_legend = legend.figure
# fig_legend.canvas.draw()
# bbox = legend.get_window_extent().transformed(fig_legend.dpi_scale_trans.inverted())
# fig_legend.savefig("legend.png", dpi=300, bbox_inches=bbox)
# plt.close()
#
#
################################### COMMENTED OUT FOR SPEED #####################################


# GOOD LINE PLOTTT
plt.figure(figsize=(13, 5))
ax = sns.lineplot(
    data=df_merged.sort_values(by=["mapper", "gene_length"]),
    x="gene_length",
    y="size_MB",
    hue="mapper",
    hue_order=tool_order_a,
    palette=tool_colors,  # use your fixed color map
    linewidth=1.8,
)

ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
ax.xaxis.grid(False)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Gene Length [log(length)]")
plt.ylabel("Size [log(MB)]")
plt.legend(title="Mapper")
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "gene_length_vs_size_log_line.png"), dpi=300)
plt.close()
#

df_total_size = df_gb.groupby("mapper", as_index=False)["size_GB"].sum()

plt.figure(figsize=(7, 5))
ax = sns.barplot(
    data=df_total_size,
    x="mapper",
    y="size_GB",
    legend=False,
    edgecolor="black",
    linewidth=1.3,
    hue_order=tool_order_a,
    palette=tool_colors,
    order=tool_order_a,
    hue="mapper",
)

ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
ax.xaxis.grid(False)  # optional: turn off vertical grid lines

# --- Add labels above bars ---
for p in ax.patches:
    height = p.get_height()
    ax.text(
        p.get_x() + p.get_width() / 2.0,
        height + 0.1,  # slightly above bar
        f"{height:.1f}",  # format value (1 decimal)
        ha="center",
        va="bottom",
        weight="bold",
    )
plt.title("Combined Index Size of all 17,575 \n Genes by Mapper", weight="bold", pad=15)
plt.xlabel("Mapper")
plt.ylabel("Size in GB")
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "total_size_by_mapper.png"), dpi=300)
plt.close()
