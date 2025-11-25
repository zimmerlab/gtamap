import seaborn as sns
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# sns.set_theme(style="ticks", context="paper", palette="Set2")
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
        "xtick.labelsize": 15,  # x tick labels
        "ytick.labelsize": 15,  # y tick labels
        "legend.fontsize": 12,  # legend
        "lines.linewidth": 1.5,  # line width
        "lines.markersize": 6,  # marker size
    },
)

plt.rcParams["figure.figsize"] = (6.5, 5)  # width x height in inches

df = pd.read_csv("./combined_logs_filter.tsv", sep="\t")
df = df.dropna()
df = df[df["readsProcessed"] != 0]

df["gene_id"] = df["gene_id"].str.replace(r"^(binned_|global_)", "", regex=True)


df_time = pd.read_csv("./combined_times.tsv", sep="\t")

# df cols
# timestamp	readsProcessed	readsAfterFiltering	readsMapped	confidentMappings	mappingLocations	percentConfidentOfMapped	percentFiltered	percentMappedTotal	percentMappedFilter	meanLocationsPerRead	allocKB	heapAllocKB	stackAllocKB	totalFileSize	bytesProcessed	percentFileProcessed	gene_id	mutation_rate type

from matplotlib.ticker import LogFormatterMathtext

lengths = pd.read_csv(
    "./../../data/genome/genes_meta_repaired.tsv",
    names=["gene_id", "start", "length", "strand", "name", "chr"],
    sep="\t",
)

# Merge on gene_id / id (adjust if column names differ)
merged = df.merge(lengths, left_on="gene_id", right_on="gene_id")

# Plot
plt.figure()
ax = sns.scatterplot(
    data=merged,
    x="length",  # Gene length
    y="percentFiltered",  # % filtered
    hue="type",
    edgecolor="black",
    s=30,
    linewidths=0.5,
)

ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
ax.xaxis.grid(False)

ax.set_xlabel("Gene Length")
ax.set_ylabel("% Surviving Reads")
ax.set_title("Filter Efficiency vs Gene Length")
plt.legend(title="Filter Strategy")
# ax.xaxis.set_major_formatter(LogFormatterMathtext())
# ax.set_xscale("log")
sns.despine()
plt.tight_layout()
plt.savefig("filter_impact.png", dpi=300)
plt.close()


# Merge on gene_id
merged_times = df_time.merge(lengths, left_on="gene", right_on="gene_id")

pivoted = merged_times.pivot(
    index="gene_id", columns="type", values="wall_time_sec"
).reset_index()


pivoted["ratio"] = pivoted["Global"] / pivoted["Binned"]


# ratio = len(pivoted[pivoted["delta"] > 0]) / len(pivoted["delta"])
pivoted = pivoted.merge(lengths[["gene_id", "length"]], on="gene_id", how="left")
ratio_counts = {
    "r >= 1": (pivoted["ratio"] >= 1).sum(),
    "r < 1": (pivoted["ratio"] < 1).sum(),
}

# Mean ratio for each group
mean_ge1 = pivoted.loc[pivoted["ratio"] >= 1, "ratio"].mean()
mean_lt1 = pivoted.loc[pivoted["ratio"] < 1, "ratio"].mean()

print(f"Average ratio for r >= 1 (Binned faster or equal): {mean_ge1:.2f}")
print(f"Average ratio for r < 1 (Global faster): {mean_lt1:.2f}")

outlier = pivoted.loc[pivoted["ratio"].idxmax()]


# Create side-by-side figure
fig, axes = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [1, 2]})

# Left: Bar plot of counts
axes[0].bar(
    ratio_counts.keys(),
    ratio_counts.values(),
    edgecolor="black",
    linewidth=2,  # Border width
)
axes[0].set_ylabel("Number of Genes")
axes[0].set_title("Frequency of Ratio r")
axes[0].tick_params(axis="x")

# Split data into two groups
green_data = pivoted[pivoted["ratio"] >= 1]
red_data = pivoted[pivoted["ratio"] < 1]
print(len(red_data))
print(len(green_data))

# Plot each group separately
axes[1].scatter(
    green_data["length"],
    green_data["ratio"],
    edgecolor="black",
    s=30,
    label="r ≥ 1",
    linewidths=0.5,
)
axes[1].scatter(
    red_data["length"],
    red_data["ratio"],
    edgecolor="black",
    s=30,
    label="r < 1",
    linewidths=0.5,
)

# # Right: Scatter plot of ratio vs gene length
# sns.scatterplot(
#     data=pivoted,
#     x="length",
#     y="ratio",
#     edgecolor="black",
#     s=20,
#     ax=axes[1],
#     color=fill_color,
# )
axes[1].set_xlabel("Gene Length")
axes[1].set_ylabel("Ratio r = (Global / Binned)")
axes[1].set_title("Difference in Wall Time vs Gene Length")
axes[1].set_xscale("log")
axes[1].legend()
axes[1].annotate(
    outlier["gene_id"],  # Text to display
    xy=(outlier["length"], outlier["ratio"]),  # Point coordinates
    xytext=(5, 5),  # Offset for the text
    textcoords="offset points",
    fontsize=10,
    color="black",
)
sns.despine(ax=axes[1])
sns.despine(ax=axes[0])

plt.tight_layout()
plt.savefig("ratio_with_bar.png", dpi=300)
plt.close()


df = pd.read_csv("./combined_eval.tsv", sep="\t")

# Split data
green_data = pivoted[pivoted["ratio"] >= 1]
red_data = pivoted[pivoted["ratio"] < 1]

# --- Scatter plot ---
plt.figure(figsize=(5.5, 3.5))
plt.scatter(
    green_data["length"],
    np.log2(green_data["ratio"]),
    edgecolor="black",
    s=25,
    label="r ≥ 1",
    linewidths=0.5,
)
plt.scatter(
    red_data["length"],
    np.log2(red_data["ratio"]),
    edgecolor="black",
    s=25,
    label="r < 1",
    linewidths=0.5,
)


ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
ax.xaxis.grid(False)
plt.xlabel("Gene Length")
plt.ylabel("r = log2(Global / Binned)")
plt.title("Difference in Wall Time \n vs Gene Length")
plt.xscale("log")
plt.legend()

sns.despine()
plt.tight_layout()
plt.savefig("ratio_scatter.png", dpi=300)
plt.close()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- Define bin size and compute averages ---
bin_size = 10  # linear bins; for log bins see note below

bins = np.arange(pivoted["length"].min(), pivoted["length"].max() + bin_size, bin_size)

bins = np.logspace(
    np.log10(pivoted["length"].min()), np.log10(pivoted["length"].max()), num=20
)
labels = bins[:-1] + bin_size / 2

pivoted["length_bin"] = pd.cut(
    pivoted["length"], bins=bins, labels=labels, include_lowest=True
)

binned = (
    pivoted.groupby("length_bin", observed=True)["ratio"].mean().reset_index().dropna()
)
binned["length_bin"] = binned["length_bin"].astype(float)

# --- Plot line plot ---
fig, ax = plt.subplots(figsize=(6, 5))
ax.plot(binned["length_bin"], np.log2(binned["ratio"]), marker="o", lw=1)
ax.scatter(
    binned["length_bin"],
    np.log2(binned["ratio"]),
    s=35,  # marker size
    color="C0",  # fill color (matches line)
    edgecolor="black",  # black outline
    linewidth=0.5,
    zorder=3,
)

ax.set_xscale("log")  # optional
ax.set_xlabel("Gene Length (log-spaced groups)")
ax.set_ylabel("Average log2(Global / Binned)")
ax.set_title("Average Ratio per Gene Length Group")

# --- Grid style ---
ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
ax.xaxis.grid(False)

ax.axhline(0, color="black", linewidth=1, linestyle="--")

sns.despine()
plt.tight_layout()
plt.savefig("ratio_new.png", dpi=300)


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- Create logarithmic bins ---
bins = np.logspace(
    np.log10(pivoted["length"].min()), np.log10(pivoted["length"].max()), num=20
)

# Use geometric mean for bin midpoints (better for log spacing)
labels = (bins[:-1] * bins[1:]) ** 0.5

pivoted["length_bin"] = pd.cut(
    pivoted["length"], bins=bins, labels=labels, include_lowest=True
)

# --- Compute mean ratio per length bin ---
binned = (
    pivoted.groupby("length_bin", observed=True)["ratio"].mean().reset_index().dropna()
)
binned["length_bin"] = binned["length_bin"].astype(float)

# --- Define colors from Seaborn husl palette ---
palette = sns.color_palette("husl")
color_low = palette[1]  # for ratio < 1
color_high = palette[0]  # for ratio >= 1

# --- Split data ---
below = binned[binned["ratio"] < 1]
above = binned[binned["ratio"] >= 1]

# --- Plot ---
fig, ax = plt.subplots()

# Line connecting all points (neutral color)
ax.plot(binned["length_bin"], np.log2(binned["ratio"]), color="gray", lw=1, zorder=1)

# Points below 1
ax.scatter(
    below["length_bin"],
    np.log2(below["ratio"]),
    s=35,
    color=color_low,
    edgecolor="black",
    linewidth=0.5,
    label="Global < Binned",
    zorder=3,
)

# Points ≥ 1
ax.scatter(
    above["length_bin"],
    np.log2(above["ratio"]),
    s=35,
    color=color_high,
    edgecolor="black",
    linewidth=0.5,
    label="Global > Binned",
    zorder=3,
)

# --- Axes, labels, and styling ---
ax.set_xscale("log")
ax.set_xlabel("Gene Length (log-spaced groups)")
ax.set_ylabel("Average log2(Global / Binned)")
ax.set_title("Average Ratio per Gene Length Group")

# Grid and reference line
ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
ax.xaxis.grid(False)
ax.axhline(0, color="black", linewidth=1, linestyle="--")

# Legend and clean up
ax.legend(title="Ratio Category")
sns.despine()
plt.tight_layout()
plt.savefig("ratio_new_colored.png", dpi=300)
plt.close()

# Filter bins where ratio > 1
above_zero = binned[binned["ratio"] > 1]

# Get the first bin (smallest gene length) where log2(ratio) > 0
first_bin = above_zero.sort_values("length_bin")["length_bin"].iloc[0]

print(f"First length bin where log2(ratio) > 0: {first_bin:.1f}")
