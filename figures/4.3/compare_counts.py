# --- Settings ---
count_files = {
    "HISAT2": "/mnt/studtemp/weyrich/results_genome_wide/counts_HISAT2.txt",
    "STAR": "/mnt/studtemp/weyrich/results_genome_wide/counts_STAR.txt",
    "Minimap2": "/mnt/studtemp/weyrich/results_genome_wide/counts_Minimap2.txt",
}

#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

# --- Styling ---
sns.set_theme(style="ticks", context="paper", palette="husl")

tool_order = ["Minimap2", "HISAT2", "STAR"]
base_colors = sns.color_palette("husl", n_colors=len(tool_order))
tool_colors = dict(zip(tool_order, base_colors))

sns.set_context(
    "paper",
    rc={
        "font.size": 15,
        "axes.titlesize": 17,
        "axes.labelsize": 15,
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
        "legend.fontsize": 12,
        "lines.linewidth": 1.5,
        "lines.markersize": 6,
    },
)
plt.rcParams["figure.figsize"] = (6.5, 5)


def style_plot(ax):
    """Apply despine and clean ticks"""
    # sns.despine(ax=ax, trim=True)
    ax.tick_params(axis="both", which="major", length=5, width=1)


# --- Paths ---
output_dir = "counts_comparison"
os.makedirs(output_dir, exist_ok=True)

# --- Load counts ---
dfs = []
for mapper, path in count_files.items():
    df = pd.read_csv(path, sep="\t", comment="#", index_col=0)
    df = df.iloc[:, -1].rename(mapper)
    dfs.append(df)

counts_df = pd.concat(dfs, axis=1)
counts_df = counts_df[counts_df.sum(axis=1) > 0]

# --- Log2 transform for plots ---
counts_log2 = np.log2(counts_df + 1)

# --- Correlation heatmap ---
plt.figure(figsize=(6.5, 5))
ax = sns.heatmap(counts_log2.corr(), annot=True, cmap="coolwarm", vmin=0.5, vmax=1)
ax.set_title("Correlation of gene counts across mappers")
style_plot(ax)
plt.tight_layout()
plt.savefig(f"{output_dir}/mapper_counts_correlation.png", dpi=300)
plt.close()

# --- Pairwise scatterplots ---
for i, m1 in enumerate(tool_order):
    for m2 in tool_order[i + 1 :]:
        plt.figure()
        ax = sns.scatterplot(
            x=counts_df[m1],
            y=counts_df[m2],
            color=tool_colors[m1],
            edgecolor="black",
            linewidth=0.5,
        )
        ax.set_xlabel(f"{m1} counts")
        ax.set_ylabel(f"{m2} counts")
        ax.set_title(f"{m1} vs {m2} Gene Counts")
        style_plot(ax)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{m1}_vs_{m2}_scatter.png", dpi=300)
        plt.close()

# --- Save merged counts ---
counts_df.to_csv(f"{output_dir}/merged_gene_counts.csv")
print("Done! Plots and merged counts saved in:", output_dir)


# --- Pairwise scatterplots with enhancements ---
scatter_color = sns.color_palette("husl")[0]  # use first husl color for all points

for i, m1 in enumerate(tool_order):
    for m2 in tool_order[i + 1 :]:
        plt.figure(figsize=(6.5, 6.5))  # square figure
        ax = sns.scatterplot(
            x=counts_df[m1],
            y=counts_df[m2],
            color=scatter_color,
            edgecolor="black",
            linewidth=0.3,
            s=35,
        )

        # Diagonal line y=x
        max_val = max(counts_df[[m1, m2]].max()) * 1.05
        ax.plot([0, max_val], [0, max_val], ls="--", color="gray", lw=1)

        # Square axes
        ax.set_xlim(0, max_val)
        ax.set_ylim(0, max_val)
        ax.set_aspect("equal", adjustable="box")

        # Labels, title
        ax.set_xlabel(f"{m1} Counts")
        ax.set_ylabel(f"{m2} Counts")
        ax.tick_params(axis="x", labelsize=12)
        ax.tick_params(axis="y", labelsize=12)
        ax.set_title(f"{m1} vs {m2} Gene Counts", size=17, weight="bold")

        # Grid
        ax.grid(True, linestyle="--", alpha=0.5)

        # Styling
        style_plot(ax)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{m1}_vs_{m2}_scatter.png", dpi=300)
        plt.close()
