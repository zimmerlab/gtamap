import pandas as pd
from scipy import stats
import seaborn as sns
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import os
import sys


def get_whiskers(series):
    q1, q3 = series.quantile([0.25, 0.75])
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    # Ensure whiskers are within actual min/max if no outliers
    lower = max(lower, series.min()) upper = min(upper, series.max())
    return lower, upper


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
        "lines.markersize": 10,  # marker size
    },
)

plt.rcParams["figure.figsize"] = (6.5, 4.5)  # width x height in inches


sns.set_theme(style="ticks", context="paper", palette="husl")
base_colors = sns.color_palette("husl", n_colors=6)

# tool_order_a = ["GTAMap", "Minimap2", "HISAT2", "STAR"]
tool_order_a = ["GTAMap", "Minimap2", "HISAT2", "STAR", "BWA", "Bowtie2"]
tool_colors = dict(zip(tool_order_a, base_colors))
tool_colors["GTAMap"] = sns.color_palette("Reds", 6)[4]  # bright red tone

tool_order_rna = ["GTAMap", "Minimap2", "HISAT2", "STAR"]
tool_order_dna = ["GTAMap", "BWA", "Bowtie2"]

tool_colors_rna = {k: tool_colors[k] for k in tool_order_rna}
tool_colors_dna = {k: tool_colors[k] for k in tool_order_dna}


ty = "dna"

if ty == "rna":  # rna analysis

    tool_map = {
        "hisat2": "HISAT2",
        "star": "STAR",
        "minimap2": "Minimap2",
        "bwa": "BWA",
        "bowtie": "Bowtie2",
        "gtamap": "GTAMap",
    }

    dfs = []
    for i in [0, 1, 2, 3, 4, 5]:
        # for i in [0, 1, 2, 3, 4]:
        df_sub = pd.read_csv(
            # f"./pipeline_all_rna_default/pipeline_all_{i}/benchmarking_dataset_{i}.all_results.tsv",
            f"./pipeline_all_rna_fair_final/benchmarking_dataset_{i}.all_results.tsv",
            sep="\t",
        )
        df_sub["mutation_rate"] = i
        dfs.append(df_sub)

    # Concatenate all together
    df = pd.concat(dfs, ignore_index=True)

    df["overlap"] = (df["fw_overlap_strict"] + df["rv_overlap_strict"]) / 2

    df["recall"] = df["level_1_recall"]

    df["recall_pos"] = (
        +df["fw_level_2_recall"]
        + df["rv_level_2_recall"]
        + df["fw_level_3_recall"]
        + df["rv_level_3_recall"]
    ) / 4

    df["tool"] = df["tool"].replace(tool_map)

    alpha, beta, gamma = 0.3, 0.5, 0.2
    df["GMS"] = alpha * df["overlap"] + beta * df["recall"] + gamma * df["recall_pos"]

    df_0 = df[df["mutation_rate"] == 0]
    worst10 = (
        df_0.sort_values(["tool", "GMS"])
        .groupby("tool")
        .head(100)
        .reset_index(drop=True)
    )

    print(worst10)
    worst10.to_csv("worst10.tsv", sep="\t")
    exit(1)

    # mean score per tool
    final_scores = df.groupby("tool")["GMS"].mean()

    # --- 1. Mean GeneScore per mutation rate ---
    df_mean = df.groupby(["mutation_rate", "tool"], as_index=False).agg(
        GMS=("GMS", "mean")
    )

    plt.figure()
    ax = sns.lineplot(
        data=df_mean,
        x="mutation_rate",
        y="GMS",
        hue="tool",
        hue_order=tool_order_rna,
        palette=tool_colors_rna,
        marker="o",
        linewidth=2,  # slightly slimmer, still bold
    )
    plt.title("Mean GeneMapperScore (GMS) per Mutation Rate - RNA-seq", weight="bold")
    plt.xlabel("Mutation Rate")
    plt.ylabel("Mean GMS")
    ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    ax.xaxis.grid(False)
    sns.despine()
    plt.tight_layout()
    plt.savefig("mean_line_rna.png", dpi=300)
    # plt.savefig("mean_line_rna_fair_final.png", dpi=300)
    plt.close()

    # --- 2. Boxplot with Zoom ---
    plt.figure()
    ax = sns.boxplot(
        data=df,
        x="mutation_rate",
        y="GMS",
        hue="tool",
        hue_order=tool_order_rna,
        palette=tool_colors_rna,
        fliersize=0,  # Hide outliers
        showfliers=False,  # Explicitly hide outliers
        boxprops=dict(edgecolor="black", linewidth=1.5),
        whiskerprops=dict(color="black", linewidth=1.5),
        capprops=dict(color="black", linewidth=1.5),
        medianprops=dict(color="black", linewidth=1.5),
    )

    # Calculate whisker bounds for zoom
    whisker_lows, whisker_highs = [], []
    for (mut_rate, tool), group in df.groupby(["mutation_rate", "tool"]):
        low, high = get_whiskers(group["GMS"])
        whisker_lows.append(low)
        whisker_highs.append(high)

    global_low = min(whisker_lows)
    global_high = max(whisker_highs)
    padding = 0.02 * (global_high - global_low)  # 2% padding
    ax.set_ylim(global_low - padding, global_high + padding)

    ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    ax.xaxis.grid(False)

    plt.title(
        "GeneMapperScore (GMS) Distribution per \n Mutation Rate and Tool",
        weight="bold",
    )
    plt.xlabel("Mutation Rate")
    plt.ylabel("GMS")
    plt.legend(title="Tool")
    sns.despine()
    plt.tight_layout()
    plt.savefig("box_rna.png", dpi=300)
    # plt.savefig("box_rna_fair_final.png", dpi=300)
    plt.close()

    import matplotlib.pyplot as plt
    import seaborn as sns

    # --- Define subsets ---
    df_low = df[df["mutation_rate"].isin([0, 1, 2])]
    df_high = df[df["mutation_rate"].isin([3, 4, 5])]

    # --- Create subplots ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    fig.suptitle(
        "GeneMapperScore (GMS) Distribution in RNA Dataset \n Across Mutation Rates (default)",
        weight="bold",
        fontsize=17,
    )

    # --- Helper function for whisker range ---
    def get_whiskers(series):
        q1 = series.quantile(0.25)
        q3 = series.quantile(0.75)
        iqr = q3 - q1
        low = max(series.min(), q1 - 1.5 * iqr)
        high = min(series.max(), q3 + 1.5 * iqr)
        return low, high

    # --- Plot for 0–2 ---
    sns.boxplot(
        data=df_low,
        x="mutation_rate",
        y="GMS",
        hue="tool",
        hue_order=tool_order_rna,
        palette=tool_colors_rna,
        fliersize=0,
        showfliers=False,
        ax=axes[0],
        boxprops=dict(edgecolor="black", linewidth=1.5),
        whiskerprops=dict(color="black", linewidth=1.5),
        capprops=dict(color="black", linewidth=1.5),
        medianprops=dict(color="black", linewidth=1.5),
    )
    axes[0].set_title("")
    axes[0].set_xlabel("Mutation Rate (%)", fontsize=16)
    axes[0].set_ylabel("GMS", fontsize=16)
    axes[0].legend_.remove()  # remove duplicate legend
    axes[0].yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    axes[0].xaxis.grid(False)

    # Increase tick label font size
    axes[0].tick_params(axis="x", labelsize=14)
    axes[0].tick_params(axis="y", labelsize=14)

    # --- Plot for 3–5 ---
    sns.boxplot(
        data=df_high,
        x="mutation_rate",
        y="GMS",
        hue="tool",
        hue_order=tool_order_rna,
        palette=tool_colors_rna,
        fliersize=0,
        showfliers=False,
        ax=axes[1],
        boxprops=dict(edgecolor="black", linewidth=1.5),
        whiskerprops=dict(color="black", linewidth=1.5),
        capprops=dict(color="black", linewidth=1.5),
        medianprops=dict(color="black", linewidth=1.5),
    )
    axes[1].set_title("")
    axes[1].set_xlabel("Mutation Rate (%)", fontsize=16)
    axes[1].set_ylabel("")  # shared y-axis label
    axes[1].legend_.remove()  # remove duplicate legend
    axes[1].yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    axes[1].xaxis.grid(False)
    axes[1].tick_params(axis="x", labelsize=14)
    axes[1].tick_params(axis="y", labelsize=14)

    # --- Place a single legend centered beneath both plots ---
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        title="Mapper",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.2),  # centered below plots
        ncol=len(labels),  # horizontal layout
        frameon=True,
        fontsize=12,
        title_fontsize=14,  # 🔹 title font size
    )

    # --- Compute whisker ranges separately for zoom ---
    def get_axis_limits(df_subset):
        lows, highs = [], []
        for (mut_rate, tool), group in df_subset.groupby(["mutation_rate", "tool"]):
            low, high = get_whiskers(group["GMS"])
            lows.append(low)
            highs.append(high)
        return min(lows), max(highs)

    # Left panel (0–2)
    low_1, high_1 = get_axis_limits(df_low)
    pad_1 = 0.02 * (high_1 - low_1)
    axes[0].set_ylim(low_1 - pad_1, high_1 + pad_1)

    # Right panel (3–5)
    low_2, high_2 = get_axis_limits(df_high)
    pad_2 = 0.02 * (high_2 - low_2)
    axes[1].set_ylim(low_2 - pad_2, high_2 + pad_2)

    sns.despine()
    plt.tight_layout()  # leave space for legend below
    plt.savefig(
        "box_rna_split.png",
        dpi=300,
        bbox_inches="tight",
        # "box_rna_split_fair_final.png",
        # dpi=300,
        # bbox_inches="tight",
    )  # ensure legend is included
    plt.close()
else:
    tool_map = {
        "hisat2": "HISAT2",
        "star": "STAR",
        "minimap2": "Minimap2",
        "bwa": "BWA",
        "bowtie2": "Bowtie2",
        "gtamap": "GTAMap",
    }

    dfs = []
    for i in [0, 3]:
        df_sub = pd.read_csv(
            f"./pipeline_all_dna/pipeline_all_{i}/benchmarking_dataset_{i}_dna.all_results.tsv",
            sep="\t",
        )
        df_sub["mutation_rate"] = i
        dfs.append(df_sub)

    # Concatenate all together
    df = pd.concat(dfs, ignore_index=True)

    df["overlap"] = (df["fw_overlap_strict"] + df["rv_overlap_strict"]) / 2

    df["recall"] = df["level_1_recall"]

    df["recall_pos"] = (
        +df["fw_level_2_recall"]
        + df["rv_level_2_recall"]
        + df["fw_level_3_recall"]
        + df["rv_level_3_recall"]
    ) / 4

    df["tool"] = df["tool"].replace(tool_map)

    alpha, beta, gamma = 0.3, 0.5, 0.2
    df["GMS"] = alpha * df["overlap"] + beta * df["recall"] + gamma * df["recall_pos"]

    # mean score per tool
    final_scores = df.groupby("tool")["GMS"].mean()

    # --- 1. Mean GeneScore per mutation rate ---
    df_mean = df.groupby(["mutation_rate", "tool"], as_index=False).agg(
        GMS=("GMS", "mean")
    )

    plt.figure()
    ax = sns.lineplot(
        data=df_mean,
        x="mutation_rate",
        y="GMS",
        hue="tool",
        hue_order=tool_order_dna,
        palette=tool_colors_dna,
        marker="o",
        linewidth=2,  # slightly slimmer, still bold
    )
    plt.title("Mean GeneMapperScore (GMS) per Mutation Rate - DNA-seq", weight="bold")

    plt.xlabel("Mutation Rate")
    plt.ylabel("Mean GMS")
    ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    ax.xaxis.grid(False)
    sns.despine()
    plt.tight_layout()
    plt.savefig("mean_line_dna.png", dpi=300)
    plt.close()

    # --- 2. Boxplot with Zoom ---
    plt.figure()
    ax = sns.boxplot(
        data=df,
        x="mutation_rate",
        y="GMS",
        hue="tool",
        hue_order=tool_order_dna,
        palette=tool_colors_dna,
        fliersize=0,  # Hide outliers
        showfliers=False,  # Explicitly hide outliers
        boxprops=dict(edgecolor="black", linewidth=1.5),
        whiskerprops=dict(color="black", linewidth=1.5),
        capprops=dict(color="black", linewidth=1.5),
        medianprops=dict(color="black", linewidth=1.5),
    )

    # Calculate whisker bounds for zoom
    whisker_lows, whisker_highs = [], []
    for (mut_rate, tool), group in df.groupby(["mutation_rate", "tool"]):
        low, high = get_whiskers(group["GMS"])
        whisker_lows.append(low)
        whisker_highs.append(high)

    global_low = min(whisker_lows)
    global_high = max(whisker_highs)
    padding = 0.02 * (global_high - global_low)  # 2% padding
    ax.set_ylim(global_low - padding, global_high + padding)

    ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    ax.xaxis.grid(False)

    plt.title(
        "GeneMapperScore Distribution per \n Mutation Rate and Tool", weight="bold"
    )
    plt.xlabel("Mutation Rate")
    plt.ylabel("GMS")
    plt.legend(title="Tool")
    sns.despine()
    plt.tight_layout()
    plt.savefig("box_dna.png", dpi=300)
    plt.close()

    import matplotlib.pyplot as plt
    import seaborn as sns

    # --- Define subsets ---
    df_low = df[df["mutation_rate"].isin([0, 1, 2])]
    df_high = df[df["mutation_rate"].isin([3, 4, 5])]

    # --- Create subplots ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    fig.suptitle(
        "GeneMapperScore Distribution in RNA Dataset Across Mutation Rates",
        fontsize=16,
        weight="bold",
    )

    # --- Helper function for whisker range ---
    def get_whiskers(series):
        q1 = series.quantile(0.25)
        q3 = series.quantile(0.75)
        iqr = q3 - q1
        low = max(series.min(), q1 - 1.5 * iqr)
        high = min(series.max(), q3 + 1.5 * iqr)
        return low, high

    # --- Plot for 0–2 ---
    sns.boxplot(
        data=df_low,
        x="mutation_rate",
        y="GMS",
        hue="tool",
        hue_order=tool_order_dna,
        palette=tool_colors_dna,
        fliersize=0,
        showfliers=False,
        ax=axes[0],
        boxprops=dict(edgecolor="black", linewidth=1.5),
        whiskerprops=dict(color="black", linewidth=1.5),
        capprops=dict(color="black", linewidth=1.5),
        medianprops=dict(color="black", linewidth=1.5),
    )
    axes[0].set_title("")
    axes[0].set_xlabel("Mutation Rate (%)", fontsize=16)
    axes[0].set_ylabel("GMS")
    axes[0].legend_.remove()  # remove duplicate legend
    axes[0].yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    axes[0].xaxis.grid(False)
    # Increase tick label font size
    axes[0].tick_params(axis="x", labelsize=14)
    axes[0].tick_params(axis="y", labelsize=14)

    # --- Plot for 3–5 ---
    sns.boxplot(
        data=df_high,
        x="mutation_rate",
        y="GMS",
        hue="tool",
        hue_order=tool_order_dna,
        palette=tool_colors_dna,
        fliersize=0,
        showfliers=False,
        ax=axes[1],
        boxprops=dict(edgecolor="black", linewidth=1.5),
        whiskerprops=dict(color="black", linewidth=1.5),
        capprops=dict(color="black", linewidth=1.5),
        medianprops=dict(color="black", linewidth=1.5),
    )
    axes[1].set_title("")
    axes[1].set_xlabel("Mutation Rate (%)", fontsize=16)
    axes[1].set_ylabel("")  # shared y-axis label
    axes[1].legend_.remove()  # remove duplicate legend
    axes[1].yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    axes[1].xaxis.grid(False)
    # Increase tick label font size
    axes[1].tick_params(axis="x", labelsize=14)
    axes[1].tick_params(axis="y", labelsize=14)

    # --- Place a single legend centered beneath both plots ---
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        title="Mapper",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.1),  # centered below plots
        ncol=len(labels),  # horizontal layout
        frameon=True,
        fontsize=12,
        title_fontsize=14,  # 🔹 title font size
    )

    # --- Compute whisker ranges separately for zoom ---
    def get_axis_limits(df_subset):
        lows, highs = [], []
        for (mut_rate, tool), group in df_subset.groupby(["mutation_rate", "tool"]):
            low, high = get_whiskers(group["GMS"])
            lows.append(low)
            highs.append(high)
        return min(lows), max(highs)

    # Left panel (0–2)
    low_1, high_1 = get_axis_limits(df_low)
    pad_1 = 0.02 * (high_1 - low_1)
    axes[0].set_ylim(low_1 - pad_1, high_1 + pad_1)

    # Right panel (3–5)
    low_2, high_2 = get_axis_limits(df_high)
    pad_2 = 0.02 * (high_2 - low_2)
    axes[1].set_ylim(low_2 - pad_2, high_2 + pad_2)

    sns.despine()
    plt.tight_layout()  # leave space for legend below
    plt.savefig(
        "box_dna_split.png", dpi=300, bbox_inches="tight"
    )  # ensure legend is included
    plt.close()
