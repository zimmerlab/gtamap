import pandas as pd
from scipy import stats
import seaborn as sns

from matplotlib.patches import Patch
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import os
import sys


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


tool_map = {
    "hisat2": "HISAT2",
    "star": "STAR",
    "minimap2": "Minimap2",
    "bwa": "BWA",
    "bowtie": "Bowtie2",
    "gtamap": "GTAMap",
}


dfs = []
for i in range(6):
    df_sub = pd.read_csv(
        f"./pipeline_all_rna_fair_final/benchmarking_dataset_{i}.all_results.tsv",
        sep="\t",
    )
    df_sub["mutation_rate"] = i
    df_sub["scenario"] = "sensitive"
    dfs.append(df_sub)
    df_sub = pd.read_csv(
        f"./pipeline_all_rna_default/pipeline_all_{i}/benchmarking_dataset_{i}.all_results.tsv",
        sep="\t",
    )
    df_sub["mutation_rate"] = i
    df_sub["scenario"] = "default"
    dfs.append(df_sub)

    df_sub = pd.read_csv(
        f"./genome_wide_simulation_sim/dataset_{i}.all_results.tsv",
        sep="\t",
    )
    df_sub["mutation_rate"] = i
    df_sub["scenario"] = "sensitive_genome"
    dfs.append(df_sub)


df = pd.concat(dfs, ignore_index=True)


# gtamap was called the same way across both scenarios, only keep first (default call)
df = df[~((df["tool"] == "gtamap") & (df["scenario"] == "sensitive"))]

gtamap_default = df[(df["tool"] == "gtamap") & (df["scenario"] == "default")].copy()
gtamap_default["scenario"] = "sensitive"
df = pd.concat([df, gtamap_default], ignore_index=True)

df["overlap"] = (df["fw_overlap_strict"] + df["rv_overlap_strict"]) / 2
df["recall"] = df["level_1_recall"]
df["recall_pos"] = (
    df["fw_level_2_recall"]
    + df["rv_level_2_recall"]
    + df["fw_level_3_recall"]
    + df["rv_level_3_recall"]
) / 4
df["tool"] = df["tool"].replace(tool_map)

alpha, beta, gamma = 0.3, 0.5, 0.2
df["GMS"] = alpha * df["overlap"] + beta * df["recall"] + gamma * df["recall_pos"]


df["tool_scenario"] = df["tool"] + " (" + df["scenario"] + ")"


summary = (
    df.groupby(["tool", "mutation_rate", "scenario"])["GMS"]
    .agg(["mean", "median"])
    .reset_index()
    .sort_values(["tool", "mutation_rate"])
)

print(summary)

plt.figure(figsize=(10, 6))
sns.lineplot(
    data=summary,
    x="mutation_rate",
    y="median",
    hue="tool",
    style="scenario",
    markers=True,
    dashes=False,
)

plt.title("Mean GMS per Tool and Mutation Rate")
plt.xlabel("Mutation Rate")
plt.ylabel("Median GMS")
plt.legend(title="Tool / Scenario", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.show()


fig, axes = plt.subplots(3, 2, figsize=(14, 12), sharex=False)
axes = axes.flatten()

fig.suptitle(
    "GeneMapperScore (GMS) per Mutation Rate and Scenario (RNA dataset)",
    fontsize=23,
    weight="bold",
)


def plot_panel(ax, df_full, title, mutation_rates, add_x):
    df_subset = df_full[df_full["mutation_rate"].isin(mutation_rates)].copy()

    tool_scenario_order = []
    for t in tool_order_rna:
        if t == "GTAMap":
            tool_scenario_order.append(f"{t} (default)")
        else:
            tool_scenario_order.extend(
                [f"{t} (default)", f"{t} (sensitive)", f"{t} (sensitive_genome)"]
            )

    palette = {}
    for t in tool_order_rna:
        palette[f"{t} (default)"] = tool_colors_rna[t]
        if t != "GTAMap":
            palette[f"{t} (sensitive)"] = tool_colors_rna[t]
            palette[f"{t} (sensitive_genome)"] = tool_colors_rna[t]

    sns.boxplot(
        data=df_subset,
        x="mutation_rate",
        y="GMS",
        hue="tool_scenario",
        hue_order=tool_scenario_order,
        palette=palette,
        fliersize=0,
        showfliers=False,
        ax=ax,
        boxprops=dict(edgecolor="black", linewidth=1.5),
        whiskerprops=dict(color="black", linewidth=1.5),
        capprops=dict(color="black", linewidth=1.5),
        medianprops=dict(color="black", linewidth=1.5),
    )

    num_mutation_rates = len(mutation_rates)
    num_tool_scenarios = len(tool_scenario_order)

    for i, patch in enumerate(ax.patches):
        tool_idx = i // num_mutation_rates % 10
        label = tool_scenario_order[tool_idx]

        if "GTAMap" in label:
            patch.set_hatch(None)
        elif "(default)" in label:
            patch.set_hatch(None)
        elif "(sensitive_genome)" in label:
            patch.set_hatch("xx")
        else:
            patch.set_hatch("////")

        patch.set_edgecolor("black")
        patch.set_linewidth(1.5)

    ax.set_title(title, fontsize=18, weight="bold")
    ax.set_xlabel("Mutation Rate (%)" if add_x else "", fontsize=16)
    ax.set_ylabel("GMS", fontsize=16)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.5)
    ax.xaxis.grid(False)
    if ax.get_legend() is not None:
        ax.get_legend().remove()
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.set_xticks(range(len(mutation_rates)))
    ax.set_xticklabels(mutation_rates)
    ax.set_xlim(-0.5, len(mutation_rates) - 0.5)


plot_panel(axes[0], df, "Mutation Rate = 0", [0], False)
plot_panel(axes[1], df, "Mutation Rate = 1", [1], False)
plot_panel(axes[2], df, "Mutation Rate = 2", [2], False)
plot_panel(axes[3], df, "Mutation Rate = 3", [3], False)
plot_panel(axes[4], df, "Mutation Rate = 4", [4], True)
plot_panel(axes[5], df, "Mutation Rate = 5", [5], True)


legend_elements = []
for tool in tool_order_rna:
    if tool == "GTAMap":
        legend_elements.append(
            Patch(
                facecolor=tool_colors_rna[tool],
                edgecolor="black",
                label=f"{tool}",
            )
        )
    else:
        legend_elements.append(
            Patch(
                facecolor=tool_colors_rna[tool],
                edgecolor="black",
                label=f"{tool} (default)",
            )
        )
        legend_elements.append(
            Patch(
                facecolor=tool_colors_rna[tool],
                edgecolor="black",
                hatch="////",
                label=f"{tool} (sensitive)",
            )
        )
        legend_elements.append(
            Patch(
                facecolor=tool_colors_rna[tool],
                edgecolor="black",
                hatch="xx",
                label=f"{tool} (sensitive_genome)",
            )
        )

fig.legend(
    handles=legend_elements,
    ncol=5,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.02),
    frameon=True,
    fontsize=12,
    title="Mapper / Scenario",
    title_fontsize=14,
)

plt.tight_layout(rect=[0, 0.05, 0.97, 0.95])
plt.savefig(
    "box_rna_split_6panels_default_sensitive_extended.png", dpi=300, bbox_inches="tight"
)
plt.close()
