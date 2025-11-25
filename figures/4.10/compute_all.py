import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os


res = sys.argv[1]
out_dir = sys.argv[2]
read_ty = sys.argv[3]

tool_order_a = ["GTAMap", "Minimap2", "HISAT2", "STAR", "BWA", "Bowtie2", "STAR_SP"]
base_colors = sns.color_palette("husl", n_colors=7)


tool_colors = dict(zip(tool_order_a, base_colors))
tool_colors["GTAMap"] = sns.color_palette("Reds", 6)[4]  # bright red tone

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

df = pd.read_csv(res, sep="\t")
# for col in df.columns.difference(["geneID", "tool"]):
#     df[col] = pd.to_numeric(df[col], errors="coerce")
#
# print(df.dtypes)
# na_rows = df[df.isna().any(axis=1)]
# print(na_rows)
# df = pd.read_csv("./0309_01/benchmarking_dataset_2_baseline.all_results.tsv", sep="\t")
if read_ty == "rna":
    print(f"gta {len(df[df["tool"] == "gtamap"])}")
    print(f"star {len(df[df["tool"] == "star"])}")
    print(f"hisat2  {len(df[df["tool"] == "hisat2"])}")
    print(f"minimap2 {len(df[df["tool"] == "minimap2"])}")
else:
    print(f"gta {len(df[df["tool"] == "gtamap"])}")
    print(f"bwa {len(df[df["tool"] == "bwa"])}")
    print(f"bowtie2  {len(df[df["tool"] == "bowtie2"])}")


sns.set_theme(style="ticks", context="paper", palette="colorblind")

fill_color = "#CFD5E8"

sns.set_context(
    "paper",
    rc={
        "font.size": 15,  # base font size
        "axes.titlesize": 17,  # title
        "axes.labelsize": 15,  # axis labels
        "xtick.labelsize": 14,  # x tick labels
        "ytick.labelsize": 14,  # y tick labels
        "legend.fontsize": 12,  # legend
        "legend.titlesize": 15,  # legend
        "lines.linewidth": 1.5,  # line width
        "lines.markersize": 10,  # marker size
    },
)

plt.rcParams["figure.figsize"] = (6.5, 4.5)  # width x height in inches

name_conv = {
    "perfect_matches_fw": "FW",
    "perfect_matches_rv": "RV",
    "true_positives": "TP",
    "false_positives": "FP",
    "false_negatives": "FN",
    "precision": "Precision",
    "level_1_recall": "Recall",
    "f1_score": "F1",
    "fw_overlap": "FW Overlap",
    "rv_overlap": "RV Overlap",
    "fw_level_2_recall": "FW Level 2",
    "rv_level_2_recall": "RV Level 2",
    "fw_level_3_recall": "FW Level 3",
    "rv_level_3_recall": "RV Level 3",
    "avg_missaligned_positions_fw": "FW Positions",
    "avg_missaligned_positions_rv": "RV Positions",
    "mm_delta_fw": "FW Delta",
    "mm_delta_rv": "RV Delta",
    "avg_missaligned_positions_fw_strict": "FW Positions",
    "avg_missaligned_positions_rv_strict": "RV Positions",
    "fw_overlap_strict": "FW Overlap",
    "rv_overlap_strict": "RV Overlap",
    "avg_entries_tp_fw": "FW TP",
    "avg_entries_tp_rv": "RV TP",
    "avg_entries_fp_fw": "FW FP",
    "avg_entries_fp_rv": "RV FP",
}


# print(df[df["avg_missaligned_positions_fw"] > 100])


def get_whiskers(series):
    q1, q3 = series.quantile([0.25, 0.75])
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    # Ensure whiskers are within actual min/max if no outliers
    lower = max(lower, series.min())
    upper = min(upper, series.max())
    return lower, upper


def plot_full_gene_metrics(row, thresholds=None, name_conv=None, save_path=None):
    """
    Plots all relevant metrics for a gene row, highlighting metrics that violate thresholds.

    Parameters:
    - row: pandas Series (one row of the original df)
    - thresholds: dict of thresholds per metric (optional)
    - name_conv: dict mapping metric names to short names (optional)
    - save_path: file path to save figure (optional)
    """
    if thresholds is None:
        thresholds = {
            "level_1_recall": 0.9,
            "fw_level_2_recall": 0.7,
            "rv_level_2_recall": 0.7,
            "fw_level_3_recall": 0.5,
            "rv_level_3_recall": 0.5,
            "fw_overlap": 0.8,
            "rv_overlap": 0.8,
        }

    metrics = list(thresholds.keys())
    values = [row[m] for m in metrics]

    # Determine colors: highlight violations
    colors = [
        (
            "#E74C3C" if row[m] < thresholds[m] else "#CFD5E8"
        )  # red for violation, default fill color otherwise
        for m in metrics
    ]

    labels = [name_conv[m] if name_conv and m in name_conv else m for m in metrics]

    plt.figure(figsize=(8, 4))
    bars = plt.bar(labels, values, color=colors, edgecolor="#324068")

    # annotate values
    for bar, val in zip(bars, values):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            val + 0.02,
            f"{val:.2f}",
            ha="center",
            fontsize=10,
        )

    plt.ylim(0, 1.05)
    plt.ylabel("Value")
    plt.title(f"Gene: {row['geneID']}", weight="bold")
    plt.xticks(rotation=25)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, facecolor="white")

    plt.show()


def get_all_outliers(df, metrics=None):
    thresholds = {
        "level_1_recall": 0.96,
        "fw_level_2_recall": 0.8,
        "rv_level_2_recall": 0.8,
        "fw_level_3_recall": 0.8,
        "rv_level_3_recall": 0.8,
        "fw_overlap": 0.8,
        "rv_overlap": 0.8,
    }

    outlier_rows = []
    for metric, thresh in thresholds.items():
        violating = df[df[metric] < thresh]
        if not violating.empty:
            temp = violating.copy()
            temp["Metric"] = metric
            temp["Value"] = temp[metric]
            outlier_rows.append(temp)

    if outlier_rows:
        return pd.concat(outlier_rows, ignore_index=True)
    else:
        return pd.DataFrame(columns=list(df.columns) + ["Metric", "Value"])


# def fancy_boxplot(
#     df, metrics, title, ylabel, filename, out_dir, figsize=(8, 6), ylim=None
# ):
#     df_long = df.melt(
#         id_vars="tool",  # include tool as grouping variable
#         value_vars=metrics,
#         var_name="Metric",
#         value_name="Value",
#     )
#     df_long["Metric_label"] = df_long["Metric"].map(name_conv)
#
#     plt.figure(figsize=figsize)
#
#     if len(metrics) > 1:
#         ax = sns.boxplot(
#             data=df_long,
#             x="Metric_label",
#             y="Value",
#             hue="tool",
#             width=0.6,
#             showfliers=True,
#             boxprops=dict(edgecolor="black", linewidth=1.5),
#             whiskerprops=dict(color="black", linewidth=1.5),
#             capprops=dict(color="black", linewidth=1.5),
#             medianprops=dict(color="black", linewidth=1.5),
#         )
#         plt.legend(title="Tool", fontsize=10, title_fontsize=11)
#         ax.set_xlabel("Metric", fontsize=13)
#
#     else:
#         # single metric: x = tool, no hue
#         metric_label = df_long["Metric_label"].iloc[0]
#         ax = sns.boxplot(
#             data=df_long,
#             x="tool",
#             y="Value",
#             width=0.6,
#             showfliers=True,
#             boxprops=dict(facecolor=fill_color, edgecolor=edge_color, linewidth=1.5),
#             whiskerprops=dict(color=edge_color, linewidth=1.5),
#             capprops=dict(color=edge_color, linewidth=1.5),
#             medianprops=dict(color=edge_color, linewidth=1.5),
#         )
#         ax.set_xlabel("Tool", fontsize=13)
#
#     ax.set_title(f"{title}", fontsize=16, weight="bold", pad=15)
#     ax.set_ylabel(ylabel, fontsize=13)
#     if ylim:
#         ax.set_ylim(ylim)
#     plt.xticks(fontsize=12)
#     plt.yticks(fontsize=12)
#     sns.despine()
#
#     out = os.path.join(out_dir, filename)
#     plt.tight_layout()
#     plt.savefig(out, dpi=300)
#     plt.close()
#


def fancy_boxplot(
    df,
    metrics,
    title,
    ylabel,
    filename,
    out_dir,
    figsize=(8, 6),
    ylim=None,
    tool_order=None,
    metric_order=None,
    tool_colors=None,
):
    df_long = df.melt(
        id_vars="tool",
        value_vars=metrics,
        var_name="Metric",
        value_name="Value",
    )
    df_long["Metric_label"] = df_long["Metric"].map(name_conv)

    plt.figure(figsize=figsize)

    if len(metrics) > 1:
        ax = sns.boxplot(
            data=df_long,
            x="Metric_label",
            y="Value",
            hue="tool",
            order=metric_order,  # <-- order of metrics
            hue_order=tool_order,  # <-- order of tools in legend + colors
            width=0.6,
            showfliers=True,
            boxprops=dict(edgecolor="black", linewidth=1.5),
            whiskerprops=dict(color="black", linewidth=1.5),
            capprops=dict(color="black", linewidth=1.5),
            medianprops=dict(color="black", linewidth=1.5),
            palette=tool_colors,
        )
        ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
        ax.xaxis.grid(False)  # optional: turn off vertical grid lines

        plt.legend(title="Tool")
        ax.set_xlabel("")

    else:
        metric_label = df_long["Metric_label"].iloc[0]
        ax = sns.boxplot(
            data=df_long,
            x="tool",
            y="Value",
            order=tool_order,  # <-- order of tools on x-axis
            hue="tool",
            width=0.6,
            showfliers=True,
            boxprops=dict(edgecolor=edge_color, linewidth=1.5),
            whiskerprops=dict(color=edge_color, linewidth=1.5),
            capprops=dict(color=edge_color, linewidth=1.5),
            medianprops=dict(color=edge_color, linewidth=1.5),
            palette=tool_colors,
        )
        ax.set_xlabel("Tool")

    ax.set_title(f"{title}", weight="bold")
    ax.set_ylabel(ylabel)
    ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    ax.xaxis.grid(False)  # optional: turn off vertical grid lines
    if ylim:
        ax.set_ylim(ylim)
    sns.despine()

    out = os.path.join(out_dir, filename)
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()


def fancy_boxplot_zoom(
    df,
    metrics,
    title,
    ylabel,
    filename,
    out_dir,
    figsize=(8, 6),
    tool_order=None,
    tool_colors=None,
):
    df_long = df.melt(
        id_vars="tool",  # include tool as grouping variable
        value_vars=metrics,
        var_name="Metric",
        value_name="Value",
    )
    df_long["Metric_label"] = df_long["Metric"].map(name_conv)

    plt.figure(figsize=figsize)

    if len(metrics) > 1:
        ax = sns.boxplot(
            data=df_long,
            x="Metric_label",
            y="Value",
            hue="tool",
            hue_order=tool_order,  # <-- order of tools in legend + colors
            width=0.6,
            showfliers=False,  # Hides outliers
            boxprops=dict(edgecolor="black", linewidth=1.5),
            whiskerprops=dict(color="black", linewidth=1.5),
            capprops=dict(color="black", linewidth=1.5),
            palette=tool_colors,
            medianprops=dict(color="black", linewidth=1.5),
        )
        plt.legend(title="Tool")
        ax.set_xlabel("")
        ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
        ax.xaxis.grid(False)  # optional: turn off vertical grid lines

    else:
        # single metric: x = tool, no hue
        metric_label = df_long["Metric_label"].iloc[0]
        ax = sns.boxplot(
            data=df_long,
            x="tool",
            y="Value",
            width=0.6,
            order=tool_order,  # <-- order of tools on x-axis
            showfliers=False,  # Hides outliers
            boxprops=dict(edgecolor=edge_color, linewidth=1.5),
            whiskerprops=dict(color=edge_color, linewidth=1.5),
            capprops=dict(color=edge_color, linewidth=1.5),
            palette=tool_colors,
            medianprops=dict(color=edge_color, linewidth=1.5),
        )
        ax.set_xlabel("Tool", fontsize=13)

    # Title and labels
    ax.set_title(f"{title}", weight="bold")
    ax.set_ylabel(ylabel)
    ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    ax.xaxis.grid(False)  # optional: turn off vertical grid lines

    if len(metrics) > 1:
        # Group by Metric and Tool
        whisker_lows, whisker_highs = [], []
        for (metric, tool), group in df_long.groupby(["Metric_label", "tool"]):
            low, high = get_whiskers(group["Value"])
            whisker_lows.append(low)
            whisker_highs.append(high)
    else:
        # Single metric: group by tool
        whisker_lows, whisker_highs = [], []
        for tool, group in df_long.groupby("tool"):
            low, high = get_whiskers(group["Value"])
            whisker_lows.append(low)
            whisker_highs.append(high)

    global_low = min(whisker_lows)
    global_high = max(whisker_highs)
    padding = 0.02 * (global_high - global_low)  # 2% padding
    ax.set_ylim(global_low - padding, global_high + padding)

    # Styling
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    sns.despine()

    # Save
    out = os.path.join(out_dir, filename)
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()


df["overlap"] = (df["fw_overlap_strict"] + df["rv_overlap_strict"]) / 2

df["recall"] = df["level_1_recall"]

df["recall_pos"] = (
    +df["fw_level_2_recall"]
    + df["rv_level_2_recall"]
    + df["fw_level_3_recall"]
    + df["rv_level_3_recall"]
) / 4


alpha, beta, gamma = 0.3, 0.5, 0.2
df["GeneScore"] = alpha * df["overlap"] + beta * df["recall"] + gamma * df["recall_pos"]

# mean score per tool
final_scores = df.groupby("tool")["GeneScore"].mean()
print("Final TSV Score per Tool:\n", final_scores)

fill_color = "#CFD5E8"
edge_color = "#324068"

tool_order_a = []
if read_ty == "rna":
    tool_order_a = ["GTAMap", "Minimap2", "HISAT2", "STAR"]
    name_conv_tools = {
        "gtamap": "GTAMap",
        "hisat2": "HISAT2",
        "minimap2": "Minimap2",
        "star": "STAR",
    }
    # tool_order_a = ["GTAMap", "STAR", "STAR_SP"]
    # name_conv_tools = {
    #     "gtamap": "GTAMap",
    #     "star": "STAR",
    #     "star_SP": "STAR_SP",
    # }
    df["tool"] = df["tool"].map(name_conv_tools)
else:
    tool_order_a = ["GTAMap", "BWA", "Bowtie2"]
    name_conv_tools = {
        "gtamap": "GTAMap",
        "bwa": "BWA",
        "bowtie2": "Bowtie2",
    }
    df["tool"] = df["tool"].map(name_conv_tools)


fancy_boxplot(
    df,
    ["avg_missaligned_positions_fw", "avg_missaligned_positions_rv"],
    "Avg Misaligned Positions",
    "Positions",
    "mis_pos.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

fancy_boxplot_zoom(
    df,
    ["avg_missaligned_positions_fw", "avg_missaligned_positions_rv"],
    "Avg Misaligned Positions",
    "Positions",
    "mis_pos_zoom.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

fancy_boxplot(
    df,
    ["avg_missaligned_positions_fw_strict", "avg_missaligned_positions_rv_strict"],
    "Strict Avg. Misaligned Positions",
    "Positions",
    "strict_mis_pos.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

fancy_boxplot_zoom(
    df,
    ["avg_missaligned_positions_fw_strict", "avg_missaligned_positions_rv_strict"],
    "Strict Avg. Misaligned Positions",
    "Positions",
    "strict_mis_pos_zoom.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

fancy_boxplot(
    df,
    ["level_1_recall"],
    "Level 1 Recall [QNAME] (Gene-wise)",
    "Score",
    "recall.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

fancy_boxplot_zoom(
    df,
    ["level_1_recall"],
    "Level 1 Recall [QNAME]",
    "Score",
    "recall_zoom.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

fancy_boxplot(
    df,
    ["precision"],
    "Precision [QNAME]",
    "Score",
    "precision.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)
fancy_boxplot(
    df,
    ["f1_score"],
    "F1 [QNAME]",
    "Score",
    "f1.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

fancy_boxplot(
    df,
    ["fw_overlap", "rv_overlap"],
    "FW / RV Overlap with Ground Truth Intervals",
    "Percent",
    "overlap.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

fancy_boxplot_zoom(
    df,
    ["fw_overlap", "rv_overlap"],
    "FW / RV Overlap with Ground Truth Intervals",
    "Percent",
    "overlap_zoom.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

fancy_boxplot(
    df,
    ["fw_overlap_strict", "rv_overlap_strict"],
    "Strict FW / RV Overlap with Ground Truth Intervals",
    "Percent",
    "strict_overlap.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

fancy_boxplot_zoom(
    df,
    ["fw_overlap_strict", "rv_overlap_strict"],
    "Strict FW / RV Overlap with Ground Truth Intervals",
    "Percent",
    "strict_overlap_zoom.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)


fancy_boxplot(
    df,
    [
        "fw_level_2_recall",
        "rv_level_2_recall",
        "fw_level_3_recall",
        "rv_level_3_recall",
    ],
    "Level 2 Recall [QNAME, START] and \n Level 3 Recall [QNAME, START, STOP]",
    "Score",
    "recall_levels.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

fancy_boxplot_zoom(
    df,
    [
        "fw_level_2_recall",
        "rv_level_2_recall",
        "fw_level_3_recall",
        "rv_level_3_recall",
    ],
    "Level 2 Recall [QNAME, START] and \n Level 3 Recall [QNAME, START, STOP]",
    "Score",
    "recall_levels_zoom.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

fancy_boxplot(
    df,
    ["mm_delta_fw", "mm_delta_rv"],
    "FW / RV mm delta (mapped-true)",
    "Delta",
    "mm_delta.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

fancy_boxplot_zoom(
    df,
    ["mm_delta_fw", "mm_delta_rv"],
    "FW / RV mm delta (mapped-true)",
    "Delta",
    "mm_delta_zoom.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)


fancy_boxplot(
    df,
    [
        "avg_entries_tp_fw",
        "avg_entries_tp_rv",
        "avg_entries_fp_fw",
        "avg_entries_fp_rv",
    ],
    "Avg. Amount of SAM Entries per TP and FP",
    "Avg. Amount",
    "entries.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

fancy_boxplot_zoom(
    df,
    [
        "avg_entries_tp_fw",
        "avg_entries_tp_rv",
        "avg_entries_fp_fw",
        "avg_entries_fp_rv",
    ],
    "Avg. Amount of SAM Entries per TP and FP",
    "Amount",
    "entries_fancy.png",
    out_dir,
    tool_order=tool_order_a,
    tool_colors=tool_colors,
)

metrics_to_check = [
    "avg_missaligned_positions_fw",
    "avg_missaligned_positions_rv",
    "precision",
    "level_1_recall",
    "f1_score",
    "fw_overlap",
    "rv_overlap",
    "fw_level_2_recall",
    "rv_level_2_recall",
    "fw_level_3_recall",
    "rv_level_3_recall",
    "mm_delta_fw",
    "mm_delta_rv",
]

outliers_df = get_all_outliers(df, metrics_to_check)
print(outliers_df)
outliers_df.to_csv(f"{out_dir}/outliers.tsv", sep="\t")

print(outliers_df["Metric"].value_counts())


for tool in outliers_df["tool"].unique():  # <-- missing column selection
    tool_df = outliers_df[outliers_df["tool"] == tool]
    outliers = len(tool_df["geneID"].unique())
    print(f"Number of outliers in {tool}: {outliers}")
    print(tool_df["Metric"].value_counts())

# print(f"Total outliers found: {len(outliers_df)}")

print(outliers_df[outliers_df["tool"] == "gtamap"])

print(outliers_df[outliers_df["geneID"] == "ENSG00000162825"])


idx = outliers_df.groupby("tool")["GeneScore"].idxmin()

# Use those indices to get the full rows, including gene_name
min_sizes_with_gene = outliers_df.loc[idx, ["tool", "GeneScore", "geneID"]].reset_index(
    drop=True
)

print(min_sizes_with_gene)

# # Filter for STAR
# star_df = df[df["tool"] == "STAR"]
#
# # Sort by GeneScore ascending (worst = lowest score)
# star_sorted = star_df.sort_values(by="GeneScore", ascending=True)
#
# # Select top 200 worst-performing genes
# top200_worst_star = star_sorted.head(200)
#
# # Display the relevant columns
# output_path = "top200_worst_star_genes.txt"
# top200_worst_star["geneID"].to_csv(output_path, index=False, header=False)
#
# print(top200_worst_star[top200_worst_star["geneID"] == "ENSG00000214300"])


##### STAR SP ###########
# # Pivot the table so each gene has one STAR and one STAR_SP score
# pivot_df = df.pivot_table(index="geneID", columns="tool", values="GeneScore")
#
# # Filter to only rows where both STAR and STAR_SP exist
# pivot_df = pivot_df.dropna(subset=["STAR", "STAR_SP"])
#
#
# # Optional: get a red tone from Seaborn's HUSL palette
# red_color = sns.color_palette("husl", 8)[0]  # first color is a vivid red/pink hue
#
# # Scatter plot
# plt.figure(figsize=(6, 6))
# plt.scatter(
#     pivot_df["STAR_SP"],
#     pivot_df["STAR"],
#     color=red_color,
#     edgecolor="black",
#     s=40,
# )
#
# # Add diagonal y = x line
# plt.plot([0, 1], [0, 1], color="gray", linestyle="--", linewidth=1)
#
# # Axis labels, title, limits
# plt.ylabel("STAR")
# plt.xlabel("STAR_SP")
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.title(
#     "Comparison of GMS between \n STAR and STAR_SP (Second Pass)",
#     weight="bold",
# )
#
# plt.grid(True, linestyle=":", alpha=0.6)
# plt.tight_layout()
# plt.savefig("star_vs_starsp.png", dpi=300)
#
#
# # Pivot the table so each gene has one STAR and one STAR_SP score
# pivot_df = df.pivot_table(index="geneID", columns="tool", values="GeneScore")
#
# # Filter to only rows where both STAR and STAR_SP exist
# pivot_df = pivot_df.dropna(subset=["STAR", "GTAMap"])
#
#
# # Optional: get a red tone from Seaborn's HUSL palette
# red_color = sns.color_palette("husl", 8)[0]  # first color is a vivid red/pink hue
#
# # Scatter plot
# plt.figure(figsize=(6, 6))
# plt.scatter(
#     pivot_df["STAR"],
#     pivot_df["GTAMap"],
#     color=red_color,
#     edgecolor="black",
#     s=40,
# )
#
# # Add diagonal y = x line
# plt.plot([0, 1], [0, 1], color="gray", linestyle="--", linewidth=1)
#
# # Axis labels, title, limits
# plt.xlabel("STAR")
# plt.ylabel("GTAMap")
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.title("Comparison of GMS between \n STAR and GTAMap", weight="bold")
#
# plt.grid(True, linestyle=":", alpha=0.6)
# plt.tight_layout()
# plt.savefig("star_vs_gtamap.png", dpi=300)
#
#
# # Pivot the table so each gene has one STAR and one STAR_SP score
# pivot_df = df.pivot_table(index="geneID", columns="tool", values="GeneScore")
#
# # Filter to only rows where both STAR and STAR_SP exist
# pivot_df = pivot_df.dropna(subset=["STAR_SP", "GTAMap"])
#
#
# # Optional: get a red tone from Seaborn's HUSL palette
# red_color = sns.color_palette("husl", 8)[0]  # first color is a vivid red/pink hue
#
# # Scatter plot
# plt.figure(figsize=(6, 6))
# plt.scatter(
#     pivot_df["STAR_SP"],
#     pivot_df["GTAMap"],
#     color=red_color,
#     edgecolor="black",
#     s=40,
# )
#
# # Add diagonal y = x line
# plt.plot([0, 1], [0, 1], color="gray", linestyle="--", linewidth=1)
#
# # Axis labels, title, limits
# plt.xlabel("STAR_SP")
# plt.ylabel("GTAMap")
# plt.xlim(0, 1)
# plt.ylim(0, 1)
# plt.title(
#     "Comparison of GMS between \n STAR_SP (Second Pass) and GTAMap ",
#     weight="bold",
# )
#
# plt.grid(True, linestyle=":", alpha=0.6)
# plt.tight_layout()
# plt.savefig("starsp_vs_gtamap.png", dpi=300)
