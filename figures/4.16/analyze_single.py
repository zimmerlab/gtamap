import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import os

# sns.set_theme(style="ticks", context="paper", palette="Set2")
sns.set_theme(style="ticks", context="paper", palette="colorblind")

tool_order_a = ["GTAMap", "Minimap2", "HISAT2", "STAR", "BWA", "Bowtie2"]

sns.set_context(
    "paper",
    rc={
        "font.size": 15,  # base font size
        "axes.titlesize": 17,  # title
        "axes.labelsize": 15,  # axis labels
        "xtick.labelsize": 14,  # x tick labels
        "ytick.labelsize": 14,  # y tick labels
        "legend.fontsize": 8,  # legend
        "lines.linewidth": 1.5,  # line width
        "lines.markersize": 6,  # marker size
    },
)

plt.rcParams["figure.figsize"] = (6.5, 4.5)  # width x height in inches

# Read data into DataFrame
df = pd.read_csv(sys.argv[1], sep="\t", engine="python")
df_stages = pd.read_csv(sys.argv[2], sep="\t", engine="python")


gene_name = sys.argv[1].split("/")[-1].split(".")[0]
print(gene_name)

if not os.path.exists(gene_name):
    os.mkdir(gene_name)


# Convert timestamp column to datetime
df["timestamp"] = pd.to_datetime(df["timestamp"])

# Replace problematic values (+Inf, NaN)
df.replace(["+Inf", "NaN"], [float("inf"), float("nan")], inplace=True)

# Convert stage start times to datetime (assuming same date as df timestamps)
date_str = df["timestamp"].iloc[0].strftime("%Y-%m-%d")
stage_start_cols = [
    col for col in df_stages.columns if col.endswith("End") or col.endswith("Start")
]

stage_starts = {}
for col in stage_start_cols:
    time_str = df_stages[col].iloc[0]
    if col == "mappingTaskProducerStart":
        stage_starts["Stage-1 Start"] = pd.to_datetime(f"{date_str} {time_str}")
        # if col == "confidentWorkerStart":
        #     stage_starts["Stage-2 Start"] = pd.to_datetime(f"{date_str} {time_str}")
        # if col == "confidentWorkerEnd":
        #     stage_starts["Stage-2 End"] = pd.to_datetime(f"{date_str} {time_str}")
    if col == "outputWorkerStart":
        stage_starts["Stage-3 Start"] = pd.to_datetime(f"{date_str} {time_str}")
    if col == "outputWorkerEnd":
        stage_starts["Stage-3 End"] = pd.to_datetime(f"{date_str} {time_str}")

print(stage_starts)
# Convert ALL stage times to datetime
all_stage_times = []
for col in stage_start_cols:
    time_str = df_stages[col].iloc[0]
    all_stage_times.append(pd.to_datetime(f"{date_str} {time_str}"))

# Find the earliest time across both metrics and stages
time_zero = min(df["timestamp"].min(), min(all_stage_times))

# Normalize time - start at 0
df["time_seconds"] = (df["timestamp"] - time_zero).dt.total_seconds()

# Normalize stage times
stage_starts_normalized = {}
for stage_name, start_time in stage_starts.items():
    stage_starts_normalized[stage_name] = (start_time - time_zero).total_seconds()

# Select numeric columns to plot
metrics = [col for col in df.columns if col not in ["timestamp", "time_seconds"]]

print(f"Time zero: {time_zero}")
print(stage_starts_normalized)

# Create separate plots
for metric in metrics:
    if metric != "allocKB":
        continue
    plt.figure(figsize=(10, 5))
    plt.plot(df["time_seconds"], df[metric], marker="o")

    # Add vertical lines for each stage start
    colors = ["red", "green", "blue", "orange", "purple"]
    for idx, (stage_name, start_time_seconds) in enumerate(
        stage_starts_normalized.items()
    ):
        color = colors[idx % len(colors)]
        plt.axvline(
            x=start_time_seconds,
            color=color,
            linestyle="--",
            alpha=1,
            label=stage_name,
        )

    plt.title(f"{metric} of TTN run")
    plt.xlabel("Time (seconds)")
    plt.ylabel(metric)
    plt.legend(loc="best", fontsize="small")
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{gene_name}/{metric}.png", dpi=300)
    plt.close()
