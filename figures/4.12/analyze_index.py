import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys

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
        "legend.fontsize": 8,  # legend
        "lines.linewidth": 1.5,  # line width
        "lines.markersize": 6,  # marker size
    },
)

plt.rcParams["figure.figsize"] = (6.5, 4.5)  # width x height in inches


# --- Input TSV ---
tsv_file = sys.argv[1]
bowtie = sys.argv[2]
bwa = sys.argv[3]

out_dir = sys.argv[4] if len(sys.argv) > 3 else "plots"
os.makedirs(out_dir, exist_ok=True)

# --- Colors ---
fill_color = "#CFD5E8"
edge_color = "black"

# --- Read data ---
df_main = pd.read_csv(tsv_file, sep="\t")

df_extra1 = pd.read_csv(bwa, sep="\t")
df_extra2 = pd.read_csv(bowtie, sep="\t")

df_extra1["exit_status"] = "0"
df_extra2["exit_status"] = "0"

# --- Combine all dataframes ---
df = pd.concat([df_main, df_extra1, df_extra2], ignore_index=True)

tool_map = {
    "hisat2": "HISAT2",
    "star": "STAR",
    "minimap": "Minimap2",
    "bwa": "BWA",
    "bowtie": "Bowtie2",
    "gtamap": "GTAMap",
}

# Replace tool names in the 'tool' column
df["tool"] = df["tool"].replace(tool_map)

# --- Convert numeric columns ---
numeric_cols = [
    "user_time",
    "system_time",
    "cpu_percent",
    "wall_time_sec",
    "max_rss_MB",
]
for col in numeric_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# --- Compute total CPU time ---
df["cpu_time_total"] = df[["user_time", "system_time"]].sum(axis=1, skipna=True)

# --- Compute efficiency ratio ---
df["efficiency_ratio"] = df["cpu_time_total"] / df["wall_time_sec"]


# --- Plotting function ---
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
        medianprops=dict(color=edge_color, linewidth=2),
        hue="tool",
        palette=tool_colors,
        hue_order=tool_order_a,
    )

    # Compute whiskers for zoom
    whisker_lows, whisker_highs = [], []
    for _, group in data.groupby(x):
        low, high = get_whiskers(group[y].dropna())
        whisker_lows.append(low)
        whisker_highs.append(high)
    global_low = min(whisker_lows)
    global_high = max(whisker_highs)
    padding = 0.002 * (global_high - global_low)  # 2% padding
    y_start = min(0, global_low - padding)
    ax.set_ylim(y_start, global_high + padding)
    if y == "wall_time_sec":
        ax.set_ylim(-0.1, global_high + padding)
    if y == "efficiency_ratio":
        ax.set_ylim(-0.01, global_high + padding)
    if y == "cpu_time_total":
        ax.set_ylim(-0.01, global_high + padding)

    # Styling
    ax.set_title(title, weight="bold", pad=15)
    ax.set_xlabel("")
    ax.legend_.remove()

    ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    ax.xaxis.grid(False)  # optional: turn off vertical grid lines
    ax.set_ylabel(ylabel)
    sns.despine()
    plt.tight_layout()

    # Save plot
    plt.savefig(os.path.join(out_dir, filename), dpi=300)
    plt.close()


# --- Metrics to plot ---
metrics = {
    "user_time": ("User CPU Time", "Seconds"),
    "system_time": ("System CPU Time", "Seconds"),
    "cpu_time_total": ("Total CPU Time", "Seconds"),
    "wall_time_sec": ("Wall Clock Time", "Seconds"),
    "cpu_percent": ("CPU Percent", "%"),
    "max_rss_MB": ("Memory Usage", "MB"),
    "efficiency_ratio": (
        "CPU Efficiency Ratio",
        "Ratio (CPU Time / Wall Time)",
    ),  # Added
}

for col, (title, ylabel) in metrics.items():
    if col in df.columns:
        fname = f"{col}_box.png"
        boxplot_metric_zoom(df, "tool", col, title, ylabel, fname, out_dir)

# --- Additional interesting plots ---

# Scatter: CPU time vs Wall time
plt.figure(figsize=(7, 5))
sns.scatterplot(
    data=df, x="wall_time_sec", y="cpu_time_total", hue="tool", palette="Set2", s=80
)
max_val = max(df["wall_time_sec"].max(), df["cpu_time_total"].max())
plt.plot(
    [0, max_val],
    [0, max_val],
    "k--",
    alpha=0.7,
    linewidth=1.5,
    label="Ideal \n (CPU Time = Wall Time)",
)
plt.title("CPU Time vs Wall Clock Time", weight="bold", pad=15)
plt.xlabel("Wall Clock Time (s)")
plt.ylabel("Total CPU Time (s)")
plt.legend(title="Tool")
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "cpu_vs_wall_scatter.png"), dpi=300)
plt.close()

# Histogram: Memory usage
plt.figure(figsize=(7, 5))
sns.histplot(
    data=df,
    x="max_rss_MB",
    hue="tool",
    multiple="stack",
    palette="Set2",
    edgecolor=edge_color,
)
plt.title("Memory Usage Distribution", weight="bold", pad=15)
plt.xlabel("Memory Usage (MB)")
plt.ylabel("Count")
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "memory_usage_hist.png"), dpi=300)
plt.close()

# NEW: Print summary statistics for efficiency ratio
print("\n=== Efficiency Ratio Summary ===")
print("Interpretation:")
print("- ~1.0: Highly CPU-bound (little waiting)")
print("- ~0.1-0.3: Mostly I/O bound (waiting for disk/network)")
print("- >1.0: Multi-threaded/parallel (using multiple cores)")
print("- ~0.0: Mostly sleeping/waiting")

print(f"\nOverall efficiency ratio statistics:")
print(f"Mean: {df['efficiency_ratio'].mean():.3f}")
print(f"Median: {df['efficiency_ratio'].median():.3f}")
print(f"Min: {df['efficiency_ratio'].min():.3f}")
print(f"Max: {df['efficiency_ratio'].max():.3f}")

print(f"\nBy tool:")
for tool in df["tool"].unique():
    tool_data = df[df["tool"] == tool]["efficiency_ratio"]
    print(f"{tool}: mean={tool_data.mean():.3f}, median={tool_data.median():.3f}")

print(f"\nPlots saved in: {out_dir}")
