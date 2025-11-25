import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
import glob


base_colors = sns.color_palette("husl", n_colors=6)

# tool_order_a = ["GTAMap", "Minimap2", "HISAT2", "STAR"]
tool_order_a = ["GTAMap", "Minimap2", "HISAT2", "STAR", "BWA", "Bowtie2"]
tool_colors = dict(zip(tool_order_a, base_colors))
tool_colors["GTAMap"] = sns.color_palette("Reds", 6)[4]  # bright red tone

tool_order_rna = ["GTAMap", "Minimap2", "HISAT2", "STAR"]
tool_order_dna = ["GTAMap", "BWA", "Bowtie2"]

tool_colors_rna = {k: tool_colors[k] for k in tool_order_rna}
tool_colors_dna = {k: tool_colors[k] for k in tool_order_dna}


# --- Input TSV ---
ty = "rna"  # use dna for rna res

if ty == "rna":
    out_dir = "results_all_fair"
    os.makedirs(out_dir, exist_ok=True)

    # --- Colors ---
    edge_color = "black"

    # tsv_files = glob.glob(os.path.join("./data/", "*.tsv"))
    # tsv_files = glob.glob(os.path.join("./rna_fair/", "*.tsv"))
    # tsv_files = glob.glob(os.path.join("./rna_fair_subset/", "*.tsv"))
    # tsv_files = glob.glob(os.path.join("./rna_def_serial/", "*.tsv"))
    tsv_files = glob.glob(os.path.join("./rna_fair_serial/", "*.tsv"))

    df_main = pd.DataFrame()
    # Read and concatenate all TSVs
    for f in tsv_files:
        df = pd.read_csv(f, sep="\t")

        # Add a column for the file name (or any label)
        df["mutation_rate"] = os.path.basename(f).split("_")[2]

        # Concatenate into the main DataFrame
        df_main = pd.concat([df_main, df], ignore_index=True)

    tool_map = {
        "hisat2": "HISAT2",
        "star": "STAR",
        "minimap2": "Minimap2",
        "bwa": "BWA",
        "bowtie": "Bowtie2",
        "gtamap": "GTAMap",
    }

    df = df_main

    # Replace tool names in the 'tool' column
    df["tool"] = df["tool"].replace(tool_map)

    # only look at mutation rate == 0
    df = df[df["mutation_rate"] == "0"].copy()

    print(df_main.shape)
    print(df_main.head())

    summary = (
        df.groupby("tool")
        .agg(
            avg_wall_time_sec=("wall_time_sec", "mean"),
            total_wall_time_sec=("wall_time_sec", "sum"),
        )
        .reset_index()
    )

    print(summary)

    # --- Convert numeric columns ---
    numeric_cols = [
        "user_time",
        "system_time",
        "cpu_percent",
        "wall_time_sec",
        "max_rss_MB",
        "mutation_rate",
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
        plt.rcParams["figure.figsize"] = (6.5, 4.5)  # width x height in inches
        sns.set_theme(style="ticks", context="paper", palette="husl")
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

        plt.figure()

        ax = sns.boxplot(
            data=data,
            x=x,
            y=y,
            showfliers=False,  # Hide outliers
            width=0.6,
            hue="tool",
            order=tool_order_rna,
            hue_order=tool_order_rna,
            palette=tool_colors_rna,
            boxprops=dict(edgecolor=edge_color, linewidth=1.5),
            whiskerprops=dict(color=edge_color, linewidth=1.5),
            capprops=dict(color=edge_color, linewidth=1.5),
            medianprops=dict(color=edge_color, linewidth=1.5),
        )

        ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
        ax.xaxis.grid(False)  # optional: turn off vertical grid lines

        # Compute whiskers per x+hue
        whisker_lows, whisker_highs = [], []
        for (_, hue_val), group in data.groupby([x, "mutation_rate"]):
            low, high = get_whiskers(group[y].dropna())
            whisker_lows.append(low)
            whisker_highs.append(high)

        global_low = min(whisker_lows)
        global_high = max(whisker_highs)
        padding = 0.02 * (global_high - global_low)
        y_start = min(0, global_low - padding)
        ax.set_ylim(0, global_high + padding)

        # Styling
        # ax.set_title(title, weight="bold", fontsize=17)
        ax.set_title(title, weight="bold")
        ax.set_xlabel("")
        ax.set_ylabel(ylabel)

        # Create legend on the main axes (optional, will be removed later)
        legend = ax.legend(
            title="Gene Mutation Rate in %",
            fontsize=12,
            title_fontsize=14,
            ncol=len(data["mutation_rate"].unique()),  # horizontal layout
        )

        # --- Save legend as a separate figure ---
        fig_legend = plt.figure(
            figsize=(4, 1)
        )  # wider than tall to fit horizontal legend
        ax_leg = fig_legend.add_subplot(111)
        ax_leg.axis("off")

        # Draw the legend in the new figure with horizontal layout
        ax_leg.legend(
            handles=legend.legend_handles,
            labels=[t.get_text() for t in legend.get_texts()],
            title=legend.get_title().get_text(),
            fontsize=12,
            title_fontsize=14,
            ncol=len(data["mutation_rate"].unique()),  # <- makes it horizontal
            loc="center",
        )
        fig_legend.canvas.draw()
        fig_legend.savefig(
            os.path.join(out_dir, "legend_" + filename), dpi=300, bbox_inches="tight"
        )
        plt.close(fig_legend)

        # Remove legend from main plot
        ax.get_legend().remove()
        sns.despine()
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, filename), dpi=300)
        plt.close()

    def boxplot_metric(
        data, x, y, title, ylabel, filename, out_dir, figsize=(6.5, 4.5)
    ):
        plt.rcParams["figure.figsize"] = (6.5, 4.5)  # width x height in inches
        sns.set_theme(style="ticks", context="paper", palette="husl")
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

        plt.figure()

        ax = sns.boxplot(
            data=data,
            x=x,
            y=y,
            width=0.6,
            hue="tool",
            order=tool_order_rna,
            hue_order=tool_order_rna,
            palette=tool_colors_rna,
            boxprops=dict(edgecolor=edge_color, linewidth=1.5),
            whiskerprops=dict(color=edge_color, linewidth=1.5),
            capprops=dict(color=edge_color, linewidth=1.5),
            medianprops=dict(color=edge_color, linewidth=1.5),
        )

        ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
        ax.xaxis.grid(False)  # optional: turn off vertical grid lines

        # Styling
        # ax.set_title(title, weight="bold", fontsize=17)
        ax.set_title(title, weight="bold")
        ax.set_xlabel("")
        ax.set_ylabel(ylabel)

        # Create legend on the main axes (optional, will be removed later)
        legend = ax.legend(
            title="Gene Mutation Rate in %",
            fontsize=12,
            title_fontsize=14,
            ncol=len(data["mutation_rate"].unique()),  # horizontal layout
        )

        # --- Save legend as a separate figure ---
        fig_legend = plt.figure(
            figsize=(4, 1)
        )  # wider than tall to fit horizontal legend
        ax_leg = fig_legend.add_subplot(111)
        ax_leg.axis("off")

        # Draw the legend in the new figure with horizontal layout
        ax_leg.legend(
            handles=legend.legend_handles,
            labels=[t.get_text() for t in legend.get_texts()],
            title=legend.get_title().get_text(),
            fontsize=12,
            title_fontsize=14,
            ncol=len(data["mutation_rate"].unique()),  # <- makes it horizontal
            loc="center",
        )
        fig_legend.canvas.draw()
        fig_legend.savefig(
            os.path.join(out_dir, "legend_" + filename), dpi=300, bbox_inches="tight"
        )
        plt.close(fig_legend)

        # Remove legend from main plot
        ax.get_legend().remove()
        sns.despine()
        plt.tight_layout()
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
            # title = f"{title} - RNA Mappers (default)"
            title = f"{title} - RNA Mappers (sensitive)"
            boxplot_metric_zoom(df, "tool", col, title, ylabel, fname, out_dir)
            fname = f"{col}_box_fliers.png"
            boxplot_metric(df, "tool", col, title, ylabel, fname, out_dir)

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
else:
    out_dir = "results_dna"
    os.makedirs(out_dir, exist_ok=True)

    # --- Colors ---
    edge_color = "black"

    tsv_files = glob.glob(os.path.join("./data_dna/", "*.tsv"))

    df_main = pd.DataFrame()
    # Read and concatenate all TSVs
    for f in tsv_files:
        df = pd.read_csv(f, sep="\t")

        # Add a column for the file name (or any label)
        df["mutation_rate"] = os.path.basename(f).split("_")[2]

        # Concatenate into the main DataFrame
        df_main = pd.concat([df_main, df], ignore_index=True)

    tool_map = {
        "hisat2": "HISAT2",
        "star": "STAR",
        "minimap2": "Minimap2",
        "bwa": "BWA",
        "bowtie2": "Bowtie2",
        "gtamap": "GTAMap",
    }

    df = df_main

    # Replace tool names in the 'tool' column
    df["tool"] = df["tool"].replace(tool_map)

    # only look at mutation rate == 0
    df = df[df["mutation_rate"] == "0"].copy()

    print(df_main.shape)
    print(df_main.head())

    # --- Convert numeric columns ---
    numeric_cols = [
        "user_time",
        "system_time",
        "cpu_percent",
        "wall_time_sec",
        "max_rss_MB",
        "mutation_rate",
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
        plt.rcParams["figure.figsize"] = (6.5, 4.5)  # width x height in inches
        sns.set_theme(style="ticks", context="paper", palette="husl")
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

        plt.figure(figsize=figsize)

        ax = sns.boxplot(
            data=data,
            x=x,
            y=y,
            showfliers=False,  # Hide outliers
            width=0.6,
            hue="tool",
            order=tool_order_dna,
            hue_order=tool_order_dna,
            palette=tool_colors_dna,
            boxprops=dict(edgecolor=edge_color, linewidth=1.5),
            whiskerprops=dict(color=edge_color, linewidth=1.5),
            capprops=dict(color=edge_color, linewidth=1.5),
            medianprops=dict(color=edge_color, linewidth=1.5),
        )

        ax.set_ylabel(ylabel)

        ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
        ax.xaxis.grid(False)  # optional: turn off vertical grid lines

        # Compute whiskers per x+hue
        whisker_lows, whisker_highs = [], []
        for (_, hue_val), group in data.groupby([x, "mutation_rate"]):
            low, high = get_whiskers(group[y].dropna())
            whisker_lows.append(low)
            whisker_highs.append(high)

        global_low = min(whisker_lows)
        global_high = max(whisker_highs)
        padding = 0.02 * (global_high - global_low)
        y_start = min(0, global_low - padding)
        ax.set_ylim(0, global_high + padding)

        # Styling
        # ax.set_title(title, weight="bold", fontsize=17)
        ax.set_title(title, weight="bold")
        ax.set_xlabel("")

        sns.despine()
        plt.tight_layout()
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
            title = f"{title} - DNA Mappers"
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

    def boxplot_wall_clock_split_dna(df, out_dir, filename="wall_clock_split.png"):
        """
        Split wall clock time into two panels:
        - Left: GTAMap + Bowtie2
        - Right: BWA
        """
        import matplotlib.pyplot as plt
        import seaborn as sns

        # --- Subsets ---
        df_left = df[df["tool"].isin(["GTAMap", "Bowtie2"])]
        df_right = df[df["tool"] == "BWA"]

        fig, axes = plt.subplots(
            1,
            2,
            figsize=(6.5, 4.5),
            gridspec_kw={"width_ratios": [2, 1]},  # left wider, right narrower
        )
        fig.suptitle(
            "Wall Clock Time - DNA Mappers",
            weight="bold",
            fontsize=17,
            x=0.55,  # center horizontally (default)
        )

        sns.set_theme(style="ticks", context="paper")

        # Helper for Tukey whiskers
        def get_whiskers(series):
            q1 = series.quantile(0.25)
            q3 = series.quantile(0.75)
            iqr = q3 - q1
            low = q1 - 1.5 * iqr
            high = q3 + 1.5 * iqr
            return low, high

        def plot_box(df_subset, ax):
            sns.boxplot(
                data=df_subset,
                x="tool",
                y="wall_time_sec",
                hue="tool",
                hue_order=tool_order_dna,
                palette=tool_colors_dna,
                fliersize=0,
                showfliers=False,
                ax=ax,
                boxprops=dict(edgecolor="black", linewidth=1.5),
                whiskerprops=dict(color="black", linewidth=1.5),
                capprops=dict(color="black", linewidth=1.5),
                medianprops=dict(color="black", linewidth=1.5),
            )
            ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
            ax.xaxis.grid(False)
            ax.set_xlabel("")
            return ax

        # --- Left panel ---
        plot_box(df_left, axes[0])
        axes[0].set_ylabel("Wall Clock Time (s)")
        axes[0].set_xticklabels(df_left["tool"].unique())

        # --- Right panel ---
        plot_box(df_right, axes[1])
        axes[1].set_ylabel("")
        axes[1].set_xticklabels(df_right["tool"].unique())

        # --- Compute whisker ranges ---
        def get_axis_limits(df_subset):
            lows, highs = [], []
            for tool, group in df_subset.groupby("tool"):
                low, high = get_whiskers(group["wall_time_sec"])
                lows.append(low)
                highs.append(high)
            return min(lows), max(highs)

        low0, high0 = get_axis_limits(df_left)
        pad0 = 0.02 * (high0 - low0)
        # axes[0].set_ylim(0, high0 + pad0)
        # axes[0].set_yscale("log")

        low1, high1 = get_axis_limits(df_right)
        pad1 = 0.02 * (high1 - low1)
        axes[1].set_ylim(0, high1 + pad1)

        sns.despine()
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, filename), dpi=300, bbox_inches="tight")
        plt.close()

    summary = (
        df.groupby("tool")
        .agg(
            avg_wall_time_sec=("wall_time_sec", "mean"),
            total_wall_time_sec=("wall_time_sec", "sum"),
        )
        .reset_index()
    )

    print(summary)

    boxplot_wall_clock_split_dna(df, out_dir)
