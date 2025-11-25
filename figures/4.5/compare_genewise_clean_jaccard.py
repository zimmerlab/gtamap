#!/usr/bin/env python3
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns


dir_gtamap = Path("/mnt/studtemp/weyrich/results_gene_wise_gtamap/")
dir_others = Path("/mnt/studtemp/weyrich/real_rna_benchmark_genes/")
output_dir = Path("genewise_comparison_results")
output_dir.mkdir(exist_ok=True)

mappers = ["Minimap2", "HISAT2", "STAR"]


def plot_mapper_agreement_heatmap(df, strand="fw"):
    metric_col = f"mapper_agreement_{strand}"
    if "mapper1" not in df.columns:
        df[["mapper1", "mapper2"]] = df["mapper_pair"].str.split("_vs_", expand=True)

    pivot_df = df.pivot_table(
        values=metric_col,
        index="mapper1",
        columns="mapper2",
        aggfunc="mean",
    )

    plt.figure(figsize=(6, 5))
    sns.heatmap(
        pivot_df,
        annot=True,
        fmt=".3f",
        cmap="viridis",
        square=True,
        vmin=0,
        vmax=1,
        cbar_kws={"label": f"Coverage of column by row"},
    )
    plt.title(f"Mapper Coverage ({strand.upper()})", weight="bold")
    plt.xlabel("Reference Mapper (covered)")
    plt.ylabel("Test Mapper (covering)")
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir / f"mapper_coverage_asymmetric_new_{strand}.png", dpi=300)


def mapper_agreement_coverage(mapper1_dict, mapper2_dict, ma1, ma2):
    """
    Compute coverage/agreement of mapper1 vs mapper2 across reads.

    For each read shared between mapper1 and mapper2:
      - Compute total bases of mapper2
      - Compute how many mapper2 bases are covered by mapper1 blocks
      - Compute coverage fraction = overlap / total_mapper2_bases

    Args:
        mapper1_dict: Dict mapping read_id -> list of alignment blocks
        mapper2_dict: Dict mapping read_id -> list of alignment blocks

    Returns:
        tuple: (avg_coverage, num_fully_covered_reads, avg_uncovered_bases)
    """
    coverages = []
    uncovered_bases_list = []
    fully_covered = 0

    for read_id, m2_aligns in mapper2_dict.items():
        if read_id not in mapper1_dict:
            continue

        m1_aligns = mapper1_dict[read_id]
        if len(m1_aligns) == 0 or len(m2_aligns) == 0:
            continue

        m2_blocks = [b for aln in m2_aligns for b in aln]
        m1_blocks = [b for aln in m1_aligns for b in aln]

        if not m2_blocks:
            continue

        total_m2_bases = sum(e - s for s, e in m2_blocks)

        covered_positions = set()
        for m_start, m_end in m2_blocks:
            for g_start, g_end in m1_blocks:
                if g_end <= m_start:
                    continue
                if g_start >= m_end:
                    break
                overlap_start = max(m_start, g_start)
                overlap_end = min(m_end, g_end)
                if overlap_start < overlap_end:
                    covered_positions.update(range(overlap_start, overlap_end))

        overlap = len(covered_positions)
        frac = overlap / total_m2_bases if total_m2_bases > 0 else 0.0
        coverages.append(frac)

        if overlap == total_m2_bases:
            fully_covered += 1
        else:
            uncovered_bases_list.append(total_m2_bases - overlap)

    avg_cov = np.mean(coverages) if coverages else 0.0
    avg_uncovered = np.mean(uncovered_bases_list) if uncovered_bases_list else 0.0

    return avg_cov, fully_covered, avg_uncovered


def merge_blocks(regions):
    """
    Merge overlapping or adjacent regions into blocks.

    Args:
        regions: List of (start, end) tuples.

    Returns:
        List of merged (start, end) blocks.
    """
    if not regions:
        return []

    merged = []
    current_start, current_end = regions[0]

    for start, end in regions[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end

    merged.append((current_start, current_end))
    return merged


def positional_overlap_accuracy(
    mapper1_intervals_dict,
    mapper2_intervals_dict,
    mapper1,
    mapper2,
    gene,
):
    """
    Compute positional agreement between two mappers (SAM files).
    For each read, checks if ANY alignment from mapper2 captures the mapper1 alignment.

    Each mapper dictionary should be structured as:
    {
        read_id: [
            { "chr": str, "blocks": [(start, end), ...] },
            ...
        ]
    }
    :param mapper1_intervals_dict: Parsed SAM data for mapper 1 (e.g., HISAT2)
    :param mapper2_intervals_dict: Parsed SAM data for mapper 2 (your mapper with multimaps)
    :return: (average accuracy, number of perfectly matching reads, average missed bases)
    """
    accuracies = []
    missed_bases_list = []
    number_of_perfect_reads = 0
    total_reads = 0

    for read_id, mapper1_aligns in mapper1_intervals_dict.items():
        if read_id not in mapper2_intervals_dict:
            continue

        mapper2_aligns = mapper2_intervals_dict[read_id]

        if len(mapper1_aligns) == 0:
            continue

        total_reads += 1

        best_overlap = 0
        total_bases_mapper1 = 0

        for aln1 in mapper1_aligns:
            aln1_total = sum(e - s for s, e in aln1)
            total_bases_mapper1 = max(total_bases_mapper1, aln1_total)

            for aln2 in mapper2_aligns:
                overlap_bases = 0
                for m1_start, m1_end in aln1:
                    for m2_start, m2_end in aln2:
                        if m2_end <= m1_start:
                            continue
                        if m2_start >= m1_end:
                            continue
                        overlap_start = max(m1_start, m2_start)
                        overlap_end = min(m1_end, m2_end)
                        if overlap_start < overlap_end:
                            overlap_bases += overlap_end - overlap_start
                best_overlap = max(best_overlap, overlap_bases)

        # if best_overlap < 50 and (mapper1 == "GTAMap" or mapper2 == "GTAMap"):
        #    print("-----")
        #    print(gene)
        #    print(read_id)
        #    print(mapper1)
        #    print(mapper1_aligns)
        #    print(mapper2_aligns)
        #    print(mapper2)
        #    print("-----")

        if total_bases_mapper1 > 0:
            acc = best_overlap / total_bases_mapper1
            accuracies.append(acc)

            if best_overlap == total_bases_mapper1:
                number_of_perfect_reads += 1
            else:
                missed_bases_list.append(total_bases_mapper1 - best_overlap)

    avg_accuracy = np.mean(accuracies) if accuracies else 0
    avg_missed_bases = np.mean(missed_bases_list) if missed_bases_list else 0

    return avg_accuracy, number_of_perfect_reads, avg_missed_bases


def convert_to_global(mapper1_data, offset):
    c_dict = dict()
    for category, read_dict in mapper1_data.items():
        if category in ["fw_mm", "rv_mm"]:
            c_dict[category] = read_dict
            continue
        converted_regions = defaultdict(list)
        for read_id, regions in read_dict.items():
            new_region_list = []
            for region in regions:
                converted_region = []
                for start, stop in region:
                    converted_region.append((start + offset - 1, stop + offset - 1))
                new_region_list.append(converted_region)
            converted_regions[read_id] = new_region_list
        c_dict[category] = converted_regions
    return c_dict


def load_gene_meta():
    # contains gene id, start, length, name, strand
    file = "/mnt/raidbio/biocluster/projekte/BA_weyrich/data/genome/genes_meta_repaired.tsv"
    df = pd.read_csv(
        file,
        sep="\t",
        names=["Id", "Start", "Length", "Strand", "GeneName", "Chr"],
        index_col=None,
    )
    return df


def parse_sam_file_filtered(sam_path, mapper, max_clip_fraction=0.01):
    """
    Parse a SAM file, ignoring reads where either mate has more than max_clip_fraction soft clipping.
    """
    samfile = pysam.AlignmentFile(sam_path, "r")
    fw_dict = defaultdict(list)
    rv_dict = defaultdict(list)

    for read in samfile:
        if not read.cigarstring:
            continue

        if mapper != "GTAMap":
            if not read.is_paired or read.mate_is_unmapped:
                continue

        softclip_bases = sum(length for op, length in read.cigartuples if op == 4)
        if (
            read.query_length > 0
            and (softclip_bases / read.query_length) > max_clip_fraction
        ):
            continue  # skip this read

        # classify forward/reverse
        if read.is_reverse:
            rv_dict[read.query_name].append(merge_blocks(read.get_blocks()))
        else:
            fw_dict[read.query_name].append(merge_blocks(read.get_blocks()))

    samfile.close()
    return {"fw": fw_dict, "rv": rv_dict}


def parse_sam_file(sam_path):

    samfile = pysam.AlignmentFile(sam_path, "r")
    fw_dict = defaultdict(list)
    rv_dict = defaultdict(list)
    rv_dict_mm = defaultdict(list)
    fw_dict_mm = defaultdict(list)

    for read in samfile:
        if read.cigarstring:
            should_ignore = False
            for mate in mates:
                if not mate.cigartuples or mate.query_length == 0:
                    continue
                softclip_bases = sum(
                    length for op, length in mate.cigartuples if op == 4
                )
                softclip_fraction = softclip_bases / mate.query_length
                if softclip_bases > 0:
                    should_ignore = True
                    break  # one mate bad → both discarded

            if should_ignore:
                continue

            if read.is_reverse and read.is_read2:
                rv_dict[read.query_name].append(merge_blocks(read.get_blocks()))
            elif read.is_reverse and read.is_read1:
                fw_dict[read.query_name].append(merge_blocks(read.get_blocks()))
            elif read.is_forward and read.is_read1:
                fw_dict[read.query_name].append(merge_blocks(read.get_blocks()))
            else:
                rv_dict[read.query_name].append(merge_blocks(read.get_blocks()))

    samfile.close()

    return {"fw": fw_dict, "rv": rv_dict}


def compare_gene(gene_id, g_path, mapper_paths, meta, i):
    """
    Compare one gene between all mappers pairwise (including GTAMap),
    computing both positional overlap accuracy and general mapper agreement.
    """
    results = {}

    # Parse GTAMap
    parsed = {"GTAMap": parse_sam_file_filtered(g_path, "GTAMap")}
    gene_start = meta.loc[meta["Id"] == gene_id, "Start"].values[0]

    # Parse all other mappers and convert coordinates
    for mapper, path in mapper_paths.items():
        if not path.exists():
            continue
        parsed[mapper] = parse_sam_file_filtered(path, mapper)
        parsed[mapper] = convert_to_global(parsed[mapper], gene_start)

    # Get all mappers that were actually parsed
    available_mappers = [
        m for m in parsed.keys() if len(parsed[m]["fw"]) > 0 or len(parsed[m]["rv"]) > 0
    ]

    # Pairwise comparisons
    for m1, m2 in combinations(available_mappers, 2):
        # Positional overlap accuracy (like before)
        fw_acc, _, _ = positional_overlap_accuracy(
            parsed[m1]["fw"], parsed[m2]["fw"], m1, m2, gene_id
        )
        rv_acc, _, _ = positional_overlap_accuracy(
            parsed[m1]["rv"], parsed[m2]["rv"], m1, m2, gene_id
        )

        pair_key = f"{m1}_vs_{m2}"
        results[pair_key] = {
            "mapper1": m1,
            "mapper2": m2,
            "fw_acc": fw_acc,
            "rv_acc": rv_acc,
            "QNAMES1": len(parsed[m1]["fw"]),
            "QNAMES2": len(parsed[m2]["fw"]),
        }

        # General mapper agreement (coverage-based)
        fw_cov, fw_full, fw_uncov = mapper_agreement_coverage(
            parsed[m1]["fw"], parsed[m2]["fw"], m1, m2
        )
        rv_cov, rv_full, rv_uncov = mapper_agreement_coverage(
            parsed[m1]["rv"], parsed[m2]["rv"], m1, m2
        )
        results[pair_key].update(
            {
                "mapper_agreement_fw": fw_cov,
                "mapper_agreement_rv": rv_cov,
                "mapper_fully_covered_fw": fw_full,
                "mapper_fully_covered_rv": rv_full,
                "mapper_uncovered_bases_fw": fw_uncov,
                "mapper_uncovered_bases_rv": rv_uncov,
            }
        )

    print(f"[{i}] {gene_id} — compared {len(results)} pairs")

    return gene_id, results


def main():
    """Main analysis pipeline for comparing mapper alignments."""
    gtamap_files = list(dir_gtamap.glob("gtamap_*.sam"))
    print(f"Found {len(gtamap_files)} GTAMap SAM files")
    meta = load_gene_meta()

    # Process genes in parallel
    data = []
    with ProcessPoolExecutor(max_workers=30) as exe:
        futures = []
        for i, f in enumerate(gtamap_files):
            gene_id = f.stem.split("_")[1]
            mapper_paths = {
                "Minimap2": dir_others / f"minimap2_{gene_id}.sam",
                "HISAT2": dir_others / f"hisat2_{gene_id}.sam",
                "STAR": dir_others / f"star_{gene_id}.sam",
            }
            futures.append(exe.submit(compare_gene, gene_id, f, mapper_paths, meta, i))
            if i == 100:  # Remove this to run all genes
                break

        for fut in futures:
            gene_id, res = fut.result()
            if res is None:
                continue
            for pair, vals in res.items():
                data.append({"gene": gene_id, "mapper_pair": pair, **vals})

    # Save results
    df = pd.DataFrame(data)
    df.to_csv(output_dir / "gene_wise_all_vs_all_metrics.csv", index=False)
    print(f"Saved metrics for {len(df)} gene-pair comparisons.")

    _plot_boxplots(df)
    _plot_heatmaps(df)
    _plot_qname_distribution(df)
    plot_mapper_agreement_heatmap(df)
    plot_mapper_agreement_heatmap(df, "rv")


def _plot_boxplots(df):
    """Generate boxplots for key metrics."""
    sample = df.sample(min(100, len(df)), random_state=42)

    for metric in ["fw_acc", "rv_acc"]:
        if metric not in df.columns:
            continue

        plt.figure()
        ax = sns.boxplot(data=sample, x="mapper_pair", y=metric, palette="husl")
        sns.despine(ax=ax)
        plt.xticks(rotation=30, ha="right")
        plt.title(f"Distribution of {metric} across 100 sampled genes")
        plt.tight_layout()
        plt.savefig(output_dir / f"{metric}_boxplot_all_vs_all.png", dpi=300)
        plt.close()


def _plot_heatmaps(df):
    """Generate heatmaps for all comparison metrics."""
    if "mapper_pair" not in df.columns:
        print("⚠️  'mapper_pair' column missing. Skipping heatmaps.")
        return

    # Split mapper pairs
    df[["mapper1", "mapper2"]] = df["mapper_pair"].str.split("_vs_", expand=True)

    metrics = ["fw_acc", "rv_acc"]

    # Add coverage metrics if available
    coverage_metrics = [
        "gtamap_cov_mapper2_fw",
        "gtamap_cov_mapper2_rv",
        "mapper2_cov_gtamap_fw",
        "mapper2_cov_gtamap_rv",
        "gtamap_fully_covered_fw",
        "gtamap_fully_covered_rv",
        "mapper2_fully_covered_fw",
        "mapper2_fully_covered_rv",
        "gtamap_uncovered_bases_fw",
        "gtamap_uncovered_bases_rv",
        "mapper2_uncovered_bases_fw",
        "mapper2_uncovered_bases_rv",
    ]
    available_coverage = [m for m in coverage_metrics if m in df.columns]
    metrics.extend(available_coverage)

    for metric in metrics:
        if metric not in df.columns:
            continue

        pivot_df = _create_symmetric_pivot(df, metric)

        plt.figure(figsize=(6, 5))
        sns.heatmap(
            pivot_df,
            annot=True,
            fmt=".3f",
            cmap="viridis",
            square=True,
            cbar_kws={"label": metric},
        )
        plt.title("Average Best Overlap per QNAME Across All Genes")
        sns.despine()
        plt.tight_layout()
        plt.xlabel("")
        plt.ylabel("")
        plt.savefig(output_dir / f"{metric}_heatmap.png", dpi=300)
        plt.close()


def _plot_qname_distribution(df):
    tool_order_a = ["GTAMap", "Minimap2", "HISAT2", "STAR", "BWA", "Bowtie2"]
    base_colors = sns.color_palette("husl", n_colors=6)

    tool_colors = dict(zip(tool_order_a, base_colors))
    tool_colors["GTAMap"] = sns.color_palette("Reds", 6)[4]  # bright red tone

    if "mapper1" not in df.columns:
        df[["mapper1", "mapper2"]] = df["mapper_pair"].str.split("_vs_", expand=True)

    qname_data = []
    for _, row in df.iterrows():
        if pd.notna(row.get("QNAMES1")):
            qname_data.append(
                {
                    "mapper": row["mapper1"],
                    "gene": row["gene"],
                    "QNAME_count": row["QNAMES1"],
                }
            )
        if pd.notna(row.get("QNAMES2")):
            qname_data.append(
                {
                    "mapper": row["mapper2"],
                    "gene": row["gene"],
                    "QNAME_count": row["QNAMES2"],
                }
            )

    qnames_df = pd.DataFrame(qname_data).dropna()

    plt.figure(figsize=(6.5, 4))
    ax = sns.boxplot(
        data=qnames_df,
        x="mapper",
        y="QNAME_count",
        palette=tool_colors,
        hue_order=tool_order_a,
        showfliers=False,
    )

    ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    ax.xaxis.grid(False)
    plt.title("Distribution of Mapped Reads per Gene", weight="bold")
    ax.set_ylabel("Mapped Reads [QNAME]")
    ax.set_xlabel("")
    ax.set_ylim(0, 2500)
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir / "qnames_per_mapper_boxplot.png", dpi=300)
    plt.close()


def _plot_gtamap_coverage(df):

    if "mapper_pair" not in df.columns:
        return

    if "mapper1" not in df.columns:
        df[["mapper1", "mapper2"]] = df["mapper_pair"].str.split("_vs_", expand=True)

    coverage_data = []
    for _, row in df.iterrows():
        if row["mapper1"] == "GTAMap":
            if "gtamap_cov_mapper2_fw" in df.columns:
                coverage_data.append(
                    {
                        "other_mapper": row["mapper2"],
                        "gene": row["gene"],
                        "strand": "forward",
                        "coverage_pct": row["gtamap_cov_mapper2_fw"] * 100,
                        "fully_covered": row.get("gtamap_fully_covered_fw", 0),
                        "uncovered_bases": row.get("gtamap_uncovered_bases_fw", 0),
                    }
                )
            if "gtamap_cov_mapper2_rv" in df.columns:
                coverage_data.append(
                    {
                        "other_mapper": row["mapper2"],
                        "gene": row["gene"],
                        "strand": "reverse",
                        "coverage_pct": row["gtamap_cov_mapper2_rv"] * 100,
                        "fully_covered": row.get("gtamap_fully_covered_rv", 0),
                        "uncovered_bases": row.get("gtamap_uncovered_bases_rv", 0),
                    }
                )
        elif row["mapper2"] == "GTAMap":
            if "mapper2_cov_gtamap_fw" in df.columns:
                coverage_data.append(
                    {
                        "other_mapper": row["mapper1"],
                        "gene": row["gene"],
                        "strand": "forward",
                        "coverage_pct": row["mapper2_cov_gtamap_fw"] * 100,
                        "fully_covered": row.get("mapper2_fully_covered_fw", 0),
                        "uncovered_bases": row.get("mapper2_uncovered_bases_fw", 0),
                    }
                )
            if "mapper2_cov_gtamap_rv" in df.columns:
                coverage_data.append(
                    {
                        "other_mapper": row["mapper1"],
                        "gene": row["gene"],
                        "strand": "reverse",
                        "coverage_pct": row["mapper2_cov_gtamap_rv"] * 100,
                        "fully_covered": row.get("mapper2_fully_covered_rv", 0),
                        "uncovered_bases": row.get("mapper2_uncovered_bases_rv", 0),
                    }
                )

    if not coverage_data:
        print(
            "⚠️  No GTAMap coverage data found. Make sure compare_gene returns coverage metrics."
        )
        return

    cov_df = pd.DataFrame(coverage_data)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for idx, strand in enumerate(["forward", "reverse"]):
        strand_data = cov_df[cov_df["strand"] == strand]
        ax = axes[idx]
        sns.boxplot(
            data=strand_data,
            x="other_mapper",
            y="coverage_pct",
            palette="husl",
            showfliers=False,
            ax=ax,
        )
        sns.stripplot(
            data=strand_data,
            x="other_mapper",
            y="coverage_pct",
            color="black",
            alpha=0.3,
            size=2,
            ax=ax,
        )
        ax.axhline(y=100, color="red", linestyle="--", alpha=0.5, label="100% coverage")
        ax.set_title(f"GTAMap Coverage - {strand.capitalize()} Strand")
        ax.set_ylabel("Coverage (%)")
        ax.set_xlabel("Other Mapper")
        ax.set_ylim(0, 105)
        ax.legend()
        sns.despine(ax=ax)

    plt.tight_layout()
    plt.savefig(output_dir / "gtamap_coverage_by_strand.png", dpi=300)
    plt.close()

    plt.figure(figsize=(8, 5))
    ax = sns.boxplot(
        data=cov_df,
        x="other_mapper",
        y="coverage_pct",
        hue="strand",
        palette="husl",
        showfliers=False,
    )
    plt.axhline(y=100, color="red", linestyle="--", alpha=0.5, label="100% coverage")
    plt.title("GTAMap Coverage of Other Mappers' Aligned Regions")
    plt.ylabel("Coverage (%)")
    plt.xlabel("Other Mapper")
    plt.ylim(0, 105)
    plt.legend(title="Strand")
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir / "gtamap_coverage_percentage.png", dpi=300)
    plt.close()

    full_coverage_stats = (
        cov_df.groupby(["other_mapper", "strand"])
        .apply(
            lambda x: (x["coverage_pct"] == 100).sum() / len(x) * 100,
            include_groups=False,
        )
        .reset_index(name="pct_genes_fully_covered")
    )

    plt.figure(figsize=(8, 5))
    ax = sns.barplot(
        data=full_coverage_stats,
        x="other_mapper",
        y="pct_genes_fully_covered",
        hue="strand",
        palette="husl",
    )
    plt.axhline(y=100, color="red", linestyle="--", alpha=0.5, label="100%")
    plt.title("Percentage of Genes with 100% GTAMap Coverage")
    plt.ylabel("% of Genes Fully Covered")
    plt.xlabel("Other Mapper")
    plt.ylim(0, 105)
    plt.legend(title="Strand")
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir / "gtamap_full_coverage_percentage.png", dpi=300)
    plt.close()

    avg_coverage = (
        cov_df.groupby(["other_mapper", "strand"])["coverage_pct"]
        .agg(["mean", "std", "count"])
        .reset_index()
    )

    print(avg_coverage.to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    for idx, strand in enumerate(["forward", "reverse"]):
        ax = axes[idx]
        strand_data = cov_df[cov_df["strand"] == strand]

        for mapper in strand_data["other_mapper"].unique():
            mapper_data = strand_data[strand_data["other_mapper"] == mapper][
                "coverage_pct"
            ]
            ax.hist(mapper_data, bins=20, alpha=0.6, label=mapper, edgecolor="black")

        ax.axvline(x=100, color="red", linestyle="--", alpha=0.5, label="100% coverage")
        ax.set_xlabel("Coverage (%)")
        ax.set_ylabel("Number of Genes")
        ax.set_title(f"{strand.capitalize()} Strand")
        ax.legend()
        sns.despine(ax=ax)

    fig.suptitle("Distribution of GTAMap Coverage Across Genes", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / "gtamap_coverage_histogram.png", dpi=300)
    plt.close()
    print("✅ Saved: gtamap_coverage_histogram.png")


def _plot_qname_delta_heatmap(df):

    if "QNAME_delta" not in df.columns:
        if "QNAMES1" in df.columns and "QNAMES2" in df.columns:
            df["QNAME_delta"] = df["QNAMES1"] - df["QNAMES2"]
        else:
            return

    if "mapper1" not in df.columns:
        df[["mapper1", "mapper2"]] = df["mapper_pair"].str.split("_vs_", expand=True)

    qname_pivot = _create_symmetric_pivot(df, "QNAME_delta", negate_mirror=True)

    plt.figure(figsize=(6, 5))
    sns.heatmap(
        qname_pivot,
        annot=True,
        fmt=".0f",
        cmap="coolwarm",
        center=0,
        square=True,
        cbar_kws={"label": "Δ QNAME count (mapper1 - mapper2)"},
    )
    plt.title("Average difference in mapped read count between mappers")
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir / "qname_delta_heatmap.png", dpi=300)
    plt.close()
    print("✅ Saved: qname_delta_heatmap.png")


def _create_symmetric_pivot(df, metric, negate_mirror=False):
    """
    Create a symmetric pivot table from asymmetric pairwise data.

    Args:
        df: DataFrame with mapper1, mapper2, and metric columns
        metric: Column name to pivot
        negate_mirror: If True, negate values when mirroring (for delta metrics)

    Returns:
        Symmetric pivot DataFrame
    """
    pivot_df = (
        df.groupby(["mapper1", "mapper2"])[metric].mean().unstack(fill_value=np.nan)
    )

    all_mappers = sorted(set(pivot_df.index).union(set(pivot_df.columns)))
    pivot_df = pivot_df.reindex(index=all_mappers, columns=all_mappers)

    for m1 in all_mappers:
        for m2 in all_mappers:
            val_12 = pivot_df.loc[m1, m2]
            val_21 = pivot_df.loc[m2, m1]

            if pd.isna(val_12) and not pd.isna(val_21):
                pivot_df.loc[m1, m2] = -val_21 if negate_mirror else val_21
            elif pd.isna(val_21) and not pd.isna(val_12):
                pivot_df.loc[m2, m1] = -val_12 if negate_mirror else val_12

    return pivot_df


def _gtamap_coverage_against_mapper(gtamap_dict, mapper_dict):
    """
    Check how well GTAMap covers another mapper's alignments.

    For each read shared between GTAMap and the mapper:
      - Compute total mapper bases
      - Compute how many mapper bases overlap GTAMap's blocks
      - Compute coverage fraction = overlap / total_mapper_bases

    Args:
        gtamap_dict: Dict mapping read_id -> list of alignment blocks
        mapper_dict: Dict mapping read_id -> list of alignment blocks

    Returns:
        tuple: (avg_coverage, num_fully_covered_reads, avg_uncovered_bases)
    """
    coverages = []
    uncovered_bases_list = []
    fully_covered = 0

    for read_id, mapper_aligns in mapper_dict.items():
        if read_id not in gtamap_dict:
            continue

        gtamap_aligns = gtamap_dict[read_id]

        mapper_blocks = [b for aln in mapper_aligns for b in aln]
        gtamap_blocks = [b for aln in gtamap_aligns for b in aln]

        if not mapper_blocks:
            continue

        total_mapper_bases = sum(e - s for s, e in mapper_blocks)
        overlap = 0

        for m_start, m_end in mapper_blocks:
            for g_start, g_end in gtamap_blocks:
                if g_end <= m_start:
                    continue
                if g_start >= m_end:
                    break
                overlap_start = max(m_start, g_start)
                overlap_end = min(m_end, g_end)
                if overlap_start < overlap_end:
                    overlap += overlap_end - overlap_start

        frac = overlap / total_mapper_bases if total_mapper_bases > 0 else 0.0
        coverages.append(frac)

        if overlap == total_mapper_bases:
            fully_covered += 1
        else:
            uncovered_bases_list.append(total_mapper_bases - overlap)

    avg_cov = np.mean(coverages) if coverages else 0.0
    avg_uncovered = np.mean(uncovered_bases_list) if uncovered_bases_list else 0.0

    return avg_cov, fully_covered, avg_uncovered


def gtamap_coverage_against_mapper(gtamap_dict, mapper_dict):
    """
    Check how well GTAMap covers another mapper's alignments.

    For each read shared between GTAMap and the mapper:
      - Compute total mapper bases
      - Compute how many mapper bases overlap GTAMap's blocks (without double-counting)
      - Compute coverage fraction = overlap / total_mapper_bases

    Args:
        gtamap_dict: Dict mapping read_id -> list of alignment blocks
        mapper_dict: Dict mapping read_id -> list of alignment blocks

    Returns:
        tuple: (avg_coverage, num_fully_covered_reads, avg_uncovered_bases)
    """
    coverages = []
    uncovered_bases_list = []
    fully_covered = 0

    for read_id, mapper_aligns in mapper_dict.items():
        if read_id not in gtamap_dict:
            continue

        gtamap_aligns = gtamap_dict[read_id]
        if len(gtamap_aligns) == 0:
            continue

        mapper_blocks = [b for aln in mapper_aligns for b in aln]
        gtamap_blocks = [b for aln in gtamap_aligns for b in aln]

        if not mapper_blocks:
            continue

        total_mapper_bases = sum(e - s for s, e in mapper_blocks)

        covered_positions = set()

        for m_start, m_end in mapper_blocks:
            for g_start, g_end in gtamap_blocks:
                if g_end <= m_start:
                    continue
                if g_start >= m_end:
                    break
                overlap_start = max(m_start, g_start)
                overlap_end = min(m_end, g_end)
                if overlap_start < overlap_end:
                    covered_positions.update(range(overlap_start, overlap_end))

        overlap = len(covered_positions)
        frac = overlap / total_mapper_bases if total_mapper_bases > 0 else 0.0
        coverages.append(frac)

        if overlap == total_mapper_bases:
            fully_covered += 1
        else:
            uncovered_bases_list.append(total_mapper_bases - overlap)

    avg_cov = np.mean(coverages) if coverages else 0.0
    avg_uncovered = np.mean(uncovered_bases_list) if uncovered_bases_list else 0.0

    return avg_cov, fully_covered, avg_uncovered


def jaccard_similarity_mappers(
    mapper1_intervals_dict, mapper2_intervals_dict, verbose=False
):
    """
    Compute Jaccard similarity between two mappers at the base-pair level.

    Jaccard = |Intersection| / |Union|

    For each read, we compute:
    - Intersection: bases covered by BOTH mappers (on the same chromosome)
    - Union: bases covered by EITHER mapper (on the same chromosome)

    Each mapper dictionary should be structured as:
        { read_id: [ { "chr": str, "blocks": [(start, end), ...] }, ... ] }

    :param mapper1_intervals_dict: Parsed SAM data for mapper 1
    :param mapper2_intervals_dict: Parsed SAM data for mapper 2
    :param verbose: If True, print detailed information for each read
    :return: (average_jaccard, perfect_matches, reads_on_diff_chr, total_reads)
    """
    jaccard_scores = []
    perfect_matches = 0
    reads_on_different_chromosomes = 0
    total_reads = 0

    all_read_ids = set(mapper1_intervals_dict.keys()) | set(
        mapper2_intervals_dict.keys()
    )

    for read_id in all_read_ids:
        total_reads += 1

        mapper1_aligns = mapper1_intervals_dict.get(read_id, [])
        mapper2_aligns = mapper2_intervals_dict.get(read_id, [])

        if not mapper1_aligns or not mapper2_aligns:
            jaccard_scores.append(0.0)
            if verbose:
                print(
                    f"Read {read_id}: Only in {'mapper1' if mapper1_aligns else 'mapper2'}"
                )
            continue

        best_jaccard = 0.0
        best_intersection = 0
        best_union = 0
        best_chr_match = False

        for aln1 in mapper1_aligns:
            chr1 = aln1["chr"]
            blocks1 = aln1["blocks"]

            for aln2 in mapper2_aligns:
                chr2 = aln2["chr"]
                blocks2 = aln2["blocks"]

                if chr1 != chr2:
                    continue

                best_chr_match = True

                intersection = calculate_blocks_intersection(blocks1, blocks2)
                union = calculate_blocks_union(blocks1, blocks2)

                if union > 0:
                    jaccard = intersection / union
                    if jaccard > best_jaccard:
                        best_jaccard = jaccard
                        best_intersection = intersection
                        best_union = union

        if not best_chr_match:
            reads_on_different_chromosomes += 1
            jaccard_scores.append(0.0)
            if verbose:
                chrs1 = {aln["chr"] for aln in mapper1_aligns}
                chrs2 = {aln["chr"] for aln in mapper2_aligns}
                print(f"Read {read_id}: Mapped to different chromosomes")
                print(f"  Mapper1: {chrs1}, Mapper2: {chrs2}")
            continue

        jaccard_scores.append(best_jaccard)

        if best_jaccard >= 0.999:
            perfect_matches += 1

        if verbose and best_jaccard < 1.0:
            print(f"Read {read_id}: Jaccard = {best_jaccard:.4f}")
            print(f"  Intersection: {best_intersection} bp, Union: {best_union} bp")

    avg_jaccard = sum(jaccard_scores) / len(jaccard_scores) if jaccard_scores else 0.0

    return {
        "avg_jaccard": avg_jaccard,
        "perfect_matches": perfect_matches,
        "reads_on_diff_chr": reads_on_different_chromosomes,
        "total_reads": total_reads,
        "jaccard_scores": jaccard_scores,
    }


def calculate_blocks_intersection(blocks1, blocks2):
    """
    Calculate the total number of overlapping bases between two sets of blocks.

    :param blocks1: List of (start, end) tuples
    :param blocks2: List of (start, end) tuples
    :return: Total overlapping bases
    """
    total_overlap = 0

    for start1, end1 in blocks1:
        for start2, end2 in blocks2:
            overlap_start = max(start1, start2)
            overlap_end = min(end1, end2)

            if overlap_start < overlap_end:
                total_overlap += overlap_end - overlap_start

    return total_overlap


def calculate_blocks_union(blocks1, blocks2):
    """
    Calculate the total number of bases covered by either set of blocks.

    Union = Total_bases_in_blocks1 + Total_bases_in_blocks2 - Intersection

    :param blocks1: List of (start, end) tuples
    :param blocks2: List of (start, end) tuples
    :return: Total bases in union
    """
    total_bases1 = sum(end - start for start, end in blocks1)
    total_bases2 = sum(end - start for start, end in blocks2)
    intersection = calculate_blocks_intersection(blocks1, blocks2)

    return total_bases1 + total_bases2 - intersection


def calculate_blocks_union_precise(blocks1, blocks2):
    """
    Alternative implementation that merges blocks before calculating union.
    More accurate but slightly slower.

    :param blocks1: List of (start, end) tuples
    :param blocks2: List of (start, end) tuples
    :return: Total bases in union
    """
    all_blocks = sorted(blocks1 + blocks2)

    if not all_blocks:
        return 0

    merged = []
    current_start, current_end = all_blocks[0]

    for start, end in all_blocks[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end

    merged.append((current_start, current_end))

    return sum(end - start for start, end in merged)


def compare_all_vs_all_jaccard(samfiles, output_dir):
    """
    Perform all-vs-all Jaccard similarity comparisons between mappers.
    """
    import os
    import itertools
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    mapper_data = {}
    mapper_names = []
    for sam_path in samfiles:
        mapper_name = os.path.basename(sam_path).split("_")[0]
        mapper_names.append(mapper_name)
        print(f"Parsing {mapper_name}...")
        mapper_data[mapper_name] = parse_global_sam_file(sam_path)

    jaccard_fw = pd.DataFrame(index=mapper_names, columns=mapper_names, dtype=float)
    jaccard_rv = pd.DataFrame(index=mapper_names, columns=mapper_names, dtype=float)

    perfect_matches_fw = pd.DataFrame(
        index=mapper_names, columns=mapper_names, dtype=int
    )
    perfect_matches_rv = pd.DataFrame(
        index=mapper_names, columns=mapper_names, dtype=int
    )

    for m1, m2 in itertools.product(mapper_names, mapper_names):
        print(f"Comparing {m1} vs {m2}")

        result_fw = jaccard_similarity_mappers(
            mapper_data[m1]["fw"], mapper_data[m2]["fw"]
        )
        jaccard_fw.loc[m1, m2] = result_fw["avg_jaccard"]
        perfect_matches_fw.loc[m1, m2] = result_fw["perfect_matches"]

        result_rv = jaccard_similarity_mappers(
            mapper_data[m1]["rv"], mapper_data[m2]["rv"]
        )
        jaccard_rv.loc[m1, m2] = result_rv["avg_jaccard"]
        perfect_matches_rv.loc[m1, m2] = result_rv["perfect_matches"]

    os.makedirs(output_dir, exist_ok=True)
    jaccard_fw.to_csv(os.path.join(output_dir, "jaccard_fw.csv"))
    jaccard_rv.to_csv(os.path.join(output_dir, "jaccard_rv.csv"))

    def plot_heatmap(matrix, title, fname):
        plt.figure(figsize=(8, 6))
        sns.heatmap(
            matrix.astype(float),
            annot=True,
            cmap="viridis",
            fmt=".3f",
            cbar=True,
            vmin=0,
            vmax=1,
        )
        plt.title(title)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, fname))
        plt.close()

    plot_heatmap(jaccard_fw, "Jaccard Similarity (FW)", "jaccard_fw_heatmap.png")
    plot_heatmap(jaccard_rv, "Jaccard Similarity (RV)", "jaccard_rv_heatmap.png")

    print(f"\nResults saved to {output_dir}/")
    return jaccard_fw, jaccard_rv


def jaccard_similarity_reads(mapper1_dict, mapper2_dict, verbose=False):
    """
    Compute Jaccard similarity between two mappers at the base-pair level.

    Jaccard = |Intersection| / |Union|

    For each read, we compute:
    - Intersection: bases covered by BOTH mappers
    - Union: bases covered by EITHER mapper

    Each mapper dictionary should be structured as:
        { read_id: [[(start, end), ...], [(start, end), ...]] }
        i.e., list of alignments, each with list of blocks

    :param mapper1_dict: Parsed alignment data for mapper 1
    :param mapper2_dict: Parsed alignment data for mapper 2
    :param verbose: If True, print detailed information for each read
    :return: dict with metrics
    """
    jaccard_scores = []
    perfect_matches = 0
    total_reads = 0

    all_read_ids = set(mapper1_dict.keys()) | set(mapper2_dict.keys())

    for read_id in all_read_ids:
        total_reads += 1

        mapper1_aligns = mapper1_dict.get(read_id, [])
        mapper2_aligns = mapper2_dict.get(read_id, [])

        if not mapper1_aligns or not mapper2_aligns:
            jaccard_scores.append(0.0)
            if verbose:
                print(
                    f"Read {read_id}: Only in {'mapper1' if mapper1_aligns else 'mapper2'}"
                )
            continue

        best_jaccard = 0.0
        best_intersection = 0
        best_union = 0

        for blocks1 in mapper1_aligns:
            for blocks2 in mapper2_aligns:
                intersection = calculate_blocks_intersection(blocks1, blocks2)
                union = calculate_blocks_union(blocks1, blocks2)

                if union > 0:
                    jaccard = intersection / union
                    if jaccard > best_jaccard:
                        best_jaccard = jaccard
                        best_intersection = intersection
                        best_union = union

        jaccard_scores.append(best_jaccard)

        if best_jaccard >= 0.999:
            perfect_matches += 1

        if verbose and best_jaccard < 1.0:
            print(f"Read {read_id}: Jaccard = {best_jaccard:.4f}")
            print(f"  Intersection: {best_intersection} bp, Union: {best_union} bp")

    avg_jaccard = sum(jaccard_scores) / len(jaccard_scores) if jaccard_scores else 0.0

    return {
        "avg_jaccard": avg_jaccard,
        "perfect_matches": perfect_matches,
        "total_reads": total_reads,
        "jaccard_scores": jaccard_scores,
    }


def calculate_blocks_intersection(blocks1, blocks2):
    positions1 = set()
    for start, end in blocks1:
        positions1.update(range(start, end))

    positions2 = set()
    for start, end in blocks2:
        positions2.update(range(start, end))

    return len(positions1 & positions2)


def calculate_blocks_union(blocks1, blocks2):
    """
    calculate the total number of bases covered by either set of blocks.
    uses a set-based approach to avoid double-counting.

    :param blocks1: List of (start, end) tuples
    :param blocks2: List of (start, end) tuples
    :return: Total bases in union
    """
    positions1 = set()
    for start, end in blocks1:
        positions1.update(range(start, end))

    positions2 = set()
    for start, end in blocks2:
        positions2.update(range(start, end))

    return len(positions1 | positions2)


def calculate_blocks_union_fast(blocks1, blocks2):
    """
    Fast approximation using the formula: Union = Total1 + Total2 - Intersection
    More efficient for large alignments but may be slightly less accurate.

    :param blocks1: List of (start, end) tuples
    :param blocks2: List of (start, end) tuples
    :return: Total bases in union
    """
    total_bases1 = sum(end - start for start, end in blocks1)
    total_bases2 = sum(end - start for start, end in blocks2)
    intersection = calculate_blocks_intersection(blocks1, blocks2)

    return total_bases1 + total_bases2 - intersection


def compare_gene_with_jaccard(gene_id, g_path, mapper_paths, meta, i):
    results = {}

    parsed = {"GTAMap": parse_sam_file_filtered(g_path, "GTAMap")}
    gene_start = meta.loc[meta["Id"] == gene_id, "Start"].values[0]

    for mapper, path in mapper_paths.items():
        if not path.exists():
            continue
        parsed[mapper] = parse_sam_file_filtered(path, mapper)
        parsed[mapper] = convert_to_global(parsed[mapper], gene_start)

    available_mappers = [
        m for m in parsed.keys() if len(parsed[m]["fw"]) > 0 or len(parsed[m]["rv"]) > 0
    ]

    for m1, m2 in combinations(available_mappers, 2):
        fw_cov_12, fw_full_12, fw_uncov_12 = mapper_agreement_coverage(
            parsed[m1]["fw"], parsed[m2]["fw"], m1, m2
        )
        rv_cov_12, rv_full_12, rv_uncov_12 = mapper_agreement_coverage(
            parsed[m1]["rv"], parsed[m2]["rv"], m1, m2
        )

        # direction 2: m2 vs m1 (REVERSED!)
        fw_cov_21, fw_full_21, fw_uncov_21 = mapper_agreement_coverage(
            parsed[m2]["fw"], parsed[m1]["fw"], m2, m1  # ← Swapped arguments
        )
        rv_cov_21, rv_full_21, rv_uncov_21 = mapper_agreement_coverage(
            parsed[m2]["rv"], parsed[m1]["rv"], m2, m1  # ← Swapped arguments
        )

        # jaccard is symmetric - compute once
        fw_jaccard = jaccard_similarity_reads(parsed[m1]["fw"], parsed[m2]["fw"])
        rv_jaccard = jaccard_similarity_reads(parsed[m1]["rv"], parsed[m2]["rv"])

        # positional overlap (check if this is also asymmetric)
        fw_acc_12, _, _ = positional_overlap_accuracy(
            parsed[m1]["fw"], parsed[m2]["fw"], m1, m2, gene_id
        )
        rv_acc_12, _, _ = positional_overlap_accuracy(
            parsed[m1]["rv"], parsed[m2]["rv"], m1, m2, gene_id
        )

        fw_acc_21, _, _ = positional_overlap_accuracy(
            parsed[m2]["fw"], parsed[m1]["fw"], m2, m1, gene_id
        )
        rv_acc_21, _, _ = positional_overlap_accuracy(
            parsed[m2]["rv"], parsed[m1]["rv"], m2, m1, gene_id
        )

        results[f"{m1}_vs_{m2}"] = {
            "mapper1": m1,
            "mapper2": m2,
            "fw_acc": fw_acc_12,
            "rv_acc": rv_acc_12,
            "QNAMES1": len(parsed[m1]["fw"]),
            "QNAMES2": len(parsed[m2]["fw"]),
            "mapper_agreement_fw": fw_cov_12,
            "mapper_agreement_rv": rv_cov_12,
            "mapper_fully_covered_fw": fw_full_12,
            "mapper_fully_covered_rv": rv_full_12,
            "mapper_uncovered_bases_fw": fw_uncov_12,
            "mapper_uncovered_bases_rv": rv_uncov_12,
            "jaccard_fw": fw_jaccard["avg_jaccard"],
            "jaccard_rv": rv_jaccard["avg_jaccard"],
            "jaccard_perfect_fw": fw_jaccard["perfect_matches"],
            "jaccard_perfect_rv": rv_jaccard["perfect_matches"],
        }

        results[f"{m2}_vs_{m1}"] = {
            "mapper1": m2,
            "mapper2": m1,
            "fw_acc": fw_acc_21,
            "rv_acc": rv_acc_21,
            "QNAMES1": len(parsed[m2]["fw"]),
            "QNAMES2": len(parsed[m1]["fw"]),
            "mapper_agreement_fw": fw_cov_21,
            "mapper_agreement_rv": rv_cov_21,
            "mapper_fully_covered_fw": fw_full_21,
            "mapper_fully_covered_rv": rv_full_21,
            "mapper_uncovered_bases_fw": fw_uncov_21,
            "mapper_uncovered_bases_rv": rv_uncov_21,
            "jaccard_fw": fw_jaccard["avg_jaccard"],  # Same as above
            "jaccard_rv": rv_jaccard["avg_jaccard"],  # Same as above
            "jaccard_perfect_fw": fw_jaccard["perfect_matches"],
            "jaccard_perfect_rv": rv_jaccard["perfect_matches"],
        }

    print(f"[{i}] {gene_id} — compared {len(results)} pairs")

    return gene_id, results


def plot_jaccard_heatmap(df, strand="fw", output_dir=Path(".")):
    metric_col = f"jaccard_{strand}"
    if metric_col not in df.columns:
        return

    if "mapper1" not in df.columns:
        df[["mapper1", "mapper2"]] = df["mapper_pair"].str.split("_vs_", expand=True)

    pivot_df = _create_symmetric_pivot(df, metric_col)

    plt.figure(figsize=(6, 5))
    sns.heatmap(
        pivot_df,
        annot=True,
        fmt=".3f",
        cmap="viridis",
        square=True,
        vmin=0,
        vmax=1,
    )
    plt.title(f"Jaccard Similarity ({strand.upper()} strand)", weight="bold")
    sns.despine()
    plt.tight_layout()
    plt.xlabel("")
    plt.ylabel("")
    plt.savefig(output_dir / f"jaccard_{strand}_heatmap.png", dpi=300)
    plt.close()
    print(f"✅ Saved: jaccard_{strand}_heatmap.png")


def plot_jaccard_vs_overlap(df, output_dir=Path(".")):
    """
    Create scatter plots comparing Jaccard similarity with overlap accuracy.
    This shows whether the metrics are correlated or capture different aspects.
    """
    if "mapper1" not in df.columns:
        df[["mapper1", "mapper2"]] = df["mapper_pair"].str.split("_vs_", expand=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for idx, strand in enumerate(["fw", "rv"]):
        ax = axes[idx]
        jaccard_col = f"jaccard_{strand}"
        overlap_col = f"{strand}_acc"

        if jaccard_col not in df.columns or overlap_col not in df.columns:
            continue

        plot_df = df[[jaccard_col, overlap_col, "mapper_pair"]].dropna()

        for mapper_pair in plot_df["mapper_pair"].unique():
            pair_data = plot_df[plot_df["mapper_pair"] == mapper_pair]
            ax.scatter(
                pair_data[overlap_col],
                pair_data[jaccard_col],
                alpha=0.5,
                label=mapper_pair,
                s=20,
            )

        ax.plot([0, 1], [0, 1], "k--", alpha=0.3, label="y=x")

        ax.set_xlabel(f"Overlap Accuracy ({strand.upper()})")
        ax.set_ylabel(f"Jaccard Similarity ({strand.upper()})")
        ax.set_title(f"{strand.upper()} Strand: Jaccard vs Overlap")
        ax.set_xlim(0, 1.05)
        ax.set_ylim(0, 1.05)
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
        ax.grid(True, alpha=0.3)
        sns.despine(ax=ax)

    plt.tight_layout()
    plt.savefig(
        output_dir / "jaccard_vs_overlap_scatter.png", dpi=300, bbox_inches="tight"
    )
    plt.close()
    print("✅ Saved: jaccard_vs_overlap_scatter.png")


def plot_jaccard_boxplot(df, output_dir=Path(".")):
    """
    Create boxplots showing distribution of Jaccard scores per mapper pair.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for idx, strand in enumerate(["fw", "rv"]):
        ax = axes[idx]
        jaccard_col = f"jaccard_{strand}"

        if jaccard_col not in df.columns:
            continue

        sample = df.sample(min(100, len(df)), random_state=42)

        sns.boxplot(data=sample, x="mapper_pair", y=jaccard_col, palette="husl", ax=ax)

        ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
        ax.set_ylabel(f"Jaccard Similarity")
        ax.set_xlabel("")
        ax.set_title(f"{strand.upper()} Strand Distribution")
        ax.set_ylim(0, 1.05)
        sns.despine(ax=ax)

    fig.suptitle(
        "Distribution of Jaccard Similarity Across Genes (100 samples)",
        fontsize=14,
        y=1.02,
    )
    plt.tight_layout()
    plt.savefig(output_dir / "jaccard_boxplot.png", dpi=300)
    plt.close()
    print("✅ Saved: jaccard_boxplot.png")


def plot_perfect_match_comparison(df, output_dir=Path(".")):
    """
    Compare the number of perfect matches across different metrics.
    """
    if "mapper1" not in df.columns:
        df[["mapper1", "mapper2"]] = df["mapper_pair"].str.split("_vs_", expand=True)

    perfect_data = []

    for mapper_pair in df["mapper_pair"].unique():
        pair_df = df[df["mapper_pair"] == mapper_pair]

        for strand in ["fw", "rv"]:
            jaccard_perfect_col = f"jaccard_perfect_{strand}"
            fully_covered_col = f"mapper_fully_covered_{strand}"

            if jaccard_perfect_col in pair_df.columns:
                avg_jaccard_perfect = pair_df[jaccard_perfect_col].mean()

                avg_fully_covered = (
                    pair_df[fully_covered_col].mean()
                    if fully_covered_col in pair_df.columns
                    else 0
                )

                perfect_data.append(
                    {
                        "mapper_pair": mapper_pair,
                        "strand": strand,
                        "metric": "Jaccard Perfect",
                        "avg_perfect": avg_jaccard_perfect,
                    }
                )

                perfect_data.append(
                    {
                        "mapper_pair": mapper_pair,
                        "strand": strand,
                        "metric": "Fully Covered",
                        "avg_perfect": avg_fully_covered,
                    }
                )

    if not perfect_data:
        print("⚠️ No perfect match data available")
        return

    perfect_df = pd.DataFrame(perfect_data)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for idx, strand in enumerate(["fw", "rv"]):
        ax = axes[idx]
        strand_data = perfect_df[perfect_df["strand"] == strand]

        sns.barplot(
            data=strand_data,
            x="mapper_pair",
            y="avg_perfect",
            hue="metric",
            palette="Set2",
            ax=ax,
        )

        ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
        ax.set_ylabel("Average Perfect Matches per Gene")
        ax.set_xlabel("")
        ax.set_title(f"{strand.upper()} Strand")
        ax.legend(title="Metric")
        sns.despine(ax=ax)

    fig.suptitle("Comparison of Perfect Match Metrics", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / "perfect_matches_comparison.png", dpi=300)
    plt.close()
    print("✅ Saved: perfect_matches_comparison.png")


def main_with_jaccard():
    gtamap_files = list(dir_gtamap.glob("gtamap_*.sam"))
    print(f"Found {len(gtamap_files)} GTAMap SAM files")
    meta = load_gene_meta()

    data = []
    with ProcessPoolExecutor(max_workers=30) as exe:
        futures = []
        for i, f in enumerate(gtamap_files):
            gene_id = f.stem.split("_")[1]
            mapper_paths = {
                "Minimap2": dir_others / f"minimap2_{gene_id}.sam",
                "HISAT2": dir_others / f"hisat2_{gene_id}.sam",
                "STAR": dir_others / f"star_{gene_id}.sam",
            }
            futures.append(
                exe.submit(compare_gene_with_jaccard, gene_id, f, mapper_paths, meta, i)
            )

        for fut in futures:
            gene_id, res = fut.result()
            if res is None:
                continue
            for pair, vals in res.items():
                data.append({"gene": gene_id, "mapper_pair": pair, **vals})

    df = pd.DataFrame(data)
    df.to_csv(output_dir / "gene_wise_all_vs_all_with_jaccard.csv", index=False)
    print(f"Saved metrics for {len(df)} gene-pair comparisons.")

    plot_mapper_agreement_heatmap(df, "fw")
    plot_mapper_agreement_heatmap(df, "rv")

    plot_jaccard_heatmap(df, "fw", output_dir)
    plot_jaccard_heatmap(df, "rv", output_dir)


def _create_symmetric_pivot(df, metric, negate_mirror=False):
    """
    Create a symmetric pivot table from asymmetric pairwise data.
    (This should already exist in your code, but included here for completeness)
    """
    pivot_df = (
        df.groupby(["mapper1", "mapper2"])[metric].mean().unstack(fill_value=np.nan)
    )

    all_mappers = sorted(set(pivot_df.index).union(set(pivot_df.columns)))
    pivot_df = pivot_df.reindex(index=all_mappers, columns=all_mappers)

    for m1 in all_mappers:
        for m2 in all_mappers:
            val_12 = pivot_df.loc[m1, m2]
            val_21 = pivot_df.loc[m2, m1]

            if pd.isna(val_12) and not pd.isna(val_21):
                pivot_df.loc[m1, m2] = -val_21 if negate_mirror else val_21
            elif pd.isna(val_21) and not pd.isna(val_12):
                pivot_df.loc[m2, m1] = -val_12 if negate_mirror else val_12

    return pivot_df


if __name__ == "__main__":
    main_with_jaccard()  # creates genewise_comparison_results/gene_wise_all_vs_all_with_jaccard.csv, so if not exist use this
    df = pd.read_csv(
        "genewise_comparison_results/gene_wise_all_vs_all_with_jaccard.csv"
    )

    # coverage metric C
    plot_mapper_agreement_heatmap(df, "fw")
    plot_mapper_agreement_heatmap(df, "rv")

    plot_jaccard_heatmap(df, "fw", output_dir)
    plot_jaccard_heatmap(df, "rv", output_dir)
