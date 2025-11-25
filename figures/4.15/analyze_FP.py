# Version 2.0
from pathlib import Path
from itertools import chain
import argparse
import polars as pl
import pysam
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import json
from collections import defaultdict
import sys
from pathlib import Path
import pysam
from collections import defaultdict
import sys
import re
from intervaltree import Interval, IntervalTree

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


def get_whiskers(series):
    q1, q3 = series.quantile([0.25, 0.75])
    iqr = q3 - q1
    lower = max(series.min(), q1 - 1.5 * iqr)
    upper = min(series.max(), q3 + 1.5 * iqr)
    return lower, upper


def boxplot_four_lists(
    fw_tp_avgs,
    rv_tp_avgs,
    fw_fp_avgs,
    rv_fp_avgs,
    title,
    ylabel,
    filename,
    out_dir,
    figsize=(8, 6),
):

    # Convert lists into DataFrame for Seaborn
    df = pd.DataFrame(
        {
            "FW TP": fw_tp_avgs,
            "RV TP": rv_tp_avgs,
            "FW FP": fw_fp_avgs,
            "RV FP": rv_fp_avgs,
        }
    )

    df_long = df.melt(var_name="Metric", value_name="Value")

    plt.figure(figsize=figsize)
    ax = sns.boxplot(
        data=df_long,
        x="Metric",
        y="Value",
        width=0.6,
        showfliers=False,
        boxprops=dict(edgecolor="black", linewidth=1.5, facecolor="#CFD5E8"),
        whiskerprops=dict(color="black", linewidth=1.5),
        capprops=dict(color="black", linewidth=1.5),
        medianprops=dict(color="black", linewidth=1.5),
    )

    # Compute zoom limits (like fancy_boxplot_zoom)
    whisker_lows, whisker_highs = [], []
    for metric, group in df_long.groupby("Metric"):
        low, high = get_whiskers(group["Value"])
        whisker_lows.append(low)
        whisker_highs.append(high)

    global_low = min(whisker_lows)
    global_high = max(whisker_highs)
    padding = 0.02 * (global_high - global_low)
    ax.set_ylim(global_low - padding, global_high + padding)

    # Titles and styling
    ax.set_title(title, fontsize=16, weight="bold", pad=15)
    ax.set_ylabel(ylabel, fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    sns.despine()

    # Save
    out_path = os.path.join(out_dir, filename)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved boxplot → {out_path}")


# Example usage


def count_mismatches(cigarstring):
    """Count total mismatched bases (X) in a CIGAR string."""
    if cigarstring is None:
        return 100000
    # Match numbers before X
    return sum(map(int, re.findall(r"(\d+)X", cigarstring)))


sns.set_theme(style="ticks", context="paper", palette="colorblind")

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

out_dir = "results_all"
os.makedirs(out_dir, exist_ok=True)

# --- Colors ---
fill_color = "#CFD5E8"
edge_color = "black"


def cigar_to_regvec(cigar, start):
    """
    Convert CIGAR string to regvec format.
    """
    pattern = re.compile(r"(\d+)([X=DNI])")  # Added D and I
    pos = start
    blocks = []

    for length, op in pattern.findall(cigar):
        length = int(length)
        if op in ("=", "X"):  # Treat matches and mismatches the same
            blocks.append(f"{pos}-{pos+length}")
            pos += length
        elif op in ("N", "D"):  # Both skip reference positions
            pos += length
        elif op == "I":  # Insertion doesn't advance reference position
            pass  # Do nothing, insertions don't consume reference

    return "|".join(blocks)

    # for f in indir.glob("*.sam"):


def parse_gtf(gtf_file):
    """
    Build interval trees per chromosome from GTF gene annotations,
    storing (gene_id, start, length) in the interval data.
    """
    chr_genes = defaultdict(IntervalTree)
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if cols[2] != "exon":
                continue
            chrom = cols[0].replace("chr", "")
            start = int(cols[3]) - 1  # GTF is 1-based
            end = int(cols[4])  # inclusive end
            length = end - start
            strand = cols[6]
            # Extract gene_id from attribute column
            attr_match = re.search(r'gene_id "([^"]+)"', cols[8])
            if not attr_match:
                continue
            gene_id = attr_match.group(1)
            if "." in gene_id:
                gene_id = gene_id.split(".")[0]

            # Store (gene_id, start, length) in the interval
            chr_genes[chrom].add(Interval(start, end, (gene_id, start, length, strand)))

    return chr_genes


def convert_to_global(mapper1_data, offset):
    c_dict = dict()
    for category, read_dict in mapper1_data.items():
        if category in ["fw_mm", "rv_mm"]:
            c_dict[category] = read_dict
            continue

        converted_regions = defaultdict(list)
        for read_id, regions in read_dict.items():
            new_region_list = []
            for start, stop in regions[0]:
                new_region_list.append((start + offset - 1, stop + offset - 1))
            converted_regions[read_id].append(new_region_list)

        c_dict[category] = converted_regions

    return c_dict


def analyze(mapper_data, ground_truth_data, verbose):
    """
    Compare mapper data with ground truth data.

    Args:
        mapper_data: Dictionary with 'fw' and 'rv' keys containing alignment blocks from pysam (each id maps to several possible maps (list of list of tuples))
        ground_truth_data: Dictionary with 'fw' and 'rv' keys containing regions

    Returns:
        Dictionary with comparison metrics
    """
    # Get the sets of read IDs
    mapper_fw_reads = set(mapper_data["fw"].keys())
    mapper_rv_reads = set(mapper_data["rv"].keys())
    gt_fw_reads = set(ground_truth_data["readid"].to_list())
    gt_rv_reads = set(ground_truth_data["readid"].to_list())

    # All mapped reads and ground truth reads
    mapped_reads = mapper_fw_reads.intersection(mapper_rv_reads)
    gt_reads = gt_fw_reads.union(gt_rv_reads)

    true_positives = mapped_reads.intersection(gt_reads)
    false_positives = mapped_reads - gt_reads

    mm_dict = {
        "fw": dict(),
        "rv": dict(),
    }

    fw_mm_tp = []
    rv_mm_tp = []
    for id in true_positives:
        mm = min(mapper_data["fw_mm"][id], key=len)
        fw_mm_tp.append(mm)
        mm_dict["fw"][id] = len(mm)
        mm = min(mapper_data["rv_mm"][id], key=len)
        rv_mm_tp.append(mm)
        mm_dict["rv"][id] = len(mm)

    fw_mm_fp = []
    rv_mm_fp = []
    for id in false_positives:
        if len(mapper_data["fw"][id]) == 0:
            continue
        mm = min(mapper_data["fw_mm"][id], key=len)
        mm_dict["fw"][id] = len(mm)
        fw_mm_fp.append(len(mm))
        mm = min(mapper_data["fw_mm"][id], key=len)
        rv_mm_fp.append(len(mm))
        mm_dict["rv"][id] = len(mm)

    data = []

    for mm in fw_mm_tp:
        data.append({"length": len(mm), "type": "TP", "read": "fw"})
    for mm in fw_mm_fp:
        data.append({"length": mm, "type": "FP", "read": "fw"})

    for mm in rv_mm_tp:
        data.append({"length": len(mm), "type": "TP", "read": "rv"})
    for mm in rv_mm_fp:
        data.append({"length": mm, "type": "FP", "read": "rv"})

    entries_per_read_fw = [
        len(alis)  # number of times this read appears in the SAM
        for read_id, alis in mapper_data["fw"].items()
        if read_id in gt_fw_reads
    ]
    avg_entries_per_read_fw = (
        sum(entries_per_read_fw) / len(entries_per_read_fw)
        if len(entries_per_read_fw) != 0
        else 0
    )
    entries_per_read_rv = [
        len(alis)
        for read_id, alis in mapper_data["rv"].items()
        if read_id in gt_rv_reads
    ]
    avg_entries_per_read_rv = (
        sum(entries_per_read_rv) / len(entries_per_read_rv)
        if len(entries_per_read_rv) != 0
        else 0
    )

    entries_per_read_fw_fp = [
        len(alis)
        for read_id, alis in mapper_data["fw"].items()
        if read_id not in gt_fw_reads
    ]
    avg_entries_per_read_fw_fp = (
        sum(entries_per_read_fw_fp) / len(entries_per_read_fw_fp)
        if len(entries_per_read_fw_fp) != 0
        else 0
    )

    entries_per_read_rv_fp = [
        len(alis)
        for read_id, alis in mapper_data["rv"].items()
        if read_id not in gt_rv_reads
    ]
    avg_entries_per_read_rv_fp = (
        sum(entries_per_read_rv_fp) / len(entries_per_read_rv_fp)
        if len(entries_per_read_rv_fp) != 0
        else 0
    )

    return (
        data,
        false_positives,
        mm_dict,
        avg_entries_per_read_fw,
        avg_entries_per_read_fw_fp,
        avg_entries_per_read_rv,
        avg_entries_per_read_rv_fp,
    )


def parse_sam_file(sam_path):
    samfile = pysam.AlignmentFile(sam_path, "r")
    fw_dict = defaultdict(list)
    rv_dict = defaultdict(list)
    rv_dict_mm = defaultdict(list)
    fw_dict_mm = defaultdict(list)

    for read in samfile:
        if read.cigarstring:
            if read.is_reverse and read.is_read2:
                rv_dict[read.query_name].append(merge_blocks(read.get_blocks()))
                mm = get_mismatches(read, True)
                rv_dict_mm[read.query_name].append(mm)
            elif read.is_reverse and read.is_read1:
                fw_dict[read.query_name].append(merge_blocks(read.get_blocks()))
                mm = get_mismatches(read, True)
                fw_dict_mm[read.query_name].append(mm)
            elif read.is_forward and read.is_read1:
                fw_dict[read.query_name].append(merge_blocks(read.get_blocks()))
                mm = get_mismatches(read, False)
                fw_dict_mm[read.query_name].append(mm)
            else:
                rv_dict[read.query_name].append(merge_blocks(read.get_blocks()))
                mm = get_mismatches(read, False)
                rv_dict_mm[read.query_name].append(mm)

    samfile.close()

    return {"fw": fw_dict, "rv": rv_dict, "fw_mm": fw_dict_mm, "rv_mm": rv_dict_mm}


def parse_ground_truth(gt_file: str):
    """Load the full ground truth TSV file and parse mutation-related columns."""

    columns_to_read = [
        "readid",
        "gene",
        "fw_mut",
        "rw_mut",
        "fw_mut_combined",
        "rw_mut_combined",
        "fw_seq_err",
        "rw_seq_err",
        "gene_length",
        "gene_start",
        "strand",
    ]

    df = pl.read_csv(gt_file, separator="\t", columns=columns_to_read).with_columns(
        pl.col("readid").cast(pl.Utf8)
    )

    # Helper to split comma-separated strings into list[Int64], handling empty strings
    def split_int_list(col):
        return (
            pl.when(pl.col(col) == "")
            .then(pl.lit([]))
            .otherwise(pl.col(col).str.split(",").cast(pl.List(pl.Int64)))
            .alias(col)
        )

    df = df.with_columns(
        split_int_list("fw_mut"),
        split_int_list("rw_mut"),
        split_int_list("fw_mut_combined"),
        split_int_list("rw_mut_combined"),
        split_int_list("fw_seq_err"),
        split_int_list("rw_seq_err"),
    )

    return df


def parse_region(region: str, is_rev: bool):
    """Convert region string '13-19|1000-10000' to [(13,19), (1000,10000)]."""
    if not region:
        return []

    regions = region.split("|")
    if not is_rev:
        return [
            # (int(start) - gene_start, int(end) - gene_start)
            (int(start) - 1, int(end) - 1)
            for start, end in (r.split("-") for r in regions)
        ]
    else:
        return [
            # (gene_length - int(end), gene_length - int(start))
            (int(start) - 1, int(end) - 1)
            for start, end in (r.split("-") for r in regions)
        ][::-1]


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
        if start <= current_end:  # Overlapping or touching
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end

    merged.append((current_start, current_end))
    return merged


def get_mismatches(read, flip_coords=False):
    if not read.cigartuples:
        # print(f"Had to skip extracting mm from {read.query_name} due to missing cigar")
        return ["*"]
    cigar_tuples = read.cigartuples
    mm = []
    read_pos = 0

    for op, length in cigar_tuples:
        if op in [0, 7, 8]:
            if op == 8:  # mismatch
                mm.extend(range(read_pos, read_pos + length))
            read_pos += length
        elif op == 1:  # insertion
            read_pos += length
        # D, N, H, P don't consume read bases — ignore

    if flip_coords:
        read_len = read.query_length
        mm = [read_len - 1 - pos for pos in mm]
        mm.reverse()

    return mm


def subset_ground_truth_by_gene(df: pl.DataFrame, gene_id: str) -> pl.DataFrame:
    return df.filter(pl.col("gene") == gene_id)


if __name__ == "__main__":
    main_df = pd.DataFrame()

    chr_genes = parse_gtf(
        "/home/malte/projects/bachelor_arbeit/data/genome/gencode.v48.annotation.gtf"
    )

    all_data = []
    indir = Path("./results_mixed_test/benchmarking_dataset_2/")
    i = -1
    gt_df = parse_ground_truth(
        "/home/malte/projects/bachelor_arbeit/data/simulated_reads/benchmarking_dataset_2/read.mappinginfo",
    )

    # counts for fps
    # TODO: ensure only counting each read once
    tp_intronic = 0
    tp_exonic = 0
    fp_intronic = 0
    fp_exonic = 0

    reads_in_analysis = 0

    fw_fp_avgs = []
    fw_tp_avgs = []
    rv_fp_avgs = []
    rv_tp_avgs = []
    for f in indir.glob("*.sam"):
        i += 1
        input_sam = f
        g_id = ""
        if g_id == "":
            g_id = f.name.split("/")[-1].split(".")[0]
            if "_" in g_id:
                g_id = g_id.split("_")[1]

        mapper1_data = parse_sam_file(input_sam.absolute())

        print(f"analyzing {g_id} {i}")

        gt_gene = subset_ground_truth_by_gene(gt_df, g_id)
        data, fps, reads_all, fw_tp_avg, fw_fp_avg, rv_tp_avg, rv_fp_avg = analyze(
            mapper1_data, gt_gene, False
        )
        fw_fp_avgs.append(fw_fp_avg)
        fw_tp_avgs.append(fw_tp_avg)
        rv_fp_avgs.append(rv_fp_avg)
        rv_tp_avgs.append(rv_tp_avg)

        mapper_reads = set(mapper1_data["fw"].keys()).intersection(
            set(mapper1_data["rv"].keys())
        )
        for d in data:
            d["g_id"] = g_id

        all_data.extend(data)

        samfile = pysam.AlignmentFile(str(input_sam.absolute()), "r")
        # Step 1: Collect min mismatch alignments for this gene
        gene_reads = {"fw": dict(), "rv": dict()}  # store relevant info per read
        for read in samfile:
            rid = read.query_name
            # if rid not in reads_all["fw"].keys() or rid not in reads_all["rv"].keys():
            #     continue

            if rid not in reads_all["fw"].keys() and rid not in reads_all["rv"].keys():
                continue

            is_rv = read.is_read2
            chrom = (
                samfile.get_reference_name(read.reference_id)
                if read.reference_id >= 0
                else "*"
            )
            if not read.cigarstring or chrom not in chr_genes:
                continue

            mm = len(get_mismatches(read, True))
            blocks = cigar_to_regvec(read.cigarstring, read.reference_start)

            # Only store the alignment if it's the minimal mm seen so far for this read/orientation
            key = "rv" if is_rv else "fw"
            if rid not in gene_reads[key] or mm < gene_reads[key][rid]["mm"]:
                gene_reads[key][rid] = {"mm": mm, "blocks": blocks, "chrom": chrom}

        samfile.close()

        # Step 2: Iterate over the selected minimal reads and count exon/intron features
        for orientation in ["fw", "rv"]:
            for rid, info in gene_reads[orientation].items():
                chrom = info["chrom"]
                blocks = info["blocks"]

                # check exon/intron overlap
                all_blocks_overlap = True
                j = 0
                # for block in blocks.split("|"):
                #     s, e = map(int, block.split("-"))
                #     if not chr_genes[chrom][s:e]:
                #         if rid not in fps:
                #             print(e - s)
                #         all_blocks_overlap = False
                #         break

                for block in blocks.split("|"):
                    s, e = map(int, block.split("-"))
                    if not chr_genes[chrom][s:e]:
                        continue
                    j += 1

                # if not all_blocks_overlap:
                if j == 0:
                    if rid in fps:
                        fp_intronic += 1
                    else:
                        tp_intronic += 1
                    reads_in_analysis += 1
                    continue

                # exonic reads
                if rid in fps:
                    fp_exonic += 1
                else:
                    tp_exonic += 1
                reads_in_analysis += 1

    main_df = pd.DataFrame(all_data)
    # --- Summary counts ---
    tp_count = (main_df["type"] == "TP").sum()
    fp_count = (main_df["type"] == "FP").sum()

    plot_df = main_df

    print(f"\nTotal True Positives (TP): {tp_count}")
    print(f"Total False Positives (FP): {fp_count}")
    print(f"Total entries: {len(main_df)}\n")
    print(f"Intronic FPs: {fp_intronic}")
    print(f"Exonic FPs: {fp_exonic}")
    print(f"Intronic TPs: {tp_intronic}")
    print(f"Exonic TPs: {tp_exonic}")
    print(f"Reads in Exon Intron Analysis: {reads_in_analysis}")

    # --- Forward reads HIST ---
    plt.figure()
    sns.histplot(
        data=plot_df[plot_df["read"] == "fw"],
        x="length",
        hue="type",
        bins=50,
        multiple="dodge",  # overlay instead of stacking
        alpha=0.6,
        palette={"TP": "green", "FP": "red"},
    )
    plt.title("Amount of Mismatches of FW Mappings", weight="bold")
    plt.xlabel("Mismatches in Mapping")
    plt.ylabel("Count")
    sns.despine()
    plt.tight_layout()
    plt.savefig("fw_hist.png", dpi=300)
    plt.close()

    # --- Reverse reads HIST ---
    plt.figure()
    sns.histplot(
        data=plot_df[plot_df["read"] == "rv"],
        x="length",
        hue="type",
        bins=50,
        multiple="dodge",  # overlay instead of stacking
        alpha=0.6,
        palette={"TP": "green", "FP": "red"},
    )
    plt.title("Amount of Mismatches of RV Mappings", weight="bold")
    plt.xlabel("Mismatches in Mapping")
    plt.ylabel("Count")
    sns.despine()
    plt.tight_layout()
    plt.savefig("rv_hist.png", dpi=300)
    plt.close()

    boxplot_four_lists(
        fw_tp_avgs,
        rv_tp_avgs,
        fw_fp_avgs,
        rv_fp_avgs,
        title="Amount of SAM Entries on Average per Gene",
        ylabel="Average Value",
        filename="gtamap_tp_fp_boxplot.png",
        out_dir="./",
    )
