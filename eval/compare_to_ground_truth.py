# Version 2.0
from itertools import chain
import argparse
import polars as pl
import pysam
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import json
from collections import defaultdict
import sys


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


def compare_to_ground_truth(mapper_data, ground_truth_data, verbose):
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
    gt_fw_reads = set(ground_truth_data["fw"].keys())
    gt_rv_reads = set(ground_truth_data["rv"].keys())

    # All mapped reads and ground truth reads
    mapped_reads = mapper_fw_reads.union(mapper_rv_reads)
    gt_reads = gt_fw_reads.union(gt_rv_reads)

    true_positives = mapped_reads.intersection(gt_reads)
    false_positives = mapped_reads - gt_reads
    false_negatives = gt_reads - mapped_reads

    precision = len(true_positives) / len(mapped_reads) if mapped_reads else 0
    recall = len(true_positives) / len(gt_reads) if gt_reads else 0
    f1_score = (
        2 * (precision * recall) / (precision + recall)
        if (precision + recall) > 0
        else 0
    )

    fw_accuracy, complient_fw, avg_missed_bases_fw, mm_delta_fw = positional_accuracy(
        mapped_intervals_dict=mapper_data["fw"],
        true_intervals_dict=ground_truth_data["fw"],
        true_mm=ground_truth_data["fw_mut_combined"],
        mapped_mm=mapper_data["fw_mm"],
        read_type="fw",
        verbose=verbose,
    )
    rw_accuracy, complient_rv, avg_missed_bases_rv, mm_delta_rv = positional_accuracy(
        mapped_intervals_dict=mapper_data["rv"],
        true_intervals_dict=ground_truth_data["rv"],
        true_mm=ground_truth_data["rv_mut_combined"],
        mapped_mm=mapper_data["rv_mm"],
        read_type="rv",
        verbose=verbose,
    )

    (
        fw_accuracy_strict,
        complient_fw_strict,
        avg_missed_bases_fw_strict,
        mm_delta_fw_strict,
    ) = positional_accuracy_strict(
        mapped_intervals_dict=mapper_data["fw"],
        true_intervals_dict=ground_truth_data["fw"],
        true_mm=ground_truth_data["fw_mut_combined"],
        mapped_mm=mapper_data["fw_mm"],
        read_type="fw",
        verbose=verbose,
    )
    (
        rw_accuracy_strict,
        complient_rv_strict,
        avg_missed_bases_rv_strict,
        mm_delta_rv_strict,
    ) = positional_accuracy_strict(
        mapped_intervals_dict=mapper_data["rv"],
        true_intervals_dict=ground_truth_data["rv"],
        true_mm=ground_truth_data["rv_mut_combined"],
        mapped_mm=mapper_data["rv_mm"],
        read_type="rv",
        verbose=verbose,
    )

    rv_level_two_recall = level_two_recall(
        pred_intervals_dict=mapper_data["rv"],
        true_intervals_dict=ground_truth_data["rv"],
    )

    fw_level_two_recall = level_two_recall(
        pred_intervals_dict=mapper_data["fw"],
        true_intervals_dict=ground_truth_data["fw"],
    )

    rv_level_three_recall = level_three_recall(
        pred_intervals_dict=mapper_data["rv"],
        true_intervals_dict=ground_truth_data["rv"],
    )

    fw_level_three_recall = level_three_recall(
        pred_intervals_dict=mapper_data["fw"],
        true_intervals_dict=ground_truth_data["fw"],
    )

    # fw_qual = calculate_avg_quality(mapper_data["fw_q"])
    # rw_qual = calculate_avg_quality(mapper_data["rv_q"])

    return {
        "true_positives": len(true_positives),
        "false_positives": len(false_positives),
        "false_negatives": len(false_negatives),
        "false_negative_ids": false_negatives,
        "precision": precision,
        "level_1_recall": recall,
        "f1_score": f1_score,
        "fw_overlap": fw_accuracy,
        "rv_overlap": rw_accuracy,
        "fw_level_2_recall": fw_level_two_recall,
        "rv_level_2_recall": rv_level_two_recall,
        "fw_level_3_recall": fw_level_three_recall,
        "rv_level_3_recall": rv_level_three_recall,
        "perfect_matches_fw": complient_fw,
        "perfect_matches_rv": complient_rv,
        "avg_missaligned_positions_fw": avg_missed_bases_fw,
        "avg_missaligned_positions_rv": avg_missed_bases_rv,
        "mm_delta_fw": mm_delta_fw,
        "mm_delta_rv": mm_delta_rv,
        "fw_overlap_strict": fw_accuracy_strict,
        "rv_overlap_strict": rw_accuracy_strict,
        "avg_missaligned_positions_fw_strict": avg_missed_bases_fw_strict,
        "avg_missaligned_positions_rv_strict": avg_missed_bases_rv_strict,
    }


def calculate_avg_quality(quality_dict):
    if len(quality_dict) == 0:
        return 0  # Avoid division by zero if the dictionary is empty
    total_quality = sum(quality_dict.values())  # Sum up the quality values
    avg_quality = total_quality / len(quality_dict)  # Calculate the average
    return avg_quality


def check_blocks_overlap(mapper_blocks, gt_regions):
    """
    Check if mapping blocks sufficiently overlap with ground truth regions.

    Args:
        mapper_blocks: List of tuples (start, end) from pysam.get_blocks()
        gt_regions: List of tuples (start, end) from ground truth

    Returns:
        Boolean indicating if there's sufficient overlap
    """
    if not mapper_blocks or not gt_regions:
        return False

    # Calculate total overlap between blocks and regions
    total_overlap = 0

    # Calculate total ground truth length
    gt_length = sum(end - start for start, end in gt_regions)

    for m_start, m_end in mapper_blocks:
        for g_start, g_end in gt_regions:
            # Calculate overlap between this block and region
            overlap_start = max(m_start, g_start)
            overlap_end = min(m_end, g_end)

            if overlap_end > overlap_start:
                total_overlap += overlap_end - overlap_start

    # Consider position correct if at least 80% of ground truth is covered
    min_required_overlap = 0.8 * gt_length

    return total_overlap >= min_required_overlap if gt_length > 0 else False


def parse_sam_file(sam_path):

    samfile = pysam.AlignmentFile(sam_path, "r")
    fw_dict = defaultdict(list)
    rv_dict = defaultdict(list)
    rv_dict_mm = defaultdict(list)
    fw_dict_mm = defaultdict(list)

    for read in samfile:
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


def parse_ground_truth(gt_file: str, gene_id: str):
    """Parse the read.mappinginfo file containing ground truth."""

    columns_to_read = [
        "readid",
        "gene",
        "fw_regvec",
        "rw_regvec",
        "fw_seq_err",
        "rw_seq_err",
        "gene_length",
        "gene_start",
        "strand",
        "fw_mut",
        "rw_mut",
        "fw_mut_combined",
        "rw_mut_combined",
    ]

    # df = pl.read_csv(gt_file, separator="\t", columns=columns_to_read)
    df = pl.read_csv(gt_file, separator="\t", columns=columns_to_read).with_columns(
        pl.col("readid").cast(pl.Utf8)
    )

    filtered_df = df.filter(df["gene"] == gene_id)

    filtered_df = filtered_df.with_columns(
        filtered_df["fw_mut"].str.split(",").cast(pl.List(pl.Int64)).alias("fw_mut"),
        filtered_df["rw_mut"].str.split(",").cast(pl.List(pl.Int64)).alias("rw_mut"),
        filtered_df["fw_seq_err"]
        .str.split(",")
        .cast(pl.List(pl.Int64))
        .alias("fw_seq_err"),
        filtered_df["rw_seq_err"]
        .str.split(",")
        .cast(pl.List(pl.Int64))
        .alias("rw_seq_err"),
        filtered_df["fw_mut_combined"]
        .str.split(",")
        .cast(pl.List(pl.Int64))
        .alias("fw_mut_combined"),
        filtered_df["rw_mut_combined"]
        .str.split(",")
        .cast(pl.List(pl.Int64))
        .alias("rw_mut_combined"),
    )
    gene_strand = filtered_df["strand"].item(0)
    is_rev = gene_strand == "-"
    gene_start = filtered_df["gene_start"].item(0)
    gene_length = filtered_df["gene_length"].item(0)

    # Create separate dictionaries for forward and reverse regions
    # is it +?
    if not is_rev:
        fw_dict = {
            row["readid"]: parse_region(row["fw_regvec"], False)
            for row in filtered_df.iter_rows(named=True)
        }

        rv_dict = {
            row["readid"]: parse_region(row["rw_regvec"], False)
            for row in filtered_df.iter_rows(named=True)
        }
    else:
        fw_dict = {
            row["readid"]: parse_region(row["fw_regvec"], True)[::-1]
            for row in filtered_df.iter_rows(named=True)
        }

        rv_dict = {
            row["readid"]: parse_region(row["rw_regvec"], True)[::-1]
            for row in filtered_df.iter_rows(named=True)
        }

    # get mutations in gene
    fw_mut = {row["readid"]: row["fw_mut"] for row in filtered_df.iter_rows(named=True)}
    rv_mut = {row["readid"]: row["rw_mut"] for row in filtered_df.iter_rows(named=True)}

    # get seq errors
    fw_seq_err = {
        row["readid"]: row["fw_seq_err"] for row in filtered_df.iter_rows(named=True)
    }
    rv_seq_err = {
        row["readid"]: row["rw_seq_err"] for row in filtered_df.iter_rows(named=True)
    }

    # get all mutations as combined set
    fw_mut_combined = {
        row["readid"]: row["fw_mut_combined"]
        for row in filtered_df.iter_rows(named=True)
    }
    rv_mut_combined = {
        row["readid"]: row["rw_mut_combined"]
        for row in filtered_df.iter_rows(named=True)
    }

    return {
        "fw": fw_dict,
        "rv": rv_dict,
        "fw_mut": fw_mut,
        "rv_mut": rv_mut,
        "fw_seq_err": fw_seq_err,
        "rv_seq_err": rv_seq_err,
        "fw_mut_combined": fw_mut_combined,
        "rv_mut_combined": rv_mut_combined,
        "gene_start": gene_start,
        "gene_length": gene_length,
    }, df.height


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


def positional_accuracy_strict(
    mapped_intervals_dict,
    true_intervals_dict,
    true_mm,
    mapped_mm,
    read_type,
    verbose,
):
    """
    Compute positional accuracy for each read and return the average accuracy.

    :param pred_intervals_dict: Dictionary of predicted intervals {read_id: list([(start, stop), ...], [(start, stop), ...)]}
    :param true_intervals_dict: Dictionary of true intervals {read_id: [(start, stop), ...]}
    :return: Average positional accuracy across all reads.
    """
    accuracies = []
    missed_bases_list = []
    mm_delta = []

    number_of_complient_intervals = 0
    total = 0
    for read_id, true_intervals in true_intervals_dict.items():
        total += 1
        total_true_bases = sum(end - start for start, end in true_intervals)

        if read_id not in mapped_intervals_dict:
            # FN case: nothing mapped for this read
            accuracies.append(0.0)
            missed_bases_list.append(total_true_bases)
            mm_delta.append(len(true_mm[read_id]) if true_mm[read_id] else 0)
            if verbose:
                print(f"[\033[31mUNMAPPED READ\033[0m] {read_id} ({read_type})")
                print(f"> [   \033[32mREF\033[0m]: {true_intervals}")
                print(f"> [\033[32mMM REF\033[0m]: {true_mm[read_id]}")
            continue

        mapped_intervals = mapped_intervals_dict[
            read_id
        ]  # is a list of lists of tuples

        total_true_bases = sum(end - start for start, end in true_intervals)

        # evaluate all mapped intervals
        overlapping_bases_list = []  # store num overlaps per map in list
        for mappedInterval in mapped_intervals:
            overlapping_bases = 0

            for true_start, true_end in true_intervals:
                for pred_start, pred_end in mappedInterval:
                    if pred_end < true_start:
                        continue  # Skip if predicted interval is before the true interval
                    if pred_start > true_end:
                        break  # No need to check further (sorted order)

                    overlap_start = max(true_start, pred_start)
                    overlap_end = min(true_end, pred_end)

                    if overlap_start <= overlap_end:
                        overlapping_bases += overlap_end - overlap_start
            overlapping_bases_list.append(overlapping_bases)

        best_overlap = max(overlapping_bases_list)
        best_overlap_index = np.argmax(overlapping_bases_list)
        if best_overlap == total_true_bases:
            accuracies.append(best_overlap / total_true_bases)
            number_of_complient_intervals += 1
            mm_delta.append(0)
        else:
            # first get rv which best matches ground truth
            if mapped_mm[read_id][best_overlap_index] == ["*"]:
                continue  # skip reads which were not mapped
            accuracies.append(best_overlap / total_true_bases)
            missed_bases = total_true_bases - best_overlap
            missed_bases_list.append(missed_bases)
            if not true_mm[read_id]:
                true_mm[read_id] = []
            delta = len(mapped_mm[read_id][best_overlap_index]) - len(true_mm[read_id])
            mm_delta.append(delta)
            seen = set()
            if verbose:

                print(f"[\033[31mINCOMPLETE MAP\033[0m] of {read_id} ({read_type}):")

                for i, interval in enumerate(mapped_intervals):
                    key = "|".join(str(x) for x in chain.from_iterable(interval))
                    if key not in seen:
                        seen.add(key)
                        missed_bases = total_true_bases - overlapping_bases_list[i]
                        print(
                            f"> [\033[34mALTERNATIVE\033[0m] Missaligned positions: ({missed_bases})"
                        )
                        if len(mapped_mm[read_id][i]) <= len(true_mm[read_id]):
                            print(f"> [   \033[35mMAP\033[0m]: {mapped_intervals[i]}")
                        else:
                            print(f"> [   \033[33mMAP\033[0m]: {mapped_intervals[i]}")
                        print(f"> [   \033[32mREF\033[0m]: {true_intervals}")
                        if len(mapped_mm[read_id][i]) <= len(true_mm[read_id]):
                            print(f"> [\033[35mMM MAP\033[0m]: {mapped_mm[read_id][i]}")
                        else:
                            print(f"> [\033[33mMM MAP\033[0m]: {mapped_mm[read_id][i]}")
                        print(f"> [\033[32mMM REF\033[0m]: {true_mm[read_id]}")

    return (
        (
            sum(accuracies) / len(accuracies),
            number_of_complient_intervals,
            (
                sum(missed_bases_list) / len(missed_bases_list)
                if len(missed_bases_list) > 0
                else 0
            ),
            (sum(mm_delta) / len(mm_delta) if len(mm_delta) > 0 else 0),
        )
        if accuracies
        else (0, 0, 0, 0)
    )


def positional_accuracy(
    mapped_intervals_dict,
    true_intervals_dict,
    true_mm,
    mapped_mm,
    read_type,
    verbose,
):
    """
    Compute positional accuracy for each read and return the average accuracy.

    :param pred_intervals_dict: Dictionary of predicted intervals {read_id: list([(start, stop), ...], [(start, stop), ...)]}
    :param true_intervals_dict: Dictionary of true intervals {read_id: [(start, stop), ...]}
    :return: Average positional accuracy across all reads.
    """
    accuracies = []
    missed_bases_list = []
    mm_delta = []

    number_of_complient_intervals = 0
    total = 0
    for read_id, true_intervals in true_intervals_dict.items():
        total += 1
        if read_id not in mapped_intervals_dict:
            continue  # Skip if no prediction exists for this read

        mapped_intervals = mapped_intervals_dict[
            read_id
        ]  # is a list of lists of tuples

        total_true_bases = sum(end - start for start, end in true_intervals)

        # evaluate all mapped intervals
        overlapping_bases_list = []  # store num overlaps per map in list
        for mappedInterval in mapped_intervals:
            overlapping_bases = 0

            for true_start, true_end in true_intervals:
                for pred_start, pred_end in mappedInterval:
                    if pred_end < true_start:
                        continue  # Skip if predicted interval is before the true interval
                    if pred_start > true_end:
                        break  # No need to check further (sorted order)

                    overlap_start = max(true_start, pred_start)
                    overlap_end = min(true_end, pred_end)

                    if overlap_start <= overlap_end:
                        overlapping_bases += overlap_end - overlap_start
            overlapping_bases_list.append(overlapping_bases)

        best_overlap = max(overlapping_bases_list)
        best_overlap_index = np.argmax(overlapping_bases_list)
        if best_overlap == total_true_bases:
            accuracies.append(best_overlap / total_true_bases)
            number_of_complient_intervals += 1
            mm_delta.append(0)
        else:
            # first get rv which best matches ground truth
            if mapped_mm[read_id][best_overlap_index] == ["*"]:
                continue  # skip reads which were not mapped
            accuracies.append(best_overlap / total_true_bases)
            missed_bases = total_true_bases - best_overlap
            missed_bases_list.append(missed_bases)
            if not true_mm[read_id]:
                true_mm[read_id] = []
            delta = len(mapped_mm[read_id][best_overlap_index]) - len(true_mm[read_id])
            mm_delta.append(delta)
            seen = set()
            if verbose:

                print(f"[\033[31mINCOMPLETE MAP\033[0m] of {read_id} ({read_type}):")

                for i, interval in enumerate(mapped_intervals):
                    key = "|".join(str(x) for x in chain.from_iterable(interval))
                    if key not in seen:
                        seen.add(key)
                        missed_bases = total_true_bases - overlapping_bases_list[i]
                        print(
                            f"> [\033[34mALTERNATIVE\033[0m] Missaligned positions: ({missed_bases})"
                        )
                        if len(mapped_mm[read_id][i]) <= len(true_mm[read_id]):
                            print(f"> [   \033[35mMAP\033[0m]: {mapped_intervals[i]}")
                        else:
                            print(f"> [   \033[33mMAP\033[0m]: {mapped_intervals[i]}")
                        print(f"> [   \033[32mREF\033[0m]: {true_intervals}")
                        if len(mapped_mm[read_id][i]) <= len(true_mm[read_id]):
                            print(f"> [\033[35mMM MAP\033[0m]: {mapped_mm[read_id][i]}")
                        else:
                            print(f"> [\033[33mMM MAP\033[0m]: {mapped_mm[read_id][i]}")
                        print(f"> [\033[32mMM REF\033[0m]: {true_mm[read_id]}")

    return (
        (
            sum(accuracies) / len(accuracies),
            number_of_complient_intervals,
            (
                sum(missed_bases_list) / len(missed_bases_list)
                if len(missed_bases_list) > 0
                else 0
            ),
            (sum(mm_delta) / len(mm_delta) if len(mm_delta) > 0 else 0),
        )
        if accuracies
        else (0, 0, 0, 0)
    )


def level_three_recall(pred_intervals_dict, true_intervals_dict):
    """
    Set1: {(True QNAME, Start, Stop, Start, Stop,...), ... }
    Set2: {(Mapped QNAME, Start, Stop, Start, Stop,...), ... }
    Recall: Set1 nn Set2 / len(Set1)
    """

    set_one = set()
    set_two = set()

    for read_id, true_intervals in true_intervals_dict.items():
        for interval in true_intervals:
            set_one.add((read_id, interval[0], interval[1]))

    for read_id, pred_intervals in pred_intervals_dict.items():
        for mapping in pred_intervals:
            for interval in mapping:
                set_two.add((read_id, interval[0], interval[1]))

    TP = set_two.intersection(set_one)
    return len(TP) / len(set_one)


def plot_metrics(results, out_dir):
    sns.set(style="whitegrid")
    sns.set_palette("colorblind")
    sns.set("paper")
    sns.set_style("ticks")

    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 18,
            "axes.labelsize": 16,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 14,
        }
    )
    # Create separate figures instead of using GridSpec
    # Figure 1: Precision, Recall, and F1 Score with TP/FP/TN/FN text box
    fig1 = plt.figure(figsize=(10, 7))
    metrics = ["precision", "level_1_recall", "f1_score"]
    values = [results[m] for m in metrics]

    bars = plt.bar(metrics, values, width=0.6, edgecolor="black", linewidth=1.5)
    # plt.ylim(0.85, 1.0)  # Set y-limit to better visualize differences
    plt.title("Precision, Recall, and F1 Score", fontweight="bold", pad=20)
    plt.ylabel("Score")

    # Remove spines
    sns.despine()

    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2.0,
            height + 0.001,
            f"{height:.4f}",
            ha="center",
            va="bottom",
            fontweight="bold",
        )

    # Add a horizontal line at y=1.0 for reference
    plt.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5)

    # Extract confusion matrix values for the text box
    tp = results["true_positives"]
    fp = results["false_positives"]
    fn = results["false_negatives"]
    tn = results["true_negatives"]

    # Add text box with TP, FP, TN, FN information
    textbox_content = (
        f"Classification Results:\n"
        f"True Positives (TP): {tp}\n"
        f"False Positives (FP): {fp}\n"
        f"False Negatives (FN): {fn}\n"
        f"True Negatives (TN): {tn}"
    )

    plt.figtext(
        1,
        0.7,
        textbox_content,
        bbox=dict(facecolor="white", alpha=0.9, edgecolor="black"),
        fontsize=14,
        verticalalignment="top",
        fontweight="bold",
    )

    plt.tight_layout()
    out = f"{out_dir}/precision_recall_f1.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")

    # Figure 2: Accuracy Metrics
    fig2 = plt.figure(figsize=(10, 7))

    # Include accuracy metrics
    accuracy_metrics = [
        "fw_overlap",
        "rv_overlap",
        "fw_level_2_recall",
        "rv_level_2_recall",
        "fw_level_3_recall",
        "rv_level_3_recall",
    ]

    acc_values = [results[m] for m in accuracy_metrics]

    bars = plt.bar(
        accuracy_metrics, acc_values, width=0.6, edgecolor="black", linewidth=1.5
    )

    plt.title("Positional Accuracy Metrics", fontweight="bold", pad=20)
    plt.ylim(0.4, 1.0)
    plt.ylabel("Score")

    sns.despine()

    # Add value labels on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2.0,
            height + 0.001,
            f"{height:.4f}",
            ha="center",
            va="bottom",
            fontweight="bold",
        )

    # Add a horizontal line at y=1.0 for reference
    plt.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5)

    plt.xticks(rotation=90)
    plt.tight_layout()

    out = f"{out_dir}/accuracy_metrics.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")

    # Show all plots
    plt.show()

    # Save results as JSON
    out = f"{out_dir}/metrics_results.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=4)


def level_two_recall(pred_intervals_dict, true_intervals_dict):
    """
    Set1: {(True QNAME, Start), ... }
    Set2: {(Mapped QNAME, Start), ... }
    Recall: Set1 nn Set2 / len(Set1)
    """
    set_one = set()
    set_two = set()

    for read_id, true_intervals in true_intervals_dict.items():
        if len(true_intervals) != 0:
            set_one.add((read_id, true_intervals[0][0]))

    for read_id, pred_intervals in pred_intervals_dict.items():
        for mapping in pred_intervals:
            if len(mapping) != 0:
                set_two.add((read_id, mapping[0][0]))

    TP = set_two.intersection(set_one)
    return len(TP) / len(set_one)


def level_one_recall(pred_intervals_dict, true_intervals_dict):
    """
    Set1: {(True QNAME), ... }
    Set2: {(Mapped QNAME), ... }
    Recall: Set1 nn Set2 / len(Set1)
    """
    set_one = set()
    set_two = set()

    for read_id, _ in true_intervals_dict.items():
        set_one.add(read_id)

    for read_id, _ in pred_intervals_dict.items():
        set_two.add(read_id)

    TP = set_two.intersection(set_one)
    return len(TP) / len(set_one)


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare RNA mapper outputs")
    parser.add_argument("--sam1", required=True, help="First SAM file (your mapper)")
    parser.add_argument("--mapper", required=True, help="Name of mapper")
    parser.add_argument("--ground_truth", help="Ground truth file (optional)")
    parser.add_argument("--output_dir", default="comparison_results", help="Output dir")
    parser.add_argument(
        "--gene_id", default="", help="Gene ID to filter reads (e.g., ENSG00000005073)"
    )
    parser.add_argument("--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    plot_prefix = args.output_dir

    g_id = args.gene_id
    if g_id == "":
        g_id = args.sam1.split("/")[-1].split(".")[0]
        if "_" in g_id:
            g_id = g_id.split("_")[1]

    ground_truth_data, N = parse_ground_truth(
        args.ground_truth,
        g_id,
    )
    print("parsed gt")

    gene_start = ground_truth_data["gene_start"]

    mapper1_data = parse_sam_file(args.sam1)
    print("parsed sam")
    if args.mapper == "hisat2" or args.mapper == "star" or args.mapper == "minimap2":
        mapper1_data = convert_to_global(mapper1_data, gene_start)
        print("converted sam")

    out_summary = os.path.join(args.output_dir, f"{args.mapper}_{g_id}.stat")
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    print("analyzing ...")
    results = compare_to_ground_truth(mapper1_data, ground_truth_data, args.verbose)
    results["true_negatives"] = (
        N
        - results["true_positives"]
        - results["false_positives"]
        - results["false_negatives"]
    )

    print_order = [
        "perfect_matches_fw",
        "perfect_matches_rv",
        "true_positives",
        "false_positives",
        "false_negatives",
        "false_negative_ids",
        "precision",
        "level_1_recall",
        "f1_score",
        "fw_overlap",
        "rv_overlap",
        "fw_level_2_recall",
        "rv_level_2_recall",
        "fw_level_3_recall",
        "rv_level_3_recall",
        "avg_missaligned_positions_fw",
        "avg_missaligned_positions_rv",
        "mm_delta_fw",
        "mm_delta_rv",
        "fw_overlap_strict",
        "rv_overlap_strict",
        "avg_missaligned_positions_fw_strict",
        "avg_missaligned_positions_rv_strict",
    ]

    for key in print_order:
        val = results[key]
        if key == "false_negative_ids":
            print(
                f"[\033[32mMETRIC\033[0m] {key} (showing 30 examples of {len(val)}) = {list(val)[0:30]}"
            )
        else:
            print(f"[\033[32mMETRIC\033[0m] {key} = {val}")

    with open(out_summary, "w") as f:
        header_fields = []
        header_fields.append("geneID")
        fields = []
        fields.append(g_id)
        header_fields.append("tool")
        fields.append(args.mapper)
        for key in print_order:
            val = results[key]
            if key == "false_negative_ids":
                continue
            else:
                header_fields.append(key)
                fields.append(str(val))
        f.write("\t".join(header_fields))
        f.write("\n")
        f.write("\t".join(fields))

    # del results["false_negative_ids"]
    # plot_metrics(results, plot_prefix)
