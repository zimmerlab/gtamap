from itertools import chain
import pandas as pd
import polars as pl
import pysam
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict
import itertools


def plot_mapper_agreement_heatmap(df, output_dir, strand="fw"):
    metric_col = f"mapper_agreement_{strand}"
    if "mapper1" not in df.columns:
        df[["mapper1", "mapper2"]] = df["mapper_pair"].str.split("_vs_", expand=True)

    # Just pivot - don't force symmetry!
    print(df)
    pivot_df = df.pivot_table(
        values=metric_col,
        index="mapper1",  # Row mapper (being tested)
        columns="mapper2",  # Column mapper (reference)
        aggfunc="mean",
    )
    print(pivot_df)

    plt.figure(figsize=(6, 5))
    sns.heatmap(
        pivot_df,
        annot=True,
        fmt=".3f",
        cmap="viridis",
        square=True,
        vmin=0,
        vmax=1,
        annot_kws={"size": 12},  # Increase numbers inside heatmap
    )
    # Increase x and y ticks size
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    # Increase colorbar tick size
    cbar = plt.gcf().axes[-1]  # Get colorbar axis
    cbar.tick_params(labelsize=12)

    plt.title(f"Mapper Coverage ({strand.upper()})", size=17, weight="bold")
    plt.xlabel("Reference Mapper (covered)", fontsize=14)
    plt.ylabel("Test Mapper (covering)", fontsize=14)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"mapper_coverage_asymmetric_new_{strand}.png", dpi=300)


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

plt.rcParams["figure.figsize"] = (6.5, 4.5)  # width x height in inches

sns.set_theme(style="ticks", context="paper", palette="husl")


def parse_global_sam_file(sam_path):
    """
    Parse a SAM file into forward/reverse read alignments.
    For each read, store:
      - reference name (chromosome)
      - merged alignment blocks (regions)
    """

    samfile = pysam.AlignmentFile(sam_path, "r")

    fw_dict = defaultdict(list)
    rv_dict = defaultdict(list)

    for read in samfile:
        # if read.cigarstring:
        if read.cigarstring and not read.is_secondary and not read.is_supplementary:

            # Extract aligned regions and chromosome
            ref_name = samfile.get_reference_name(read.reference_id)
            merged_blocks = merge_blocks(read.get_blocks())

            # Store in correct direction dict
            if read.is_reverse and read.is_read2:
                rv_dict[read.query_name].append(
                    {"chr": ref_name, "blocks": merged_blocks}
                )
            elif read.is_reverse and read.is_read1:
                fw_dict[read.query_name].append(
                    {"chr": ref_name, "blocks": merged_blocks}
                )
            elif read.is_forward and read.is_read1:
                fw_dict[read.query_name].append(
                    {"chr": ref_name, "blocks": merged_blocks}
                )
            else:
                rv_dict[read.query_name].append(
                    {"chr": ref_name, "blocks": merged_blocks}
                )

    samfile.close()

    return {"fw": fw_dict, "rv": rv_dict}


def mapper_agreement_coverage_genomewide(mapper1_dict, mapper2_dict, ma1, ma2):
    """
    Compute coverage/agreement of mapper1 vs mapper2 across reads (genome-wide).
    For genome-wide data, only compare blocks on the SAME chromosome.

    For each read shared between mapper1 and mapper2:
      - Compute total bases of mapper2 (across all chromosomes)
      - Compute how many mapper2 bases are covered by mapper1 blocks (per chromosome)
      - Compute coverage fraction = overlap / total_mapper2_bases

    Args:
        mapper1_dict: Dict mapping read_id -> list of alignment dicts
                      Each dict has: {"chr": str, "blocks": [(start, end), ...]}
        mapper2_dict: Dict mapping read_id -> list of alignment dicts
        ma1, ma2: Mapper names (for debugging)

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

        # Group blocks by chromosome for both mappers
        m2_by_chr = defaultdict(list)
        m1_by_chr = defaultdict(list)

        for aln in m2_aligns:
            m2_by_chr[aln["chr"]].extend(aln["blocks"])

        for aln in m1_aligns:
            m1_by_chr[aln["chr"]].extend(aln["blocks"])

        # Calculate total mapper2 bases across all chromosomes
        total_m2_bases = sum(
            sum(e - s for s, e in blocks) for blocks in m2_by_chr.values()
        )

        if total_m2_bases == 0:
            continue

        # Track covered positions per chromosome
        total_overlap = 0

        for chr_name, m2_blocks in m2_by_chr.items():
            # Only compare blocks on the SAME chromosome
            if chr_name not in m1_by_chr:
                continue  # mapper1 didn't align to this chromosome

            m1_blocks = m1_by_chr[chr_name]

            # Sort blocks for efficient overlap detection
            m2_blocks_sorted = sorted(m2_blocks)
            m1_blocks_sorted = sorted(m1_blocks)

            # Track covered positions on this chromosome
            covered_positions = set()

            for m_start, m_end in m2_blocks_sorted:
                for g_start, g_end in m1_blocks_sorted:
                    if g_end <= m_start:
                        continue
                    if g_start >= m_end:
                        break

                    overlap_start = max(m_start, g_start)
                    overlap_end = min(m_end, g_end)

                    if overlap_start < overlap_end:
                        covered_positions.update(range(overlap_start, overlap_end))

            total_overlap += len(covered_positions)

        # Compute coverage fraction
        frac = total_overlap / total_m2_bases if total_m2_bases > 0 else 0.0
        coverages.append(frac)

        if total_overlap == total_m2_bases:
            fully_covered += 1
        else:
            uncovered_bases_list.append(total_m2_bases - total_overlap)

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
        if start <= current_end:  # Overlapping or touching
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end

    merged.append((current_start, current_end))
    return merged


def jaccard_similarity_mappers(
    mapper1_intervals_dict, mapper2_intervals_dict, verbose=False
):
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
    total_overlap = 0

    for start1, end1 in blocks1:
        for start2, end2 in blocks2:
            overlap_start = max(start1, start2)
            overlap_end = min(end1, end2)

            if overlap_start < overlap_end:
                total_overlap += overlap_end - overlap_start

    return total_overlap


def calculate_blocks_union(blocks1, blocks2):
    total_bases1 = sum(end - start for start, end in blocks1)
    total_bases2 = sum(end - start for start, end in blocks2)
    intersection = calculate_blocks_intersection(blocks1, blocks2)

    return total_bases1 + total_bases2 - intersection


def plot_jaccard_heatmap(matrix, strand="fw"):

    plt.figure(figsize=(6, 5))

    sns.heatmap(
        matrix.astype(float),
        annot=True,
        fmt=".3f",
        cmap="viridis",
        cbar=True,
        vmin=0,
        vmax=1,
        annot_kws={"size": 12},
    )

    plt.title(f"Jaccard Similarity ({strand.upper()} strand)", weight="bold", size=17)
    plt.xlabel("")
    plt.ylabel("")

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    cbar = plt.gcf().axes[-1]
    cbar.tick_params(labelsize=12)

    sns.despine()
    plt.tight_layout()
    plt.savefig(f"jaccard_{strand}_heatmap.png", dpi=300)
    plt.close()


def compare_all_vs_all_jaccard(samfiles, output_dir):
    """
    Perform all-vs-all Jaccard similarity comparisons between mappers.
    """

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

    results = dict()
    for m1, m2 in itertools.product(mapper_names, mapper_names):
        print(f"Comparing {m1} vs {m2}")
        if m1 == m2:
            continue

        # forward reads
        result_fw = jaccard_similarity_mappers(
            mapper_data[m1]["fw"], mapper_data[m2]["fw"]
        )
        jaccard_fw.loc[m1, m2] = result_fw["avg_jaccard"]
        perfect_matches_fw.loc[m1, m2] = result_fw["perfect_matches"]

        # reverse reads
        result_rv = jaccard_similarity_mappers(
            mapper_data[m1]["rv"], mapper_data[m2]["rv"]
        )
        jaccard_rv.loc[m1, m2] = result_rv["avg_jaccard"]
        perfect_matches_rv.loc[m1, m2] = result_rv["perfect_matches"]

        fw_cov_12, fw_full_12, fw_uncov_12 = mapper_agreement_coverage_genomewide(
            mapper_data[m1]["fw"], mapper_data[m2]["fw"], m1, m2
        )
        rv_cov_12, rv_full_12, rv_uncov_12 = mapper_agreement_coverage_genomewide(
            mapper_data[m1]["rv"], mapper_data[m2]["rv"], m1, m2
        )

        # m2 vs m1 (REVERSED!)
        fw_cov_21, fw_full_21, fw_uncov_21 = mapper_agreement_coverage_genomewide(
            mapper_data[m2]["fw"], mapper_data[m1]["fw"], m2, m1  # ← Swapped arguments
        )
        rv_cov_21, rv_full_21, rv_uncov_21 = mapper_agreement_coverage_genomewide(
            mapper_data[m2]["rv"], mapper_data[m1]["rv"], m2, m1  # ← Swapped arguments
        )

        results[f"{m1}_vs_{m2}"] = {
            "mapper1": m1,
            "mapper2": m2,
            "mapper_agreement_fw": fw_cov_12,
            "mapper_agreement_rv": rv_cov_12,
            "mapper_fully_covered_fw": fw_full_12,
            "mapper_fully_covered_rv": rv_full_12,
            "mapper_uncovered_bases_fw": fw_uncov_12,
            "mapper_uncovered_bases_rv": rv_uncov_12,
        }

        results[f"{m2}_vs_{m1}"] = {
            "mapper1": m2,
            "mapper2": m1,
            "mapper_agreement_fw": fw_cov_21,
            "mapper_agreement_rv": rv_cov_21,
            "mapper_fully_covered_fw": fw_full_21,
            "mapper_fully_covered_rv": rv_full_21,
            "mapper_uncovered_bases_fw": fw_uncov_21,
            "mapper_uncovered_bases_rv": rv_uncov_21,
        }

    data = []
    for pair_key, metrics in results.items():
        row = {"mapper_pair": pair_key}
        row.update(metrics)
        data.append(row)

    df = pd.DataFrame(data)

    plot_mapper_agreement_heatmap(df, output_dir, "fw")
    plot_mapper_agreement_heatmap(df, output_dir, "rv")

    # save csv
    os.makedirs(output_dir, exist_ok=True)
    jaccard_fw.to_csv(os.path.join(output_dir, "jaccard_fw.csv"))
    jaccard_rv.to_csv(os.path.join(output_dir, "jaccard_rv.csv"))

    plot_jaccard_heatmap(jaccard_fw, strand="fw")
    plot_jaccard_heatmap(jaccard_rv, strand="rv")


if __name__ == "__main__":
    samfiles = [
        "/mnt/studtemp/weyrich/results_genome_wide/HISAT2_subset_all.sam",
        "/mnt/studtemp/weyrich/results_genome_wide/Minimap2_subset_all.sam",
        "/mnt/studtemp/weyrich/results_genome_wide/STAR_subset_all.sam",
    ]
    compare_all_vs_all_jaccard(samfiles, "res_jaccard")
