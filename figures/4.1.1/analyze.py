import sys
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from pathlib import Path


def process_sam_file(input_sam):
    """Process a single SAM file and return mapper type, sets of read categories."""
    mapper = input_sam.name.split(".")[0].split("_")[0]
    gene = input_sam.name.split(".")[0].split("_")[1]
    if mapper == "gtamap":
        (
            all_reads_in_sam,
            unique_properly_paired,
            unique_singleton,
            multi_properly_paired,
            multi_singleton,
            unmapped_flagged,
        ) = parse_sam_ids_count_multimap_paired(input_sam.absolute())
    else:
        (
            all_reads_in_sam,
            unique_properly_paired,
            unique_singleton,
            multi_properly_paired,
            multi_singleton,
            unmapped_flagged,
        ) = parse_sam_ids(input_sam.absolute())
    return mapper, {
        "comb_all_reads_in_sam": all_reads_in_sam,
        "comb_unique_properly_paired": unique_properly_paired,
        "comb_unique_singleton": unique_singleton,
        "comb_multi_properly_paired": multi_properly_paired,
        "comb_multi_singleton": multi_singleton,
        "comb_unmapped_flagged": unmapped_flagged,
    }


def parse_sam(filename, total_reads=500000):
    """Parse SAM file to track unique/multimapped/unmapped reads, and split mapped reads into paired/singletons."""
    unique_properly_paired = set()
    unique_singleton = set()
    multi_properly_paired = set()
    multi_singleton = set()
    unmapped_flagged = set()

    with open(filename, "r") as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            read_id = fields[0]
            flag = int(fields[1])

            # Skip secondary/supplementary alignments
            if flag & 256 or flag & 2048:
                continue

            # Unmapped read
            if flag & 4:
                unmapped_flagged.add(read_id)
                continue

            # This read is mapped - determine if unique or multi
            nh_tag = None
            for tag in fields[11:]:
                if tag.startswith("NH:i:"):
                    nh_tag = int(tag.split(":")[-1])
                    break

            is_multi = nh_tag and nh_tag > 1

            # Determine pairing status (only for paired-end reads)
            if flag & 1:  # paired in sequencing
                if flag & 2:  # properly paired
                    if is_multi:
                        multi_properly_paired.add(read_id)
                    else:
                        unique_properly_paired.add(read_id)
                elif flag & 8:  # mate unmapped - this is a singleton
                    if is_multi:
                        multi_singleton.add(read_id)
                    else:
                        unique_singleton.add(read_id)
                else:  # mate mapped but not properly paired - still a singleton
                    if is_multi:
                        multi_singleton.add(read_id)
                    else:
                        unique_singleton.add(read_id)

    all_reads_in_sam = (
        unique_properly_paired
        | unique_singleton
        | multi_properly_paired
        | multi_singleton
    )

    unmapped_absent = total_reads - len(all_reads_in_sam)
    unmapped_flagged = unmapped_flagged - all_reads_in_sam

    return {
        "Total mapped": len(all_reads_in_sam),
        "Unique properly paired": len(unique_properly_paired),
        "Unique singletons": len(unique_singleton),
        "Multi properly paired": len(multi_properly_paired),
        "Multi singletons": len(multi_singleton),
        "Unmapped": unmapped_absent,
    }


def parse_sam_ids(filename, total_reads=500000):
    """Parse SAM file to track unique/multimapped/unmapped reads, and split mapped reads into paired/singletons."""
    unique_properly_paired = set()
    unique_singleton = set()
    multi_properly_paired = set()
    multi_singleton = set()
    unmapped_flagged = set()

    with open(filename, "r") as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            read_id = fields[0]
            flag = int(fields[1])

            # Skip secondary/supplementary alignments
            # if flag & 256 or flag & 2048:
            #    continue

            # Unmapped read
            if flag & 4:
                unmapped_flagged.add(read_id)
                continue

            # This read is mapped - determine if unique or multi
            nh_tag = None
            for tag in fields[11:]:
                if tag.startswith("NH:i:"):
                    nh_tag = int(tag.split(":")[-1])
                    break

            is_multi = nh_tag and nh_tag > 1

            # Determine pairing status (only for paired-end reads)
            if flag & 1:  # paired in sequencing
                if flag & 2:  # properly paired
                    if is_multi:
                        multi_properly_paired.add(read_id)
                    else:
                        unique_properly_paired.add(read_id)
                elif flag & 8:  # mate unmapped - this is a singleton
                    if is_multi:
                        multi_singleton.add(read_id)
                    else:
                        unique_singleton.add(read_id)
                else:  # mate mapped but not properly paired - still a singleton
                    if is_multi:
                        multi_singleton.add(read_id)
                    else:
                        unique_singleton.add(read_id)

    all_reads_in_sam = (
        unique_properly_paired
        | unique_singleton
        | multi_properly_paired
        | multi_singleton
    )

    unmapped_flagged = unmapped_flagged - all_reads_in_sam
    unique_properly_paired = unique_properly_paired - unmapped_flagged
    unique_singleton = unique_singleton - unmapped_flagged
    multi_properly_paired = multi_properly_paired - unmapped_flagged
    multi_singleton = multi_singleton - unmapped_flagged

    unique_properly_paired = unique_properly_paired - multi_properly_paired
    unique_singleton = unique_singleton - multi_singleton
    multi_properly_paired = multi_properly_paired - unique_properly_paired
    multi_singleton = multi_singleton - unique_singleton

    unique_singleton = unique_singleton - unique_properly_paired
    unique_properly_paired = unique_properly_paired - unique_singleton

    return (
        all_reads_in_sam,
        unique_properly_paired,
        unique_singleton,
        multi_properly_paired,
        multi_singleton,
        unmapped_flagged,
    )


from collections import defaultdict


def parse_sam_ids_count_multimap_paired(filename, total_reads=500000):
    """
    Parse a SAM file (paired-end aware) and classify reads.
    Multimapped reads are determined by how many times a QNAME appears.
    """
    unique_properly_paired = set()
    unique_singleton = set()
    multi_properly_paired = set()
    multi_singleton = set()
    unmapped_flagged = set()

    # Count occurrences of each read ID
    read_counts = defaultdict(int)
    paired_flags = {}  # store whether a read is paired

    # First pass: count mapped occurrences and store pairing info
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            read_id = fields[0]
            flag = int(fields[1])

            # Skip secondary/supplementary alignments
            if flag & 256 or flag & 2048:
                continue

            # Track if read is paired
            paired_flags[read_id] = bool(flag & 1)

            # Only count mapped reads
            if not (flag & 4):
                read_counts[read_id] += 1

    # Second pass: classify reads
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            read_id = fields[0]
            flag = int(fields[1])

            if flag & 256 or flag & 2048:
                continue

            if flag & 4:  # unmapped
                unmapped_flagged.add(read_id)
                continue

            # Determine expected count for paired vs single
            expected_count = 2 if paired_flags.get(read_id, False) else 1
            is_multi = read_counts[read_id] > expected_count

            # Determine pairing status
            if flag & 1:  # paired
                if flag & 2:  # properly paired
                    if is_multi:
                        multi_properly_paired.add(read_id)
                    else:
                        unique_properly_paired.add(read_id)
                elif flag & 8:  # mate unmapped -> singleton
                    if is_multi:
                        multi_singleton.add(read_id)
                    else:
                        unique_singleton.add(read_id)
                else:  # mate mapped but not properly paired -> singleton
                    if is_multi:
                        multi_singleton.add(read_id)
                    else:
                        unique_singleton.add(read_id)
            else:  # single-end
                if is_multi:
                    multi_singleton.add(read_id)
                else:
                    unique_singleton.add(read_id)

    all_reads_in_sam = (
        unique_properly_paired
        | unique_singleton
        | multi_properly_paired
        | multi_singleton
    )
    unmapped_flagged -= all_reads_in_sam

    unmapped_flagged = unmapped_flagged - all_reads_in_sam
    unique_properly_paired = unique_properly_paired - unmapped_flagged
    unique_singleton = unique_singleton - unmapped_flagged
    multi_properly_paired = multi_properly_paired - unmapped_flagged
    multi_singleton = multi_singleton - unmapped_flagged

    unique_properly_paired = unique_properly_paired - multi_properly_paired
    unique_singleton = unique_singleton - multi_singleton
    multi_properly_paired = multi_properly_paired - unique_properly_paired
    multi_singleton = multi_singleton - unique_singleton

    return (
        all_reads_in_sam,
        unique_properly_paired,
        unique_singleton,
        multi_properly_paired,
        multi_singleton,
        unmapped_flagged,
    )


def genome_wide():
    # these are the
    sam_files = [
        "/mnt/studtemp/weyrich/results_genome_wide/HISAT2_subset_all.sam",
        "/mnt/studtemp/weyrich/results_genome_wide/STAR_subset_all.sam",
        "/mnt/studtemp/weyrich/results_genome_wide/Minimap2_subset_all.sam",
    ]

    results = {}
    results_u = {}
    for sam_file in sam_files:
        if not os.path.exists(sam_file):
            print(f"Warning: {sam_file} not found")
            continue

        # data = parse_sam(sam_file)
        mapper = sam_file.split("/")[-1].split("_")[0]

        (
            all_reads_in_sam,
            unique_properly_paired,
            unique_singleton,
            multi_properly_paired,
            multi_singleton,
            unmapped_flagged,
        ) = parse_sam_ids(sam_file)

        results_u[mapper] = {
            "comb_all_reads_in_sam": all_reads_in_sam,
            "comb_unique_properly_paired": unique_properly_paired,
            "comb_unique_singleton": unique_singleton,
            "comb_multi_properly_paired": multi_properly_paired,
            "comb_multi_singleton": multi_singleton,
            "comb_unmapped_flagged": unmapped_flagged,
        }

        results[mapper] = {
            "Total mapped": len(all_reads_in_sam),
            "Properly paired": len(unique_properly_paired),
            "Singletons": len(unique_singleton),
            "Multi properly paired": len(multi_properly_paired),
            "Multi singletons": len(multi_singleton),
            "Unmapped": len(unmapped_flagged),
        }

    # --- Plotting ---
    fig, axes = plt.subplots(1, 3, figsize=(12, 7))
    axes = axes.flatten()

    categories = [
        "Total mapped",
        "Properly paired",
        "Singletons",
        "Multi properly paired",
        "Multi singletons",
    ]
    colors = sns.color_palette("husl", n_colors=len(categories) + 7)[7:]

    i = 0
    for idx, (sam_file, data) in enumerate(results.items()):
        ax = axes[idx]
        main_values = [data[c] for c in categories]

        # Plot main read categories
        bars1 = ax.bar(
            categories, main_values, color=colors, edgecolor="black", linewidth=1.1
        )

        # Add value labels
        for bars in [bars1]:
            for bar in bars:
                height = bar.get_height()
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    height,
                    f"{int(height)}",
                    ha="center",
                    va="bottom",
                    fontsize=13,
                    fontweight="bold",
                    rotation=90,
                )

        ax.set_title(
            os.path.basename(sam_file).split("_")[0], fontsize=20, fontweight="bold"
        )
        if i == 0:
            ax.set_ylabel("Number of Mapped QNAMES", fontsize=18)
        else:
            ax.set_ylabel("")
        ax.grid(axis="y", alpha=0.3, linestyle="--")
        ax.set_xticklabels(categories, rotation=90)
        ax.set_ylim(0, 500000)
        ax.tick_params(axis="x", labelsize=15)
        ax.tick_params(axis="y", labelsize=15)
        i += 1

    plt.suptitle(
        "Genome-wide Read Mapping Statistics by Mapper \n on Subset of 500,000 Read-pairs",
        fontsize=23,
        fontweight="bold",
        y=1.02,
    )
    sns.despine()
    plt.tight_layout()
    plt.savefig("mapping_statistics_paired.png", dpi=300, bbox_inches="tight")
    plt.show()

    from upsetplot import UpSet, from_contents

    # --- UpSet Plot for overlap of total mapped read IDs ---
    mapped_sets = {
        "HISAT2": results_u["HISAT2"]["comb_all_reads_in_sam"],
        "STAR": results_u["STAR"]["comb_all_reads_in_sam"],
        "Minimap2": results_u["Minimap2"]["comb_all_reads_in_sam"],
    }

    # Build the UpSet input data
    upset_data = from_contents(mapped_sets)

    # Plot
    fig = plt.figure(figsize=(16, 8))
    upset = UpSet(
        upset_data,
        subset_size="count",
        show_counts=True,
        sort_by="degree",  # group intersections by how many sets overlap
    )
    upset.plot()
    # Reduce text size of built-in labels
    for text in plt.gcf().findobj(match=plt.Text):
        if text.get_text().isdigit():
            text.set_fontsize(8)

    plt.suptitle(
        "Overlap of Total Mapped QNAMEs Across Mappers Genome-wide",
        fontsize=14,
        fontweight="bold",
    )
    plt.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1)
    plt.savefig(
        "upset_total_mapped_qnames_genome_wide.png", dpi=300, bbox_inches="tight"
    )
    plt.show()


def gene_wise():
    d = Path("/mnt/studtemp/weyrich/real_rna_benchmark_genes")
    i = 0
    data_minimap2 = {
        "comb_all_reads_in_sam": set(),
        "comb_unique_properly_paired": set(),
        "comb_unique_singleton": set(),
        "comb_multi_properly_paired": set(),
        "comb_multi_singleton": set(),
        "comb_unmapped_flagged": set(),
    }

    data_hisat2 = {
        "comb_all_reads_in_sam": set(),
        "comb_unique_properly_paired": set(),
        "comb_unique_singleton": set(),
        "comb_multi_properly_paired": set(),
        "comb_multi_singleton": set(),
        "comb_unmapped_flagged": set(),
    }

    data_star = {
        "comb_all_reads_in_sam": set(),
        "comb_unique_properly_paired": set(),
        "comb_unique_singleton": set(),
        "comb_multi_properly_paired": set(),
        "comb_multi_singleton": set(),
        "comb_unmapped_flagged": set(),
    }

    for f in d.glob("*.sam"):
        input_sam = f
        gene = input_sam.name.split("/")[-1].split(".")[0].split("_")[1]
        mapper = input_sam.name.split("/")[-1].split(".")[0].split("_")[0]
        print(f"{gene} {i}/17,575")
        i += 1
        (
            all_reads_in_sam,
            unique_properly_paired,
            unique_singleton,
            multi_properly_paired,
            multi_singleton,
            unmapped_flagged,
        ) = parse_sam_ids(input_sam.absolute())

        if mapper == "hisat2":
            data_hisat2["comb_all_reads_in_sam"] = (
                data_hisat2["comb_all_reads_in_sam"] | all_reads_in_sam
            )
            data_hisat2["comb_unique_properly_paired"] = (
                data_hisat2["comb_unique_properly_paired"] | unique_properly_paired
            )
            data_hisat2["comb_unique_singleton"] = (
                data_hisat2["comb_unique_singleton"] | unique_singleton
            )
            data_hisat2["comb_multi_properly_paired"] = (
                data_hisat2["comb_multi_properly_paired"] | multi_properly_paired
            )
            data_hisat2["comb_multi_singleton"] = (
                data_hisat2["comb_multi_singleton"] | multi_singleton
            )
            data_hisat2["comb_unmapped_flagged"] = (
                data_hisat2["comb_unmapped_flagged"] | unmapped_flagged
            )
        elif mapper == "star":
            data_star["comb_all_reads_in_sam"] = (
                data_star["comb_all_reads_in_sam"] | all_reads_in_sam
            )
            data_star["comb_unique_properly_paired"] = (
                data_star["comb_unique_properly_paired"] | unique_properly_paired
            )
            data_star["comb_unique_singleton"] = (
                data_star["comb_unique_singleton"] | unique_singleton
            )
            data_star["comb_multi_properly_paired"] = (
                data_star["comb_multi_properly_paired"] | multi_properly_paired
            )
            data_star["comb_multi_singleton"] = (
                data_star["comb_multi_singleton"] | multi_singleton
            )
            data_star["comb_unmapped_flagged"] = (
                data_star["comb_unmapped_flagged"] | unmapped_flagged
            )
        else:
            data_minimap2["comb_all_reads_in_sam"] = (
                data_minimap2["comb_all_reads_in_sam"] | all_reads_in_sam
            )
            data_minimap2["comb_unique_properly_paired"] = (
                data_minimap2["comb_unique_properly_paired"] | unique_properly_paired
            )
            data_minimap2["comb_unique_singleton"] = (
                data_minimap2["comb_unique_singleton"] | unique_singleton
            )
            data_minimap2["comb_multi_properly_paired"] = (
                data_minimap2["comb_multi_properly_paired"] | multi_properly_paired
            )
            data_minimap2["comb_multi_singleton"] = (
                data_minimap2["comb_multi_singleton"] | multi_singleton
            )
            data_minimap2["comb_unmapped_flagged"] = (
                data_minimap2["comb_unmapped_flagged"] | unmapped_flagged
            )

    for key, val in data_hisat2.items():
        print(f"{key}: {len(val)}")
    for key, val in data_star.items():
        print(f"{key}: {len(val)}")
    for key, val in data_minimap2.items():
        print(f"{key}: {len(val)}")

    # Prepare results dictionary for plotting
    results = {
        "HISAT2": {
            "Total mapped": len(data_hisat2["comb_all_reads_in_sam"])
            - len(data_hisat2["comb_unmapped_flagged"]),
            "Properly paired": len(data_hisat2["comb_unique_properly_paired"]),
            "Singletons": len(data_hisat2["comb_unique_singleton"]),
            "Multi properly paired": len(data_hisat2["comb_multi_properly_paired"]),
            "Multi singletons": len(data_hisat2["comb_multi_singleton"]),
            "Unmapped": len(data_hisat2["comb_unmapped_flagged"]),
        },
        "STAR": {
            "Total mapped": len(data_star["comb_all_reads_in_sam"])
            - len(data_star["comb_unmapped_flagged"]),
            "Properly paired": len(data_star["comb_unique_properly_paired"]),
            "Singletons": len(data_star["comb_unique_singleton"]),
            "Multi properly paired": len(data_star["comb_multi_properly_paired"]),
            "Multi singletons": len(data_star["comb_multi_singleton"]),
            "Unmapped": len(data_star["comb_unmapped_flagged"]),
        },
        "Minimap2": {
            "Total mapped": len(data_minimap2["comb_all_reads_in_sam"])
            - len(data_minimap2["comb_unmapped_flagged"]),
            "Properly paired": len(data_minimap2["comb_unique_properly_paired"]),
            "Singletons": len(data_minimap2["comb_unique_singleton"]),
            "Multi properly paired": len(data_minimap2["comb_multi_properly_paired"]),
            "Multi singletons": len(data_minimap2["comb_multi_singleton"]),
            "Unmapped": len(data_minimap2["comb_unmapped_flagged"]),
        },
    }

    # --- Plotting ---
    fig, axes = plt.subplots(1, 3, figsize=(12, 7))
    axes = axes.flatten()

    categories = [
        "Total mapped",
        "Properly paired",
        "Singletons",
        "Multi properly paired",
        "Multi singletons",
        "Unmapped",
    ]
    colors = sns.color_palette("husl", n_colors=len(categories) + 7)[7:]

    i = 0
    for idx, (mapper_name, data) in enumerate(results.items()):
        ax = axes[idx]
        main_values = [data[c] for c in categories]

        # Plot main read categories
        bars1 = ax.bar(
            categories, main_values, color=colors, edgecolor="black", linewidth=1.1
        )

        # Add value labels
        for bar in bars1:
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                height,
                f"{int(height)}",
                ha="center",
                va="bottom",
                fontsize=13,
                fontweight="bold",
            )

        ax.set_title(mapper_name, fontsize=14, fontweight="bold")
        if i == 0 or i == 2:
            ax.set_ylabel("Number of Reads", fontsize=11)
        else:
            ax.set_ylabel("")
        ax.grid(axis="y", alpha=0.3, linestyle="--")
        ax.set_xticklabels(categories, rotation=45)
        # Adjust ylim based on your data - you might need a different scale
        ax.set_ylim(0, 500000)  # or set a fixed value
        i += 1

    plt.suptitle(
        "Combined Gene-wise QNAME Statistics by Mapper",
        fontsize=22,
        fontweight="bold",
        y=1.02,
    )
    sns.despine()
    plt.tight_layout()
    plt.savefig("mapping_statistics_paired_genes.png", dpi=300, bbox_inches="tight")


def gene_wise_parallel():

    # Initialize data containers
    data = {
        "hisat2": {
            k: set()
            for k in [
                "comb_all_reads_in_sam",
                "comb_unique_properly_paired",
                "comb_unique_singleton",
                "comb_multi_properly_paired",
                "comb_multi_singleton",
                "comb_unmapped_flagged",
            ]
        },
        "star": {
            k: set()
            for k in [
                "comb_all_reads_in_sam",
                "comb_unique_properly_paired",
                "comb_unique_singleton",
                "comb_multi_properly_paired",
                "comb_multi_singleton",
                "comb_unmapped_flagged",
            ]
        },
        "minimap2": {
            k: set()
            for k in [
                "comb_all_reads_in_sam",
                "comb_unique_properly_paired",
                "comb_unique_singleton",
                "comb_multi_properly_paired",
                "comb_multi_singleton",
                "comb_unmapped_flagged",
            ]
        },
        "gtamap": {
            k: set()
            for k in [
                "comb_all_reads_in_sam",
                "comb_unique_properly_paired",
                "comb_unique_singleton",
                "comb_multi_properly_paired",
                "comb_multi_singleton",
                "comb_unmapped_flagged",
            ]
        },
    }

    d = Path("/mnt/studtemp/weyrich/results_gene_wise_gtamap")
    sam_files = list(d.glob("*.sam"))
    total_files = len(sam_files)
    print(f"Processing {total_files} SAM files in parallel...")

    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_sam_file, f): f for f in sam_files}
        for i, future in enumerate(as_completed(futures), 1):
            mapper, result = future.result()
            print(f"[{i}/{total_files}] Processed {futures[future].name}")
            for key, val in result.items():
                data[mapper][key] |= val  # Union of sets

    d = Path("/mnt/studtemp/weyrich/real_rna_benchmark_genes")
    sam_files = list(d.glob("*.sam"))
    total_files = len(sam_files)
    print(f"Processing {total_files} SAM files in parallel...")
    with ProcessPoolExecutor(max_workers=30) as executor:
        futures = {executor.submit(process_sam_file, f): f for f in sam_files}
        for i, future in enumerate(as_completed(futures), 1):
            mapper, result = future.result()
            print(f"[{i}/{total_files}] Processed {futures[future].name}")
            for key, val in result.items():
                data[mapper][key] |= val  # Union of sets

    # --- Summary ---
    for mapper_name, mapper_data in data.items():
        print(f"\n--- {mapper_name.upper()} ---")
        for key, val in mapper_data.items():
            print(f"{key}: {len(val)}")

    # --- Summary ---
    for mapper_name, mapper_data in data.items():
        print(f"\n--- {mapper_name.upper()} ---")
        for key, val in mapper_data.items():
            print(f"{key}: {len(val)}")

    # --- Prepare results for plotting ---
    results = {
        "HISAT2": {
            "Total mapped": len(data["hisat2"]["comb_all_reads_in_sam"])
            - len(data["hisat2"]["comb_unmapped_flagged"]),
            "Properly paired": len(data["hisat2"]["comb_unique_properly_paired"]),
            "Singletons": len(data["hisat2"]["comb_unique_singleton"]),
            "Multi properly paired": len(data["hisat2"]["comb_multi_properly_paired"]),
            "Multi singletons": len(data["hisat2"]["comb_multi_singleton"]),
            "Unmapped": len(data["hisat2"]["comb_unmapped_flagged"]),
        },
        "STAR": {
            "Total mapped": len(data["star"]["comb_all_reads_in_sam"])
            - len(data["star"]["comb_unmapped_flagged"]),
            "Properly paired": len(data["star"]["comb_unique_properly_paired"]),
            "Singletons": len(data["star"]["comb_unique_singleton"]),
            "Multi properly paired": len(data["star"]["comb_multi_properly_paired"]),
            "Multi singletons": len(data["star"]["comb_multi_singleton"]),
            "Unmapped": len(data["star"]["comb_unmapped_flagged"]),
        },
        "Minimap2": {
            "Total mapped": len(data["minimap2"]["comb_all_reads_in_sam"])
            - len(data["minimap2"]["comb_unmapped_flagged"]),
            "Properly paired": len(data["minimap2"]["comb_unique_properly_paired"]),
            "Singletons": len(data["minimap2"]["comb_unique_singleton"]),
            "Multi properly paired": len(
                data["minimap2"]["comb_multi_properly_paired"]
            ),
            "Multi singletons": len(data["minimap2"]["comb_multi_singleton"]),
            "Unmapped": len(data["minimap2"]["comb_unmapped_flagged"]),
        },
        #'GTAMap': {
        #    'Total mapped': len(data['gtamap']['comb_all_reads_in_sam']) - len(data['gtamap']['comb_unmapped_flagged']),
        #    'Properly paired': len(data['gtamap']['comb_unique_properly_paired']),
        #    'Singletons': len(data['gtamap']['comb_unique_singleton']),
        #    'Multi properly paired': len(data['gtamap']['comb_multi_properly_paired']),
        #    'Multi singletons': len(data['gtamap']['comb_multi_singleton']),
        #    'Unmapped': len(data['gtamap']['comb_unmapped_flagged'])
        # }
        "GTAMap": {
            "Total mapped": len(data["gtamap"]["comb_all_reads_in_sam"])
            - len(data["gtamap"]["comb_unmapped_flagged"]),
            "Unique mapped": len(data["gtamap"]["comb_unique_singleton"]),
            "Multi mapped": len(data["gtamap"]["comb_multi_singleton"]),
        },
    }

    # --- Plotting ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes = axes.flatten()
    categories = [
        "Total mapped",
        "Properly paired",
        "Singletons",
        "Multi properly paired",
        "Multi singletons",
    ]
    categories_gta = ["Total mapped", "Unique mapped", "Multi mapped"]

    colors = sns.color_palette("husl", n_colors=len(categories) + 7)[7:]
    gta_colors = sns.color_palette("husl", n_colors=len(categories_gta) + 8)[8:]

    i = 0
    for idx, (mapper_name, data_) in enumerate(results.items()):
        ax = axes[idx]
        if mapper_name != "GTAMap":
            main_values = [data_[c] for c in categories]
            bars1 = ax.bar(
                categories, main_values, color=colors, edgecolor="black", linewidth=1.1
            )

            for bar in bars1:
                height = bar.get_height()
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    height,
                    f"{int(height)}",
                    ha="center",
                    va="bottom",
                    fontsize=13,
                    fontweight="bold",
                )

            ax.set_title(mapper_name, fontsize=14, fontweight="bold")
            if i == 0 or i == 2:
                ax.set_ylabel("Number of \n Mapped QNAMES", fontsize=14)
            else:
                ax.set_ylabel("")
            ax.grid(axis="y", alpha=0.3, linestyle="--")
            ax.set_xticklabels(categories, rotation=90)
            ax.set_ylim(0, 500000)
            ax.tick_params(axis="x", labelsize=15)
            ax.tick_params(axis="y", labelsize=15)
            i += 1
        else:
            main_values = [data_[c] for c in categories_gta]
            bars1 = ax.bar(
                categories_gta,
                main_values,
                color=gta_colors,
                edgecolor="black",
                linewidth=1.1,
            )

            for bar in bars1:
                height = bar.get_height()
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    height,
                    f"{int(height)}",
                    ha="center",
                    va="bottom",
                    fontsize=13,
                    fontweight="bold",
                )

            ax.set_title(mapper_name, fontsize=14, fontweight="bold")
            if i == 0 or i == 2:
                ax.set_ylabel("Number of \n Mapped QNAMES", fontsize=14)
            else:
                ax.set_ylabel("")
            ax.grid(axis="y", alpha=0.3, linestyle="--")
            ax.set_xticklabels(categories_gta, rotation=90)
            ax.set_ylim(0, 500000)
            ax.tick_params(axis="x", labelsize=15)
            ax.tick_params(axis="y", labelsize=15)
            i += 1

    plt.suptitle(
        "Combined Gene-wise Read Mapping Statistics by Mapper \n on Subset of 500,000 Read-pairs",
        fontsize=22,
        fontweight="bold",
        y=1.02,
    )
    sns.despine()
    plt.tight_layout()
    plt.savefig("mapping_statistics_paired_genes.png", dpi=300, bbox_inches="tight")
    from upsetplot import UpSet, from_contents

    # --- UpSet Plot for overlap of total mapped read IDs ---
    mapped_sets = {
        "HISAT2": data["hisat2"]["comb_all_reads_in_sam"],
        "STAR": data["star"]["comb_all_reads_in_sam"],
        "Minimap2": data["minimap2"]["comb_all_reads_in_sam"],
        "GTAMap": data["gtamap"]["comb_all_reads_in_sam"],
    }

    # Build the UpSet input data
    upset_data = from_contents(mapped_sets)

    # Plot
    fig = plt.figure(figsize=(16, 8))
    upset = UpSet(
        upset_data,
        subset_size="count",
        show_counts=True,
        sort_by="degree",  # group intersections by how many sets overlap
    )
    upset.plot()
    # Reduce text size of built-in labels
    for text in plt.gcf().findobj(match=plt.Text):
        if text.get_text().isdigit():
            text.set_fontsize(8)

    plt.suptitle(
        "Overlap of Total Mapped QNAMEs Across Mappers and Genes",
        fontsize=14,
        fontweight="bold",
    )
    plt.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1)
    plt.savefig("upset_total_mapped_qnames.png", dpi=300, bbox_inches="tight")
    plt.show()

    gtamap_set = data["gtamap"]["comb_all_reads_in_sam"]
    other_mappers = (
        data["hisat2"]["comb_all_reads_in_sam"]
        | data["star"]["comb_all_reads_in_sam"]
        | data["minimap2"]["comb_all_reads_in_sam"]
    )

    unique_gtamap_reads = gtamap_set - other_mappers  # reads only in GTAMap

    output_file = "gtamap_unique_reads.txt"
    with open(output_file, "w") as f:
        for read_id in sorted(unique_gtamap_reads):
            f.write(read_id + "\n")

    print(f"Written {len(unique_gtamap_reads)} GTAMap-unique reads to {output_file}")

    minimap_set = data["minimap2"]["comb_all_reads_in_sam"]
    other_mappers = (
        data["hisat2"]["comb_all_reads_in_sam"]
        | data["star"]["comb_all_reads_in_sam"]
        | data["gtamap"]["comb_all_reads_in_sam"]
    )

    unique_m_reads = minimap_set - other_mappers

    output_file = "minimap_unique_reads.txt"
    with open(output_file, "w") as f:
        for read_id in sorted(unique_m_reads):
            f.write(read_id + "\n")


if __name__ == "__main__":
    # gene_wise_parallel() # for gene wise comp
    genome_wide()
