#!/usr/bin/env python3

"""
Converts RepeatMasker .out files to BED format for genomic analysis.

This script parses RepeatMasker output files (which contain information about repetitive
DNA sequences) and converts them to the standard BED format. It extracts chromosome/contig
coordinates, repeat types, and provides filtering options for specific repeat types and
chromosomes. The output is tab-separated with columns: chromosome, start, end, repeat_type.

Parameters:
    -f, --file: Path to RepeatMasker .out file (required)
    -t, --type: Filter by specific repeat type (column 11). Can be used multiple times
    --replace-chr: Remove "chr" prefix from chromosome names
    --use-regular-chr: Only include standard human chromosomes (1-22, X, Y)
"""

import argparse
import gzip
from typing import Set


def parse_file(file: str, type_filters: Set, replace_chr: bool, use_regular_chr: bool):

    regular_chromosomes = {str(i) for i in range(1, 23)}.union({"X", "Y"})

    # check if file is gzipped
    if file.endswith(".gz"):
        fh = gzip.open(file, "rt")
    else:
        fh = open(file, "r")

    # skip header lines
    for _ in range(3):
        next(fh)

    for line in fh:
        line = line.strip()

        line_list = line.split()

        if len(line_list) < 11:
            continue

        contig = line_list[4]
        contig_no_chr = contig[3:] if contig.startswith("chr") else contig

        if use_regular_chr and contig_no_chr not in regular_chromosomes:
            continue

        if replace_chr:
            contig = contig_no_chr

        start = int(line_list[5]) - 1
        end = int(line_list[6])
        repeat_type = line_list[10]

        if len(type_filters) > 0 and repeat_type not in type_filters:
            continue

        print(f"{contig}\t{start}\t{end}\t{repeat_type}")

    fh.close()


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description="Convert RepeatMasker .out file to BED format"
    )

    argparser.add_argument(
        "-f",
        "--file",
        type=str,
        dest="file",
        help="RepeatMasker .out file",
        required=True,
    )
    argparser.add_argument(
        "-t",
        "--type",
        action="append",
        help="Use only this repeat type (column #11). Can be specified multiple times.",
        default=[],
    )
    argparser.add_argument(
        "--replace-chr",
        action="store_true",
        help='Replace "chr" prefix in contig names',
        default=False,
    )
    argparser.add_argument(
        "--use-regular-chr",
        action="store_true",
        help="Use only the regular human chromosomes (1-22XY)",
        default=False,
    )

    args = argparser.parse_args()

    type_filters = set(args.type) if args.type else set()

    parse_file(args.file, type_filters, args.replace_chr, args.use_regular_chr)
