#!/usr/bin/env python3

import argparse
import pysam
import numpy as np


def get_depth_array(bam, chrom, chrom_len, min_baseq=0):
    """
    Return total depth array for one chromosome.
    """
    cov = bam.count_coverage(
        chrom,
        start=0,
        end=chrom_len,
        quality_threshold=min_baseq
    )

    depth = (
        np.array(cov[0], dtype=np.int32) +
        np.array(cov[1], dtype=np.int32) +
        np.array(cov[2], dtype=np.int32) +
        np.array(cov[3], dtype=np.int32)
    )

    return depth


def merge_low_depth_regions(mask, max_gap):
    """
    Find True regions in boolean mask.
    Allow gaps <= max_gap to be merged.
    """

    regions = []
    n = len(mask)

    i = 0

    while i < n:

        if not mask[i]:
            i += 1
            continue

        start = i
        last_low = i
        gap = 0
        i += 1

        while i < n:

            if mask[i]:
                last_low = i
                gap = 0
            else:
                gap += 1

                if gap > max_gap:
                    break

            i += 1

        end = last_low + 1

        regions.append((start, end))

    return regions


def flank_mean(depth, start, end, flank_size):
    """
    Mean depth in left and right flanks.
    """

    n = len(depth)

    left_start = max(0, start - flank_size)
    left_end = start

    right_start = end
    right_end = min(n, end + flank_size)

    if left_end > left_start:
        left_mean = np.mean(depth[left_start:left_end])
    else:
        left_mean = np.nan

    if right_end > right_start:
        right_mean = np.mean(depth[right_start:right_end])
    else:
        right_mean = np.nan

    return left_mean, right_mean


def main():

    parser = argparse.ArgumentParser(
        description="Detect deletion-like low-coverage regions from bacterial BAM files."
    )

    parser.add_argument(
        "-b", "--bam",
        required=True,
        help="Input sorted and indexed BAM"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file"
    )

    parser.add_argument(
        "--min-length",
        type=int,
        default=100,
        help="Minimum deletion length in bp [default: 100]"
    )

    parser.add_argument(
        "--depth-ratio",
        type=float,
        default=0.10,
        help="Maximum depth relative to genome median depth [default: 0.10]"
    )

    parser.add_argument(
        "--max-gap",
        type=int,
        default=20,
        help="Maximum internal non-low-depth gap allowed [default: 20 bp]"
    )

    parser.add_argument(
        "--flank-size",
        type=int,
        default=200,
        help="Flanking region size used for validation [default: 200 bp]"
    )

    parser.add_argument(
        "--min-flank-ratio",
        type=float,
        default=0.50,
        help="Minimum flank depth relative to genome median [default: 0.50]"
    )

    parser.add_argument(
        "--min-median-depth",
        type=float,
        default=10,
        help="Minimum whole-genome median depth required [default: 10]"
    )

    parser.add_argument(
        "--min-baseq",
        type=int,
        default=0,
        help="Minimum base quality for depth calculation [default: 0]"
    )

    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")

    genome_depths = []

    depth_cache = {}

    print("Calculating depth...")

    for chrom, chrom_len in zip(bam.references, bam.lengths):

        depth = get_depth_array(
            bam,
            chrom,
            chrom_len,
            args.min_baseq
        )

        depth_cache[chrom] = depth

        genome_depths.append(depth)

    all_depth = np.concatenate(genome_depths)

    nonzero_depth = all_depth[all_depth > 0]

    if len(nonzero_depth) == 0:
        raise RuntimeError("No mapped coverage detected.")

    genome_median = np.median(nonzero_depth)

    print(f"Genome median depth: {genome_median:.2f}")

    if genome_median < args.min_median_depth:
        print(
            f"WARNING: median depth {genome_median:.2f} "
            f"is below recommended threshold "
            f"{args.min_median_depth}"
        )

    deletion_depth_threshold = genome_median * args.depth_ratio
    flank_threshold = genome_median * args.min_flank_ratio

    print(
        f"Deletion depth threshold: "
        f"{deletion_depth_threshold:.2f}"
    )

    print(
        f"Minimum flank depth: "
        f"{flank_threshold:.2f}"
    )

    results = []

    for chrom in bam.references:

        depth = depth_cache[chrom]

        low_mask = depth <= deletion_depth_threshold

        regions = merge_low_depth_regions(
            low_mask,
            args.max_gap
        )

        for start, end in regions:

            length = end - start

            if length < args.min_length:
                continue

            region_depth = depth[start:end]

            mean_depth = np.mean(region_depth)
            median_depth = np.median(region_depth)
            min_depth = np.min(region_depth)
            max_depth = np.max(region_depth)

            left_mean, right_mean = flank_mean(
                depth,
                start,
                end,
                args.flank_size
            )

            left_ok = (
                not np.isnan(left_mean)
                and left_mean >= flank_threshold
            )

            right_ok = (
                not np.isnan(right_mean)
                and right_mean >= flank_threshold
            )

            flank_status = "PASS" if left_ok and right_ok else "FAIL"

            relative_depth = mean_depth / genome_median

            results.append(
                (
                    chrom,
                    start + 1,
                    end,
                    length,
                    mean_depth,
                    median_depth,
                    min_depth,
                    max_depth,
                    relative_depth,
                    left_mean,
                    right_mean,
                    flank_status
                )
            )

    with open(args.output, "w") as out:

        header = [
            "chrom",
            "start",
            "end",
            "length",
            "mean_depth",
            "median_depth",
            "min_depth",
            "max_depth",
            "relative_depth",
            "left_flank_mean",
            "right_flank_mean",
            "flank_status"
        ]

        out.write("\t".join(header) + "\n")

        for r in results:

            out.write(
                f"{r[0]}\t"
                f"{r[1]}\t"
                f"{r[2]}\t"
                f"{r[3]}\t"
                f"{r[4]:.2f}\t"
                f"{r[5]:.2f}\t"
                f"{r[6]}\t"
                f"{r[7]}\t"
                f"{r[8]:.4f}\t"
                f"{r[9]:.2f}\t"
                f"{r[10]:.2f}\t"
                f"{r[11]}\n"
            )

    bam.close()

    print(f"Done.")
    print(f"Detected {len(results)} candidate regions.")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
