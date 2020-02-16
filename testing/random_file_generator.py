#!/usr/bin/env python3

import argparse
import random
import sys


def random_0125_map(
    n,
    outfile=sys.stdout,
    seed=None,
    chromosomes=[str(i) for i in range(1, 30)] + ["X", "Y"],
    min_pos=100000,
    max_pos=999999,
    start_from_id=1,
):
    r = random.Random(seed)

    outfile.write("Name\tChromosome\tPosition\n")
    outfile.writelines(
        (
            f"{'SNP' + str(start_from_id + i)}\t"
            + f"{r.choice(chromosomes)}\t"
            + f"{r.randint(min_pos, max_pos)}\n"
            for i in range(n)
        )
    )


def random_0125_samples(n, map_size, outfile=sys.stdout, seed=None, start_from_id=1):
    r = random.Random(seed)
    outfile.write("ID\tCALL...\n")
    outfile.writelines(
        (
            f"{'SAM' + str(start_from_id + i)}\t"
            + f"{''.join(r.choices('0125', k=map_size))}\n"
            for i in range(n)
        )
    )


def random_plink_map(
    n,
    outfile=sys.stdout,
    seed=None,
    chromosomes=[str(i) for i in range(1, 30)] + ["X", "Y"],
    min_pos=100000,
    max_pos=999999,
    min_dist=0,
    max_dist=0,
    start_from_id=1,
):
    r = random.Random(seed)

    outfile.writelines(
        (
            f"{r.choice(chromosomes)} "
            + f"{'SNP' + str(start_from_id + i)} "
            + f"{r.randint(min_dist, max_dist)} "
            + f"{r.randint(min_pos, max_pos)}\n"
        )
        for i in range(n)
    )


def random_plink_samples(
    n, map_size, outfile=sys.stdout, seed=None, bases="ATCG", start_from_id=1
):
    r = random.Random(seed)

    outfile.writelines(
        (
            f"0 {'SAM' + str(start_from_id + i)} 0 0 "
            + f"{r.choice('12')} 0 "
            + f"{' '.join(r.choices(bases, k=2*map_size))}\n"
            for i in range(n)
        )
    )


def random_final_report(
    n,
    map_size,
    outfile=sys.stdout,
    seed=None,
    start_snps_from_id=1,
    start_samples_from_id=1,
):
    r = random.Random(seed)
    TAB = "\t"
    outfile.write("[Header]\n[Data]\n")
    outfile.write(
        "SNP Name\tSample ID\tAllele1 - Forward\t"
        + "Allele2 - Forward\tAllele1 - Top\tAllele2 - Top\t"
        + "Allele1 - AB\tAllele2 - AB\tGC Score\tX\tY\n"
    )
    outfile.writelines(
        (
            f"{'SNP' + str(start_snps_from_id + i%map_size)}\t"
            + f"{'SAM' + str(start_samples_from_id + i//map_size)}\t"
            + f"{TAB.join(r.choices('ATCG', k=4))}\t"
            + f"{TAB.join(r.choices('AB', k=2))}\t"
            + f"{r.random():.4f}\t"
            + f"{r.random():.3f}\t"
            + f"{r.random():.3f}\n"
            for i in range(map_size * n)
        )
    )


def random_vcf(
    n,
    map_size,
    outfile=sys.stdout,
    seed=None,
    chromosomes=[str(i) for i in range(1, 30)] + ["X", "Y"],
    bases="ATCG",
    min_pos=100000,
    max_pos=999999,
    min_qual=1,
    max_qual=100,
    start_snps_from_id=1,
    start_samples_from_id=1,
):
    r = random.Random(seed)
    outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    outfile.writelines(
        [f"\t{'SAM' + str(start_samples_from_id + i)}" for i in range(n)]
    )
    outfile.write("\n")
    outfile.writelines(
        (
            f"{r.choice(chromosomes)}\t"
            + f"{r.randint(min_pos, max_pos)}\t"
            + f"{'SNP' + str(start_snps_from_id + i)}\t"
            + f"{r.choice(bases)}\t"
            + f"{r.choice(bases)}\t"
            + f"{r.randint(min_qual, max_qual)}\t"
            + f"PASS\t"
            + f"RND\t"
            + f"GT\t"
            + "\t".join((f"{r.randint(0, 1)}|{r.randint(0, 1)}" for j in range(n)))
            + "\n"
            for i in range(map_size)
        )
    )


def id_mapping(n, outfile=sys.stdout, first_sample_id=1):
    outfile.writelines(
        (
            f"{'SAM' + str(first_sample_id + i)}\t" + f"{'IND' + str(i+1)}\n"
            for i in range(n)
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "format",
        help="file format to generate",
        choices=["vcf", "fr", "0125map", "0125ped", "plmap", "plped", "ids"],
    )
    parser.add_argument(
        "-k", help="number of files to generate, default 1", type=int, default=1
    )
    parser.add_argument(
        "-n", type=int, default=0, help="number of samples in each file, if applicable"
    )
    parser.add_argument("-m", type=int, help="number of SNPs in each file")
    parser.add_argument(
        "--first-snp-id",
        type=int,
        default=1,
        help="number on the ID of the first snp "
        + "(following snps will have sequential IDs)"
        + ", default 1",
    )
    parser.add_argument(
        "--first-sample-id",
        type=int,
        default=1,
        help="number on the ID of the first sample "
        + "(following samples will have sequential IDs)"
        + ", default 1",
    )
    parser.add_argument("-o", help="Output file prefix", required=True)
    args = parser.parse_args()
    digits = len(str(args.k))

    for i in range(args.k):
        filename = args.o + "_" + str(i + 1).zfill(digits) + "." + args.format
        print(f"Generating {filename}...")
        with open(filename, "w+") as f:
            if args.format == "vcf":
                random_vcf(
                    args.n,
                    args.m,
                    outfile=f,
                    start_snps_from_id=args.first_snp_id + i * args.m,
                    start_samples_from_id=args.first_sample_id + i * args.n,
                )
            elif args.format == "fr":
                random_final_report(
                    args.n,
                    args.m,
                    outfile=f,
                    start_snps_from_id=args.first_snp_id + i * args.m,
                    start_samples_from_id=args.first_sample_id + i * args.n,
                )
            elif args.format == "0125map":
                random_0125_map(
                    args.m, outfile=f, start_from_id=args.first_snp_id + i * args.m
                )
            elif args.format == "plmap":
                random_plink_map(
                    args.m, outfile=f, start_from_id=args.first_snp_id + i * args.m
                )
            elif args.format == "0125ped":
                random_0125_samples(
                    args.n,
                    args.m,
                    outfile=f,
                    start_from_id=args.first_sample_id + i * args.n,
                )
            elif args.format == "plped":
                random_plink_samples(
                    args.n,
                    args.m,
                    outfile=f,
                    start_from_id=args.first_sample_id + i * args.n,
                )
            elif args.format == "ids":
                id_mapping(
                    args.n, outfile=f, first_sample_id=args.first_sample_id + i * args.n
                )
            else:
                raise Exception("Unknown format.")
