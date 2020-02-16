#!/usr/bin/env python3
import argparse
import snpdb
import time

from pprint import pprint
from readers import *
from writers import *


_FORMAT_CHOICES = ["0125", "pl", "fr", "vcf"]
_MAP_READERS = [Z125MapReader, PlinkMapReader, FinalReportMapReader, VcfMapReader]
_SAMPLE_READERS = [
    Z125SampleReader,
    PlinkSampleReader,
    FinalReportSampleReader,
    VcfSampleReader,
]

_MAP_WRITERS = [Z125MapWriter, PlinkMapWriter]
_SAMPLE_WRITERS = [Z125SampleWriter, PlinkSampleWriter]


def import_map(filename, fmt, mapname, **kwargs):
    reader = _MAP_READERS[_FORMAT_CHOICES.index(fmt)](filename)
    snpdb.import_map(reader, mapname, **kwargs)


def import_samples(filename, fmt, mapname, idfilename=None, **kwargs):
    id_map = {}
    if idfilename is not None:
        with open(idfilename, "r") as f:
            for line in f:
                (sample, individual) = line.split()
                id_map[sample] = individual
    reader = _SAMPLE_READERS[_FORMAT_CHOICES.index(fmt)](filename)
    snpdb.import_samples(reader, mapname, id_map=id_map, **kwargs)


def export_map(mapname, fmt, out_file_path):
    writer = _MAP_WRITERS[_FORMAT_CHOICES.index(fmt)]
    snpdb.export_map(mapname, writer, out_file_path)


def export_samples(samples, map, fmt, out_file_path):
    writer = _SAMPLE_WRITERS[_FORMAT_CHOICES.index(fmt)]
    snpdb.export_samples(samples, map, writer, out_file_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")

    # import-map
    p = subparsers.add_parser(
        "import-map", help="import map from proper " + "map or ped files"
    )

    p.add_argument("format", help="map file format", choices=_FORMAT_CHOICES)
    p.add_argument("mapfile", help="map file path")
    p.add_argument("mapname", help="id of the map to be stored")
    p.add_argument("-q", "--quiet", help="omit all output", action="store_true")
    group = p.add_mutually_exclusive_group()
    group.add_argument(
        "--force-create-new", help="always create new snps", action="store_true"
    )
    group.add_argument(
        "--force-use-existing",
        help="always use existing  "
        + "snps when available (may ask user to decide which)",
        action="store_true",
    )

    # import-sample
    p = subparsers.add_parser("import-samples", help="import samples from file")
    p.add_argument("format", help="map file format", choices=_FORMAT_CHOICES)
    p.add_argument("samplefile", help="samples file path")
    p.add_argument(
        "mapname", help="id of the map to be used " + "(must exist in the database)"
    )
    p.add_argument("-q", "--quiet", help="omit all output", action="store_true")
    p.add_argument(
        "--idfile",
        help="specify an id mapping (associates "
        + "individual ids to the samples), "
        + "file should contain sample id on the "
        + "first column and individual id onthe "
        + "second",
    )

    # find-snps
    p = subparsers.add_parser("find-snps", help="search snps in the database")
    p.add_argument("--id", help="match numeric internal id exactly", type=int)
    p.add_argument("--name", help="match name exactly")
    p.add_argument("--chr", help="match chromosome exactly")
    p.add_argument(
        "--min-chr",
        help="match only snps whose chromosome is " + "at least MIN_CHR (integer)",
        type=int,
    )
    p.add_argument(
        "--max-chr",
        help="match only snps whose chromosome is " + "at most MAX_CHR (integer)",
        type=int,
    )
    p.add_argument(
        "--min-pos",
        help="match only snps whose position is " + "at least MIN_POS",
        type=int,
    )
    p.add_argument(
        "--max-pos",
        help="match only snps whose position is " + "at most MAX_POS",
        type=int,
    )
    p.add_argument("--map", help="match snps that belong to MAP")

    # find-maps
    p = subparsers.add_parser("find-maps", help="search maps in the database")
    p.add_argument("--name", help="match name exactly")
    p.add_argument("--format", help="match format string exactly")
    p.add_argument(
        "--min-size", help="match only maps of size " + "at least MIN_SIZE", type=int
    )
    p.add_argument(
        "--max-size", help="match only maps of size " + "at most MAX_SIZE", type=int
    )

    # find-individuals
    p = subparsers.add_parser(
        "find-individuals", help="search individuals in the database"
    )
    p.add_argument(
        "--name",
        help="match at least one of the individuals' " + "names (tatoos) exactly",
    )
    p.add_argument(
        "--sample",
        help="match only individuals which have a " + "sample with this specific id",
    )
    p.add_argument(
        "--map",
        help="match only individuals which have a " + "sample with this specific map",
    )

    # find-samples
    p = subparsers.add_parser("find-samples", help="search samples in the database")
    p.add_argument("--id", help="match sample id within map")
    p.add_argument("--map", help="match map the sample belongs to")

    # get-snp-genotype
    p = subparsers.add_parser(
        "get-snp-genotype", help="retrive the genotype given a sample " + "and a SNP"
    )
    p.add_argument("map", help="map to which the sample belongs")
    p.add_argument("sample", help="id of the sample")
    p.add_argument("snp", help="internal id of the SNP", type=int)

    # put-file
    p = subparsers.add_parser("put-file", help="upload file to the database")
    p.add_argument("file", help="path of file to upload", nargs="+")
    p.add_argument(
        "--individual",
        help="internal id of an individual" + " to be associated with the file",
    )

    # find-files
    p = subparsers.add_parser("find-files", help="search files in the database")
    p.add_argument(
        "--individual",
        help="match only files associated with "
        + "individual whose internal id is INDIVIDUAL",
    )
    p.add_argument("--name", help="match file name exactly")

    # get-files
    p = subparsers.add_parser("get-files", help="download files from the database")
    p.add_argument(
        "--individual",
        help="download only files associated "
        + "with individual whose internal id is INDIVIDUAL",
    )
    p.add_argument("--name", help="download only files with the specified NAME")

    # export-map
    p = subparsers.add_parser("export-map", help="export map from database to file")
    p.add_argument("format", help="format of the output file", choices=["0125", "pl"])
    p.add_argument("map", help="name of the map to export")
    p.add_argument("outfile", help="path of the file to export to")

    # export-samples
    p = subparsers.add_parser(
        "export-samples", help="export samples from database to file"
    )
    p.add_argument("format", help="format of the output file", choices=["0125", "pl"])
    p.add_argument("map", help="name of the map to which the samples belong")
    p.add_argument("outfile", help="path of the file to export to")
    p.add_argument(
        "sample",
        help="ids of the samples to export, if none "
        + "is specified, will export every sample "
        + "within map",
        nargs="*",
    )

    # summarize
    p = subparsers.add_parser(
        "summarize",
        help="search individuals in the database and summarize all data on each match",
    )
    p.add_argument(
        "--name",
        help="match at least one of the individuals' " + "names (tatoos) exactly",
    )
    p.add_argument(
        "--sample",
        help="match only individuals which have a " + "sample with this specific id",
    )
    p.add_argument(
        "--map",
        help="match only individuals which have a " + "sample with this specific map",
    )

    args = parser.parse_args()
    if args.subcommand == "import-map":
        report = not args.quiet
        start = time.time()
        import_map(
            args.mapfile,
            args.format,
            args.mapname,
            force_create_new=args.force_create_new,
            force_use_existing=args.force_use_existing,
            report=report,
        )
        print(f"Done in {time.time() - start:.3f} s.")
    elif args.subcommand == "import-samples":
        report = not args.quiet
        start = time.time()
        import_samples(
            args.samplefile,
            args.format,
            args.mapname,
            idfilename=args.idfile,
            report=report,
        )
        print(f"Done in {time.time() - start:.3f} s.")
    elif args.subcommand == "find-snps":
        for snp in snpdb.find_snp(
            args.name,
            args.min_chr,
            args.max_chr,
            args.min_pos,
            args.max_pos,
            args.map,
            args.id,
            args.chr,
        ):
            print(snp)
    elif args.subcommand == "find-maps":
        for map in snpdb.find_maps(
            args.name, args.min_size, args.max_size, args.format
        ):
            print(map)
    elif args.subcommand == "find-individuals":
        for ind in snpdb.find_individuals(None, args.name, args.map, args.sample):
            print(ind)
    elif args.subcommand == "find-samples":
        for sample in snpdb.find_sample(args.id, args.map):
            print(sample)
    elif args.subcommand == "get-snp-genotype":
        print(snpdb.find_snp_of_sample(args.map, args.sample, args.snp))
    elif args.subcommand == "put-file":
        for fname in args.file:
            with open(fname, "rb") as f:
                snpdb.insert_file(f, args.individual)
    elif args.subcommand == "find-files":
        for file in snpdb.list_files(args.individual, args.name):
            print(file)
    elif args.subcommand == "get-files":
        snpdb.get_files(snpdb.list_files(args.individual, args.name))
    elif args.subcommand == "export-map":
        export_map(args.map, args.format, args.outfile)
    elif args.subcommand == "export-samples":
        export_samples(args.sample, args.map, args.format, args.outfile)
    elif args.subcommand == "summarize":
        for ind in snpdb.find_individuals(None, args.name, args.map, args.sample):
            pprint(snpdb.summarize(ind))
    elif args.subcommand is None:
        print("Subcommand required. Use -h for help.")
