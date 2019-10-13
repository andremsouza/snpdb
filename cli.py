#!/usr/bin/env python3
import argparse
import snpdb
import time

from readers import *

_FORMAT_CHOICES = ["0125", "pl", "fr", "vcf"]
_MAP_READERS = [Z125MapReader, PlinkMapReader,
          FinalReportMapReader, VcfMapReader]
_SAMPLE_READERS = [Z125SampleReader, PlinkSampleReader,
                  FinalReportSampleReader, VcfSampleReader]




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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")
   

    # import-map
    p = subparsers.add_parser("import-map",
                              help="import map from proper " +
                                   "map or ped files")

    p.add_argument("format", help="map file format", choices=_FORMAT_CHOICES)
    p.add_argument("mapfile", help="map file path")
    p.add_argument("mapname", help="id of the map to be stored")
    p.add_argument("-q", "--quiet", help="omit all output",
                        action="store_true")
    group = p.add_mutually_exclusive_group()
    group.add_argument("--force-create-new", help="always create new snps",
                       action="store_true")
    group.add_argument("--force-use-existing", help="always use existing  " +
                       "snps when available (may ask user to decide which)",
                       action="store_true")


    # import-sample
    p = subparsers.add_parser("import-samples",
                              help="import samples from file")
    p.add_argument("format", help="map file format", choices=_FORMAT_CHOICES)
    p.add_argument("samplefile", help="samples file path")
    p.add_argument("mapname", help="id of the map to be used " + 
                                   "(must exist in the database)")
    p.add_argument("-q", "--quiet", help="omit all output",
                   action="store_true")
    p.add_argument("--idfile", help="specify an id mapping (associates " +
                                    "individual ids to the samples), "+ 
                                    "file should contain sample id on the " +
                                    "first column and individual id onthe " +
                                    "second")


    # find-snps
    p = subparsers.add_parser("find-snps",
                             help="search snps in the database")
    p.add_argument("--name", help="match name exactly")
    p.add_argument("--min-chr", help="match only snps whose chromosome is " + 
                                     "at least MIN_CHR (lexicographically)")
    p.add_argument("--max-chr", help="match only snps whose chromosome is " + 
                                     "at most MAX_CHR (lexicographically)")
    p.add_argument("--min-pos", help="match only snps whose position is " + 
                                     "at least MIN_POS",
                   type=int)
    p.add_argument("--max-pos", help="match only snps whose position is " + 
                                     "at most MAX_POS",
                   type=int) 


    # find-maps
    p = subparsers.add_parser("find-maps",
                             help="search maps in the database")
    p.add_argument("--name", help="match name exactly")
    p.add_argument("--format", help="match format string exactly")
    p.add_argument("--min-size", help="match only maps of size " + 
                                     "at least MIN_SIZE",
                   type=int)
    p.add_argument("--max-size", help="match only maps of size " + 
                                     "at most MAX_SIZE",
                   type=int)
    
    # find-individuals
    p = subparsers.add_parser("find-individuals",
                             help="search individuals in the database")
    p.add_argument("--name", help="match at least one of the individuals' " +
                                  "names (tatoos) exactly")
    p.add_argument("--sample", help="match only individuals which have a " + 
                                    "sample with this specific id")
    p.add_argument("--map", help="match only individuals which have a " +
                                 "sample with this specific map")


    args = parser.parse_args()
    if args.subcommand == "import-map":
        report = not args.quiet
        start = time.time()
        import_map(args.mapfile, args.format, args.mapname,
                   force_create_new=args.force_create_new,
                   force_use_existing=args.force_use_existing,
                   report=report)
        print(f"Done in {time.time() - start:.3f} s.")
    elif args.subcommand == "import-samples":
        report = not args.quiet
        start = time.time()
        import_samples(args.samplefile, args.format, args.mapname,
                       idfilename=args.idfile,
                       report=report)
        print(f"Done in {time.time() - start:.3f} s.")
    elif args.subcommand == "find-snps":
        for snp in snpdb.find_snp(args.name, args.min_chr, args.max_chr,
                                  args.min_pos, args.max_pos):
            print(snp)
    elif args.subcommand == "find-maps":
        for map in snpdb.find_maps(args.name, args.min_size, args.max_size,
                                   args.format):
            print(map)
    elif args.subcommand == "find-individuals":
        for ind in snpdb.find_individuals(None, args.name, args.map, args.sample):
            print(ind)
    elif args.subcommand is None:
        print("Subcommand required. Use -h for help.")
