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
    import_map_parser = subparsers.add_parser("import-map",
                                              help="import map from proper " +
                                                   "map or ped files")

    import_map_parser.add_argument("format", help="map file format",
                        choices=_FORMAT_CHOICES)
    import_map_parser.add_argument("mapfile", help="map file path")
    import_map_parser.add_argument("mapname", help="id of the map to be stored")
    import_map_parser.add_argument("-q", "--quiet", help="omit all output",
                        action="store_true")
    group = import_map_parser.add_mutually_exclusive_group()
    group.add_argument("--force-create-new", help="always create new snps",
                       action="store_true")
    group.add_argument("--force-use-existing", help="always use existing  " +
                       "snps when available (may ask user to decide which)",
                       action="store_true")


    # import-sample
    import_sample_parser = subparsers.add_parser("import-samples",
                                                 help="import samples from file")
    import_sample_parser.add_argument("format", help="map file format",
                        choices=_FORMAT_CHOICES)
    import_sample_parser.add_argument("samplefile", help="samples file path")
    import_sample_parser.add_argument("mapname", help="id of the map to be used"+ 
                                        "(must exist in the database)")
    import_sample_parser.add_argument("-q", "--quiet", help="omit all output",
                        action="store_true")
    

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
                       report=report)     
        print(f"Done in {time.time() - start:.3f} s.")
    elif args.subcommand is None:
        print("Subcommand required. Use -h for help.")
