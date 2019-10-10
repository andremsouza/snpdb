#!/usr/bin/env python3
import argparse
import snpdb
import time

from readers import Z125MapReader, PlinkMapReader, FinalReportMapReader, VcfMapReader

Z125 = 0
PLINK = 1
ILMFR = 2
VCF = 3

READERS = [Z125MapReader, PlinkMapReader,
          FinalReportMapReader, VcfMapReader]

def import_map(filename, fileformat, mapname, **kwargs):
   reader = READERS[fileformat](filename)
   snpdb.import_map(reader, mapname, **kwargs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("format", help="map file format",
                        choices=["vcf", "fr", "z125", "pl"])
    parser.add_argument("mapfile", help="map file path")
    parser.add_argument("mapname", help="id of the map to be stored")
    parser.add_argument("-q", "--quiet", help="omit all output",
                        action="store_true")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--force-create-new", help="always create new snps",
                       action="store_true")
    group.add_argument("--force-use-existing", help="always use existing snps " + 
                       "when available (may ask user to decide which)",
                       action="store_true")
    args = parser.parse_args()

    if args.format == "vcf":
        fmt = VCF
    elif args.format == "fr":
        fmt = ILMFR
    elif args.format == "z125":
        fmt = Z125
    else:
        fmt = PLINK
    
    report = not args.quiet

    start = time.time()

    import_map(args.mapfile, fmt, args.mapname,
               force_create_new=args.force_create_new,
               force_use_existing=args.force_use_existing,
               report=report)
    
    if report:
        print(f"Done in {time.time() - start:.3f} s.")
