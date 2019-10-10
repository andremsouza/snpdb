#!/usr/bin/env python3
import argparse
import snpdb
import time

from readers import Z125SampleReader, PlinkSampleReader, FinalReportSampleReader, VcfSampleReader

Z125 = 0
PLINK = 1
ILMFR = 2
VCF = 3

READERS = [Z125SampleReader, PlinkSampleReader,
          FinalReportSampleReader, VcfSampleReader]

def import_samples(filename, fileformat, mapname, idfilename=None, **kwargs):
    id_map = {}
    if idfilename is not None:
        with open(idfilename, "r") as f:
            for line in f:
                (sample, individual) = line.split()
                id_map[sample] = individual
    reader = READERS[fileformat](filename)
    snpdb.import_samples(reader, mapname, id_map=id_map, **kwargs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("format", help="map file format",
                        choices=["vcf", "fr", "z125", "pl"])
    parser.add_argument("samplefile", help="samples file path")
    parser.add_argument("mapname", help="id of the map to be used"+ 
                                        "(must exist in the database)")
    parser.add_argument("-q", "--quiet", help="omit all output",
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

    import_samples(args.samplefile, fmt, args.mapname,
                   report=report)
    
    if report:
        print(f"Done in {time.time() - start:.3f} s.")
