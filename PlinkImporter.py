from AbstractImporter import AbstractImporter

class PlinkImporter(AbstractImporter):
    
    def read_map(self, filename):
        snps = []
        with open(filename, "r") as f:
            for line in f:
                (chrom, id, dist, pos) = line.split()
                current = {
                    self.SNPS_NAME_ATTR: id,
                    self.SNPS_CHROMOSOME_ATTR: chrom,
                    self.SNPS_POSITION_ATTR: int(pos)
                }
                if dist != "0":
                    current["dist"] = int(dist)
                snps.append(current)
        return (snps, {self.MAPS_FORMAT_ATTR: "PLINK"})

    def read_samples(self, filename):
        COLUMN_NAMES = ["fid", self.SAMPLE_DICT_ID, "pid",
                        "mid", "sex", "aff"]
        samples = []
        with open(filename, "r") as f:
            for line in f:
                sample = dict()
                tokens = line.split()
                for (attr_name, attr) in zip(COLUMN_NAMES, tokens[0:6]):
                    if attr != "0":
                        sample[attr_name] = attr
                g = [{
                    self.SNPBLOCKS_SNP_GENOTYPE_INSIDE_LIST: tokens[i] + tokens[i+1]}
                    for i in range(6, len(tokens), 2)]
                sample[self.SAMPLE_DICT_GENOTYPE] = g
                samples.append(sample)
        return samples
