from AbstractImporter import AbstractImporter

class Zero125Importer(AbstractImporter):
    
    def read_map(self, filename):
        snps = []
        with open(filename, "r") as f:
            next(f)
            for line in f:
                (id, chrom, pos) = line.split()
                current = {
                    self.SNPS_NAME_ATTR: id,
                    self.SNPS_CHROMOSOME_ATTR: chrom,
                    self.SNPS_POSITION_ATTR: pos
                }
                snps.append(current)
        return (snps, {self.MAPS_FORMAT_ATTR: "0125"})

    def read_samples(self, filename):
        samples = []
        with open(filename, "r") as f:
            next(f)
            for line in f:
                (id, genotype) = line.split()
                g = [{
                    self.SNPBLOCKS_SNP_GENOTYPE_INSIDE_LIST: genotype[i]}
                    for i in range(len(genotype))]
                samples.append(
                    {self.SAMPLE_DICT_ID: id,
                    self.SAMPLE_DICT_GENOTYPE: g})
        return samples
