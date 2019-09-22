from AbstractImporter import AbstractImporter

class FinalReportImporter(AbstractImporter):

    def read_map(self, filename):
        snps = [] 
        with open(filename, "r") as f:
            while "[Data]" not in next(f):
                pass
            next(f)
            prev_sample_id = None
            for line in f:
                values = line.split()
                (snp_id, sample_id) = (values[0], values[1])
                if prev_sample_id not in (None, sample_id):
                    break
                prev_sample_id = sample_id
                snps.append({self.SNPS_NAME_ATTR: snp_id})
        return (snps, {self.MAPS_FORMAT_ATTR: "FR"}) 

    def read_samples(self, filename):
        samples = []
        with open(filename, "r") as f:
            while "[Data]" not in next(f):
                pass
            next(f)
            prev_sample_id = None
            sample, g = None, None
            for line in f: 
                (snp_id, sample_id, a1f, a2f, a1t,
                a2t, a1ab, a2ab, gc, x, y) = line.split()
                if prev_sample_id != sample_id:
                    if prev_sample_id is not None:
                        sample[self.SAMPLE_DICT_GENOTYPE] = g
                        samples.append(sample)
                    sample = {self.SAMPLE_DICT_ID: sample_id}
                    g = []
                g.append(
                    {self.SNPBLOCKS_SNP_GENOTYPE_INSIDE_LIST:
                        {"a1f": a1f,
                        "a2f": a2f,
                        "a1t": a1t,
                        "a2t": a2t,
                        "a1ab": a1ab,
                        "a2ab": a2ab},
                     "gc": float(gc),
                     "x": float(x),
                     "y": float(y)})
                prev_sample_id = sample_id
        sample[self.SAMPLE_DICT_GENOTYPE] = g
        samples.append(sample)
        return samples
