from snpdb import MapWriter, SampleWriter


class Z125MapWriter(MapWriter):
    def write(self, out_file_path):
        with open(out_file_path, "w") as f:
            f.write("Name\tChromosome\tTEMPPOS\n")
            f.writelines(
                (
                    f"{snp[self.SNP_NAME]}\t"
                    + f"{snp[self.SNP_CHROM]}\t"
                    + f"{snp[self.SNP_POS]}\n"
                    for snp in self._snps
                )
            )
            f.close()


class PlinkMapWriter(MapWriter):
    def write(self, out_file_path):
        with open(out_file_path, "w") as f:
            f.writelines(
                (
                    f"{snp[self.SNP_CHROM]} "
                    + f"{snp[self.SNP_NAME]} "
                    + (f"{snp['dist']} " if "dist" in snp else "0 ")
                    + f"{snp[self.SNP_POS]}\n"
                    for snp in self._snps
                )
            )


class Z125SampleWriter(SampleWriter):
    def write(self, out_file_path):
        with open(out_file_path, "w") as f:
            f.write("ID CALL...\n")
            f.writelines(
                (
                    f"{sample[self.SAMPLE_ID]} "
                    + f"{''.join(sample[self.SAMPLE_GENOTYPE]['g'])}\n"
                    for sample in self._samples
                )
            )


class PlinkSampleWriter(SampleWriter):
    def write(self, out_file_path):
        COLUMN_NAMES = ["fid", self.SAMPLE_ID, "pid", "mid", "sex", "aff"]
        with open(out_file_path, "w") as f:
            f.writelines(
                (
                    (f"{sample['fid']} " if "fid" in sample else "0 ")
                    + f"{sample[self.SAMPLE_ID]} "
                    + (f"{sample['pid']} " if "pid" in sample else "0 ")
                    + (f"{sample['mid']} " if "mid" in sample else "0 ")
                    + (f"{sample['sex']} " if "sex" in sample else "0 ")
                    + (f"{sample['aff']} " if "aff" in sample else "0 ")
                    + f"{self.__fmt(sample[self.SAMPLE_GENOTYPE]['g'])}\n"
                    for sample in self._samples
                )
            )

    def __fmt(self, g):
        return " ".join((c for c in "".join(g)))
