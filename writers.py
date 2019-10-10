from snpdb import MapWriter


class Z125MapWriter(MapWriter):
    
    def write(self, out_file_path):
        with open(out_file_path, "w") as f:
            f.write("Name\tChromosome\tTEMPPOS\n")
            f.writelines((f"{snp[self.SNP_NAME]}\t" + 
                          f"{snp[self.SNP_CHROM]}\t" +
                          f"{snp[self.SNP_POS]}\n"
                          for snp in self._snps))
            f.close()




class PlinkMapWriter(MapWriter):
    
    def write(self, out_file_path):
        with open(out_file_path, "w") as f:
            f.writelines((f"{snp[self.SNP_CHROM]} " +
                          f"{snp[self.SNP_NAME]} " +
                          (f"{snp['dist']} " if "dist" in snp else "0 ") +
                          f"{snp[self.SNP_POS]}\n"
                          for snp in self._snps))
