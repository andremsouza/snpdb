from MapReader import MapReader

class PlinkMapReader(MapReader):
    
    def __iter__(self):
        with open(self._MAP_FILE, "r") as f:
            for line in f:
                (chrom, id, dist, pos) = line.split()
                snp = {self.SNP_NAME: id,
                       self.SNP_CHROM: chrom,
                       self.SNP_POS: int(pos)}
                if dist != "0":
                    snp["dist"] = int(dist)
                yield snp

    def __len__(self):
        try:
            return self.__len
        except:
            pass
        with open(self._MAP_FILE, "r") as f:
            for cnt, l in enumerate(f):
                pass
            self.__len = cnt + 1
        return self.__len
    
    def map_meta(self):
        return {self.MAP_FORMAT: "PLINK"}
