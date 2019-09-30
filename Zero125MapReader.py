from MapReader import MapReader

class Zero125MapReader(MapReader):
    
    def __iter__(self):
        with open(self._MAP_FILE, "r") as f:
            next(f)
            for line in f:
                (name, chrom, pos) = line.split()
                snp = {self.SNP_NAME: name,
                       self.SNP_CHROM: chrom,
                       self.SNP_POS: int(pos)}
                yield snp

    def __len__(self):
        try:
            return self.__len
        except:
            pass
        with open(self._MAP_FILE, "r") as f:
            for cnt, l in enumerate(f):
                pass
            self.__len = cnt
        return self.__len
    
    def map_meta(self):
        return {self.MAP_FORMAT: "0125"}
