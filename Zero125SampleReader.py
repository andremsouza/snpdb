from SampleReader import SampleReader

class Zero125SampleReader(SampleReader):
    
    def __iter__(self):
        with open(self._PED_FILE, "r") as f:
            next(f)
            for line in f:
                (id, genotype) = line.split()
                sample = {self.SAMPLE_ID: id,
                          self.SAMPLE_GENOTYPE: [{"g":genotype[i]}
                                           for i in range(len(genotype))]}
                yield sample
    
    def __len__(self):
        try:
            return self.__len
        except:
            pass
        with open(self._PED_FILE, "r") as f:
            for cnt, l in enumerate(f):
                pass
            self.__len = cnt
        return self.__len
