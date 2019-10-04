from SampleReader import SampleReader

class PlinkSampleReader(SampleReader):
    
    def __iter__(self):
        COLUMN_NAMES = ["fid", self.SAMPLE_ID, "pid",
                        "mid", "sex", "aff"]
        with open(self._PED_FILE, "r") as f:
            for line in f:
                sample = dict()
                tokens = line.split()
                for (attr_name, attr) in zip(COLUMN_NAMES, tokens[0:6]):
                    if attr != "0":
                        sample[attr_name] = attr
                g = {"g": [tokens[i] + tokens[i+1]
                          for i in range(6, len(tokens), 2)]}
                sample[self.SAMPLE_GENOTYPE] = g
                yield sample
            
    def __len__(self):
        try:
            return self.__len
        except:
            pass
        with open(self._PED_FILE, "r") as f:
            for cnt, l in enumerate(f):
                pass
            self.__len = cnt+1
        return self.__len
