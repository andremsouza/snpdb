from SampleReader import SampleReader

class FinalReportSampleReader(SampleReader):
    
    def __iter__(self):
        with open(self._PED_FILE, "r") as f:
            while "[Data]" not in next(f):
                pass
            next(f)
            prev_sample_id = None
            sample = None
            g = {"a1f": [], "a2f": [], "a1t": [],
                 "a2t": [], "a1ab": [], "a2ab": [],
                 "gc": [], "x": [], "y": []}
            for line in f: 
                (snp_id, sample_id, a1f, a2f, a1t,
                a2t, a1ab, a2ab, gc, x, y) = line.split()
                if prev_sample_id != sample_id:
                    if prev_sample_id is not None:
                        sample[self.SAMPLE_GENOTYPE] = g
                        yield sample
                    sample = {self.SAMPLE_ID: sample_id}
                    g = {"a1f": [], "a2f": [], "a1t": [],
                         "a2t": [], "a1ab": [], "a2ab": [],
                         "gc": [], "x": [], "y": []}
                g["a1f"].append(a1f)
                g["a2f"].append(a2f)
                g["a1t"].append(a1t)
                g["a2t"].append(a2t)
                g["a1ab"].append(a1ab)
                g["a2ab"].append(a2ab)
                g["gc"].append(float(gc))
                g["x"].append(float(x))
                g["y"].append(float(y))
                prev_sample_id = sample_id
        sample[self.SAMPLE_GENOTYPE] = g
        yield sample
            
    def __len__(self):
        try:
            return self.__len
        except:
            pass
        with open(self._PED_FILE, "r") as f:
            while "[Data]" not in next(f):
                pass
            next(f)
            prev_sample_id = None
            cnt = 0
            for line in f: 
                sample_id = line.split()[1]
                if prev_sample_id != sample_id:
                    cnt += 1
                prev_sample_id = sample_id
        self.__len = cnt
        return cnt
