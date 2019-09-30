from MapReader import MapReader

class FinalReportMapReader(MapReader):
    
    def __iter__(self):
        with open(self._MAP_FILE, "r") as f:
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
                yield {self.SNP_NAME: snp_id}

    def __len__(self):
        try:
            return self.__len
        except:
            pass
        with open(self._MAP_FILE, "r") as f:
            while "[Data]" not in next(f):
                pass
            next(f)
            cnt = 0
            prev_sample_id = None
            for line in f:
                values = line.split()
                (snp_id, sample_id) = (values[0], values[1])
                if prev_sample_id not in (None, sample_id):
                    break
                prev_sample_id = sample_id
                cnt += 1
            self.__len = cnt
        return self.__len
    
    def map_meta(self):
        return {self.MAP_FORMAT: "FR"}
