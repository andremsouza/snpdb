from snpdb import MapReader, SampleReader
import re

class Z125MapReader(MapReader):
    
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




class Z125SampleReader(SampleReader):
    
    def __iter__(self):
        with open(self._PED_FILE, "r") as f:
            next(f)
            for line in f:
                (id, genotype) = line.split()
                sample = {self.SAMPLE_ID: id,
                          self.SAMPLE_GENOTYPE: {"g": genotype}}
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




class VcfMapReader(MapReader):
    
    def __iter__(self):
        with open(self._MAP_FILE) as f:
            while next(f)[0:2] == "##":
                pass
            for line in f:
                values = line.split()
                (chrom, pos, snp_id, ref, alt, qual, fil, info) = values[0:8]
                snp = {
                    self.SNP_CHROM: chrom,
                    self.SNP_POS: int(pos),
                    "ref": ref,
                    "alt": alt.split(","),
                    "qual": float(qual),
                    "filter": fil,
                    "info": info
                }

                if snp_id != ".":
                    snp.update({self.SNP_NAME: snp_id})
                yield snp


    def __len__(self):
        with open(self._MAP_FILE) as f:
            while next(f)[0:2] == "##":
                pass
            for cnt, line in enumerate(f):
                pass
            return cnt+1

    def map_meta(self):
        meta = {self.MAP_FORMAT: "VCF"}
        meta_line = re.compile("##(.+?)=(.+)")
        key_value = re.compile("([^,]+?)=((?:\".*\")|(?:[^\"][^,]*)(?:,|$))")
        with open(self._MAP_FILE, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                match = meta_line.fullmatch(line)
                if match is None:
                    break
                left, right = match.group(1, 2)
                if right[0] == "<" and right[-1] == ">":
                    if left not in meta:
                        meta[left] = []
                    data = dict()
                    for (le, ri) in key_value.findall(right[1:-1]):
                      data[le] = ri.rstrip(",").rstrip("\"").lstrip("\"")
                    meta[left].append(data)
                else:
                    meta[left] = right
        return meta





class VcfSampleReader(SampleReader):

    def __iter__(self):
        with open(self._PED_FILE) as f:
            header = None
            for line in f:
                if line[0:2] != "##":
                    header = line.split()
                    break
            sample_strs = [header[i] + " " for i in range(9, len(header))]
            for line in f:
                tokens = line.split()
                if len(tokens) != len(header):
                    raise Exception("Missing column value.")
                for i in range(9, len(tokens)):
                    sample_strs[i-9] += tokens[i] + " "
        for i, sample_str in enumerate(sample_strs):
            snps = sample_str.split()
     
            gen = {"g": snps[1:]}
            sample = {self.SAMPLE_ID: snps[0],
                      self.SAMPLE_GENOTYPE: gen}
            yield sample


    def __len__(self):
        try:
            return self.__len
        except:
            pass
        with open(self._PED_FILE) as f:
            for line in f:
                if line[0:2] != "##":
                    self.__len = len(line.split())-9
                    break
        return self.__len
                

