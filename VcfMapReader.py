from MapReader import MapReader
import re

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
